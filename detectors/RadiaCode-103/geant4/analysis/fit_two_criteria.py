# -*- coding: utf-8 -*-
"""
Нуклидное разложение измеренного гамма-фона детектора RadiaCode-103 по расчётным шаблонам.
Два независимых разложения (A — критерий chi2/ndf, пуассоновские веса; B — критерий невязки формы, относительные веса),
цепочки разбиты на изотопы, шаблоны только из новой GDML-модели (run_field/output/).
"""
import math, os, sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
from scipy.optimize import nnls

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "detectors", "RadiaCode-103", "analysis"))
sys.path.insert(0, os.path.join(REPO, "common", "py"))

import rcspec
import read_rcxml
import fit_lines as fl

# G4MODELS_TEMPLATE_DIR — подмена каталога шаблонов, нужна для мутационной
# проверки (#SA-3): скрытая порча шаблона обязана менять результат.
TEMPLATE_DIR = os.environ.get(
    "G4MODELS_TEMPLATE_DIR", os.path.join(HERE, "..", "run_field", "output"))
# Дефолт — шаблоны ДЕЙСТВУЮЩЕЙ модели: поле ЕРН в реальной комнате
# (run_roomfield) -> отклик GDML-модели прибора. Прежний дефолт
# rc103_field_m1_%s_1e8.csv считался на флюенсе results/wallion/, то есть на
# старой линии; она вместе со всеми своими результатами убрана в архив по
# решению оператора 27.08.2026 и не используется.
TEMPLATE_FMT = os.environ.get("G4MODELS_TEMPLATE_FMT",
                              "rc103_field_room_%s.csv")
NUCS = ["K40", "Ra226", "Pb214", "Bi214", "Pb212", "Ac228", "Bi212", "Tl208"]
# Мюонный столбец. ЕДИНИЦЫ ДРУГИЕ: не Бк/кг, а «мюонов в секунду через диск
# источника» — шаблон нормирован на ОДИН мюон (колонка per_muon), амплитуду
# подбирает NNLS наравне с активностями. Радиус 300 мм принят рабочим по
# проверке насыщения: при 150 мм абсолютная светосила занижена на 21 %
# (обрезаются наклонные треки), форма при 300 и 600 мм совпадает.
MUON_CSV = os.environ.get(
    "G4MODELS_MUON_CSV",
    os.path.join(HERE, "..", "run_muon", "output", "rc103_muon_r300_1e8.csv"))
MUON_PDG_PER_S = 47.21813758   # 0.0167 см^-2 c^-1 * площадь диска, справочно

MEAS_DIR = os.environ.get("G4MODELS_MEASURED", r"C:\g4work\measured\RadiaCode-103")
MEAS_NAME = "\u0424\u043e\u043d 7 \u0434\u043d\u0435\u0439 \u0431\u0435\u0437 \u0434\u043e\u043c\u0438\u043a\u0430.xml"
CAL_ROOM = [-3.711311, 2.444318, 0.000321]
E_LO, E_HI = 20.0, 2830.0
BANDS = [(20,100),(100,300),(300,700),(700,1500),(1500,2000),(2000,2400),(2400,2830)]
EQUILIBRIUM = [
    ("Bi214","Pb214",1.0,"радон: поздние/ранние ДПР"),
    ("Pb214","Ra226",1.0,"радон: эманирование Ra-226"),
    ("Tl208","Ac228",0.3594,"торий: ветвление Bi-212 -> Tl-208"),
    ("Bi212","Pb212",1.0,"торий: Bi-212 / Pb-212"),
    ("Ac228","Pb212",1.0,"торий: Ac-228 / Pb-212")
]

def read_template(path):
    with open(path, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]
    
    meta = {}
    i = 0
    while i < len(lines) and not lines[i].startswith("bin_keV"):
        if lines[i] == "metric,value":
            i += 1
            while i < len(lines) and not lines[i].startswith("bin_keV"):
                key, val = lines[i].split(",", 1)
                try:
                    meta[key] = float(val)
                except ValueError:
                    meta[key] = val
                i += 1
        else:
            i += 1

    if i >= len(lines):
        raise SystemExit(f"Файл {path} не содержит данных")

    # ВАЖНО: lines[i] — это строка-ЗАГОЛОВОК "bin_keV,counts,cps", данные с i+1.
    # Размер массива берём по ЧИСЛОВОМУ максимуму бина: max() по строкам дал бы
    # лексикографический максимум ("999" > "1000") и обрезал бы гистограмму.
    rows = []
    for line in lines[i + 1:]:
        parts = line.split(",")
        if len(parts) < 3:
            continue
        # (бин, cps/per_muon, сырые отсчёты МК) — последние нужны для оценки
        # статистической ошибки самого шаблона: для невзвешенного МК
        # sumw2 = counts, дисперсия в единицах cps есть counts / T_run^2.
        rows.append((int(float(parts[0])), float(parts[2]), float(parts[1])))
    if not rows:
        raise SystemExit(f"Файл {path} не содержит гистограммы")

    nbin = max(b for b, _, _ in rows) + 1
    arr = np.zeros(nbin)
    cnt_mc = np.zeros(nbin)
    for b, cps, c in rows:
        arr[b] = cps
        cnt_mc[b] = c

    if not np.any(arr):
        raise SystemExit(f"Файл {path} пустой")

    return (meta, arr, cnt_mc)

def template_variance(cnt_mc, t_run):
    """Дисперсия шаблона в единицах cps^2, по конвенции контура.

    Для невзвешенного МК sumw2 = counts. Дисперсия сворачивается тем же ядром,
    что и значения (а не квадратом ядра) — принятое в проекте приближение,
    см. analysis/run_bg_shield.py:216. Оно КОНСЕРВАТИВНО, слегка завышает
    ошибку, поскольку sum(w) = 1 больше, чем sum(w^2).
    """
    if t_run <= 0:
        return np.zeros_like(cnt_mc)
    return rcspec.fold(cnt_mc, "103") / (t_run * t_run)


def load_templates():
    missing = []
    names = []
    cols = []
    metas = []
    varis = []          # дисперсии шаблонов, cps^2 на канал

    for nuc in NUCS:
        path = os.path.join(TEMPLATE_DIR, TEMPLATE_FMT % nuc)
        if not os.path.exists(path):
            missing.append(path)
            continue
        meta, arr, cnt_mc = read_template(path)
        folded = rcspec.fold(arr, "103")
        names.append(nuc)
        cols.append(folded)
        varis.append(template_variance(cnt_mc, float(meta.get("t_run_s", 0.0))))
        metas.append(meta)

    # Мюонный столбец — последним. Формат CSV тот же (третья колонка per_muon),
    # поэтому читается тем же ридером.
    if os.path.exists(MUON_CSV):
        meta_mu, raw_mu, cnt_mu = read_template(MUON_CSV)
        names.append("mu")
        cols.append(rcspec.fold(raw_mu, "103"))
        # У мюонов шаблон нормирован на ОДИН мюон (per_muon = counts/n_events),
        # поэтому роль T_run играет само число разыгранных мюонов.
        varis.append(template_variance(cnt_mu, float(meta_mu.get("n_events", 0.0))))
        metas.append(meta_mu)
    else:
        print("[--] мюонный шаблон НЕ найден: %s" % os.path.abspath(MUON_CSV))

    if not names:
        raise SystemExit("Не найдено ни одного шаблона. Запустите run_field по спеке run_field/_spec_run_field.md")

    for i, (name, meta) in enumerate(zip(names, metas)):
        print("[%-6s] hits %10s  cps_total %12s  T_run %12s"
              % (name, meta.get("n_hits_in_crystal", "?"),
                 meta.get("cps_total", "?"), meta.get("t_run_s", "?")))

    if missing:
        print("Отсутствующие шаблоны:")
        for path in missing:
            print(path)

    return (names, cols, metas, varis)

def degeneracy_report(A, names):
    norm = A / np.maximum(np.linalg.norm(A, axis=0), 1e-300)
    C = norm.T @ norm
    cond = np.linalg.svd(norm, compute_uv=False)[0] / np.linalg.svd(norm, compute_uv=False)[-1]
    print("Диагностика вырожденности:")
    print("=" * 68)
    print(" ".join(f"{name:>8}" for name in names))
    for i, row in enumerate(C):
        print(" ".join(f"{row[j]:>8.3f}" if j != i else f"{row[j]:>8.3f}" for j in range(len(row))))
    print(f"cond = {cond:.2f}")
    if cond > 100:
        print("Предупреждение: часть амплитуд не определится из-за высокой вырожденности")
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            if C[i,j] > 0.995:
                print(f"Предупреждение: звенья {names[i]} и {names[j]} плохо разделены (C[{i},{j}] = {C[i,j]:.3f})")
    return cond

def metrics(pred_counts, meas_counts, n_par):
    var = np.maximum(meas_counts, 1.0)
    chi2 = sum((pred_counts - meas_counts)**2 / var)
    ndf = max(len(meas_counts) - n_par, 1)
    chi2_ndf = chi2 / ndf
    pred_norm = pred_counts.sum()
    meas_norm = meas_counts.sum()
    if pred_norm <= 0 or meas_norm <= 0:
        shape = float("nan")
    else:
        shape = 0.5 * sum(abs(pred_counts/pred_norm - meas_counts/meas_norm))
    return (chi2_ndf, shape)

def fit(A_counts, meas_counts, weights, names, title, note, var_counts=None):
    """var_counts — дисперсии шаблонов [канал x звено] в отсчётах^2 на 1 Бк/кг.

    Если она передана, полная дисперсия канала = дисперсия измерения (она
    заложена в исходных weights как 1/sigma^2) ПЛЮС дисперсия модели
    sum_k a_k^2 * var_k. Последняя зависит от искомых амплитуд, поэтому веса
    пересчитываются итеративно. Без этого ошибки амплитуд занижены: подгонка
    считает шаблоны бесшумными, а в них 10^4 событий на 1341 канал.
    """
    w = weights
    amp = None
    for _ in range(6 if var_counts is not None else 1):
        Aw = A_counts * w[:, None]
        amp, _ = nnls(Aw, meas_counts * w)
        if var_counts is None:
            break
        model_var = var_counts @ (amp ** 2)
        base_var = 1.0 / np.maximum(weights, 1e-300) ** 2
        w = 1.0 / np.sqrt(base_var + model_var)
    Aw = A_counts * w[:, None]
    yw = meas_counts * w
    # Ковариация может не обратиться при вырожденном базисе — тогда амплитуды
    # ЕСТЬ, а ошибок нет. Ловить обе в одном try нельзя: amp затирался бы nan.
    try:
        cov = np.linalg.inv(Aw.T @ Aw)
        sd = np.sqrt(np.maximum(np.diag(cov), 0))
    except np.linalg.LinAlgError:
        sd = np.full(len(amp), np.nan)

    pred = A_counts @ amp
    chi2ndf, shape = metrics(pred, meas_counts, len(amp))

    print("=" * 68)
    print(f"РАЗЛОЖЕНИЕ {title}")
    print(note)
    print("=" * 68)
    for i, (name, a, s) in enumerate(zip(names, amp, sd)):
        rel_err = (s / a * 100) if a != 0 else float("inf")
        mark = ""
        if a == 0:
            mark = " <-- обнулена NNLS"
        elif s > a:
            mark = " <-- НЕ ОПРЕДЕЛЕНА (ошибка больше значения)"
        unit = "мюон/с" if name == "mu" else "Бк/кг"
        print(f"{name:>8} | {a:10.3f} | {s:10.3f} | {rel_err:6.1f}% {unit}{mark}")
        if name == "mu" and a > 0:
            print(f"{'':>8} |   ^ против PDG {MUON_PDG_PER_S:.1f} мюон/с "
                  f"-> отношение {a / MUON_PDG_PER_S:.2f} "
                  f"(справочная сверка порядка, НЕ ограничение подгонки)")

    print(f"chi2/ndf = {chi2ndf:.3f}")
    print(f"невязка формы = {shape:.4f} (0 = формы совпали)")

    return (amp, sd, pred, chi2ndf, shape)

def bands_report(pred, meas, e_meas, live):
    print("Таблица по полосам:")
    print("=" * 68)
    print("полоса,кэВ     | измерено   | модель     | м/и")
    for lo, hi in BANDS:
        sel = (e_meas >= lo) & (e_meas < hi)
        m = meas[sel].sum() / live
        p = pred[sel].sum() / live
        ratio = p/m if m != 0 else float("nan")
        print(f"{lo:4d}-{hi:4d} кэВ | {m:10.2f} | {p:10.2f} | {ratio:6.3f}")
    sel = (e_meas >= E_LO) & (e_meas < E_HI)
    m = meas[sel].sum() / live
    p = pred[sel].sum() / live
    ratio = p/m if m != 0 else float("nan")
    print(f"полный        | {m:10.2f} | {p:10.2f} | {ratio:6.3f}")

def equilibrium_report(names, amp, sd, title):
    print(f"\nПроверка равновесия ({title}):")
    print("=" * 68)
    found = False
    for a_name, b_name, expected, note in EQUILIBRIUM:
        # Звена может не быть вовсе (шаблон не прогнан) — .index() бросил бы
        # ValueError и уронил отчёт на неполном наборе шаблонов.
        if a_name not in names or b_name not in names:
            continue
        i_a = names.index(a_name)
        i_b = names.index(b_name)
        if amp[i_a] > 0 and amp[i_b] > 0:
            found = True
            r = amp[i_a] / amp[i_b]
            er = r * np.sqrt((sd[i_a]/amp[i_a])**2 + (sd[i_b]/amp[i_b])**2)
            # Отклонение в сигмах — без него «1,26 против 1,00» невозможно
            # отличить от «в пределах ошибки» и от реального нарушения.
            nsig = abs(r - expected) / er if er > 0 else float("nan")
            if not (nsig == nsig):
                verdict = "не определено"
            elif er > expected:
                verdict = "не определено (ошибка больше ожидаемого)"
            elif nsig < 2.0:
                verdict = "равновесие"
            elif nsig < 3.0:
                verdict = "на границе"
            else:
                verdict = "НАРУШЕНО"
            print("    %-6s/%-6s = %7.3f +- %6.3f   ожид %.4f   отн %5.2f   "
                  "%5.1f sigma  -> %-12s (%s)"
                  % (a_name, b_name, r, er, expected, r / expected, nsig,
                     verdict, note))
    if not found:
        print("Равновесие не проверялось — ни одна пара амплитуд положительна")

def main():
    meas_path = os.path.join(MEAS_DIR, MEAS_NAME)
    if not os.path.exists(meas_path):
        raise SystemExit(f"Файл {meas_path} не найден. Это личные измерения оператора, в репозиторий не входят. Укажите переменную окружения G4MODELS_MEASURED")

    smp = read_rcxml.read(meas_path)[0]
    cnt = smp.counts[:-1].astype(float)
    ch = np.arange(len(cnt))
    e_meas = np.asarray(sum(c * ch**i for i, c in enumerate(CAL_ROOM)))
    live = float(smp.live)
    print(f"Имя файла: {MEAS_NAME}")
    print(f"Живое время: {live/3600:.2f} ч ({live:.0f} сек)")
    print(f"Всего отсчётов: {cnt.sum():.0f}")

    names, cols, metas, varis = load_templates()

    # Перевод модели на реальную шкалу
    A = np.zeros((len(e_meas), len(cols)))
    VAR = np.zeros((len(e_meas), len(cols)))
    for k, c in enumerate(cols):
        A[:, k] = fl.rebin_model_to_meas(np.arange(len(c)) + 0.5, c, e_meas)
        # Ремешок линеен и суммирует, поэтому применим к дисперсиям напрямую.
        v = varis[k]
        VAR[:, k] = fl.rebin_model_to_meas(np.arange(len(v)) + 0.5, v, e_meas)

    sel = (e_meas >= E_LO) & (e_meas < E_HI)
    A = A[sel]
    VAR = VAR[sel]
    meas_counts = cnt[sel]
    e_meas = e_meas[sel]

    print(f"Число каналов в подгонке: {len(meas_counts)}")
    print(f"Диапазон: {E_LO:.1f} - {E_HI:.1f} кэВ")

    cond = degeneracy_report(A, names)
    A_counts = A * live

    VAR_counts = VAR * live * live   # отсчёты^2 на 1 Бк/кг

    # Разложение A
    weights_a = 1.0 / np.sqrt(np.maximum(meas_counts, 1.0))
    amp_a, sd_a, pred_a, chi2ndf_a, shape_a = fit(
        A_counts, meas_counts, weights_a, names,
        "A — критерий chi2/ndf",
        "веса пуассоновские sigma_i = sqrt(N_i) ПЛЮС дисперсия шаблонов "
        "(итеративно), мягкая область 20-300 кэВ несёт около 91 % отсчётов",
        var_counts=VAR_counts
    )
    bands_report(pred_a, meas_counts, e_meas, live)
    equilibrium_report(names, amp_a, sd_a, "A")

    # Разложение B
    weights_b = 1.0 / np.maximum(meas_counts, 1.0)
    amp_b, sd_b, pred_b, chi2ndf_b, shape_b = fit(
        A_counts, meas_counts, weights_b, names,
        "B — критерий невязки формы",
        "веса относительные sigma_i = N_i ПЛЮС дисперсия шаблонов (итеративно), "
        "отклонение в 1 % весит одинаково в мягкой и жёсткой части",
        var_counts=VAR_counts
    )
    bands_report(pred_b, meas_counts, e_meas, live)
    equilibrium_report(names, amp_b, sd_b, "B")

    # Сравнение
    print("\nСРАВНЕНИЕ ДВУХ РАЗЛОЖЕНИЙ")
    print("=" * 68)
    print("метрика      | A        | B        | лучшее")
    print("-" * 68)
    print(f"chi2/ndf     | {chi2ndf_a:.3f} | {chi2ndf_b:.3f} | {'A' if chi2ndf_a < chi2ndf_b else 'B'}")
    print(f"форма        | {shape_a:.4f} | {shape_b:.4f} | {'A' if shape_a < shape_b else 'B'}")

    print("\nАмплитуды:")
    print("звено     | A, Бк/кг   | B, Бк/кг   | B/A")
    print("-" * 68)
    for i, name in enumerate(names):
        a = amp_a[i]
        b = amp_b[i]
        ratio = b/a if a != 0 else float("nan")
        print(f"{name:>8} | {a:10.3f} | {b:10.3f} | {ratio:6.3f}")

    # Требует толкования
    print("\nТРЕБУЕТ ТОЛКОВАНИЯ")
    print("=" * 68)
    trigger = False

    for i, name in enumerate(names):
        if amp_a[i] == 0 or amp_b[i] == 0:
            print(f"Звено {name} обнулено в одном из разложений — вырожденность или физическое отсутствие?")
            trigger = True

    for i, name in enumerate(names):
        a = amp_a[i]
        b = amp_b[i]
        # «Более чем вдвое в ЛЮБУЮ сторону»: и рост, и падение.
        if a > 0 and b > 0 and (b / a > 2.0 or b / a < 0.5):
            print(f"Амплитуды звена {name} расходятся более чем вдвое "
                  f"(A={a:.2f}, B={b:.2f}) — какой канал тянет?")
            trigger = True

    if max(chi2ndf_a, chi2ndf_b) > 10:
        print("Одно из разложений имеет chi2/ndf > 10 — модель формы не описывает либо ошибки занижены")
        trigger = True

    if not trigger:
        print("(пусто — ни один триггер не сработал)")

if __name__ == "__main__":
    main()
