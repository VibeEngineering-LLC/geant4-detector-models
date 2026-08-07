# -*- coding: utf-8 -*-
"""ПРЯМАЯ задача: спектр пачки WT-20 по составу, без подгонки под измерение.

Постановка обратная той, что решает `wt20_unfold.py`. Там активности искались
из измеренного спектра; здесь они ВЫЧИСЛЯЮТСЯ из состава по этикетке и закона
радиоактивного накопления, а спектр строится вперёд. Ни один параметр под
измеренный спектр не настраивается: измерение появляется только в конце, для
сравнения.

    активность звена  ->  спектр звена  ->  сумма  ->  сравнение с замером

Ключевое допущение о состоянии ряда — торий в электродах ХИМИЧЕСКИ ЧИСТЫЙ, и
дочерние далеко не в равновесии (постановка оператора 07.08.2026). Химическая
очистка отделяет ЭЛЕМЕНТЫ, а не изотопы, и это определяет начальные условия:

  * Th-232 и Th-228 — оба торий. Разделить их химией нельзя, поэтому Th-228
    остаётся в металле полностью, в той активности, что была на момент
    очистки (в исходном сырье ряд в вековом равновесии, значит A(Th-228) = A0).
  * Ra-228, Ac-228 — радий и актиний, другие элементы: удаляются нацело.
  * Ra-224, Rn-220, Po-216, Pb-212, Bi-212, Tl-208, Po-212 — радий, радон,
    полоний, свинец, висмут, таллий: удаляются нацело.

Дальше ряд восстанавливается сам, каждое звено по своему периоду. Верхняя
ветвь (Ra-228, T½ = 5,75 года) нарастает медленно; нижняя (всё ниже Th-228)
следует за Th-228, у которого T½ = 1,91 года, поэтому в первые годы после
очистки нижняя часть ряда СИЛЬНЕЕ верхней — ровно то соотношение, которое даёт
разложение измеренного спектра.

Спектр строится ДВУМЯ независимыми способами (методики §2 и §4 приложения):

  метод 1 — понуклидные Монте-Карло-шаблоны: розыгрыш распада целого звена
            в геометрии замера, каскад и суммирование входят переносом;
  метод 2 — матрица отклика на моноэнергетический квант, свёрнутая с
            библиотечными выходами линий МАГАТЭ; суммирования каскада не даёт.

Расхождение методов — мера вклада суммирования и полноты библиотеки, а не
свободный параметр.

    python analysis/wt20_forward.py <спектр.xml> <шаблоны> <матрица> [выход]
"""
import argparse
import csv
import io
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.environ.get("SPECTRAVIBE_ROOT")
if not _ROOT:
    raise SystemExit("не задана переменная окружения SPECTRAVIBE_ROOT")
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
sys.path.insert(0, _HERE)

from wt20_unfold import (E_MAX, E_STEP, LIB, broaden,                 # noqa: E402
                         fwhm_from_file, read_atomspectra_xml,
                         read_correction, read_template, rebin_to_grid,
                         SPEC_ACT_TH232)

# Ряд тория от Th-232 до стабильного Pb-208. Порядок — по цепочке; поле
# `element` нужно для начальных условий: химическая очистка удаляет элемент
# целиком, поэтому два изотопа тория ведут себя одинаково, а всё остальное
# обнуляется.
#
# Периоды полураспада читаются из библиотеки МАГАТЭ при запуске (поле
# `half_life_sec` файлов линий), а не пишутся сюда числами.
CHAIN = [
    # ключ,   подпись,   элемент, файл линий,   ветвление от родителя
    ("Th232", "Th-232",  "Th", "232th", 1.0),
    ("Ra228", "Ra-228",  "Ra", "228ra", 1.0),
    ("Ac228", "Ac-228",  "Ac", "228ac", 1.0),
    ("Th228", "Th-228",  "Th", "228th", 1.0),
    ("Ra224", "Ra-224",  "Ra", "224ra", 1.0),
    ("Rn220", "Rn-220",  "Rn", "220rn", 1.0),
    ("Po216", "Po-216",  "Po", None,    1.0),
    ("Pb212", "Pb-212",  "Pb", "212pb", 1.0),
    ("Bi212", "Bi-212",  "Bi", "212bi", 1.0),
    ("Tl208", "Tl-208",  "Tl", "208tl", None),   # ветвление из библиотеки
    ("Po212", "Po-212",  "Po", None,    None),   # дополнение до единицы
]

# Периоды звеньев, у которых нет своего файла линий: у Po-216 и Po-212
# гамма-выхода практически нет, и выгрузка МАГАТЭ по ним пуста.
#
# Значения взяты из NUBASE2020 — оценённого свода ядерных данных
# (Kondev F.G. et al., Chinese Physics C 45 (2021) 030001, open access;
# сырой файл `references/data/nubase_3.mas20.txt` в библиотеке SpectraVibe,
# выгрузка www-nds.iaea.org/amdc/ame2020/).
#
# Прежде здесь стояли 0,148 с и 2,99e-7 с, внесённые ПО ПАМЯТИ с пометкой
# «ENSDF»: сверка с NUBASE2020 показала расхождение 2,8 % и 1,9 %. На
# результат оно не влияет (оба звена без гамма-выхода и в равновесии с
# родителем), но вспомненное числом в коде стоять не должно.
HL_EXTRA_SEC = {
    "Po216": 0.1440,     # NUBASE2020: 144,0 мс
    "Po212": 2.944e-7,   # NUBASE2020: 294,4 нс
}

YEAR_SEC = 365.25 * 86400.0


def M2_chain_files():
    """Звенья с файлами линий — вход метода 2 (свёртка матрицы с выходами)."""
    return [(k, l, e, f, b) for k, l, e, f, b in CHAIN if f is not None]


def half_lives():
    """Периоды полураспада звеньев, с. Из библиотеки МАГАТЭ, где она есть."""
    out = {}
    for key, lab, _el, fn, _br in CHAIN:
        if fn is None:
            out[key] = HL_EXTRA_SEC[key]
            continue
        p = os.path.join(LIB, "%s_gammas.csv" % fn)
        r = next(csv.DictReader(io.open(p, encoding="utf-8")), None)
        if not r or not r.get("half_life_sec"):
            raise SystemExit("в %s нет поля half_life_sec" % p)
        out[key] = float(r["half_life_sec"])
    return out


def branch_tl208():
    """Доля распадов Bi-212 по альфа-ветви (на Tl-208), из библиотеки МАГАТЭ."""
    p = os.path.join(LIB, "212bi_gammas.csv")
    for r in csv.DictReader(io.open(p, encoding="utf-8")):
        d, pc = r.get("decay"), r.get("decay_%")
        if d and d.strip().upper().startswith("A") and pc:
            try:
                return float(pc) / 100.0
            except ValueError:
                pass
    raise SystemExit("в %s нет альфа-ветви Bi-212" % p)


def activities_after_purification(t_sec, a0, hl, br_tl):
    """Активности звеньев ряда через t после химической очистки тория, Бк.

    Начальные условия (t = 0, момент очистки):
        A(Th-232) = a0      — тория ничем не тронули;
        A(Th-228) = a0      — тот же элемент, изотопы химией не делятся,
                              а в исходном сырье ряд был в равновесии;
        всё остальное = 0   — другие элементы удалены нацело.

    Решается последовательным интегрированием уравнений цепочки. Th-232 с
    периодом 1,4·10¹⁰ лет на любом мыслимом сроке хранения постоянен, поэтому
    его активность держится равной a0, а звенья ниже считаются по формуле
    двухчленного накопления от своего родителя.

    Для звеньев, у которых период много меньше родительского (все ниже
    Th-228: сутки, часы, минуты, секунды), решение выходит на вековое
    равновесие с родителем за считанные периоды, и через год после очистки
    A_k = A(родителя) с точностью много лучше процента. Это учтено точной
    формулой, а не приближением.
    """
    lam = {k: math.log(2.0) / hl[k] for k in hl}
    A = {}
    A["Th232"] = a0                                   # постоянен

    # --- верхняя ветвь: Th-232 -> Ra-228 -> Ac-228 ------------------------
    # Ra-228 нарастает из нуля к равновесию с Th-232 с периодом 5,75 года
    A["Ra228"] = a0 * (1.0 - math.exp(-lam["Ra228"] * t_sec))
    # Ac-228 (6,15 ч) следует за Ra-228 практически мгновенно
    lr, la = lam["Ra228"], lam["Ac228"]
    A["Ac228"] = a0 * (1.0 + (lr * math.exp(-la * t_sec)
                              - la * math.exp(-lr * t_sec)) / (la - lr))

    # --- Th-228: остаётся от очистки И подпитывается от Ra-228 ------------
    # Собственный распад начальной активности плюс накопление от растущего
    # Ra-228 (через Ac-228, чей период пренебрежимо мал против остальных).
    lt = lam["Th228"]
    decay_initial = a0 * math.exp(-lt * t_sec)
    # вклад от Ra-228: решение dN/dt = lam_Ra*N_Ra(t) - lam_Th*N_Th,
    # где A_Ra(t) = a0*(1 - exp(-lam_Ra*t))
    ingrowth = a0 * (1.0 - lt / (lt - lr) * math.exp(-lr * t_sec)
                     + lr / (lt - lr) * math.exp(-lt * t_sec))
    A["Th228"] = decay_initial + ingrowth

    # --- нижняя часть: всё ниже Th-228 следует за ним ---------------------
    # Периоды: Ra-224 3,63 сут, Rn-220 55,6 с, Po-216 0,148 с, Pb-212 10,6 ч,
    # Bi-212 60,6 мин, Tl-208 3,05 мин, Po-212 299 нс — все много меньше
    # периода Th-228 (1,91 года). Для такой пары решение цепочки за время,
    # много большее периода дочернего, выходит на ПОДВИЖНОЕ равновесие:
    #
    #     A_k = A_род · λ_k / (λ_k − λ_Th228),
    #
    # то есть активность дочернего превышает родительскую в λ_k/(λ_k − λ_Th)
    # раз. При λ_k ≫ λ_Th множитель равен единице с точностью λ_Th/λ_k: для
    # самого медленного из этих звеньев, Ra-224, превышение составляет
    # 0,52 %, для Pb-212 — 0,006 %. Множитель считается явно, чтобы величина
    # была верной, а не «примерно верной».
    #
    # Начальный переход (первые сутки-двое после очистки, пока нижняя часть
    # набирается от нуля) здесь не воспроизводится: возраст электродов —
    # годы, и множитель (1 − exp(−λ_k·t)) отличается от единицы на величину
    # порядка exp(−λ_k·t), то есть на десятки порядков ниже точности.
    for key in ("Ra224", "Rn220", "Po216", "Pb212", "Bi212"):
        A[key] = A["Th228"] * lam[key] / (lam[key] - lt)
    # Ветвление Bi-212: альфа-ветвь на Tl-208, бета-ветвь на Po-212
    A["Tl208"] = A["Bi212"] * br_tl
    A["Po212"] = A["Bi212"] * (1.0 - br_tl)
    return A


def load_templates(tdir, centres, fwhm_fn):
    """МК-шаблоны нуклидов, уширенные аппаратной линией. -> {ключ: спектр}."""
    out, head0 = {}, {}
    for key, lab, _el, _fn, _br in CHAIN:
        p = os.path.join(tdir, "%s.csv" % key)
        if not os.path.exists(p):
            continue
        head, e, c = read_template(p)
        head0 = head0 or head
        n = float(head["N_primaries"])
        if c.sum() <= 0:
            continue
        out[key] = broaden(e, c / n, centres, 0.0, 0.0, fwhm_fn=fwhm_fn)
    return out, head0


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("xml", help="файл замера (для сетки, времени и сравнения)")
    ap.add_argument("tdir", help="каталог МК-шаблонов нуклидов")
    ap.add_argument("mdir", nargs="?", help="каталог матрицы отклика (метод 2)")
    ap.add_argument("-o", "--out", help="каталог вывода")
    ap.add_argument("-t", "--years", default="1,3,5,10,20,50",
                    help="возраст после очистки, лет (через запятую)")
    args = ap.parse_args()
    outdir = args.out or os.path.dirname(args.xml)
    os.makedirs(outdir, exist_ok=True)

    spec = read_atomspectra_xml(args.xml)
    if isinstance(spec, (list, tuple)):
        spec = spec[0]
    bg = getattr(spec, "background_embedded", None)
    t_smp = float(spec.real_time)

    corr = read_correction(os.path.join(outdir, "calibration_fitted.csv"))
    edges = np.arange(0.0, E_MAX + E_STEP, E_STEP)
    centres = 0.5 * (edges[:-1] + edges[1:])
    y = rebin_to_grid(np.asarray(spec.counts, float), list(spec.energy_cal),
                      corr.get("sample"), edges)
    if bg is not None:
        ybg = rebin_to_grid(np.asarray(bg.counts, float), list(bg.energy_cal),
                            corr.get("background"), edges) \
            * (t_smp / float(bg.real_time))
    else:
        ybg = np.zeros_like(y)

    fwhm_fn = fwhm_from_file(args.xml)
    if fwhm_fn is None:
        raise SystemExit("в файле замера нет кривой ПШПВ прибора")

    hl = half_lives()
    br_tl = branch_tl208()
    print("периоды полураспада (МАГАТЭ, half_life_sec):")
    for key, lab, _el, fn, _br in CHAIN:
        t = hl[key]
        u = ("%.4g лет" % (t / YEAR_SEC) if t > YEAR_SEC
             else "%.4g сут" % (t / 86400.0) if t > 86400.0
             else "%.4g ч" % (t / 3600.0) if t > 3600.0
             else "%.4g с" % t)
        print("  %-8s %14s  %s" % (lab, u,
                                   "ENSDF, внесён явно" if fn is None
                                   else "выгрузка %s" % fn))
    print("ветвление Bi-212 -> Tl-208: %.2f %% (МАГАТЭ)" % (100.0 * br_tl))

    # --- активность по составу этикетки ------------------------------------
    tmap, head0 = load_templates(args.tdir, centres, fwhm_fn)
    mass_g = float(head0.get("wt20_mass_g", "0").split()[0])
    if mass_g <= 0:
        raise SystemExit("в шапке шаблона нет массы пачки (wt20_mass_g)")
    th_g = mass_g * 0.02 * 0.878809          # 2 % масс. ThO2, доля Th в ThO2
    a0 = th_g * SPEC_ACT_TH232
    print("\nсостав по этикетке: пачка %.1f г, 2 %% масс. ThO2 -> тория %.3f г"
          % (mass_g, th_g))
    print("активность Th-232: %.0f Бк (вычислена из T½ и молярной массы)" % a0)

    years = [float(v) for v in args.years.split(",") if v.strip()]
    print("\n--- активности звеньев после химической очистки, Бк ---")
    print("  %-8s %s" % ("возраст", "".join("%10s" % l
                                            for _k, l, _e, _f, _b in CHAIN)))
    acts_by_year = {}
    for yr in years:
        A = activities_after_purification(yr * YEAR_SEC, a0, hl, br_tl)
        acts_by_year[yr] = A
        print("  %5.1f лет %s" % (yr, "".join("%10.0f" % A[k]
                                              for k, _l, _e, _f, _b in CHAIN)))

    print("\n  отношения, по которым видно состояние ряда:")
    print("  %-8s %12s %12s" % ("возраст", "Ac-228/Pb-212", "Ra-228/Th-232"))
    for yr in years:
        A = acts_by_year[yr]
        print("  %5.1f лет %12.3f %12.3f"
              % (yr, A["Ac228"] / A["Pb212"] if A["Pb212"] else float("nan"),
                 A["Ra228"] / A["Th232"]))

    # --- спектры: метод 1 ---------------------------------------------------
    print("\n--- метод 1: понуклидные МК-шаблоны ---")
    print("  шаблонов загружено: %d (%s)"
          % (len(tmap), ", ".join(sorted(tmap))))
    cols = {}
    for yr in years:
        A = acts_by_year[yr]
        s = np.zeros(len(centres))
        for key, spec_k in tmap.items():
            s += A[key] * spec_k * t_smp
        cols["метод1_%gлет" % yr] = s

    # --- спектры: метод 2 (матрица отклика × выходы линий) ------------------
    # Тот же состав, тот же закон накопления — другой способ получить отклик
    # звена. Расхождение методов есть мера вклада суммирования каскада и
    # полноты библиотеки линий, а не свободный параметр.
    if args.mdir and os.path.isdir(args.mdir):
        import wt20_unfold_matrix as M2
        grid_E, mats, sigs = M2.read_matrix(args.mdir, centres)
        print("\n--- метод 2: матрица отклика × выходы линий МАГАТЭ ---")
        raw2 = {}
        for key, lab, _el, fn, _br in M2_chain_files():
            if fn is None:
                continue
            lines = M2.read_lines_ext(fn, 1.0, E_MAX)
            r, _v = M2.convolve_nuc(lines, grid_E, mats, sigs, centres,
                                    0.0, 0.0)
            raw2[key] = broaden(centres, r, centres, 0.0, 0.0,
                                fwhm_fn=fwhm_fn)
            print("  %-8s линий %d, Σ(y·ε) = %.4g" % (lab, len(lines),
                                                      float(r.sum())))
        for yr in years:
            A = acts_by_year[yr]
            s = np.zeros(len(centres))
            for key, spec_k in raw2.items():
                s += A[key] * spec_k * t_smp
            cols["метод2_%gлет" % yr] = s
    else:
        print("\n[метод 2 пропущен: каталог матрицы отклика не задан]")

    # --- запись -------------------------------------------------------------
    meas = y - ybg
    p_out = os.path.join(outdir, "forward_spectra.csv")
    names = ["E_keV", "измерено_минус_фон"] + list(cols)
    with io.open(p_out, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(names)
        for i, e in enumerate(centres):
            w.writerow(["%.1f" % e, "%.4g" % meas[i]]
                       + ["%.4g" % cols[c][i] for c in cols])
    print("\nзаписано: %s" % p_out)

    # Активности по возрастам — отдельным файлом: страница отчёта берёт числа
    # оттуда, а не повторяет их разметкой.
    p_act = os.path.join(outdir, "forward_activities.csv")
    with io.open(p_act, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["# активности звеньев ряда после химической очистки "
                    "тория, Бк; вычислены из состава этикетки и закона "
                    "накопления, под измеренный спектр не подгонялись"])
        w.writerow(["# A(Th-232) = %.0f Бк при 2 %% масс. ThO2 и массе "
                    "пачки %.1f г" % (a0, mass_g)])
        w.writerow(["возраст_лет"] + [l for _k, l, _e, _f, _b in CHAIN])
        for yr in years:
            A = acts_by_year[yr]
            w.writerow(["%g" % yr] + ["%.0f" % A[k]
                                      for k, _l, _e, _f, _b in CHAIN])
    print("записано: %s" % p_act)

    # Вклад отдельных звеньев в спектр — для заливки на рисунке наложения.
    # Берётся при первом из перечисленных возрастов.
    y0 = years[0]
    p_comp = os.path.join(outdir, "forward_components.csv")
    keys = [k for k, _l, _e, _f, _b in CHAIN if k in tmap]
    with io.open(p_comp, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["E_keV", "измерено_минус_фон", "сумма_модели"]
                   + [dict((k, l) for k, l, _e, _f, _b in CHAIN)[k]
                      for k in keys])
        A = acts_by_year[y0]
        tot = cols["метод1_%gлет" % y0]
        for i, e in enumerate(centres):
            w.writerow(["%.1f" % e, "%.4g" % meas[i], "%.4g" % tot[i]]
                       + ["%.4g" % (A[k] * tmap[k][i] * t_smp) for k in keys])
    print("записано: %s  (вклады звеньев при %g годах)" % (p_comp, y0))

    # --- сверка по полосам --------------------------------------------------
    BANDS = [(50, 72, "K-серия вольфрама"), (72, 100, "K-серия дочерних"),
             (150, 300, "полоса 238,63"), (500, 650, "полоса 583,19"),
             (850, 1000, "полоса 911,20"), (2500, 2700, "полоса 2614,51")]
    band_out = []
    for meth in ("метод1", "метод2"):
        if "%s_%gлет" % (meth, years[0]) not in cols:
            continue
        print("\n--- модель/измерено по полосам, %s ---" % meth)
        print("  %-28s %s" % ("полоса", "".join("%9s" % ("%g лет" % y_)
                                                for y_ in years)))
        for lo, hi, lab in BANDS:
            sel = (centres >= lo) & (centres < hi)
            sm = float(meas[sel].sum())
            vals = [float(cols["%s_%gлет" % (meth, y_)][sel].sum()) / sm
                    if sm > 0 else float("nan") for y_ in years]
            print("  %-28s %s" % ("%d-%d кэВ, %s" % (lo, hi, lab),
                                  "".join("%9.3f" % v for v in vals)))
            band_out.append([meth, lo, hi, lab, "%.0f" % sm]
                            + ["%.4f" % v for v in vals])
    if band_out:
        p_band = os.path.join(outdir, "forward_bands.csv")
        with io.open(p_band, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter=";")
            w.writerow(["# отношение модель/измерено по полосам энергии;"])
            w.writerow(["# модель считана вперёд по составу, не подгонялась"])
            w.writerow(["метод", "E_нач_кэВ", "E_кон_кэВ", "полоса",
                        "измерено"] + ["%g_лет" % y_ for y_ in years])
            w.writerows(band_out)
        print("\nзаписано: %s" % p_band)

    # --- согласие двух методов между собой ---------------------------------
    # Сравнение не с измерением, а методов друг с другом: оно не зависит от
    # того, верны ли активности, и потому проверяет ОТКЛИК. Расхождение вне
    # сумм-областей означает расхождение самих способов построения отклика;
    # в областях сумм-пиков каскада метод 2 занижен по построению.
    if "метод2_%gлет" % years[0] in cols:
        y0 = years[0]
        print("\n--- согласие методов между собой (метод2/метод1) ---")
        print("  отношение не зависит от активностей: они одни и те же")
        for lo, hi, lab in BANDS + [(3100, 3300, "сумм-пик 2614+583")]:
            sel = (centres >= lo) & (centres < hi)
            s1 = float(cols["метод1_%gлет" % y0][sel].sum())
            s2 = float(cols["метод2_%gлет" % y0][sel].sum())
            print("  %-28s %8.3f" % ("%d-%d кэВ, %s" % (lo, hi, lab),
                                     s2 / s1 if s1 > 0 else float("nan")))

    # --- какой возраст ближе всего -----------------------------------------
    # Не подгонка: активности при каждом возрасте вычислены заранее и не
    # трогаются. Здесь только называется, при каком из посчитанных возрастов
    # невязка по полосам наименьшая, и с какой она величиной.
    print("\n--- невязка по полосам (среднеквадратичное log(мод/изм)) ---")
    for meth in ("метод1", "метод2"):
        if "%s_%gлет" % (meth, years[0]) not in cols:
            continue
        best = None
        for y_ in years:
            ls = []
            for lo, hi, lab in BANDS:
                sel = (centres >= lo) & (centres < hi)
                sm = float(meas[sel].sum())
                sd = float(cols["%s_%gлет" % (meth, y_)][sel].sum())
                if sm > 0 and sd > 0:
                    ls.append(math.log(sd / sm) ** 2)
            q = math.sqrt(sum(ls) / len(ls)) if ls else float("nan")
            print("  %s, %5.1f лет: %.3f" % (meth, y_, q))
            if best is None or q < best[0]:
                best = (q, y_)
        if best:
            print("  -> наименьшая невязка %s при %g годах: %.3f"
                  % (meth, best[1], best[0]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
