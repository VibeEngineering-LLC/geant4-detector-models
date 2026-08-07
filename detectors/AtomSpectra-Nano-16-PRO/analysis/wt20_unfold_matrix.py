# -*- coding: utf-8 -*-
"""Разложение спектра WT-20 МЕТОДОМ 2 (матрица отклика × выходы линий).

Собирает шаблоны нуклидов свёрткой матрицы отклика (моноэнергетические γ,
разыгранные в объёме стержней, каталог `wt20_resp/`) с абсолютными выходами
линий из библиотеки МАГАТЭ (гамма + K/L-рентген атомной релаксации). Дальше
задача решается тем же способом, что и в методике 1 (`wt20_unfold.py`):
нормальная система взвешенного МНК с полной дисперсией весов, оценка
δa < 1 по гл. 14; ширина линии — ExpGaussExp с хвостами файла замера,
ПШПВ² = f₀ + f₁·E берётся из выходов методики 1 (единая пара, чтобы
сравнение форм шло при одинаковом уширении).

Что этот способ УМЕЕТ и ЧЕГО НЕ УМЕЕТ — см. `docs/wt20-methods-compare.md`.
Кратко: одиночный отклик каждой линии каскада воспроизводится линия-в-линию,
электронный канал источника (бета, тормозное, конверсионные электроны) и
истинное совпадающее суммирование каскада отсутствуют — вводятся отдельно.

    python analysis/wt20_unfold_matrix.py <спектр.xml> <каталог матрицы> [каталог вывода]
"""
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

from wt20_unfold import (E_MAX, E_STEP, LIB, N_SIG_MIN,               # noqa: E402
                         TAIL_L, TAIL_R,
                         branch_to_tl208, broaden, broaden_var,
                         lsrm_solve, read_atomspectra_xml,
                         read_correction, read_template, rebin_to_grid)

# Полоса подгонки — та же, что у методики 1, иначе методы сравнивать нельзя.
# Низ 150 кэВ: в полосе 72-100 кэВ модель даёт 0,24 измеренного (задача №9),
# и подгонка, вынужденная описывать эту дыру, ломает остальное. Прогон от
# 30 кэВ доступен переменной WT20_FIT_LO.
E_FIT = (float(os.environ.get("WT20_FIT_LO", "150")),
         float(os.environ.get("WT20_FIT_HI", "3000")))

# Нуклиды разложения — тот же набор, что в методике 1, включая слабые
# звенья (Th-232, Ra-228, Th-228, Rn-220): по гамма-выходу их вклад мал,
# но на рисунке нуклидная заливка должна быть одинаковой в обоих методах.
# Po-216 и Po-212 файлов линий не имеют — в NUCS не включены.
NUCS = [
    ("Th232", "Th-232", "232th", "#4a4034"),
    ("Ra228", "Ra-228", "228ra", "#6b5b45"),
    ("Ac228", "Ac-228", "228ac", "#b8347f"),
    ("Th228", "Th-228", "228th", "#8a6d3b"),
    ("Ra224", "Ra-224", "224ra", "#c98b1e"),
    ("Rn220", "Rn-220", "220rn", "#5f7d3a"),
    ("Pb212", "Pb-212", "212pb", "#c07b2a"),
    ("Bi212", "Bi-212", "212bi", "#3d7048"),
    ("Tl208", "Tl-208", "208tl", "#c9a227"),
]

# Ветви ряда — как в методике 1: A1 = Th-232 + Ra-228 + Ac-228 (равновесие
# внутри ветви — предположение, оговорено в тексте отчёта; замечание
# оператора «после очистки на заводе равновесие могло не установиться»
# закрывает задача №45); A2 = Th-228 + Ra-224 + Rn-220 + Pb-212 + Bi-212
# + Tl-208·br. Ветвление Bi-212 → Tl-208 читается из библиотеки МАГАТЭ.
def make_groups(br):
    """Компоненты разложения — ПО НУКЛИДАМ, по одной переменной на звено.

    Группировка в подцепочки A1/A2 снята (задача №68, директива оператора
    «а надо по нуклидам»): она держалась на предположении о вековом
    равновесии внутри ветви, которое этим замером не проверяется и после
    заводской очистки тория выполняться не обязано. Ветвление Bi-212 в
    матрицу не входит — Tl-208 свободен, а библиотечное значение служит
    величиной для сверки.

    Аргумент `br` сохранён для совместимости вызова и не используется.
    """
    return {lab: {key: 1.0} for key, lab, _fn, _c in NUCS}


def usable_nuclides(tdir=None):
    """Нуклиды, чей МК-шаблон набран достаточной статистикой (методика 1).

    Порог тот же, что в analysis/wt20_unfold.py (N_SIG_MIN): при малом числе
    событий шаблон есть реализация шума, а не спектр, и величина замером не
    определяется. Возвращает set ключей либо None, если каталог шаблонов не
    задан (переменная WT20_TEMPLATES).
    """
    if tdir is None:
        tdir = os.environ.get("WT20_TEMPLATES", "")
    if not tdir or not os.path.isdir(tdir):
        return None
    ok = set()
    for key, _lab, _fn, _c in NUCS:
        p = os.path.join(tdir, "%s.csv" % key)
        if not os.path.exists(p):
            continue
        n_sig = 0.0
        for ln in io.open(p, encoding="utf-8"):
            if not ln.startswith("#"):
                break
            if "N_with_signal" in ln and "=" in ln:
                try:
                    n_sig = float(ln.split("=", 1)[1].split()[0])
                except (ValueError, IndexError):
                    pass
        if n_sig >= N_SIG_MIN:
            ok.add(key)
    return ok or None


def read_lines_ext(nuc_key, e_lo, e_hi):
    """Все линии нуклида (гамма + K/L-рентген), доля на распад.

    Один источник — `*_gammas.csv`: выкачка МАГАТЭ по `rad_types=g` уже
    содержит рентген атомной релаксации в хвосте файла (побайтово совпадает
    с `*_xrays.csv`, проверено). Прежняя редакция читала оба файла — двойной
    счёт рентгена на +10…+48 %/распад (аудитор, 07.08.2026).

    Второй дубль: в `*_gammas.csv` для Kβ лежат ТРИ строки — сумма (`KB` в
    рентгеновском файле) и её компоненты (`KpB1` + `KpB2`). Метки `shell`
    в гамма-файле нет, поэтому дедупликация — по ЭНЕРГИЯМ компонент, взятым
    из рентгеновского файла (там shell есть).
    """
    drop_e = set()
    xp = os.path.join(LIB, "%s_xrays.csv" % nuc_key)
    if os.path.exists(xp):
        for r in csv.DictReader(io.open(xp, encoding="utf-8")):
            if r.get("shell") in ("KpB1", "KpB2"):
                try:
                    drop_e.add(round(float(r["energy"]), 3))
                except (TypeError, ValueError):
                    pass
    out = []
    p = os.path.join(LIB, "%s_gammas.csv" % nuc_key)
    if not os.path.exists(p):
        return out
    for r in csv.DictReader(io.open(p, encoding="utf-8")):
        try:
            e = float(r["energy"])
            i = float(r["intensity"])
        except (TypeError, ValueError, KeyError):
            continue
        if round(e, 3) in drop_e:
            continue
        if i > 0 and e_lo <= e <= e_hi:
            out.append((e, i / 100.0))
    return out


def read_matrix(mdir, e_grid):
    """Читает матрицу R(E_l → i) с сетки моноэнергетических прогонов.

    Возвращает (grid_E — сетка входных энергий линий, кэВ; M — массив
    [len(grid_E), len(e_grid)], нормированный на распад, то есть на квант
    входа). Сырые файлы имеют шаг 1 кэВ, ребиннятся к сетке разложения.
    """
    files = sorted(f for f in os.listdir(mdir) if f.startswith("rod_E")
                   and f.endswith(".csv") and "_" not in f[5:-4])
    if not files:
        raise SystemExit("в %s нет файлов rod_E*.csv" % mdir)
    Es, mats, sigs = [], [], []
    stamps = set()
    for f in files:
        head, e, c = read_template(os.path.join(mdir, f))
        n = float(head["N_primaries"])
        stamps.add(head.get("src_sha1", "?"))
        # энергия входа снята из имени файла (rod_E<XXXX.X>.csv)
        Ein = float(f[5:-4])
        # ребиннинг гистограммы 1 кэВ к сетке разложения
        edges = np.concatenate(([e_grid[0] - 0.5*(e_grid[1]-e_grid[0])],
                                0.5*(e_grid[:-1]+e_grid[1:]),
                                [e_grid[-1] + 0.5*(e_grid[-1]-e_grid[-2])]))
        hist = np.zeros(len(e_grid))
        # e — центры бинов 1 кэВ, c — счёт в этих бинах
        idx = np.digitize(e, edges) - 1
        for j, cnt in zip(idx, c):
            if 0 <= j < len(e_grid):
                hist[j] += cnt
        Es.append(Ein)
        mats.append(hist / n)      # на квант входа
        sigs.append(np.sqrt(hist) / n)
    if len(stamps) > 1:
        raise SystemExit("матрица набрана разными ревизиями: %s"
                         % ", ".join(sorted(stamps)))
    print("матрица: %d точек %.0f-%.0f кэВ, штамп %s"
          % (len(Es), min(Es), max(Es), stamps.pop() if stamps else "?"))
    return np.array(Es), np.array(mats), np.array(sigs)


def convolve_nuc(lines, grid_E, mats, sigs, e_grid, f0, f1):
    """Собирает шаблон нуклида: (r_i, σr_i) на распад.

    Для каждой линии выбирается ДВА ближайших узла матрицы; отклик линии —
    линейная интерполяция откликов узлов (интерполяция ФОРМЫ и МАСШТАБА
    вместе). Сдвиг ППП на разность энергий узла и линии делается сдвигом
    гистограммы на целое число бинов сетки разложения.

    Всё собранное — сумма ПОДПИК-ПОДКОНТИНУУМОВ, без приборного уширения.
    Уширение навешивается СНАРУЖИ единой свёрткой ExpGaussExp, чтобы
    результат сопоставлялся с методикой 1 при одинаковой ширине.
    """
    step = e_grid[1] - e_grid[0]
    raw = np.zeros(len(e_grid))
    var = np.zeros(len(e_grid))
    for E_l, y in lines:
        if y <= 0:
            continue
        # два ближайших узла
        j = int(np.searchsorted(grid_E, E_l))
        j_lo = max(0, min(len(grid_E) - 2, j - 1))
        j_hi = j_lo + 1
        E_lo, E_hi = grid_E[j_lo], grid_E[j_hi]
        w_hi = 0.0 if E_hi == E_lo else (E_l - E_lo) / (E_hi - E_lo)
        w_hi = max(0.0, min(1.0, w_hi))
        w_lo = 1.0 - w_hi
        for w, jk, En in ((w_lo, j_lo, E_lo), (w_hi, j_hi, E_hi)):
            if w <= 0:
                continue
            shift = int(round((E_l - En) / step))
            m = mats[jk]
            s = sigs[jk]
            if shift == 0:
                raw += y * w * m
                var += (y * w) ** 2 * s * s
            elif shift > 0:
                raw[shift:] += y * w * m[:-shift]
                var[shift:] += (y * w) ** 2 * s[:-shift] * s[:-shift]
            else:
                raw[:shift] += y * w * m[-shift:]
                var[:shift] += (y * w) ** 2 * s[-shift:] * s[-shift:]
    # приборное уширение — той же формой, что в методике 1
    return raw, var


def broaden_dispersion(raw_var, e_grid, f0, f1):
    """Уширение дисперсии сырого шаблона тем же ядром в квадрате.

    Использует broaden_var из wt20_unfold, но здесь сырой массив уже
    гистограмма (шаг сетки), а не список одиночных линий. Разложение: точку
    сетки трактуем как отдельную линию с интенсивностью raw_var[k] * N,
    N=1 — усреднение раз просто вернёт то же значение. Проще применить ту
    же процедуру broaden к массиву σ² и результат оставить как есть — это
    сумма Kᵢⱼ² · c_j — дисперсия свёртки пуассоновских счётчиков.
    """
    lin = broaden(e_grid, raw_var, e_grid, f0, f1)
    return lin


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    src, mdir = sys.argv[1], sys.argv[2]
    outdir = sys.argv[3] if len(sys.argv) > 3 else os.path.dirname(src)
    os.makedirs(outdir, exist_ok=True)

    spec = read_atomspectra_xml(src)
    if isinstance(spec, (list, tuple)):
        spec = spec[0]
    bg = getattr(spec, "background_embedded", None)
    t_smp = float(spec.real_time)

    corr = read_correction(os.path.join(outdir, "calibration_fitted.csv"))
    if not corr:
        print("ВНИМАНИЕ: поправок калибровки нет — работаем по заводской шкале")

    edges = np.arange(0.0, E_MAX + E_STEP, E_STEP)
    centres = 0.5 * (edges[:-1] + edges[1:])

    y = rebin_to_grid(np.asarray(spec.counts, float), list(spec.energy_cal),
                      corr.get("sample"), edges)
    if bg is not None:
        ybg = rebin_to_grid(np.asarray(bg.counts, float),
                            list(bg.energy_cal), corr.get("background"),
                            edges) * (t_smp / float(bg.real_time))
    else:
        ybg = np.zeros_like(y)

    grid_E, mats, sigs = read_matrix(mdir, centres)

    # Пара (f₀, f₁) — из выходов методики 1, чтобы уширение было одинаковым
    # и сравнение форм не примешивало разницу разрешений.
    f0, f1 = -400.0, 2.8
    p_act = os.path.join(outdir, "unfold_activities.csv")
    if os.path.exists(p_act):
        for row in csv.reader(io.open(p_act, encoding="utf-8")):
            if not row or row[0].startswith("#"):
                continue
            if row[0] == "ПШПВ2_f0":
                f0 = float(row[1])
            if row[0] == "ПШПВ2_f1":
                f1 = float(row[1])
    print("ПШПВ² = %.1f + %.3f·E  (взято из методики 1)" % (f0, f1))

    # Набор нуклидов берётся ТОТ ЖЕ, что у методики 1 — иначе методы
    # сравнивать нельзя. Отбраковка там идёт по статистике прогона распада
    # (N_with_signal в шапке шаблона): если прогон не набрал событий, у звена
    # нет заметного гамма-выхода, и его вклад неопределим в обоих методах
    # одинаково. Список читается из тех же файлов шаблонов; при их отсутствии
    # берутся все нуклиды, а расхождение наборов печатается.
    usable = usable_nuclides()
    if usable is not None:
        skipped = [lab for key, lab, _f, _c in NUCS if key not in usable]
        if skipped:
            print("исключены (нет статистики прогона у методики 1): %s"
                  % ", ".join(skipped))

    # Понуклидные шаблоны
    br_tl = branch_to_tl208()
    raw = {}
    varm = {}
    for key, lab, fn, _c in NUCS:
        if usable is not None and key not in usable:
            continue
        lines = read_lines_ext(fn, 1.0, E_MAX)
        r, v = convolve_nuc(lines, grid_E, mats, sigs, centres, f0, f1)
        raw[key] = r
        varm[key] = v
        print("  %-7s линий %d, Σ(y·ε) = %.4g" %
              (lab, len(lines), float(r.sum())))

    # Уширение — тем же broaden, что в методике 1
    tmap = {}
    vmap = {}
    for key in raw:
        # broaden ждёт линии (E, c) — передаём сетку и содержимое; получится
        # свёртка ядром ExpGaussExp той же ширины
        tmap[key] = broaden(centres, raw[key], centres, f0, f1) * t_smp
        vmap[key] = broaden_dispersion(varm[key], centres, f0, f1) * t_smp**2

    # Компоненты — только те, чьи шаблоны построены (набор согласован с
    # методикой 1 выше). Отброшенные нуклиды в GROUP не входят: иначе столбец
    # матрицы подгонки собирается из пустой суммы и ломает сборку.
    GROUP = {lab: mem for lab, mem in make_groups(br_tl).items()
             if all(k in tmap for k in mem)}
    print("ветвление Bi-212 -> Tl-208: %.2f %% (МАГАТЭ)" % (100.0 * br_tl))

    m = (centres >= E_FIT[0]) & (centres <= E_FIT[1])
    yy = y[m]

    def build():
        A = np.column_stack([sum(w * tmap[k][m] for k, w in members.items()
                                 if k in tmap)
                             for members in GROUP.values()])
        V = np.column_stack([sum(w * w * vmap[k][m] for k, w in members.items()
                                 if k in vmap)
                             for members in GROUP.values()])
        return A, V

    A, V = build()
    x, sig, chi2, act = lsrm_solve(A, V, yy, ybg[m])
    print("\n[метод 2, окно %.0f-%.0f кэВ, решение по ЛСРМ гл. 12]" % E_FIT)
    print("хи²/n = %.1f" % chi2)
    print("--- активности (A ± σ, δa) ---")
    rows = []
    for i, gname in enumerate(GROUP):
        if not act[i]:
            print("  %-38s не идентифицирован" % gname)
            rows.append((gname, 0.0, 0.0))
            continue
        da = sig[i] / x[i] if x[i] > 0 else float("inf")
        print("  %-38s %8.0f ± %5.0f Бк  δa=%.3f" %
              (gname, x[i], sig[i], da))
        rows.append((gname, float(x[i]), float(sig[i])))

    # --- запись выходов, совместимых с draw_wt20_unfold и wt20_report -----
    fw662 = math.sqrt(max(f0 + f1 * 661.657, 1.0))
    with io.open(os.path.join(outdir, "unfold_matrix_activities.csv"), "w",
                 encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["компонента", "активность_Бк", "сигма_Бк"])
        for lab, a, s in rows:
            w.writerow([lab, "%.5g" % a, "%.4g" % s])
        w.writerow(["ПШПВ_662_кэВ", "%.3g" % fw662, ""])
        w.writerow(["ПШПВ2_f0", "%.4g" % f0, ""])
        w.writerow(["ПШПВ2_f1", "%.4g" % f1, ""])
        w.writerow(["хи2_на_канал", "%.4g" % chi2, ""])
        w.writerow(["метод", "матрица отклика × выходы линий МАГАТЭ", ""])

    # Пофрагментный спектр — как у методики 1: заливка ПО ОТДЕЛЬНЫМ
    # НУКЛИДАМ, а не по подцепочкам (директива оператора 07.08.2026).
    # Подгонка сама остаётся двухпараметрической (A1, A2) — иначе на
    # слабых звеньях свобода теряет смысл. Вклад нуклида k в подцепочке g
    # с внутренним весом w — это x[g]·w·T_k(E), где T_k(E) — уширенный
    # понуклидный шаблон.
    model = A @ x + ybg[m]
    group_list = list(GROUP)
    names_out = ["E_keV", "измерено", "модель"]
    parts = []
    for key, lab, _fn, _c in NUCS:
        if key not in tmap:
            continue
        col = np.zeros(int(m.sum()))
        for gi, gname in enumerate(group_list):
            w = GROUP[gname].get(key, 0.0)
            if w > 0 and act[gi]:
                col += x[gi] * w * tmap[key][m]
        if col.sum() <= 0:
            continue
        parts.append(col)
        names_out.append(lab)
    parts.append(ybg[m])
    names_out.append("фон")

    with io.open(os.path.join(outdir, "unfold_matrix_spectrum.csv"), "w",
                 encoding="utf-8", newline="") as f:
        f.write(",".join(names_out) + "\n")
        for i in range(int(m.sum())):
            row = [centres[m][i], yy[i], model[i]] + [p[i] for p in parts]
            f.write(",".join("%.6g" % v for v in row) + "\n")
    print("\nзаписано: unfold_matrix_activities.csv, unfold_matrix_spectrum.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
