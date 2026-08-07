# -*- coding: utf-8 -*-
"""Разложение спектра WT-20 КВАЗИШАБЛОННЫМ методом ЛСРМ (гл. 13) — для сравнения.

Задача №39, директива оператора 07.08.2026: «для сравнения сделай разложение
по методике ЛСРМ». Метод — по «Алгоритмическим основам», гл. 13 дословно:

  «эталонный спектр каждого радионуклида представляется в виде суммы пиков
   полного поглощения всех линий этого радионуклида. Разложение спектра
   анализируемого образца осуществляется по спектрам полного поглощения
   отдельных радионуклидов и общей комптоновской подложке, которая
   рассчитывается для спектра анализируемого образца. Спектры полного
   поглощения рассчитываются из библиотечных линий радионуклида с учетом
   эффективности регистрации и формы пика … библиотека должна содержать все
   линии радионуклида в используемом … диапазоне.»

Отличия от шаблонного разложения §2 (wt20_unfold.py) — ровно два:
  1) эталон нуклида — НЕ полный Монте-Карло-отклик, а сумма ППП всех линий
     из библиотеки МАГАТЭ, умноженных на эффективность регистрации ε(E);
  2) континуум — НЕ из переноса, а «общая комптоновская подложка»,
     оцениваемая из самого измеренного спектра.
Решатель тот же — нормальная система гл. 12 (lsrm_solve), идентификация
по δa < 1 (гл. 14). Разложение ПОНУКЛИДНОЕ, как в SpectraLine: каждый
нуклид со своей активностью; ветви ряда собираются после решения.

Происхождение входов:
  * ε(E) — снята с тех же Монте-Карло-шаблонов: нетто ППП на распад в СЫРОМ
    (неуширенном) шаблоне, делённое на библиотечный выход линии; точки
    сглажены полиномом в лог-лог (стандартная форма кривых ЛСРМ). Другого
    источника эффективности для этой геометрии нет; так сравнение изолирует
    сам МЕТОД (ППП + подложка против полного отклика) при одной физике.
    Каскадное суммирование при этом сидит В ТОЧКАХ ε (это эффективные ППП
    данной геометрии), а квазишаблон переносит его на все линии нуклида
    гладкой кривой — в отличие от полного шаблона, где оно понуклидное.
  * Форма пика — ExpGaussExp с хвостами из файла замера, ПШПВ² = f₀ + f₁·E,
    пара перебирается по сетке вместе с активностями (директива оператора).
  * Подложка — алгоритм в гл. 13 НЕ расшифрован (страницы про подложку в
    книге описывают полином под пиком в зонах). Принята стандартная оценка
    континуума SNIP (Ryan et al., NIM B 9 (1985) 396) с энергозависимым
    полуокном 1,5·ПШПВ(E); это ВЫБОР РЕАЛИЗАЦИИ, помечен явно.
  * Веса гл. 12: дисперсия точки = счёт образца + вычтенный фон + подложка
    (подложка оценена из тех же данных; её включение в дисперсию —
    консервативный запас, σ выходят слегка шире).

    python analysis/wt20_unfold_lsrm.py <спектр.xml> <каталог шаблонов> [вывод]
"""
import csv
import io
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from wt20_unfold import (E_MAX, E_STEP, LIB, TAIL_L, TAIL_R,          # noqa: E402
                         branch_to_tl208, broaden, line_shape,
                         lsrm_solve, poly, read_correction,
                         read_template, rebin_to_grid)
from wt20_unfold import read_atomspectra_xml                          # noqa: E402

E_FIT = (150.0, 3000.0)   # та же полоса, что у отката §2 — сравнение честное

# Нуклиды разложения и их файлы линий МАГАТЭ. Th-232, Ra-228, Th-228 и прочие
# звенья почти без гамма-выхода в полосе не раскладываются (их квазишаблоны
# пусты) — это то же ограничение, что и в §2, только здесь оно видно прямо
# по библиотеке.
NUCS = [
    ("Ac228", "Ac-228", "228ac_gammas.csv", "#d81b8c"),
    ("Pb212", "Pb-212", "212pb_gammas.csv", "#b07d2a"),
    ("Bi212", "Bi-212", "212bi_gammas.csv", "#2f6b34"),
    ("Tl208", "Tl-208", "208tl_gammas.csv", "#8a8f2a"),
    ("Ra224", "Ra-224", "224ra_gammas.csv", "#c98b1e"),
    ("Rn220", "Rn-220", "220rn_gammas.csv", "#6b8f3a"),
]
# Файлы линий для съёма ε — всегда полный набор: исключение нуклида из
# РАЗЛОЖЕНИЯ не отменяет его опорных точек эффективности.
LINE_FILES = {k: f for k, _, f, _ in NUCS}

# Чистка библиотеки, как в SpectraLine при неоднозначности (гл. 14.3-14.4):
# WT20_LSRM_DROP="Ra224,Rn220" исключает нуклиды из набора разложения.
_drop = {s.strip() for s in os.environ.get("WT20_LSRM_DROP", "").split(",")
         if s.strip()}
NUCS = [t for t in NUCS if t[0] not in _drop]

# Опорные линии съёма ε(E): сильные, отделённые от соседних линий СВОЕГО же
# нуклида (чужие в шаблоне одного нуклида не присутствуют). Выход берётся из
# файла МАГАТЭ суммой по линиям, попавшим в окно пика, — пара 964,77 + 968,97
# так снимается одним окном автоматически.
EFF_POINTS = [
    ("Ac228", 153.98), ("Ac228", 209.25), ("Ac228", 270.25),
    ("Ac228", 338.32), ("Ac228", 911.20), ("Ac228", 968.97),
    ("Ac228", 1588.19),
    ("Pb212", 238.63), ("Pb212", 300.09),
    ("Ra224", 240.99),
    ("Tl208", 277.37), ("Tl208", 510.77), ("Tl208", 583.19),
    ("Tl208", 860.56), ("Tl208", 2614.51),
    ("Bi212", 727.33), ("Bi212", 1620.50),
]
PEAK_HALF = 2.5    # полуокно пика в сыром шаблоне (1 кэВ/бин), кэВ
SIDE_LO, SIDE_HI = 4.0, 12.0   # боковые полосы континуума, кэВ от линии


def read_lines(fname, e_lo, e_hi):
    """Все линии нуклида из файла МАГАТЭ в полосе: [(E, выход_долей), ...]."""
    out = []
    p = os.path.join(LIB, fname)
    for r in csv.DictReader(io.open(p, encoding="utf-8")):
        try:
            e = float(r["energy"])
            i = float(r["intensity"])
        except (TypeError, ValueError, KeyError):
            continue
        if i > 0 and e_lo <= e <= e_hi:
            out.append((e, i / 100.0))
    if not out:
        raise SystemExit("в %s нет линий в полосе %.0f-%.0f" % (p, e_lo, e_hi))
    return out


def eff_from_templates(tdir):
    """Точки ε(E) из сырых шаблонов: нетто ППП на распад / выход из файла.

    Сырой шаблон — гистограмма 1 кэВ ДО уширения: полное поглощение лежит в
    одном-двух бинах на энергии линии, соседние линии того же нуклида
    разделены. Континуум под пиком снимается средним по боковым полосам.
    Погрешность — счётная, по сырым отсчётам Монте-Карло.
    """
    pts = []
    lines_cache = {}
    for nk, e0 in EFF_POINTS:
        head, e, c = read_template(os.path.join(tdir, "%s.csv" % nk))
        n_prim = float(head["N_primaries"])
        if nk not in lines_cache:
            lines_cache[nk] = read_lines(LINE_FILES[nk], 1.0, 4000.0)
        sel = (e >= e0 - PEAK_HALF) & (e <= e0 + PEAK_HALF)
        side = (((e >= e0 - SIDE_HI) & (e <= e0 - SIDE_LO)) |
                ((e >= e0 + SIDE_LO) & (e <= e0 + SIDE_HI)))
        gross = float(c[sel].sum())
        cont = float(c[side].mean()) if side.any() else 0.0
        net = gross - cont * int(sel.sum())
        y_sum = sum(y for el, y in lines_cache[nk]
                    if e0 - PEAK_HALF <= el <= e0 + PEAK_HALF)
        if net <= 0 or y_sum <= 0:
            print("  ε-точка %s %.2f кэВ пропущена (нетто %.0f, выход %.4f)"
                  % (nk, e0, net, y_sum))
            continue
        eps = net / n_prim / y_sum
        sig = math.sqrt(gross + cont * int(sel.sum())) / n_prim / y_sum
        pts.append((e0, eps, sig, nk, y_sum))
    return pts


def fit_eff(pts):
    """Кривая ε(E): полином 3-й степени в лог-лог, веса по σ точек."""
    le = np.log([p[0] for p in pts])
    lf = np.log([p[1] for p in pts])
    w = np.array([p[1] / max(p[2], 1e-12) for p in pts])   # 1/σ_ln
    cf = np.polyfit(le, lf, 3, w=w)

    def eps(e):
        return np.exp(np.polyval(cf, np.log(np.asarray(e, float))))
    return eps, cf


def snip_baseline(y, half_win_bins):
    """Оценка континуума SNIP (LLS-вариант) с энергозависимым полуокном.

    Алгоритм гл. 13 не расшифрован в книге — это выбор реализации: Ryan et
    al. 1985, возрастающее окно, клип полуокна по 1,5·ПШПВ(E) в бинах.
    """
    v = np.log(np.log(np.sqrt(np.maximum(y, 0.0) + 1.0) + 1.0) + 1.0)
    n = len(v)
    hw = np.asarray(half_win_bins, int)
    idx = np.arange(n)
    for p in range(1, int(hw.max()) + 1):
        pe = np.minimum(p, hw)
        lo = np.maximum(idx - pe, 0)
        hi = np.minimum(idx + pe, n - 1)
        v = np.minimum(v, 0.5 * (v[lo] + v[hi]))
    b = (np.exp(np.exp(v) - 1.0) - 1.0) ** 2 - 1.0
    return np.maximum(b, 0.0)


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    src, tdir = sys.argv[1], sys.argv[2]
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
    ybg = np.zeros_like(y)
    if bg is not None:
        ybg = rebin_to_grid(np.asarray(bg.counts, float),
                            list(bg.energy_cal), corr.get("background"),
                            edges) * (t_smp / float(bg.real_time))

    m = (centres >= E_FIT[0]) & (centres <= E_FIT[1])
    yy = y[m]

    # --- эффективность регистрации из сырых шаблонов -------------------------
    print("[квазишаблонный метод ЛСРМ гл. 13; ε(E) — с МК-шаблонов, задача №39]")
    pts = eff_from_templates(tdir)
    eps_fn, cf = fit_eff(pts)
    print("ε(E): лог-лог полином 3-й ст. по %d точкам; коэф. %s" %
          (len(pts), ", ".join("%.5g" % v for v in cf)))
    print("  %8s %7s %10s %10s %7s" % ("E, кэВ", "нуклид", "ε_МК", "ε_кривая",
                                       "Δ, %"))
    for e0, ev, sg, nk, ys in sorted(pts):
        ef = float(eps_fn(e0))
        print("  %8.2f %7s %10.3e %10.3e %+7.1f" %
              (e0, nk, ev, ef, 100.0 * (ef - ev) / ev))

    # --- квазишаблоны: все линии МАГАТЭ в полосе × ε(E) ----------------------
    nuc_lines = {}
    for nk, lab, fname, _col in NUCS:
        ls = read_lines(fname, E_FIT[0], E_FIT[1])
        e_l = np.array([e for e, _ in ls])
        s_l = np.array([yl for _, yl in ls]) * eps_fn(e_l)   # отсч/распад в ППП
        nuc_lines[nk] = (e_l, s_l)
        print("  %s: линий в полосе %d, Σ(выход·ε) = %.4g" %
              (lab, len(ls), float(s_l.sum())))

    _qcache = {}

    def qtemplate(nk, f0, f1):
        key = (nk, round(f0, 6), round(f1, 6))
        if key not in _qcache:
            e_l, s_l = nuc_lines[nk]
            _qcache[key] = broaden(e_l, s_l, centres, f0, f1) * t_smp
        return _qcache[key]

    def fwhm_bins(f0, f1):
        fw = np.sqrt(np.maximum(f0 + f1 * centres, 4.0))
        return np.clip(np.round(1.5 * fw / E_STEP), 2, 80)[m]

    # --- форма пика: по гл. 13 она ЗАДАНА ЗАРАНЕЕ (файл пика-образа), а не
    # подгоняется. Свободный перебор здесь вырожден с подложкой: узкая ширина
    # сжимает SNIP-окно, подложка забирает основание пиков, и хи² выбирает
    # нефизичные 27,6 кэВ на 662 (перебор оставлен за WT20_LSRM_SCAN=1).
    # Умолчание — пара из шаблонного отката §2 (−400 + 2,80·E, 38,1 кэВ на
    # 662), согласованная с промером одиночных пиков этого спектра.
    y_net = np.maximum(y - ybg, 0.0)[m]

    def solve_at(kk, pp):
        # WT20_LSRM_NOBASE=1 — демонстрационный режим «только пики, без
        # подложки»: показывает, что достаётся активностям, если комптоновский
        # континуум не описан ничем.
        if os.environ.get("WT20_LSRM_NOBASE") == "1":
            base = np.zeros_like(y_net)
        else:
            base = snip_baseline(y_net, fwhm_bins(kk, pp))
        A = np.column_stack([qtemplate(nk, kk, pp)[m]
                             for nk, _, _, _ in NUCS])
        VarA = np.zeros_like(A)
        x, sig, chi2, act = lsrm_solve(A, VarA, yy, ybg[m] + base)
        return chi2, kk, pp, x, sig, act, A, base

    if os.environ.get("WT20_LSRM_SCAN") == "1":
        best = None
        for kk in np.arange(-600.0, 601.0, 100.0):
            for pp in np.arange(1.0, 5.01, 0.2):
                if kk + pp * E_FIT[0] <= 4.0:
                    continue
                r = solve_at(kk, pp)
                if r[5].any() and (best is None or r[0] < best[0]):
                    best = r
    else:
        f0f1 = os.environ.get("WT20_LSRM_FW", "-400,2.8")
        f0v, f1v = (float(v) for v in f0f1.split(","))
        print("форма пика фиксирована: ПШПВ² = %.0f + %.2f·E "
              "(гл. 13: пик-образ задаётся заранее)" % (f0v, f1v))
        best = solve_at(f0v, f1v)
    chi2, k_fw, p_fw, x, sig, act, A, base = best
    fw662 = math.sqrt(max(k_fw + p_fw * 661.657, 1.0))
    print("\nПШПВ² = %.0f + %.2f·E  (ПШПВ 662 = %.1f кэВ, %.2f %%), "
          "хи²/n = %.1f, полоса %.0f-%.0f кэВ"
          % (k_fw, p_fw, fw662, 100.0 * fw662 / 661.657, chi2, *E_FIT))
    print("подложка SNIP: %.0f отсчётов в полосе (%.1f %% измеренного)"
          % (base.sum(), 100.0 * base.sum() / max(yy.sum(), 1.0)))

    # --- идентификация и активности (гл. 14) --------------------------------
    br = branch_to_tl208()
    print("\n--- активности понуклидно (A ± σ, δa = σ/A) ---")
    rows = []
    for i, (nk, lab, _f, _c) in enumerate(NUCS):
        if not act[i]:
            print("  %-8s не идентифицирован (a < 0 или пустой шаблон)" % lab)
            rows.append((lab, 0.0, 0.0))
            continue
        da = sig[i] / x[i] if x[i] > 0 else float("inf")
        verdict = "идентифицирован" if da < 1.0 else "не идентифицирован (δa≥1)"
        print("  %-8s %10.0f ± %6.0f Бк  δa=%.3f  %s"
              % (lab, x[i], sig[i], da, verdict))
        rows.append((lab, float(x[i]), float(sig[i])))
    a_by = {lab: (a, s) for lab, a, s in rows}

    print("\n--- сборка в ветви ряда (для сравнения с §2) ---")
    a1, s1 = a_by.get("Ac-228", (0, 0))
    print("  A1 (по Ac-228)            = %6.0f ± %4.0f Бк" % (a1, s1))
    for lab in ("Pb-212", "Bi-212"):
        a, s = a_by.get(lab, (0, 0))
        print("  A2 (по %-7s)          = %6.0f ± %4.0f Бк" % (lab, a, s))
    a_tl, s_tl = a_by.get("Tl-208", (0, 0))
    print("  A2 (по Tl-208 / %.4f)   = %6.0f ± %4.0f Бк" %
          (br, a_tl / br, s_tl / br))

    # --- CSV -----------------------------------------------------------------
    with io.open(os.path.join(outdir, "unfold_lsrm_activities.csv"), "w",
                 encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["нуклид", "активность_Бк", "сигма_Бк"])
        for lab, a, s in rows:
            w.writerow([lab, "%.5g" % a, "%.4g" % s])
        w.writerow(["ПШПВ2_f0", "%.4g" % k_fw, ""])
        w.writerow(["ПШПВ2_f1", "%.4g" % p_fw, ""])
        w.writerow(["хи2_на_канал", "%.4g" % chi2, ""])
        w.writerow(["подложка_отсч", "%.6g" % base.sum(), ""])

    # --- рисунок в духе рис. 13.1 ЛСРМ ---------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    gx = centres[m]
    model = A @ x
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(14, 9), sharex=True,
        gridspec_kw={"height_ratios": [3.2, 1.0], "hspace": 0.06})
    ax.fill_between(gx, 1e-2, ybg[m] + base, color="#9aa4ad", alpha=0.55,
                    lw=0, label="подложка (SNIP) + фон")
    for i, (nk, lab, _f, col) in enumerate(NUCS):
        if act[i] and x[i] > 0:
            ax.plot(gx, ybg[m] + base + A[:, i] * x[i], color=col, lw=1.1,
                    label="%s — %.0f Бк" % (lab, x[i]))
    ax.plot(gx, ybg[m] + base + model, color="#c22", lw=1.3,
            label="огибающая (пики+подложка+фон)")
    ax.plot(gx, yy, color="k", lw=0.7, drawstyle="steps-mid", label="измерено")
    ax.set_yscale("log")
    ax.set_ylim(1.0, 2.5 * yy.max())
    ax.set_xlim(*E_FIT)
    ax.set_ylabel("отсчётов в канале %.0f кэВ" % E_STEP)
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=9)
    ax.set_title("WT-20: квазишаблонный метод ЛСРМ (гл. 13) — ППП всех линий "
                 "МАГАТЭ × ε(E) + подложка SNIP\nε(E) снята с МК-шаблонов; "
                 "ПШПВ² = %.0f + %.2f·E; хи²/n = %.1f" % (k_fw, p_fw, chi2))
    sd = np.sqrt(np.maximum(yy + ybg[m] + base, 1.0))
    axr.plot(gx, (yy - (ybg[m] + base + model)) / sd, color="#802",
             lw=0.7, drawstyle="steps-mid")
    axr.axhline(0, color="k", lw=0.6)
    axr.set_ylim(-8, 8)
    axr.set_ylabel("(изм. − мод.)/σ")
    axr.set_xlabel("Энергия, кэВ")
    axr.grid(alpha=0.25)
    out_png = os.path.join(outdir, "wt20_lsrm_quasitemplate.png")
    fig.savefig(out_png, dpi=130, bbox_inches="tight")
    print("\nзаписано: %s" % out_png)
    return 0


if __name__ == "__main__":
    sys.exit(main())
