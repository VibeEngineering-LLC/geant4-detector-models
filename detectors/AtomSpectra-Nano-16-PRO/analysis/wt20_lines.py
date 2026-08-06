# -*- coding: utf-8 -*-
"""Активность ПО ОТДЕЛЬНЫМ ЛИНИЯМ: оконный съём измеренного и модельного пика.

Зачем при наличии шаблонной подгонки. Подгонка по всему спектру опирается на
форму континуума, а континуум в контактной геометрии зависит от того, что лежит
ВОКРУГ (стол, стена, рука) — то есть от вещей, которых в модели нет. Оконный
съём этого не касается: и в измерении, и в модели берётся ОДНО И ТО ЖЕ окно
вокруг линии, из обоих вычитается подложка по одному правилу, и отношение даёт
активность. Если два способа расходятся, расхождение — это результат, а не шум.

Величина определена только вместе с окном: окно ±1 ПШПВ, подложка — пьедестал
по ГОСТ 26874-86 (метод выбирается по симметрии), опорная плоскость — дно
корпуса, геометрия — пачка WT-20 под прибором.

    python analysis/wt20_lines.py <спектр.xml> <каталог шаблонов> [каталог вывода]
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

from gamma.peaks.centroid_gost import gost_select_pedestal_method  # noqa: E402

import wt20_unfold as U                                            # noqa: E402
import wt20_calibration as C                                       # noqa: E402

# Линия -> файл шаблона нуклида, который её даёт. Совпадение линий разных
# нуклидов внутри окна учитывается моделью само: в окно модельного спектра
# попадает всё, что туда попадает и в измерении, — но только от ЭТОГО нуклида,
# поэтому линии выбраны так, чтобы соседей из других звеньев ряда не было.
LINES = [
    ("Pb212", 238.63, "Pb-212"),
    ("Tl208", 583.19, "Tl-208"),
    ("Bi212", 727.33, "Bi-212"),
    ("Ac228", 911.20, "Ac-228"),
    ("Ac228", 968.97, "Ac-228"),
    ("Tl208", 2614.51, "Tl-208"),
]
FWHM_662 = 39.0     # кэВ, по подгонке ПШПВ в wt20_unfold (ПШПВ² = 200 + 2,0·E)


def fwhm_keV(e):
    return math.sqrt(max(200.0 + 2.0 * e, 1.0))


def window_area(counts_e, centres, e0, half_fwhm=1.0):
    """Площадь пика в окне ±half_fwhm·ПШПВ с вычетом линейной подложки.

    Подложка строится по средним в двух окнах шириной 0,5 ПШПВ, отставленных на
    1,5 ПШПВ от центра линии, — тот же приём, что у пьедестала ГОСТ, но на
    равномерной энергетической сетке, где ширина канала одна и та же.
    """
    fw = fwhm_keV(e0)
    step = centres[1] - centres[0]
    lo = e0 - half_fwhm * fw
    hi = e0 + half_fwhm * fw
    sel = (centres >= lo) & (centres <= hi)
    if sel.sum() < 3:
        return None
    bl = (centres >= e0 - 2.0 * fw) & (centres <= e0 - 1.5 * fw)
    br = (centres >= e0 + 1.5 * fw) & (centres <= e0 + 2.0 * fw)
    if bl.sum() < 2 or br.sum() < 2:
        return None
    b = 0.5 * (counts_e[bl].mean() + counts_e[br].mean())
    gross = float(counts_e[sel].sum())
    back = b * sel.sum()
    return dict(gross=gross, back=back, net=gross - back, n=int(sel.sum()),
                fwhm=fw)


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    src, tdir = sys.argv[1], sys.argv[2]
    outdir = sys.argv[3] if len(sys.argv) > 3 else os.path.dirname(src)

    from gamma.io.atomspectra_xml import read_atomspectra_xml
    spec = read_atomspectra_xml(src)
    if isinstance(spec, (list, tuple)):
        spec = spec[0]
    bg = getattr(spec, "background_embedded", None)
    t = float(spec.real_time)

    corr = U.read_correction(os.path.join(outdir, "calibration_fitted.csv"))
    edges = np.arange(0.0, U.E_MAX + U.E_STEP, U.E_STEP)
    centres = 0.5 * (edges[:-1] + edges[1:])
    y = U.rebin_to_grid(np.asarray(spec.counts, float),
                        list(spec.energy_cal), corr.get("sample"), edges)
    if bg is not None:
        yb = U.rebin_to_grid(np.asarray(bg.counts, float),
                             list(bg.energy_cal), corr.get("background"),
                             edges) * (t / float(bg.real_time))
        y = y - yb
        print("фон вычтен: %.0f отсчётов, приведён к %.0f с" % (yb.sum(), t))

    print("\n%-8s %9s %10s %10s %12s %10s"
          % ("нуклид", "E, кэВ", "изм. нетто", "модель/расп", "A, Бк", "±, %"))
    rows = []
    for key, e0, label in LINES:
        p = os.path.join(tdir, "%s.csv" % key)
        if not os.path.exists(p):
            continue
        head, te, tc = U.read_template(p)
        n_prim = float(head["N_primaries"])
        # модельный спектр на сетке разложения, уширенный тем же законом
        tm = U.broaden(te, tc / n_prim, centres, 200.0, 2.0)
        wm = window_area(tm, centres, e0)
        wy = window_area(y, centres, e0)
        if not wm or not wy or wm["net"] <= 0:
            print("  %-8s %9.2f — окно не строится" % (label, e0))
            continue
        a = wy["net"] / (wm["net"] * t)
        # неопределённость: счётная у измерения (брутто + подложка) и
        # статистика розыгрыша у модели
        du = math.sqrt(wy["gross"] + wy["back"]) / max(wy["net"], 1.0)
        dm = 1.0 / math.sqrt(max(wm["net"] * n_prim, 1.0))
        err = 100.0 * math.hypot(du, dm)
        print("%-8s %9.2f %10.0f %10.3e %12.0f %10.1f"
              % (label, e0, wy["net"], wm["net"], a, err))
        rows.append((label, e0, wy["net"], wm["net"], a, err))

    if rows:
        with io.open(os.path.join(outdir, "line_activities.csv"), "w",
                     encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["нуклид", "E_кэВ", "нетто_изм", "модель_на_распад",
                        "A_Бк", "неопр_%"])
            for r in rows:
                w.writerow([r[0], "%.2f" % r[1], "%.0f" % r[2], "%.4e" % r[3],
                            "%.4g" % r[4], "%.1f" % r[5]])
        print("\nзаписано: %s" % os.path.join(outdir, "line_activities.csv"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
