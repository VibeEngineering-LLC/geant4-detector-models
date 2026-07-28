"""Континуум расчётного спектра против измеренного — ВНЕ пиков.

Зачем. `deconv_balance.py` показал, что расхождение групп 583,2 и 911,2 с
одиночной 2614,5 сидит не в связке линий и не в нормировке, а в ФОРМЕ
подложки: под группой модель даёт континуума больше измеренного, под
одиночной линией меньше. Проверять после этого надо не деконволюцию, а сам
расчётный спектр.

Как ставится. Модель даёт отсчёты на N разыгранных распадов, измерение —
отсчёты за живое время при неизвестной активности A. Их отношение

    A_уч = N · (счёт измерения) / (счёт модели) / t_живое

есть активность, которую дал бы участок спектра, будь модель верна. На пике
эта величина совпадает с активностью по деконволюции — так и нормируется
сверка. На УЧАСТКАХ БЕЗ ЛИНИЙ она показывает, где расчётный континуум
расходится с измеренным и в какую сторону.

Участки выбираются автоматически: сетка по энергии, из неё выбрасывается
всё, что ближе KEEP_OUT сигм к любой линии спектра испускания ярче
MIN_YIELD. Так в сравнение попадает только подложка.

Запуск:  python detectors/Gamma-1S/analysis/continuum.py
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402
import becqmoni as bm  # noqa: E402
import deconv as dc  # noqa: E402
import kit_recalc as kr  # noqa: E402

NUC = "Th-232"
E_LO, E_HI, E_STEP, E_HALF = 120.0, 2500.0, 40.0, 12.0
# Порог выхода и запретная зона вокруг линии подобраны так, чтобы в НИЖНЕЙ
# половине спектра вообще осталось где мерить: ряд тория густой, и при
# запрете в 3σ вокруг всего ярче полупроцента свободных участков ниже
# 1800 кэВ не остаётся ни одного. Отсюда порог 2 % и зона 2,5σ — линии
# слабее этого дают под окном вклад меньше процента.
MIN_YIELD = 0.02
KEEP_OUT = 2.5             # в сигмах
ANCHOR = 2614.511          # чистая одиночная линия, по ней нормируется сверка


def _emit_lines(base):
    p = os.path.join(dc.BUILD, base + "_emit.csv")
    if not os.path.exists(p):
        return None
    emit, N = kr.load_hist(p)
    if not N:
        return None
    return [(e, c / N) for e, c in emit.items() if c / N > MIN_YIELD]


def _counts(sp, bg, arr, lo, hi):
    """Счёт измерения за вычетом фона и счёт модели в том же окне."""
    ch = np.arange(len(sp.n), dtype=float)
    en = sp.energy(ch)
    m = (en >= lo) & (en <= hi)
    if m.sum() < 4:
        return None
    gross = float(sp.n[m].sum())
    ybg = 0.0
    if bg is not None:
        bch = np.arange(len(bg.n), dtype=float)
        ybg = float(np.interp(en[m], bg.energy(bch),
                              bg.n.astype(float)).sum()) * (sp.live / bg.live)
    net = gross - ybg
    i0, i1 = max(0, int(round(lo))), min(len(arr), int(round(hi)))
    mod = float(arr[i0:i1].sum())
    if mod <= 0 or net <= 0:
        return None
    return net, math.sqrt(gross + ybg), mod


def main():
    print(__doc__.split("Запуск:")[0].strip())
    print()
    rows = []
    for geom, mask, nuc, _a, _d, _d0, _m, _v in kr.VOLUME_RECORDS:
        if nuc != NUC:
            continue
        kd = paths.kit_dir(geom)
        files = sorted(str(p) for p in kd.rglob(mask)) if kd else []
        if not files:
            continue
        sp, bg, _cal = bm.read_checked(files[0])
        _lines_of, ckey = kr.VLINES[nuc]
        base = kr.RUNBASE.get((geom, ckey))
        if not base:
            continue
        arr, N = dc._broadened(base)
        lines = _emit_lines(base)
        if not lines or not N:
            continue

        half = dc.SPAN * dc.fwhm(ANCHOR)
        ref = _counts(sp, bg, arr, ANCHOR - half, ANCHOR + half)
        if not ref:
            continue
        A_ref = N * ref[0] / ref[2] / sp.live

        print("%s — опора %.1f кэВ даёт %.1f Бк; ниже отношение к ней"
              % (geom, ANCHOR, A_ref))
        print("   %9s %10s %8s %8s" % ("окно, кэВ", "A уч., Бк", "к опоре", "±"))
        E = E_LO
        while E <= E_HI:
            lo, hi = E - E_HALF, E + E_HALF
            near = min((abs(e - E) / dc.sigma(E) for e, _I in lines),
                       default=99.0)
            if near > KEEP_OUT:
                r = _counts(sp, bg, arr, lo, hi)
                if r:
                    A = N * r[0] / r[2] / sp.live
                    print("   %4.0f..%-4.0f %10.1f %8.3f %8.3f"
                          % (lo, hi, A, A / A_ref, A / A_ref * r[1] / r[0]))
                    rows.append("%s,%.0f,%.0f,%.2f,%.4f,%.4f"
                                % (geom, lo, hi, A, A / A_ref,
                                   A / A_ref * r[1] / r[0]))
            E += E_STEP
        print()
    out = os.path.join(str(paths.results("Gamma-1S")), "continuum.csv")
    with open(out, "w", encoding="utf-8", newline="") as fh:
        fh.write("# континуум вне пиков: активность по участку против опоры\n")
        fh.write("geometry,E_lo_keV,E_hi_keV,A_Bq,ratio_to_anchor,d_ratio\n")
        fh.write("\n".join(rows) + "\n")
    print("    таблица: %s (%d строк)" % (out, len(rows)))
    print("Единица в столбце «к опоре» означает, что расчётный континуум на\n"
          "этом участке согласен с измеренным в той же шкале, в какой согласован\n"
          "пик опоры. Уклонение вверх — модель НЕДОдаёт подложки, вниз —\n"
          "ПЕРЕдаёт.")


if __name__ == "__main__":
    main()
