"""Баланс «пик против континуума» в окне деконволюции: измерение и модель.

Зачем. Связанная деконволюция (`deconv.py`) даёт по группам 583,2 и 911,2
кэВ ряда тория активность на 12-26 % выше, чем по одиночной 2614,5. Ни
веса, ни нормировка, ни аннигиляционная линия этого не объясняли: на всех
ОДИНОЧНЫХ линиях деконволюция сходится с оконным съёмом до долей процента.

Проверка ставится так. Активность в `deconv.py` есть отношение амплитуд
одной и той же формы, снятых с измерения и с уширенной модели. Отношение
точно ровно настолько, насколько СОВПАДАЮТ ФОРМЫ. Поэтому здесь в каждом
окне считается, какую долю полной площади подгонка отдала пикам, а какую
подложке и ступеньке — отдельно у измерения и у модели. Столбец `изм/мод`
и есть множитель, на который активность по этой группе разойдётся с
истинной, если формы разошлись.

Что показала (28.07.2026): расхождение не случайно и не в связке линий.
Модель ЗАВЫШАЕТ континуум под группой 583,2 и ЗАНИЖАЕТ его под одиночной
2614,5 — то есть у расчётного спектра другой наклон подложки, а не другой
уровень. Эффект растёт от маринелли к Петри, то есть с уменьшением пробы
и приближением её к кристаллу.

Запуск:  python detectors/Gamma-1S/analysis/deconv_balance.py
"""
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


def _area(M, coef, dE, cols):
    return float(np.sum((M[:, cols] @ coef[cols]) * dE))


def _shares(M, coef, dE):
    """Доли пика, подложки и ступеньки в полной площади окна."""
    pk = _area(M, coef, dE, [0])
    fl = _area(M, coef, dE, [1, 2])
    st = _area(M, coef, dE, [3])
    tot = pk + fl + st
    return (None if tot <= 0 else (pk / tot, fl / tot, st / tot))


def main():
    print(__doc__.split("Запуск:")[0].strip())
    print()
    print("%-13s %8s %5s | %7s %7s %7s | %7s %7s %7s | %7s"
          % ("геометрия", "E, кэВ", "линий",
             "пик изм", "конт", "полка", "пик мод", "конт", "полка",
             "изм/мод"))
    for geom, mask, nuc, _asp, _dp, _d0, _m, _v in kr.VOLUME_RECORDS:
        if nuc != NUC:
            continue
        kd = paths.kit_dir(geom)
        files = sorted(str(p) for p in kd.rglob(mask)) if kd else []
        if not files:
            continue
        sp, bg, _cal = bm.read_checked(files[0])
        lines_of, ckey = kr.VLINES[nuc]
        base = kr.RUNBASE.get((geom, ckey))
        if not base:
            continue
        for E0 in lines_of:
            half = dc.SPAN * dc.fwhm(E0)
            lines = dc.group_lines(base, E0, half)
            if not lines:
                continue
            lines = sorted(lines)
            lo, hi = E0 - half, E0 + half
            arr, _N = dc._broadened(base)
            fm = dc._fit_model(arr, lines, lo, hi)
            best = None
            for d in dc.SHIFT_GRID:
                r = dc._fit_measured(sp, bg, lines, lo, hi, float(d))
                if r and r[0] > 0 and (best is None or r[2] < best[2]):
                    best = r
            if best is None or not fm:
                continue
            fit = best[5]
            x, dE, shift = fit["x"], fit["dE"], fit["shift"]
            sig_step = max(dc.sigma(e) for e, _ in lines)
            mid = 0.5 * (lines[0][0] + lines[-1][0])
            sh_m = _shares(dc._design(x, dE, lines, shift, mid + shift,
                                      sig_step), best[4], dE)
            xg = np.arange(int(round(lo)), int(round(hi)), dtype=float)
            dEg = np.ones_like(xg)
            sh_g = _shares(dc._design(xg, dEg, lines, 0.0, mid, sig_step),
                           fm[4], dEg)
            if not sh_m or not sh_g:
                continue
            print("%-13s %8.1f %5d | %6.1f%% %6.1f%% %6.1f%% | "
                  "%6.1f%% %6.1f%% %6.1f%% | %7.3f"
                  % (geom, E0, len(lines),
                     100 * sh_m[0], 100 * sh_m[1], 100 * sh_m[2],
                     100 * sh_g[0], 100 * sh_g[1], 100 * sh_g[2],
                     sh_m[0] / sh_g[0]))
    print("\nСтолбец изм/мод — во сколько раз подгонка отдала пикам больше на\n"
          "измерении, чем на модели. Единица означает, что формы совпали и\n"
          "нормировка деконволюции точна. Отклонение переносится в активность\n"
          "по этой группе множителем один к одному.")


if __name__ == "__main__":
    main()
