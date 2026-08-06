# -*- coding: utf-8 -*-
"""Отношение площадей Tl-208 583,19 / 2614,51 против толщины кристалла по пучку.

Диагностика расхождения: измеренное отношение 7,84, модель при паспортных
15,00 мм даёт 3,51. Отношение не зависит ни от активности, ни от нормировки
источника — только от отклика. Здесь считается, как оно ведёт себя с толщиной
кристалла ПО ПУЧКУ (в этой постановке излучение входит через грань 18 × 60 мм).

Окна те же, что в `wt20_lines.py`: ±1 ПШПВ, подложка по двум боковым окнам.

ОГОВОРКА К СКАНУ. Начало координат модели — ЦЕНТР кристалла, а корпус строится
от его граней, поэтому при изменении cryY вместе с толщиной уезжает и плечо
«источник — кристалл» (при 15 мм оно 14,5 мм, при 9 мм — 20,5 мм). Абсолютные
эффективности скана из-за этого сравнивать нельзя. ОТНОШЕНИЕ линий от плеча
почти не зависит: расстояние меняет телесный угол для обеих линий одинаково, —
и только оно здесь и читается.

    python analysis/wt20_depth_ratio.py <каталог прогона> [измеренное отношение]
"""
import os
import re
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import wt20_lines as L                                            # noqa: E402
import wt20_unfold as U                                           # noqa: E402


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    d = sys.argv[1]
    meas = float(sys.argv[2]) if len(sys.argv) > 2 else 7.84

    edges = np.arange(0.0, U.E_MAX + U.E_STEP, U.E_STEP)
    centres = 0.5 * (edges[:-1] + edges[1:])
    print("%8s %12s %12s %10s %12s"
          % ("cryY, мм", "583 на расп", "2614 на расп", "отношение",
             "изм./модель"))
    rows = []
    for fn in sorted(os.listdir(d)):
        m = re.match(r"Tl208_cryY([0-9.]+)\.csv$", fn)
        if not m:
            continue
        head, e, c = U.read_template(os.path.join(d, fn))
        n = float(head["N_primaries"])
        sp = U.broaden(e, c / n, centres, 200.0, 2.0)
        w1 = L.window_area(sp, centres, 583.19)
        w2 = L.window_area(sp, centres, 2614.51)
        if not w1 or not w2 or w2["net"] <= 0:
            continue
        r = w1["net"] / w2["net"]
        rows.append((float(m.group(1)), r))
        print("%8.1f %12.4e %12.4e %10.2f %12.2f"
              % (float(m.group(1)), w1["net"], w2["net"], r, meas / r))
    if len(rows) >= 2:
        rows.sort()
        xs = [x for x, _ in rows]
        ys = [y for _, y in rows]
        print("\nизмеренное отношение %.2f; линейная интерполяция по сетке "
              "даёт толщину %s"
              % (meas,
                 ("%.1f мм" % float(np.interp(meas, ys[::-1], xs[::-1])))
                 if min(ys) <= meas <= max(ys)
                 else "ВНЕ диапазона сканирования — расхождение толщиной не "
                      "объясняется"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
