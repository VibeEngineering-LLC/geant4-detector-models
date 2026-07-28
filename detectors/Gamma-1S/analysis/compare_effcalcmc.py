"""Сверка нашей кривой Geant4 с расчётом EffCalcMC на ТЕХ ЖЕ слоях.

Зачем. У модели есть систематический ход по энергии против ИЗМЕРЕННОЙ кривой
ЛСРМ (REPORT §2, §5а): занижение мягкого края, завышение жёсткого. Ход может
принадлежать либо коду (розыгрыш, съём площади, телесный угол), либо описанию
слоёв (толщины/плотности против реального прибора). Расчёт EffCalcMC по нашим
же .din/.sin (nuclidemaster/) разделяет гипотезы: независимая реализация с
теми же слоями воспроизводит либо нашу кривую (дефект в описании), либо
измеренную (дефект в коде).

Итог первого прогона (маринелли, 10^8 испытаний, 28.07.2026): среднее
отношение наш/EffCalcMC = 1,000, наклон −1,2 % на декаду, перепад по
48–3000 кэВ −4,8 % — при наблюдаемых против измерения 18–52 %. Два МК
согласны; расхождение с измерением принадлежит ОПИСАНИЮ СЛОЁВ.

Остаток ниже 100 кэВ (наш выше на 9–25 % у края, ниже на 3–8 % в 85–240)
— кандидаты: раскладка засыпки, сечения в ОИСН-16 (71 % железа),
определение площади пика. Для точечной геометрии этих факторов нет.

Запуск:
    python compare_effcalcmc.py [файл.efa] [наша_кривая.csv]
По умолчанию — nuclidemaster/EffReg_G1S_Marinelli.efa против
results/eff_rho1.60.csv.
"""
import csv
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

sys.path.insert(0, str(paths.tools()))
from fetch_efr import parse_efr  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_point import zoned_fit  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DET = os.path.join(HERE, "..")


def main(argv):
    efa = argv[0] if argv else os.path.join(
        DET, "nuclidemaster", "EffReg_G1S_Marinelli.efa")
    ours = argv[1] if len(argv) > 1 else os.path.join(
        DET, "results", "eff_rho1.60.csv")
    # .efa ЛСРМ — cp1251
    pts = parse_efr(open(efa, encoding="cp1251", errors="replace").read())
    pts = [p for s in pts for p in s["points"]] if isinstance(pts, list) and \
        pts and isinstance(pts[0], dict) else pts[0]["points"]

    rows = list(csv.DictReader(open(ours, encoding="utf-8", newline="")))
    E = np.array([float(r["E_keV"]) for r in rows])
    i = np.argsort(E)
    E = E[i]
    y = np.array([float(r["eps_net"]) for r in rows])[i]
    dy = np.array([float(r["d_eps"]) for r in rows])[i]
    _fits, ev = zoned_fit(E, y, dy)

    print("сверка: %s\n против %s" % (os.path.basename(efa),
                                      os.path.basename(ours)))
    print("%8s %11s %11s %9s" % ("E, кэВ", "EffCalcMC", "наш МК", "наш/ECM"))
    lr = []
    for Ee, eff, _dp, _n in sorted(pts):
        if Ee < E[0] or Ee > E[-1]:
            continue
        m = ev(Ee)
        lr.append((Ee, m / eff))
        print("%8.1f %11.4e %11.4e %9.3f" % (Ee, eff, m, m / eff))
    if len(lr) > 2:
        x = np.log([e for e, _ in lr])
        z = np.log([r for _, r in lr])
        b = np.polyfit(x, z, 1)[0]
        print("\nсреднее наш/ECM = %.3f, наклон d lnR/d lnE = %+.4f,"
              " перепад по диапазону %.1f %%"
              % (math.exp(z.mean()), b,
                 100 * (math.exp(b * (x[-1] - x[0])) - 1)))


if __name__ == "__main__":
    main(sys.argv[1:])
