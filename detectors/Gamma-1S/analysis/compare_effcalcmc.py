"""Сверка нашей кривой Geant4 с расчётом EffCalcMC на ТЕХ ЖЕ слоях.

Зачем. У модели есть систематический ход по энергии против ИЗМЕРЕННОЙ кривой
ЛСРМ (REPORT §2, §5а): занижение мягкого края, завышение жёсткого. Ход может
принадлежать либо коду (розыгрыш, съём площади, телесный угол), либо описанию
слоёв (толщины/плотности против реального прибора). Расчёт EffCalcMC по нашим
же .din/.sin (nuclidemaster/) разделяет гипотезы.

НАПРАВЛЕНИЕ ИНТЕРПОЛЯЦИИ (урок от аудитора). Интерполируется ГУСТАЯ кривая
(ECM, 50 лог-точек) в узлы РЕДКОЙ (наша сетка, 22–24 узла) локальной
квадратикой в log-log — ошибка такой интерполяции ничтожна. Первая редакция
делала наоборот: тянула зонную аппроксимацию нашей сетки в точки ECM, и
аппроксимация сглаживала ровно тот мягкий край, где вся структура.

ОДНИМ ЧИСЛОМ НЕ СВОДИТСЯ — печатаются срезы. Итог первого прогона
(маринелли, 10^8, 28.07.2026): отношение наш/ECM НЕ плоское —
+20 % на трёх узлах ниже 70 кэВ, −7 % в 122–166, монотонный рост ~+10 %
на декаду выше 90 кэВ. Прямая по всем точкам даёт около нуля лишь потому,
что мягкое смещение гасит жёсткий рост — то же маскирующее среднее, что
и χ²/ν поверх несогласного набора. Согласие кодов ±3 % — ТОЛЬКО выше
90 кэВ; главному выводу это не мешает: расхождение кодов на порядок меньше
измеренного хода и противоположно ему по знаку.

Скрипту нужны ДВА файла, оба в репозитории; каталог расчётных спектров не
нужен и не проверяется.

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

HERE = os.path.dirname(os.path.abspath(__file__))
DET = os.path.join(HERE, "..")

# parse_efr — из tools/ репозитория; paths для этого не нужен
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "tools"))
from fetch_efr import parse_efr  # noqa: E402

sys.path.insert(0, HERE)
from curvefit import local_quad  # noqa: E402


def slices(pairs):
    """Срезы (метка, [пары]): весь набор и жёсткая часть без мягкого края."""
    return [("все узлы", pairs),
            ("E >= 90", [p for p in pairs if p[0] >= 90.0])]


def main(argv):
    efa = argv[0] if argv else os.path.join(
        DET, "nuclidemaster", "EffReg_G1S_Marinelli.efa")
    ours = argv[1] if len(argv) > 1 else os.path.join(
        DET, "results", "eff_rho1.60.csv")
    # .efa ЛСРМ — cp1251
    secs = parse_efr(open(efa, encoding="cp1251", errors="replace").read())
    pts = sorted(p for s in secs for p in s["points"])
    Ec = [p[0] for p in pts]
    yc = [p[1] for p in pts]
    ev = local_quad(Ec, yc)

    rows = list(csv.DictReader(open(ours, encoding="utf-8", newline="")))
    grid = sorted((float(r["E_keV"]), float(r["eps_net"]), float(r["d_eps"]))
                  for r in rows)

    print("сверка: %s (%d точек, интерполируется)\n против %s (узлы)"
          % (os.path.basename(efa), len(pts), os.path.basename(ours)))
    print("%8s %11s %11s %9s %7s" % ("E, кэВ", "наш узел", "ECM интерп.",
                                     "наш/ECM", "стат.%"))
    pairs = []
    for E, y, dy in grid:
        if E < Ec[0] or E > Ec[-1]:
            print("%8.1f %11.4e   вне сетки ECM" % (E, y))
            continue
        m = ev(E)
        pairs.append((E, y / m))
        print("%8.1f %11.4e %11.4e %9.3f %7.2f"
              % (E, y, m, y / m, 100 * dy / y))

    for name, sel in slices(pairs):
        if len(sel) < 3:
            continue
        x = np.log([e for e, _ in sel])
        z = np.log([r for _, r in sel])
        b = np.polyfit(x, z, 1)[0]
        sko = float(np.std(np.exp(z - z.mean()) - 1))
        print("\n   %-10s (%2d узлов): среднее %.3f, наклон %+.1f %% на "
              "декаду, перепад %+.1f %%, СКО %.1f %%"
              % (name, len(sel), math.exp(z.mean()),
                 100 * (math.exp(b * math.log(10)) - 1),
                 100 * (math.exp(b * (x[-1] - x[0])) - 1), 100 * sko))


if __name__ == "__main__":
    main(sys.argv[1:])
