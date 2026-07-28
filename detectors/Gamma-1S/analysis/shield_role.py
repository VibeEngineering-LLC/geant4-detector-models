# -*- coding: utf-8 -*-
"""Сколько защита добавляет в пик полного поглощения.

ЗАЧЕМ. Два разных вопроса упираются в одно число.

1. В точечных прогонах источник разыгрывается КОНУСОМ на детектор, а результат
   делится на долю телесного угла. Кванты вне конуса не рождаются, значит путь
   «вылетел в сторону — рассеялся на свинце — вернулся в кристалл» не
   воспроизводится. Насколько это занижает ППП?
2. Опорные расчёты сторонних кодов защиту могут не моделировать вовсе. Тогда
   сравнивать наш расчёт с их расчётом можно только зная этот вклад.

КАК. Один и тот же точечный источник на 5 см, полный 4π, два прогона:
режим shield (защита собрана) и bare (устройство детектирования в воздухе).
Геометрия детектора, расстояние и статистика одинаковы, отличается только
наличие защиты. Разница площадей пика — целиком её вклад.

    python detectors/Gamma-1S/analysis/shield_role.py

Прогоны делаются макросами scat_test.mac (shield) и scat_bare.mac (bare) —
см. REPORT, раздел про роль защиты.
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))
WIN, BG0, BG1 = 6.0, 30.0, 10.0
CASES = [(122.1, "E122.1"), (661.657, "E661.7"), (2614.511, "E2614.5")]


def area(path, E):
    N, hist = None, {}
    for ln in open(path, encoding="utf-8"):
        if ln.startswith("#"):
            if "N_primaries" in ln:
                N = int(ln.split("=")[1])
            continue
        if not ln[:1].isdigit():
            continue
        e, c = ln.split(",")
        hist[float(e)] = int(c)
    n = math.floor(E + WIN - 0.5) - math.ceil(E - WIN - 0.5) + 1
    ns = math.floor(E - BG1 - 0.5) - math.ceil(E - BG0 - 0.5) + 1
    gross = sum(c for e, c in hist.items() if abs(e - E) <= WIN)
    side = sum(c for e, c in hist.items() if E - BG0 <= e <= E - BG1)
    bg = side / ns * n
    net = gross - bg
    d = math.sqrt(max(gross + (n / ns) * bg, 1.0))
    tot = sum(hist.values())
    return N, net, d, tot


def main():
    print("Вклад защиты в пик полного поглощения; точечный источник 5 см, "
          "полный 4π\n")
    print("%9s %12s %12s %9s %12s" %
          ("E, кэВ", "ППП с защитой", "ППП без неё", "отношение", "полный счёт"))
    bad = 0
    for E, tag in CASES:
        ps = os.path.join(BUILD, "scat_p5_full_%s.csv" % tag)
        pb = os.path.join(BUILD, "scat_p5_bare_%s.csv" % tag)
        if not (os.path.exists(ps) and os.path.exists(pb)):
            print("%9.1f  нет пары прогонов (%s)" % (E, tag))
            bad += 1
            continue
        Ns, ns_, dns, ts = area(ps, E)
        Nb, nb_, dnb, tb = area(pb, E)
        if Ns != Nb:
            print("   ВНИМАНИЕ: статистика разная, %d против %d" % (Ns, Nb))
        r = ns_ / nb_ if nb_ else float("nan")
        dr = r * math.hypot(dns / max(ns_, 1), dnb / max(nb_, 1))
        print("%9.1f %12.0f %12.0f  %.4f±%.4f %12s"
              % (E, ns_, nb_, r, dr, "%.4f / %.4f" % (ts / Ns, tb / Nb)))
    if bad:
        print("\nПрогоны: g1s.exe scat_test.mac shield и "
              "g1s.exe scat_bare.mac bare")
        return 1
    print("\nОтношение больше единицы — защита ДОБАВЛЯЕТ в пик обратно\n"
          "рассеянные кванты; меньше — отнимает больше, чем добавляет.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
