"""Устойчива ли аппроксимация на жёстком крае: две зоны против трёх.

ЗАЧЕМ. Прежде чем объяснять расхождение на 2614,5 кэВ физикой (каскадное
суммирование в самой аттестации), надо убедиться, что само модельное значение
на краю не создано аппроксимацией. Сейчас сшивка двухзонная, и вторая зона
тянется от 338 до 3552 кэВ — декада с лишним на полиноме второй степени. На
её дальнем краю полином ничем не удерживается, кроме собственной формы.

Основание для подозрения не умозрительное: на МЯГКОМ крае зонная сшивка уже
однажды сгладила структуру — дала −1,2 % на декаду там, где сырые узлы
показывали −5,3 %. Тот же механизм на жёстком крае дал бы вклад в
«расхождение», которое затем лечили бы физической поправкой. Поправка,
наложенная на артефакт аппроксимации, увела бы вывод дальше от истины, чем
его отсутствие.

КРИТЕРИЙ — тот же, что при выборе нынешней раскладки: RMS ошибки предсказания
на скользящем исключении узла (leave-one-out). χ² подгонки не годится, его
занижают лишние степени свободы: полином, проведённый через все узлы, имеет
нулевую невязку и нулевую предсказательную силу.

ЧТО СЧИТАЕТСЯ РЕШАЮЩИМ. Сравнивается не качество подгонки само по себе, а
СДВИГ аппроксимированного значения на 2614,5 кэВ при переходе от двух зон к
трём. Сдвиг больше 1…2 % означает, что часть «жёсткого края» сидела в
аппроксимации и должна быть снята до физических поправок.

Прогонов не требует: считает по уже готовым узлам сетки.
"""
import glob
import math
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from curvefit import ZONES, zoned_fit  # noqa: E402

WIN = 6.0
EDGE = 2614.51          # линия, ради которой всё затевалось
TAGS = ["p5cm", "p25cm"]

# Кандидаты раскладок. Первая — нынешняя. Дальше трёхзонные с отдельной
# высокоэнергетической зоной; граница перебирается в диапазоне, названном
# оператором, перекрытие берётся широким, как у действующей пары.
def candidates():
    out = [("2 зоны (нынешняя)", list(ZONES))]
    for cut in (900.0, 1000.0, 1100.0):
        for deg_hi in (2, 3):
            for ov in (300.0, 450.0):
                z = [(None, 661.7, 5),
                     (338.3, cut + ov, 2),
                     (cut, None, deg_hi)]
                out.append(("3 зоны: граница %.0f, перекрытие %.0f, степень %d"
                            % (cut, ov, deg_hi), z))
    return out


def load(p):
    """Формат заголовка — N_primaries, а не N. Свой парсер тут уже дал
    молча пустую кривую и «сетка не готова» при готовой сетке, поэтому
    читается тем же кодом, что и везде."""
    hist, N = {}, None
    for line in open(p, encoding="utf-8"):
        if line.startswith("#"):
            if "N_primaries" in line:
                N = int(line.split("=")[1])
            continue
        if line and line[0].isdigit():
            e, c = line.split(",")
            hist[float(e)] = int(c)
    return hist, N


def mc_curve(tag):
    build = str(paths.build("Gamma-1S"))
    saf = os.path.join(build, "grid", "%s_solidangle.txt" % tag)
    if not os.path.exists(saf):
        return {}
    frac = float(open(saf).read().strip())
    out = {}
    for p in glob.glob(os.path.join(build, "grid", tag + "_E*.csv")):
        m = re.search(r"_E(\d+\.\d)\.csv$", p)
        if not m:
            continue
        E = float(m.group(1))
        hist, N = load(p)
        peak = sum(c for e, c in hist.items() if abs(e - E) <= WIN)
        if peak > 0 and N:
            out[round(E, 1)] = ((peak / N) * frac, math.sqrt(peak) / N * frac)
    return out


def loo_rms(Eg, yg, dyg, zones):
    """RMS относительной ошибки предсказания узла, исключённого из подгонки."""
    err = []
    for i in range(len(Eg)):
        keep = np.ones(len(Eg), bool)
        keep[i] = False
        if keep.sum() < 8:
            continue
        try:
            _f, ev = zoned_fit(Eg[keep], yg[keep], dyg[keep], zones)
            pred = ev(float(Eg[i]))
        except Exception:
            continue
        if pred > 0:
            err.append((pred - yg[i]) / yg[i])
    if not err:
        return None
    return math.sqrt(sum(e * e for e in err) / len(err))


def main():
    for tag in TAGS:
        mc = mc_curve(tag)
        if not mc:
            print("%-6s сетка не готова" % tag)
            continue
        Eg = np.array(sorted(mc))
        yg = np.array([mc[e][0] for e in Eg])
        dyg = np.array([mc[e][1] for e in Eg])
        print("\n===== %s: узлов %d, диапазон %.1f…%.1f кэВ"
              % (tag, len(Eg), Eg[0], Eg[-1]))
        print("%-52s %10s %14s %9s"
              % ("раскладка", "LOO RMS", "ε(2614,5)", "сдвиг"))

        base = None
        for name, zones in candidates():
            try:
                _fits, ev = zoned_fit(Eg, yg, dyg, zones)
                at_edge = ev(EDGE)
            except Exception as ex:
                print("%-52s   не считается: %s" % (name, ex))
                continue
            rms = loo_rms(Eg, yg, dyg, zones)
            if base is None:
                base = at_edge
                shift = 0.0
            else:
                shift = 100.0 * (at_edge - base) / base
            print("%-52s %9s %14.5e %+8.2f %%"
                  % (name, "%.2f %%" % (100 * rms) if rms else "—",
                     at_edge, shift))

        # Узел 2614,5 есть в сетке — с ним и сверяется аппроксимация.
        key = min(mc, key=lambda k: abs(k - EDGE))
        if abs(key - EDGE) < 2.0:
            print("сырой узел сетки на %.1f кэВ: %.5e" % (key, mc[key][0]))
            print("(аппроксимация обязана его воспроизводить; расхождение с"
                  " ним — мера того,")
            print(" насколько полином гнёт край ради остального диапазона)")


if __name__ == "__main__":
    main()
