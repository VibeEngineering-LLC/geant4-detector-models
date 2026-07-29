"""Сколько «жёсткого края» сидит в аппроксимации: полный перебор разбивок.

ЗАЧЕМ. Прежде чем объяснять расхождение на 2614,5 кэВ физикой, надо знать,
определено ли само модельное значение на краю данными — или выбором разбивки.
Нынешняя сшивка двухзонная, вторая зона тянется от 338 до 3552 кэВ: декада с
лишним на полиноме второй степени, и на дальнем краю его ничто не удерживает.
Основание для подозрения не умозрительное — на мягком крае зонная сшивка уже
сглаживала структуру.

ПЕРЕБИРАЕТСЯ ВСЁ, что раньше выбиралось человеком: число зон (2…5), положение
границ (по узлам), степень каждой зоны, ширина перекрытия. Решает не мнение, а
критерий.

КРИТЕРИЙ — RMS ошибки предсказания на скользящем исключении узла. Ошибка
ПОДГОНКИ для выбора негодна и печатается рядом только как индикатор: полином
степени N через N+1 точку даёт структурно нулевую невязку и нулевую
предсказательную силу. Сильное расхождение двух ошибок — прямой признак
переобучения.

ОГРАНИЧИТЕЛЬ ПЕРЕОБУЧЕНИЯ: в каждой зоне не меньше (степень + 2) узлов,
считая БЕЗ перекрытия. Иначе перебор находит разбивку, где каждая зона
проходит ровно через свои точки.

ГЛАВНОЕ — НЕ ПОБЕДИТЕЛЬ, А РАЗБРОС. Когда границы подбираются по LOO, сам
выбор границы становится подгонкой, и лучший вариант слегка переобучен этим
выбором. Поэтому решающий вывод даёт не одна лучшая раскладка, а разброс
значения на 2614,5 среди топ-вариантов, почти равных по LOO:

  * разброс мал  — край определён данными, аппроксимации можно верить, и
    физические поправки лягут на измеренную величину;
  * разброс велик — край текущей сеткой узлов не определён вовсе. Тогда спор
    о нём решается добавлением узлов выше 1500 кэВ, а не выбором разбивки, и
    любая физическая поправка на краю преждевременна.

Прогонов не требует: считает по уже готовым узлам сетки.
"""
import glob
import itertools
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
EDGE = 2614.51
TAGS = ["p5cm", "p25cm"]
NZONES = (2, 3, 4, 5)
DEGREES = (2, 3, 4, 5)
OVERLAPS = (0.25, 0.45)     # доля декады, на которую зоны заходят друг в друга
TOPN = 10
MIN_NODES_OVER_DEG = 2      # узлов в зоне не меньше, чем степень + это


def load(p):
    """Заголовок содержит N_primaries, а не N. Свой парсер здесь уже дал
    молча пустую кривую и «сетка не готова» при готовой сетке."""
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


def build_zones(Eg, cuts, degs, ov):
    """Границы -> список зон с перекрытием. cuts — индексы узлов-разделителей."""
    edges = [None] + [float(Eg[i]) for i in cuts] + [None]
    zones = []
    for k, deg in enumerate(degs):
        lo, hi = edges[k], edges[k + 1]
        # перекрытие: зона заползает в соседнюю на долю ov по логарифму
        if lo is not None:
            lo = lo * math.exp(-ov)
        if hi is not None:
            hi = hi * math.exp(ov)
        zones.append((lo, hi, deg))
    return zones


def zone_counts(Eg, cuts):
    """Сколько узлов в каждой зоне БЕЗ учёта перекрытия."""
    bnds = [0] + list(cuts) + [len(Eg)]
    return [bnds[i + 1] - bnds[i] for i in range(len(bnds) - 1)]


def fit_errors(Eg, yg, dyg, zones):
    """(RMS подгонки, RMS предсказания LOO, значение на краю)."""
    try:
        _f, ev = zoned_fit(Eg, yg, dyg, zones)
        at_edge = ev(EDGE)
        fit_err = math.sqrt(sum(((ev(float(e)) - y) / y) ** 2
                                for e, y in zip(Eg, yg)) / len(Eg))
    except Exception:
        return None
    err = []
    for i in range(len(Eg)):
        keep = np.ones(len(Eg), bool)
        keep[i] = False
        try:
            _f2, ev2 = zoned_fit(Eg[keep], yg[keep], dyg[keep], zones)
            pred = ev2(float(Eg[i]))
        except Exception:
            return None
        if pred <= 0 or not np.isfinite(pred):
            return None
        err.append((pred - yg[i]) / yg[i])
    loo = math.sqrt(sum(e * e for e in err) / len(err))
    return fit_err, loo, at_edge


def main():
    for tag in TAGS:
        mc = mc_curve(tag)
        if not mc:
            print("%-6s сетка не готова" % tag)
            continue
        Eg = np.array(sorted(mc))
        yg = np.array([mc[e][0] for e in Eg])
        dyg = np.array([mc[e][1] for e in Eg])
        n = len(Eg)
        print("\n" + "=" * 78)
        print("%s: узлов %d, диапазон %.1f…%.1f кэВ" % (tag, n, Eg[0], Eg[-1]))

        results, skipped = [], 0
        for nz in NZONES:
            # границы — между узлами; индекс cut означает «зона кончается перед
            # узлом cut». Комбинаторика ограничена требованием узлов на зону.
            for cuts in itertools.combinations(range(1, n), nz - 1):
                for degs in itertools.product(DEGREES, repeat=nz):
                    cnt = zone_counts(Eg, cuts)
                    if any(c < d + MIN_NODES_OVER_DEG
                           for c, d in zip(cnt, degs)):
                        skipped += 1
                        continue
                    for ov in OVERLAPS:
                        z = build_zones(Eg, cuts, degs, ov)
                        r = fit_errors(Eg, yg, dyg, z)
                        if r is None:
                            skipped += 1
                            continue
                        fit_err, loo, edge = r
                        results.append((loo, fit_err, edge, nz, cuts, degs, ov))

        if not results:
            print("  ни один вариант не прошёл ограничитель узлов")
            continue
        results.sort(key=lambda t: t[0])
        print("  вариантов посчитано %d, отброшено ограничителем %d"
              % (len(results), skipped))

        # Нынешняя раскладка — для сравнения.
        cur = fit_errors(Eg, yg, dyg, list(ZONES))
        if cur:
            print("  нынешняя (2 зоны, 5/2): LOO %.2f %%, подгонка %.2f %%,"
                  " ε(2614,5) %.5e" % (100 * cur[1], 100 * cur[0], cur[2]))

        print("\n  ТОП-%d по предсказанию:" % TOPN)
        print("  %-4s %-22s %-14s %6s %8s %8s %13s"
              % ("зон", "границы, кэВ", "степени", "перекр", "LOO", "подгонка",
                 "ε(2614,5)"))
        top = results[:TOPN]
        for loo, fit_err, edge, nz, cuts, degs, ov in top:
            b = "/".join("%.0f" % Eg[i] for i in cuts) or "—"
            d = "/".join(str(x) for x in degs)
            print("  %-4d %-22s %-14s %6.2f %7.2f%% %7.2f%% %13.5e"
                  % (nz, b, d, ov, 100 * loo, 100 * fit_err, edge))

        # РЕШАЮЩЕЕ: разброс края среди почти равных по LOO.
        edges = [t[2] for t in top]
        spread = 100.0 * (max(edges) - min(edges)) / np.mean(edges)
        best_loo, worst_loo = top[0][0], top[-1][0]
        print("\n  разброс ε(2614,5) по топ-%d: %.2f %% "
              "(LOO при этом от %.2f %% до %.2f %%)"
              % (len(top), spread, 100 * best_loo, 100 * worst_loo))

        key = min(mc, key=lambda k: abs(k - EDGE))
        if abs(key - EDGE) < 2.0:
            raw = mc[key][0]
            print("  сырой узел сетки на %.1f кэВ: %.5e" % (key, raw))
            print("  топ-варианты отклоняются от сырого узла на %.2f … %.2f %%"
                  % (100 * (min(edges) - raw) / raw,
                     100 * (max(edges) - raw) / raw))

        if spread < 1.0:
            print("  ВЫВОД: край определён данными — почти равные по"
                  " предсказанию разбивки")
            print("  дают одно значение. Физические поправки лягут на"
                  " измеренную величину.")
        else:
            print("  ВЫВОД: край разбивкой НЕ определён — равные по"
                  " предсказанию варианты")
            print("  расходятся на %.1f %%. Спор о жёстком крае решается"
                  " добавлением узлов" % spread)
            print("  сетки выше 1500 кэВ, а не выбором разбивки; поправки на"
                  " краю преждевременны.")


if __name__ == "__main__":
    main()
