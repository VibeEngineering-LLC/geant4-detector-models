"""Разбор целевых прогонов: глубина колодца и плотность MgO.

Печатает, как отношение МК/эксперимент едет при изменении каждого параметра,
и говорит, сводится ли расхождение — или параметр не при чём.
"""
import glob
import json
import math
import os
import re
import sys

# Корни путей — из переменных окружения (common/py/paths.py), чтобы в коде не
# было ни одного пути, привязанного к конкретной машине.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

# parse_efr живёт в инструментах репозитория, а не среди данных
sys.path.insert(0, str(paths.tools()))
from fetch_efr import parse_efr  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))

if not os.path.isdir(BUILD):
    raise SystemExit(
        "Нет каталога расчётных спектров %s.\n"
        "Они не коммитятся (сотни файлов), а воспроизводятся драйверами:\n"
        "    python detectors/Gamma-1S/drivers/run_grid.py\n"
        "    python detectors/Gamma-1S/drivers/run_all_grids.py\n"
        "Либо укажите G4MODELS_BUILD_GAMMA_1S на готовый каталог."
        % BUILD)
OUT = os.path.join(BUILD, "probe2")
REF = str(paths.ref("Gamma-1S"))
WIN = 6.0
FRAC_POINT = (1 - math.cos(math.radians(60.0))) / 2


def load(p):
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


def eff(tag, E, frac=1.0):
    p = os.path.join(OUT, "%s_E%07.1f.csv" % (tag, E))
    if not os.path.exists(p):
        return None
    hist, N = load(p)
    peak = sum(c for e, c in hist.items() if abs(e - E) <= WIN)
    return (peak / N) * frac if (N and peak) else None


def lsrm(fn):
    p = paths.find_data(fn)
    if p is None:
        return None
    p = str(p)
    if not os.path.exists(p):
        return {}
    return {round(x[0], 1): x[1]
            for s in parse_efr(paths.read_text(p))
            for x in s["points"]}


if __name__ == "__main__":
    mar = lsrm("УДС-ГЦ-63х63-USB__SN-01_-_Маринелли.efr")
    pnt = lsrm("УДС-ГЦ-63х63-USB__SN-01_-_Точечная-5см.efr")

    print("=== 1. Глубина колодца маринелли ===")
    print("Опорное: точечная геометрия даёт МК/эксп = 0,971, то есть детектор")
    print("верен; в маринельке при колодце 74 мм было 1,171.\n")
    print("%10s %12s %12s %12s %12s" %
          ("колодец", "662 МК", "662 отн.", "1461 МК", "1461 отн."))
    for w in (74.0, 65.0, 55.0, 45.0):
        tag = "well%.0f" % w
        row = [tag]
        vals = []
        for E in (661.7, 1460.8):
            m = eff(tag, 661.657 if E < 1000 else 1460.822)
            e = mar.get(E)
            vals += [m, (m / e) if (m and e) else None]
        if vals[0] is None:
            continue
        print("%10.0f %12.4e %12.3f %12.4e %12.3f"
              % (w, vals[0], vals[1], vals[2], vals[3]))
    print("\nЕсли отношение садится к 1,0 при правдоподобной глубине —")
    print("расхождение объясняется геометрией сосуда. Если нет — искать дальше.")

    print("\n=== 2. Плотность отражателя MgO, точечная 5 см ===")
    print("Критерий: мягкий край должен выправиться, НЕ ломая середину и")
    print("жёсткий край. Иначе это подгонка, а не уточнение.\n")
    ES = [59.5, 88.0, 122.1, 661.7, 2614.5]
    EF = [59.5, 88.0, 122.1, 661.657, 2614.511]
    print("%8s" % "MgO" + "".join("%12s" % ("%.0f кэВ" % e) for e in ES))
    for m in (1.30, 1.50, 2.00):
        tag = "mgo%.2f" % m
        row = "%8.2f" % m
        ok = False
        for E, Ef in zip(ES, EF):
            v = eff(tag, Ef, FRAC_POINT)
            ref = pnt.get(E)
            row += "%12s" % ("%.3f" % (v / ref) if (v and ref) else "-")
            ok = ok or bool(v)
        if ok:
            print(row)
    print("\n(в клетках — отношение МК/эксперимент)")
