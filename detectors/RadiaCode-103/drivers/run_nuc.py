# -*- coding: utf-8 -*-
"""Прогоны реальных нуклидов в пробе: полный распад со всеми продуктами.

Даёт спектр на один распад родителя при равновесии цепочки — с гамма-линиями,
бета-континуумом, конверсионными электронами и тормозным излучением сразу.
Это второй, независимый продукт рядом с кривыми эффективности: свёртка линий по
кривой против полного распада разделяет гамма- и бета-вклад и проверяет кривые.
"""
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor

import sys

# Корни путей — из переменных окружения (common/py/paths.py).
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402
# Модули прибора лежат в двух каталогах: разбор в analysis/, запуск прогонов
# в drivers/. Импорт через каталог-сосед иначе не находится: python кладёт
# в sys.path только каталог запускаемого файла.
for _d in ("analysis", "drivers"):
    _p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import rcspec
from run_grid import VESSELS

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = str(paths.build("RadiaCode-103"))

if not os.path.isdir(BUILD):
    raise SystemExit(
        "Нет каталога расчётных спектров %s.\n"
        "Они не коммитятся (сотни файлов), а воспроизводятся драйверами:\n"
        "    python detectors/RadiaCode-103/drivers/run_grid.py\n"
        "Либо укажите G4MODELS_BUILD_RADIACODE_103 на готовый каталог."
        % BUILD)
EXE = os.path.join(BUILD, "rc_curves.exe")
RESULTS = rcspec.RESULTS

# (имя, Z, A, границы цепочки Amin Amax Zmin Zmax)
# Границы нужны, чтобы цепочка останавливалась там, где кончается равновесие:
# радиевая подцепочка обрывается на Pb-210 (A=210 исключён), урановая голова —
# на U-234 (A<234 исключены).
NUCLIDES = [
    ("K40",   19, 40,  (40, 40, 19, 19)),
    ("Cs137", 55, 137, (137, 137, 55, 56)),
    ("Ra226", 88, 226, (211, 226, 82, 88)),
    ("Th232", 90, 232, (208, 232, 81, 90)),
    ("U238",  92, 238, (234, 238, 90, 92)),
]

CONFIGS = [("full", "water", 1.00), ("full", "soil", 1.60),
           ("full", "organic", 0.49)]
NEV = 1_500_000


def macro(path, outdir, vessel):
    # Область розыгрыша — своя для каждого сосуда, см. пояснение к VESSELS в
    # run_grid.py: цилиндр обязан охватывать ВСЁ тело пробы, иначе распады
    # окажутся только в центральной части и эффективность выйдет завышенной.
    ves = VESSELS[vessel]
    src = ["/run/verbose 0", "/event/verbose 0", "/tracking/verbose 0",
           "/run/printProgress 0",
           # иначе долгоживущее первичное ядро (Ra-226, 1600 лет) будет убито
           # по порогу времени и не распадётся вовсе
           "/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+30 ns",
           "/gps/particle ion", "/gps/energy 0 keV",
           "/gps/pos/type Volume", "/gps/pos/shape Cylinder",
           "/gps/pos/radius %.2f mm" % ves["radius"],
           "/gps/pos/halfz %.2f mm" % ves["halfz"],
           "/gps/pos/centre 0 0 %.2f mm" % ves["centre"],
           "/gps/pos/confine sample",
           "/gps/ang/type iso", ""]
    todo = 0
    for name, z, a, lim in NUCLIDES:
        out = os.path.join(outdir, "nuc_%s.csv" % name)
        if os.path.exists(out):
            continue
        todo += 1
        src.append("/process/had/rdm/nucleusLimits %d %d %d %d" % lim)
        src.append("/gps/ion %d %d 0 0" % (z, a))
        src.append("/rc/outFile %s" % out.replace("\\", "/"))
        src.append("/run/beamOn %d" % NEV)
    open(path, "w", encoding="utf-8").write("\n".join(src) + "\n")
    return todo


def run_one(cfg, vessel):
    mode, matrix, rho = cfg
    name = "%s_%.2f" % (matrix, rho)
    outdir = rcspec.rdir("nuclides", name, v=vessel)
    os.makedirs(outdir, exist_ok=True)
    mac = os.path.join(BUILD, "nuc_%s_%s.mac" % (vessel, name))
    if macro(mac, outdir, vessel) == 0:
        print("[--] %s/%s уже посчитано" % (vessel, name), flush=True)
        return 0
    t0 = time.time()
    with open(os.path.join(outdir, "run.log"), "w", encoding="utf-8") as lf:
        p = subprocess.run([EXE, mac, mode, matrix, "%.6f" % rho, vessel],
                           cwd=BUILD, stdout=lf, stderr=subprocess.STDOUT)
    print("[%s] %s/%s  %.1f мин" % ("ok" if p.returncode == 0 else "СБОЙ",
                                    vessel, name, (time.time() - t0) / 60),
          flush=True)
    return p.returncode


def main():
    v = rcspec.vessel()
    print("сосуд %s, область розыгрыша: R %.2f, halfz %.2f, центр %.2f мм"
          % (v, VESSELS[v]["radius"], VESSELS[v]["halfz"], VESSELS[v]["centre"]))
    with ThreadPoolExecutor(max_workers=len(CONFIGS)) as ex:
        rc = list(ex.map(lambda c: run_one(c, v), CONFIGS))
    print("сбоев:", sum(1 for r in rc if r != 0))


if __name__ == "__main__":
    main()
