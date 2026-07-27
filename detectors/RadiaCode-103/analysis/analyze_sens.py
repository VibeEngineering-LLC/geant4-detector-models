# -*- coding: utf-8 -*-
"""Сверка модели с паспортной чувствительностью: 30 cps на 1 мкЗв/ч по Cs-137.

Пучок 662 кэВ, диск радиусом 40 мм (площадь 50.265 см²). Скорость счёта
    N = (счёт/первичных) * Ф * A,
где Ф — плотность потока, дающая H*(10) = 1 мкЗв/ч. По ICRP 74 на 662 кэВ
h*(10) = 4.13 пЗв·см², откуда Ф = 1e6/3600/4.13 = 67.3 см⁻²·с⁻¹.
"""
import os

import numpy as np

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

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = str(paths.build("RadiaCode-103"))

if not os.path.isdir(BUILD):
    raise SystemExit(
        "Нет каталога расчётных спектров %s.\n"
        "Они не коммитятся (сотни файлов), а воспроизводятся драйверами:\n"
        "    python detectors/RadiaCode-103/drivers/run_grid.py\n"
        "Либо укажите G4MODELS_BUILD_RADIACODE_103 на готовый каталог."
        % BUILD)

AREA = np.pi * 4.0 ** 2          # см², диск радиусом 40 мм
H10_662 = 4.13                   # пЗв·см², ICRP 74 с интерполяцией
PHI_1USVH = 1e6 / 3600.0 / H10_662

DIRS = [("sens_cs137_front.csv", "фронт (дисплей, кристалл ближе)"),
        ("sens_cs137_back.csv", "тыл (задняя крышка)"),
        ("sens_cs137_nose.csv", "торец (нос)")]


def main():
    print("плотность потока для 1 мкЗв/ч на 662 кэВ: %.1f см⁻²·с⁻¹" % PHI_1USVH)
    print("площадь пучка: %.2f см²\n" % AREA)
    print("%-34s %10s %12s %10s" % ("направление", "порог 0", "порог 20 кэВ",
                                    "порог 50"))
    for fn, title in DIRS:
        p = os.path.join(BUILD, fn)
        if not os.path.exists(p):
            print("%-34s нет файла" % title)
            continue
        meta, hist = rcspec.read_spec(p)
        n = float(meta["N_primaries"])
        row = []
        for thr in (0, 20, 50):
            cps = hist[thr:].sum() / n * PHI_1USVH * AREA
            row.append(cps)
        print("%-34s %10.1f %12.1f %10.1f" % (title, *row))
    print("\nпаспорт RadiaCode: 30 cps на 1 мкЗв/ч (Cs-137)")


if __name__ == "__main__":
    main()
