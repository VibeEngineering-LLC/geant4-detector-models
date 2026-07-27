# -*- coding: utf-8 -*-
"""Разбор спектра поля ЕРН: проверка по мощности дозы, доля рассеянных,
генерация макроса-источника для прогонов фона.

Проверка: для обычного помещения H*(10) должна получиться 0.1..0.2 мкЗв/ч.
Если модель поля даёт это без всякой подгонки, значит и активности, и геометрия
стены, и нормировка сходятся.
"""
import os
import sys

import numpy as np

# Корни путей — из переменных окружения (common/py/paths.py).
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.abspath(os.path.join(HERE, "..", "results"))

# ICRP 74, монолинейные коэффициенты перехода флюенс -> H*(10), пЗв·см²
ICRP74 = np.array([
    [0.010, 0.061], [0.015, 0.83], [0.020, 1.05], [0.030, 0.81],
    [0.040, 0.64], [0.050, 0.55], [0.060, 0.51], [0.080, 0.61],
    [0.100, 0.83], [0.150, 1.20], [0.200, 1.63], [0.300, 2.34],
    [0.400, 3.01], [0.500, 3.58], [0.600, 4.00], [0.800, 4.79],
    [1.000, 5.44], [1.500, 6.99], [2.000, 8.31], [3.000, 10.5],
    [4.000, 12.5], [5.000, 14.4], [6.000, 16.3], [8.000, 20.4],
    [10.00, 24.5],
])

# Линии, попадающие в поле (для разделения первичных и рассеянных)
LINES_KEV = [186.2, 238.6, 241.9, 295.2, 338.3, 351.9, 463.0, 510.7, 583.2,
             609.3, 727.3, 794.9, 860.6, 911.2, 964.8, 968.9, 1120.3, 1238.1,
             1377.7, 1408.0, 1460.8, 1509.2, 1588.2, 1729.6, 1764.5, 1847.4,
             2204.2, 2614.5]


def h10(e_keV):
    """Интерполяция коэффициента h*(10) в логарифмических координатах."""
    x = np.log(np.clip(e_keV / 1000.0, ICRP74[0, 0], ICRP74[-1, 0]))
    return np.exp(np.interp(x, np.log(ICRP74[:, 0]), np.log(ICRP74[:, 1])))


def read(path):
    meta, e, f = {}, [], []
    for line in open(path, encoding="utf-8"):
        if line.startswith("#"):
            if "=" in line:
                k, v = line[1:].split("=", 1)
                meta[k.strip()] = v.strip()
        elif line[0].isdigit():
            a, b = line.split(",")
            e.append(float(a))
            f.append(float(b))
    return meta, np.array(e), np.array(f)


def main():
    # Карта поля — РЕЗУЛЬТАТ (её считает wallfield.exe), поэтому она
    # коммитится в results/. В каталоге расчётов ищем только как запасной
    # вариант, для свежепосчитанной.
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        RESULTS, "wallfield.csv")
    if not os.path.exists(path):
        path = os.path.join(str(paths.build("RadiaCode-103")), "wallfield.csv")
    meta, e, flu = read(path)
    tot = flu.sum()

    # мощность дозы
    dose = (flu * h10(e)).sum() * 1e-6 * 3600      # пЗв/с -> мкЗв/ч
    print("полный флюенс:      %.2f см⁻²·с⁻¹" % tot)
    print("H*(10):             %.3f мкЗв/ч   (обычное помещение 0.10..0.20)" % dose)
    print("средняя энергия:    %.0f кэВ" % ((flu * e).sum() / tot))

    # первичные (в каналах с линиями) против рассеянных
    bw = e[1] - e[0] if len(e) > 1 else 10.0
    mask = np.zeros_like(e, dtype=bool)
    for L in LINES_KEV:
        mask |= np.abs(e - L) < bw
    print("в каналах с линиями: %.1f %%   вне линий (рассеяние): %.1f %%"
          % (100 * flu[mask].sum() / tot, 100 * flu[~mask].sum() / tot))
    lo = e < 400
    print("доля флюенса < 400 кэВ: %.1f %%  — здесь у CsI максимум отклика"
          % (100 * flu[lo].sum() / tot))

    # Энергетический спектр источника отдельным макросом: геометрию поверхности
    # задаёт уже прогон (run_bg.py), потому что от неё зависит нормировка.
    out = os.path.join(RESULTS, "field_spectrum.mac")
    os.makedirs(RESULTS, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        f.write("# Спектр поля ЕРН помещения, посчитан wallfield.exe.\n")
        f.write("# FLUENCE_TOTAL_CM2_S = %.6f\n" % tot)
        f.write("# H10_USV_PER_H = %.6f\n" % dose)
        f.write("/gps/particle gamma\n/gps/ene/type Arb\n/gps/hist/type arb\n")
        for ei, fi in zip(e, flu):
            if fi > 0:
                f.write("/gps/hist/point %.4f %.6e\n" % (ei / 1000.0, fi))
        f.write("/gps/hist/inter Lin\n")
    print("\nспектр источника:", out)


if __name__ == "__main__":
    main()
