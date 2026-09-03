# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import csv
import argparse

sys.stdout.reconfigure(encoding="utf-8")

TEMPS = {"T_m40": -40, "T_0": 0, "T_p40": 40}
MATERIALS = ["CsI_Tl", "NaI_Tl"]

def load(material, label):
    path = os.path.join(os.path.dirname(__file__), "out", f"payne3_fig2_{material}_{label}.csv")
    data = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if row[0] == "E_keV":
                continue
            try:
                # CSV несёт 2–4 колонки; берём первые две (E, выход).
                e, y = float(row[0]), float(row[1])
                data.append((e, y))
            except (ValueError, IndexError):
                continue
    arr = np.array(data)
    idx = np.argsort(arr[:, 0])
    return arr[idx].T

def interp_at_T(y_dict, T):
    # Температура -> метка (ключ y_dict), а не строка "T_<число>": метки фиксированы.
    by_temp = {temp: label for label, temp in TEMPS.items()}
    y = {label: np.asarray(arr, dtype=float) for label, arr in y_dict.items()}
    ts = sorted(TEMPS.values())
    if T <= ts[0]:
        return y[by_temp[ts[0]]]
    if T >= ts[-1]:
        return y[by_temp[ts[-1]]]
    i = 1 if T > ts[1] else 0
    t0, t1 = ts[i], ts[i + 1]
    w = (T - t0) / (t1 - t0)
    return (1 - w) * y[by_temp[t0]] + w * y[by_temp[t1]]

def selftest():
    grid = np.logspace(np.log10(10), np.log10(100), 2)
    y = {
        "T_m40": [1.30, 1.00],
        "T_0": [1.20, 0.99],
        "T_p40": [1.10, 0.98]
    }
    # Проверка, что interp_at_T(y, 0) возвращает y[T_0]
    result = interp_at_T(y, 0)
    if not np.allclose(result, y["T_0"], atol=1e-9):
        print("SELFTEST FAIL: interp_at_T(y, 0) не равен T_0")
        return False
    # Проверка, что interp_at_T(y, 20) дает полусумму
    result = interp_at_T(y, 20)
    expected = 0.5 * (np.array(y["T_0"]) + np.array(y["T_p40"]))
    if not np.allclose(result, expected, atol=1e-9):
        print("SELFTEST FAIL: interp_at_T(y, 20) не равен полусумме")
        return False
    # Проверка, что interp_at_T(y, 26) не совпадает с T_0
    result = interp_at_T(y, 26)
    if np.allclose(result, y["T_0"], atol=1e-9):
        print("SELFTEST FAIL: interp_at_T(y, 26) совпадает с T_0")
        return False
    print("SELFTEST OK")
    return True

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--T", type=float, default=26)
    parser.add_argument("--material", choices=MATERIALS, default="CsI_Tl")
    parser.add_argument("--nodes", type=int, default=60)
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args(argv)

    if args.selftest:
        return 0 if selftest() else 1

    material = args.material
    T = args.T
    nodes = args.nodes

    # Загрузка данных
    y_dict = {}
    for label in ["T_m40", "T_0", "T_p40"]:
        data = load(material, label)
        E, Y = data[0], data[1]
        y_dict[label] = Y
        if not np.all(np.diff(E) > 0):
            print(f"Ошибка: энергии не отсортированы в {label}")
            return 1

    # Объединение диапазонов
    E_m40, Y_m40 = load(material, "T_m40")
    E_0, Y_0 = load(material, "T_0")
    E_p40, Y_p40 = load(material, "T_p40")

    lo = max(E_m40[0], E_0[0], E_p40[0])
    hi = min(E_m40[-1], E_0[-1], E_p40[-1])

    if lo >= hi:
        print("Ошибка: нет пересечения диапазонов")
        return 1

    grid = np.logspace(np.log10(lo), np.log10(hi), nodes)

    # Интерполяция по log10(E)
    Y_m40_interp = np.interp(np.log10(grid), np.log10(E_m40), Y_m40)
    Y_0_interp = np.interp(np.log10(grid), np.log10(E_0), Y_0)
    Y_p40_interp = np.interp(np.log10(grid), np.log10(E_p40), Y_p40)

    # Проверка линейности
    mid = 0.5 * (Y_m40_interp + Y_p40_interp)
    dev = np.abs(mid - Y_0_interp)
    max_dev = np.max(dev)
    mean_dev = np.mean(dev)
    threshold = 0.005

    print(f"Максимальное отклонение: {max_dev:.6f} при E = {grid[np.argmax(dev)]:.1f} кэВ")
    print(f"Среднее отклонение: {mean_dev:.6f}")
    if max_dev <= threshold:
        print("ЛИНЕЙНА")
    else:
        print("НЕ линейна")

    # Интерполяция по температуре
    y_interp = interp_at_T({
        "T_m40": Y_m40_interp,
        "T_0": Y_0_interp,
        "T_p40": Y_p40_interp
    }, T)

    # Запись результата
    out_dir = os.path.join(os.path.dirname(__file__), "out")
    os.makedirs(out_dir, exist_ok=True)
    filename = f"payne3_fig2_{material}_T{int(round(T)):+d}C.csv"
    path = os.path.join(out_dir, filename)

    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write(f"# Температура: {T} °C\n")
        f.write("# Способ: линейная интерполяция по значениям\n")
        f.write(f"# Проверка линейности: max|откл.| = {max_dev:.6f}, порог = {threshold}\n")
        f.write("# Источник: Payne S.A. et al., IEEE TNS 61 (2014) 2771, OSTI 1762905\n")
        f.write("# Оцифровка графика, не табличные данные автора\n")
        f.write("E_keV,rel_light_yield\n")
        for e, y in zip(grid, y_interp):
            f.write(f"{e:.6f},{y:.6f}\n")

    # Вывод таблицы
    print("\nТаблица значений:")
    print("E кэВ\tT=-40\tT=0\tT=+40\tT={T}°C".format(T=T))
    E_vals = [10, 20, 40, 80, 150, 300, 430]
    for e in E_vals:
        if e < lo or e > hi:
            continue
        i = np.argmin(np.abs(grid - e))
        y_m40 = Y_m40_interp[i]
        y_0 = Y_0_interp[i]
        y_p40 = Y_p40_interp[i]
        y_T = y_interp[i]
        print(f"{e}\t{y_m40:.6f}\t{y_0:.6f}\t{y_p40:.6f}\t{y_T:.6f}")

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
