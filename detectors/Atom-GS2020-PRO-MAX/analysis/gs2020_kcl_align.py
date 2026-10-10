import sys
import os
import json
import xml.etree.ElementTree as ET
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

# Входы (константы)
X1 = r"<DOSIM>\Спектры\Atom GS2020 PRO MAX\Референсы\KCl ч 1 л 1085 г (41,5 ч, 05-07.10).xml"
X2 = r"<WORKDIR>\GEANT4\results\kcl_parts\KCl ч 1 л 1085 г часть 2 (88,3 − 41,5 ч).xml"
TEMPLATE = r"<DOSIM>\Спектры\Atom GS2020 PRO MAX\Референсы\KCl ч 1 л 1085 г (88,3 ч, 05-09.10).xml"
J1 = r"C:\g4work\gs2020\kcl_parts\cal_shape_p1.json"
J2 = r"C:\g4work\gs2020\kcl_parts\cal_shape_p2.json"
OUT = r"<WORKDIR>\GEANT4\results\kcl_parts\KCl ч 1 л 1085 г сумма выровненная (88,3 ч).xml"


def read_xml(path):
    """Чтение спектра из XML. Возвращает (tree, es, counts)."""
    tree = ET.parse(path)
    es = tree.getroot().find(".//EnergySpectrum")
    counts = np.array([float(x.text) for x in es.find("Spectrum").findall("DataPoint")])
    return tree, es, counts


def energy(x, r):
    """
    Шкала части (точная копия формулы gs2020_calib.py:113-118).
    x — массив дробных каналов, r — словарь json.load(open(J))["kcl"].
    """
    c = r["file_coeffs"]
    e_file = sum(c[k] * x ** k for k in range(len(c)))
    th = np.asarray(r["theta"], float)
    if r.get("knots"):
        return e_file + np.interp(e_file, r["knots"], th)
    u = (e_file - 1000) / 1000.0
    return e_file + sum(t * u ** j for j, t in enumerate(th))


def main():
    # Чтение входных файлов
    tree1, es1, c1 = read_xml(X1)
    tree2, es2, c2 = read_xml(X2)

    # Живые времена
    lt1 = float(es1.find("LiveTime").text)
    lt2 = float(es2.find("LiveTime").text)

    # Загрузка калибровочных коэффициентов
    with open(J1, "r", encoding="utf-8") as f:
        r1 = json.load(f)["kcl"]
    with open(J2, "r", encoding="utf-8") as f:
        r2 = json.load(f)["kcl"]

    # Размерность спектра
    n = len(c1)
    if len(c2) != n:
        raise SystemExit("ОТКАЗ: разные размеры спектров")

    # Рёбра каналов (общий массив для обеих частей)
    edges = np.arange(n + 1, dtype=float) - 0.5

    # Энергии на рёбрах
    E1e = energy(edges, r1)
    E2e = energy(edges, r2)

    # Проверка монотонности
    if np.any(np.diff(E1e) <= 0):
        raise SystemExit("ОТКАЗ: шкала немонотонна (часть 1)")
    if np.any(np.diff(E2e) <= 0):
        raise SystemExit("ОТКАЗ: шкала немонотонна (часть 2)")

    # Перебиновка части 2 на шкалу части 1
    # Дробное положение рёбер части 1 в индексах рёбер части 2
    j2 = np.interp(E1e, E2e, np.arange(n + 1, dtype=float))

    # Кумулятивная сумма отсчётов части 2
    cum2 = np.concatenate(([0.0], np.cumsum(c2)))

    # Отсчёты части 2 в каналах части 1
    c2r = np.diff(np.interp(j2, np.arange(n + 1, dtype=float), cum2))

    # Сумма
    tot = c1 + c2r

    # Проверки
    # 1. Сохранение отсчётов
    rel_diff = abs(c2r.sum() - c2.sum()) / c2.sum()
    print(f"Сохранение отсчётов: c2r.sum()={c2r.sum():.6f}, c2.sum()={c2.sum():.6f}, отн. разность={rel_diff:.2e}")
    if rel_diff > 1e-3:
        raise SystemExit("ОТКАЗ: нарушение сохранения отсчётов")

    # 2. Сдвиг в каналах на опорных энергиях
    ref_energies = [80, 240, 610, 1460, 2615]
    for E in ref_energies:
        ch1 = np.interp(E, E1e, edges)
        ch2 = np.interp(E, E2e, edges)
        shift = ch2 - ch1
        print(f"E {E} кэВ: кан ч.1 {ch1:.4f}, ч.2 {ch2:.4f}, сдвиг {shift:.4f}")

    # 3. Центроид пика K-40 (1380–1540 кэВ)
    # Каналы части 1 в окне
    e1_channels = energy(np.arange(n), r1)
    mask = (e1_channels >= 1380) & (e1_channels <= 1540)
    if not np.any(mask):
        raise SystemExit("ОТКАЗ: нет каналов в окне 1380–1540 кэВ")

    # Центроид по c1
    centroid_c1 = np.sum(e1_channels[mask] * c1[mask]) / np.sum(c1[mask])

    # Центроид по c2r (перебиновка)
    centroid_c2r = np.sum(e1_channels[mask] * c2r[mask]) / np.sum(c2r[mask])

    # Центроид по c2 без перебиновки (энергии части 1 — та же ось)
    # Для c2 используем те же энергии (e1_channels), так как перебиновка уже привела к общей шкале
    # Но по условию: "по c2 без перебиновки (энергии части 1 — та же ось)"
    # Это значит, что мы берём c2 в тех же каналах, что и c1, и используем e1_channels
    centroid_c2 = np.sum(e1_channels[mask] * c2[mask]) / np.sum(c2[mask])

    print(f"Центроид K-40: c1={centroid_c1:.4f} кэВ, c2r={centroid_c2r:.4f} кэВ, c2={centroid_c2:.4f} кэВ")
    print(f"Центроид до перебиновки (c2): {centroid_c2:.4f} кэВ, после (c2r): {centroid_c2r:.4f} кэВ")

    # 3а. Гейт выравнивания: центроид K-40 части 2 после перебиновки = части 1 с точностью 0,25 кэВ (до — 2,2 кэВ)
    if abs(centroid_c2r - centroid_c1) > 0.25:
        raise SystemExit("ОТКАЗ: центроид K-40 после выравнивания отличается от части 1 на %.2f кэВ" % (centroid_c2r - centroid_c1))

    # 4. Отрицательные значения в c2r
    if np.any(c2r < -1e-9):
        raise SystemExit("ОТКАЗ: отрицательные значения в c2r")

    # Запись результата
    tree_out, es_out, _ = read_xml(TEMPLATE)

    # Обновление DataPoint
    data_points = es_out.find("Spectrum").findall("DataPoint")
    if len(data_points) != len(tot):
        raise SystemExit("ОТКАЗ: несовпадение количества DataPoint в шаблоне")

    for i, dp in enumerate(data_points):
        dp.text = "%.6f" % tot[i]

    # Обновление LiveTime
    es_out.find("LiveTime").text = "%.1f" % (lt1 + lt2)

    # Обновление MeasurementTime
    mt1 = int(float(es1.find("MeasurementTime").text))
    mt2 = int(float(es2.find("MeasurementTime").text))
    es_out.find("MeasurementTime").text = str(mt1 + mt2)

    # Запись файла
    tree_out.write(OUT, encoding="utf-8", xml_declaration=True)

    # Печать итогов
    print(f"Живое время: {lt1 + lt2:.1f} с")
    print(f"Сумма отсчётов: {tot.sum():.6f}")
    print(f"Файл сохранён: {OUT}")


if __name__ == "__main__":
    main()
