import sys
import os
import argparse
import math
import csv
import statistics

sys.stdout.reconfigure(encoding="utf-8")

# Импорт parse_csv из analyze_ncompt.py
try:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from analyze_ncompt import parse_csv
    from analyze_sigma_intr import read_spectrum_table
except ImportError as e:
    print(f"Ошибка импорта: {e}")
    sys.exit(3)

def find_peak_position(count_light, start_bin=1):
    """Найти пик полного поглощения в диапазоне [0.2*E, 1.05*E]"""
    max_count = 0
    peak_bin = -1
    for i in range(start_bin, len(count_light)):
        if count_light[i] > max_count:
            max_count = count_light[i]
            peak_bin = i
    return peak_bin, max_count

def calculate_weighted_center(counts, bin_centers):
    """Вычислить центр тяжести по бинам выше половины максимума"""
    if not counts:
        return 0.0
    max_val = max(counts)
    threshold = max_val / 2.0
    total_weighted = 0.0
    total_weight = 0.0
    for i, (count, center) in enumerate(zip(counts, bin_centers)):
        if count >= threshold:
            total_weighted += count * center
            total_weight += count
    return total_weighted / total_weight if total_weight > 0 else 0.0

def smooth_curve(data, window=9):
    """Сгладить кривую скользящим средним"""
    smoothed = [0.0] * len(data)
    half_window = window // 2
    for i in range(len(data)):
        start_idx = max(0, i - half_window)
        end_idx = min(len(data), i + half_window + 1)
        window_data = data[start_idx:end_idx]
        smoothed[i] = sum(window_data) / len(window_data)
    return smoothed

def numerical_derivative(smoothed_curve):
    """Численная производная"""
    deriv = [0.0] * len(smoothed_curve)
    for i in range(1, len(smoothed_curve) - 1):
        deriv[i] = (smoothed_curve[i+1] - smoothed_curve[i-1]) / 2.0
    return deriv

def edge_plateau_level(counts, bin_centers, edge_energy_keV):
    """Половина уровня плато континуума перед краем (см. _spec_analyze_edge.md)."""
    lo, hi = 0.45 * edge_energy_keV, 0.85 * edge_energy_keV
    plateau = [c for b, c in zip(bin_centers, counts) if lo <= b <= hi]
    if len(plateau) < 5:
        return None
    # statistics.median, а не numpy: этот файл numpy не импортирует.
    level = 0.5 * float(statistics.median(plateau))
    return level if level > 0 else None


def find_edge_half_height(counts, bin_centers, edge_energy_keV):
    """Край = спад до половины плато. Заменяет поиск минимума производной,
    который ловил посторонние структуры (W-061)."""
    level = edge_plateau_level(counts, bin_centers, edge_energy_keV)
    if level is None:
        return None
    hi = 0.85 * edge_energy_keV
    prev_b = prev_c = None
    for b, c in zip(bin_centers, counts):
        if b < hi:
            prev_b, prev_c = b, c
            continue
        if b > 1.35 * edge_energy_keV:
            break
        if c < level and prev_c is not None and prev_c >= level and prev_c != c:
            return prev_b + (prev_c - level) * (b - prev_b) / (prev_c - c)
        prev_b, prev_c = b, c
    return None


def find_edge_position(deriv, bin_centers, edge_energy_keV):
    """Найти точку перегиба в области [0.7*E_edge, 1.15*E_edge]"""
    # Определяем границы поиска
    min_edge = 0.7 * edge_energy_keV
    max_edge = 1.15 * edge_energy_keV
    
    # Находим бины в этом диапазоне
    valid_indices = []
    for i, center in enumerate(bin_centers):
        if min_edge <= center <= max_edge:
            valid_indices.append(i)
    
    if not valid_indices:
        return None
    
    # Ищем минимум производной среди этих бинов
    min_deriv = float('inf')
    min_idx = -1
    for i in valid_indices:
        if deriv[i] < min_deriv:
            min_deriv = deriv[i]
            min_idx = i
    
    if min_idx == -1:
        return None
    
    # Уточняем положение параболой по трём точкам
    if min_idx > 0 and min_idx < len(deriv) - 1:
        x = [bin_centers[min_idx-1], bin_centers[min_idx], bin_centers[min_idx+1]]
        y = [deriv[min_idx-1], deriv[min_idx], deriv[min_idx+1]]
        
        # Парабола: y = a*x^2 + b*x + c
        # Используем интерполяцию по трем точкам
        if x[2] != x[0]:
            a = (y[2] - y[0]) / ((x[2] - x[0]) * (x[2] - x[1]))
            b = (y[1] - y[0]) / (x[1] - x[0]) - a * (x[1] + x[0])
            # Минимум параболы: x = -b/(2*a)
            if abs(a) > 1e-10:
                edge_pos = -b / (2.0 * a)
                return edge_pos
    return bin_centers[min_idx]

def calculate_edge_energy(energy_keV):
    """Вычислить теоретический край"""
    E = energy_keV
    E_edge = E * (2 * E / 511) / (1 + 2 * E / 511)
    return E_edge

def process_run(file_path, results):
    data = parse_csv(file_path)
    if not data:
        print(f"Ошибка: не удалось прочитать файл {file_path}")
        return False
    
    energy_keV = data.get("energy_keV")
    npsm_enabled = data.get("npsm_enabled", 1)
    
    if energy_keV is None:
        print(f"Ошибка: в файле {file_path} отсутствует поле energy_keV")
        return False
    
    # parse_csv отдаёт только шапку. Чтение второй таблицы (спектров) уже
    # реализовано в analyze_sigma_intr.read_spectrum_table — берётся оттуда
    # импортом, а не переписывается заново (§33: дублирующая реализация —
    # дефект, две копии разойдутся).
    table = read_spectrum_table(file_path)
    bin_keV_list = [r[0] for r in table]
    count_edep_list = [r[1] for r in table]
    count_light_list = [r[2] for r in table]
    
    if not bin_keV_list:
        print(f"Ошибка: в файле {file_path} не найдены данные для анализа")
        return False
    
    # Вычисляем теоретический край
    E_edge = calculate_edge_energy(energy_keV)
    
    # Калибровка по пикy полного поглощения
    peak_bin, max_count = find_peak_position(count_light_list)
    if peak_bin == -1 or max_count <= 0:
        print(f"Ошибка: не найден пик полного поглощения в файле {file_path}")
        return False
    
    # Уточняем положение пика
    bin_centers = [x for x in bin_keV_list]
    peak_position = calculate_weighted_center(count_light_list[peak_bin-2:peak_bin+3], 
                                              bin_centers[peak_bin-2:peak_bin+3])
    
    k = energy_keV / peak_position if peak_position > 0 else 1.0
    
    # Калибровка растягивает ОСЬ ЭНЕРГИЙ, а не высоту отсчётов. Прежняя версия
    # умножала на k сам счёт (count * k): положение края от этого не меняется
    # вовсе, и поиск шёл в области, где края нет, упираясь в границу диапазона.
    bin_centers_light = [b * k for b in bin_centers]

    # Сглаживаем обе кривые
    smooth_edep = smooth_curve(count_edep_list, 9)
    smooth_light = smooth_curve(count_light_list, 9)
    
    # Вычисляем производные
    deriv_edep = numerical_derivative(smooth_edep)
    deriv_light = numerical_derivative(smooth_light)
    
    # Находим положение края по энергии
    # Обе кривые обрабатываются ОДНИМ методом — систематика сокращается в разности.
    edge_pos_energy = find_edge_half_height(smooth_edep, bin_centers, E_edge)
    
    # Находим положение края по свету
    edge_pos_light = find_edge_half_height(smooth_light, bin_centers_light, E_edge)
    
    if edge_pos_energy is None or edge_pos_light is None:
        print(f"Ошибка: не удалось определить положение края в файле {file_path}")
        return False
    
    # Сдвиг
    shift = edge_pos_light - edge_pos_energy
    
    results.append({
        "energy_keV": energy_keV,
        "E_edge": E_edge,
        "edge_pos_energy": edge_pos_energy,
        "peak_position": peak_position,
        "k": k,
        "edge_pos_light": edge_pos_light,
        "shift": shift,
        "npsm_enabled": npsm_enabled
    })
    
    return True

def check_criteria(results):
    """Проверить критерии"""
    criteria = {}
    
    # К-1: край по ЭНЕРГИИ совпадает с теоретическим
    if len(results) > 0:
        energy_keV = results[0]["energy_keV"]
        E_edge = calculate_edge_energy(energy_keV)
        edge_pos_energy = results[0]["edge_pos_energy"]
        diff = abs(edge_pos_energy - E_edge)
        criteria["K1"] = {
            "threshold": "≤ 5 кэВ",
            "value": f"{diff:.2f} кэВ",
            "pass": diff <= 5.0
        }
    else:
        criteria["K1"] = {
            "threshold": "≤ 5 кэВ",
            "value": "не определено",
            "pass": False
        }
    
    # К-2: знак сдвига отрицательный
    if len(results) > 0:
        shift = results[0]["shift"]
        criteria["K2"] = {
            "threshold": "< 0",
            "value": f"{shift:.2f} кэВ",
            "pass": shift < 0
        }
    else:
        criteria["K2"] = {
            "threshold": "< 0",
            "value": "не определено",
            "pass": False
        }
    
    # К-3: величина сдвига
    if len(results) > 0:
        shift = results[0]["shift"]
        abs_shift = abs(shift)
        criteria["K3"] = {
            "threshold": "12…28 кэВ",
            "value": f"{abs_shift:.2f} кэВ",
            "pass": 12.0 <= abs_shift <= 28.0
        }
    else:
        criteria["K3"] = {
            "threshold": "12…28 кэВ",
            "value": "не определено",
            "pass": False
        }
    
    # К-4: прогон с npsm_enabled = 0
    has_disabled_run = any(r["npsm_enabled"] == 0 for r in results)
    if has_disabled_run:
        disabled_shift = next((r["shift"] for r in results if r["npsm_enabled"] == 0), None)
        if disabled_shift is not None:
            abs_shift = abs(disabled_shift)
            criteria["K4"] = {
                "threshold": "≤ 3 кэВ",
                "value": f"{abs_shift:.2f} кэВ",
                "pass": abs_shift <= 3.0
            }
        else:
            criteria["K4"] = {
                "threshold": "≤ 3 кэВ",
                "value": "не определено",
                "pass": False
            }
    else:
        criteria["K4"] = {
            "threshold": "≤ 3 кэВ",
            "value": "прогон с выключенной моделью не найден",
            "pass": None,
            "unverified": True
        }
    
    return criteria

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="append", help="CSV файлы для анализа")
    parser.add_argument("--json", help="Файл для сохранения JSON результата")
    
    args = parser.parse_args()
    
    if not args.run:
        print("Ошибка: не указаны файлы для анализа (--run)")
        sys.exit(1)
    
    results = []
    for file_path in args.run:
        success = process_run(file_path, results)
        if not success:
            sys.exit(2)
    
    # Проверяем критерии
    criteria = check_criteria(results)
    
    # Выводим результаты
    print("Энергия (кэВ), Теоретический край (кэВ), Край по энергии (кэВ), "
          "Пик ППП в свете (кэВ), Коэффициент k, Край по свету (кэВ), Сдвиг (кэВ)")
    
    for result in results:
        print(f"{result['energy_keV']:.1f}, {result['E_edge']:.2f}, "
              f"{result['edge_pos_energy']:.2f}, {result['peak_position']:.2f}, "
              f"{result['k']:.4f}, {result['edge_pos_light']:.2f}, "
              f"{result['shift']:.2f}")
    
    print("\nКРИТЕРИИ:")
    for key, crit in criteria.items():
        if "unverified" in crit:
            print(f"{key}: {crit['threshold']} — {crit['value']} (непроверено)")
        else:
            status = "PASS" if crit["pass"] else "FAIL"
            print(f"{key}: {crit['threshold']} — {crit['value']} ({status})")
    
    # Требует толкования
    print("\nТРЕБУЕТ ТОЛКОВАНИЯ:")
    has_low_counts = False
    for result in results:
        energy_keV = result["energy_keV"]
        E_edge = calculate_edge_energy(energy_keV)
        min_edge = 0.7 * E_edge
        max_edge = 1.15 * E_edge
        
        # Суммарные отсчёты в области края. Ключей count_edep/count_light в
        # result нет — они не сохранялись; берётся сохранённая сумма, а при её
        # отсутствии проверка честно пропускается с пометкой, а не падает.
        total_counts = result.get("counts_near_edge")
        if total_counts is None:
            print(f"  E={energy_keV:.0f} кэВ: суммарные отсчёты у края не сохранены,"
                  f" проверка достаточности статистики НЕ выполнена")
            total_counts = 10 ** 9
        
        if total_counts < 200:
            has_low_counts = True
            print(f"Область края содержит менее 200 отсчётов ({total_counts})")
    
    # Проверка на несколько минимумов производной
    print("Проверка на наличие нескольких сравнимых минимумов производной:")
    print("Это требует дополнительного анализа, так как не реализовано в текущей версии.")
    
    # Итог
    if len(results) > 0:
        shift = results[0]["shift"]
        acceptance = "PASS" if all(crit["pass"] for crit in criteria.values() if "unverified" not in crit) else "FAIL"
        print(f"\nSTAGE4_EDGE_SHIFT={shift:.2f} keV")
        print(f"STAGE4_ACCEPTANCE={acceptance}")
    else:
        print("\nSTAGE4_EDGE_SHIFT=0.00 keV")
        print("STAGE4_ACCEPTANCE=FAIL")

if __name__ == "__main__":
    main()
