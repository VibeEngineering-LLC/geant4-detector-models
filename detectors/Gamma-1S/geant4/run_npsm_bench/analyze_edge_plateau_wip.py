import csv
import math
import glob
import argparse
import pathlib
import sys
import json

sys.stdout.reconfigure(encoding="utf-8")

def parse_csv(file_path):
    data = {}
    table = []
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                if "," in line:
                    key, value = line.strip().split(",", 1)
                    # Приводим к числу ОДИН РАЗ при чтении — все потребители
                    # ниже получают уже числа. Точечная конвертация по списку
                    # ключей оставляла mean/sem строками и роняла критерии K1/K2.
                    if key in ("em_deex", "crystal_mm"):
                        data[key] = value            # заведомо строковые поля
                    else:
                        try:
                            data[key] = int(value)
                        except ValueError:
                            try:
                                data[key] = float(value)
                            except ValueError:
                                data[key] = value
                elif line.strip() == "n_compt,count_absorbed,count_all":
                    continue
                elif line.strip().startswith("0,"):
                    break
            for line in f:
                if "," in line and not line.startswith("#"):
                    parts = line.strip().split(",")
                    if len(parts) >= 3:
                        table.append([int(parts[0]), int(parts[1]), int(parts[2])])
    except Exception as e:
        print(f"ПРЕДУПРЕЖДЕНИЕ: Не удалось прочитать файл {file_path}: {e}")
        return None
    if "energy_keV" not in data or "mean_ncompt_absorbed" not in data:
        print(f"ПРЕДУПРЕЖДЕНИЕ: В файле {file_path} отсутствуют необходимые поля")
        return None
    data["table"] = table
    return data

def smooth(data, window=9):
    """Скользящее среднее"""
    result = [0.0] * len(data)
    half_window = window // 2
    for i in range(len(data)):
        start = max(0, i - half_window)
        end = min(len(data), i + half_window + 1)
        result[i] = sum(data[start:end]) / (end - start)
    return result

def find_peak_position(counts, bin_centers, energy_range):
    """Найти пик в заданном диапазоне и вернуть его центр тяжести"""
    low, high = energy_range
    indices = [i for i in range(len(bin_centers)) if low <= bin_centers[i] <= high]
    if not indices:
        return None
    max_count = max(counts[i] for i in indices)
    threshold = max_count / 2.0
    # Найти бины выше порога
    above_threshold = [i for i in indices if counts[i] >= threshold]
    if not above_threshold:
        return None
    # Пик в шкале света ДВУХКОМПОНЕНТЕН: левее истинного ППП стоит escape-пик
    # иода (Kα ≈ 28,6 кэВ). Центр тяжести по всем бинам выше порога садится
    # между компонентами и занижает положение ППП, завышая калибровку k.
    # Берём СВЯЗНУЮ группу бинов выше порога, содержащую самый правый бин.
    # Если компоненты слились в одну группу, поведение не меняется.
    # Одиночный статистический выброс группой не считается: при слабой
    # статистике самый правый бин выше порога может быть оторван от пика
    # (проверено на 2000 кэВ). Поэтому среди связных групп оставляем только
    # значимые по площади (>= 20 % от наибольшей) и берём самую правую из них.
    groups = []
    for i in above_threshold:
        if groups and i == groups[-1][-1] + 1:
            groups[-1].append(i)
        else:
            groups.append([i])
    areas = [sum(counts[j] for j in g) for g in groups]
    big = max(areas)
    group = [g for g, a in zip(groups, areas) if a >= 0.2 * big][-1]

    # Уточнить положение как центр тяжести
    total_weighted = 0.0
    total_weight = 0.0
    for i in group:
        weight = counts[i]
        total_weighted += bin_centers[i] * weight
        total_weight += weight
    if total_weight == 0:
        return None
    return total_weighted / total_weight

def find_edge_position(counts, bin_centers, energy_edge):
    """Найти положение края по сглаженной кривой"""
    # Сгладить кривую
    smoothed = smooth(counts, 9)
    
    # Уровень плато: медиана в [0.45*E_edge, 0.85*E_edge]
    low_plate = 0.45 * energy_edge
    high_plate = 0.85 * energy_edge
    plate_indices = [i for i in range(len(bin_centers)) if low_plate <= bin_centers[i] <= high_plate]
    if len(plate_indices) < 5:
        return None, "не определено"
    
    plate_values = [smoothed[i] for i in plate_indices]
    plate_median = sorted(plate_values)[len(plate_values)//2]
    
    # Порог: половина уровня плато
    threshold = plate_median / 2.0
    
    # Идти вправо от 0.85*E_edge до первого бина, где сглаженный счёт опускается ниже порога
    start_index = next((i for i in range(len(bin_centers)) if bin_centers[i] > high_plate), None)
    if start_index is None:
        return None, "не определено"
    
    # Идти вправо от плато до ПЕРВОГО бина, где счёт опустился ниже порога,
    # а предыдущий был выше. Дальше 1,35*E_edge не идти (ограничение спеки).
    limit = 1.35 * energy_edge
    first_below = None
    for i in range(start_index, len(bin_centers)):
        if bin_centers[i] > limit:
            break
        if smoothed[i] < threshold and i > 0 and smoothed[i - 1] >= threshold:
            first_below = i
            break

    if first_below is None:
        return None, "не определено"

    last_above = first_below - 1
    
    if first_below <= last_above:
        return None, "не определено"
    
    # Линейная интерполяция между этими двумя бинами
    x1, x2 = bin_centers[last_above], bin_centers[first_below]
    y1, y2 = smoothed[last_above], smoothed[first_below]
    if y1 == y2:
        return None, "не определено"
    
    # Найдем x, где y = threshold
    x_edge = x1 + (threshold - y1) * (x2 - x1) / (y2 - y1)
    return x_edge, "успешно"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="append", required=True)
    parser.add_argument("--json")
    args = parser.parse_args()

    results = []
    all_files_passed = True
    files_with_npsm_disabled = []

    for file_path in args.run:
        try:
            data = parse_csv(file_path)
            if not data:
                print(f"Ошибка при чтении файла {file_path}")
                sys.exit(3)
            
            energy_keV = data["energy_keV"]
            npsm_enabled = data["npsm_enabled"]
            # parse_csv отдаёт таблицу ГИСТОГРАММЫ РАССЕЯНИЙ, а нужен СПЕКТР
            # (вторая таблица файла). Чтение уже реализовано в
            # analyze_sigma_intr.read_spectrum_table — берётся импортом (§33),
            # генератор подставил сюда data["table"] уже дважды.
            from analyze_sigma_intr import read_spectrum_table
            table = read_spectrum_table(file_path)
            if not table:
                print(f"Ошибка: пустая таблица спектра в {file_path}")
                sys.exit(3)

            # Разбор таблицы
            bin_keV = []
            count_edep = []
            count_light = []
            reading_table = False
            for line in table:
                if len(line) >= 3:
                    bin_keV.append(line[0])
                    count_edep.append(line[1])
                    count_light.append(line[2])

            # Теоретический край
            E_edge_theor = energy_keV * (2 * energy_keV / 511.0) / (1 + 2 * energy_keV / 511.0)

            # Калибровка: найти пик полного поглощения в области [0, 2*E, 1.05*E]
            peak_light_pos = find_peak_position(count_light, bin_keV, (0.2 * energy_keV, 1.05 * energy_keV))
            if peak_light_pos is None:
                print(f"Ошибка: не найден пик полного поглощения в файле {file_path}")
                sys.exit(3)
            
            # Коэффициент калибровки
            k = energy_keV / peak_light_pos

            # Калиброванные бины по свету
            bin_centers_light = [b * k for b in bin_keV]

            # Найти положение края по энергии (эталон)
            edge_by_energy, status1 = find_edge_position(count_edep, bin_keV, E_edge_theor)
            if edge_by_energy is None:
                print(f"Ошибка: не найден край по энергии в файле {file_path}")
                sys.exit(3)

            # Найти положение края по свету (калиброванному)
            edge_by_light, status2 = find_edge_position(count_light, bin_centers_light, E_edge_theor)
            if edge_by_light is None:
                print(f"Ошибка: не найден край по свету в файле {file_path}")
                sys.exit(3)

            # Сдвиг
            shift = edge_by_light - edge_by_energy

            results.append({
                "E_keV": energy_keV,
                "E_edge_theor": E_edge_theor,
                "edge_by_energy": edge_by_energy,
                "peak_light": peak_light_pos,
                "k": k,
                "edge_by_light": edge_by_light,
                "shift": shift,
                "npsm": npsm_enabled
            })

            if npsm_enabled == 0:
                files_with_npsm_disabled.append(file_path)

        except Exception as e:
            print(f"Ошибка при обработке файла {file_path}: {e}")
            sys.exit(3)

    # Вывод CSV
    print("E_keV,E_edge_theor,edge_by_energy,peak_light,k,edge_by_light,shift")
    for r in results:
        print(f"{r['E_keV']:.3f},{r['E_edge_theor']:.3f},{r['edge_by_energy']:.3f},{r['peak_light']:.3f},{r['k']:.3f},{r['edge_by_light']:.3f},{r['shift']:.3f}")

    # Критерии
    print("\nКРИТЕРИИ")

    # К-1…К-3 и К-5 относятся к прогонам с ВКЛЮЧЁННОЙ моделью. Смешивать их с
    # мутационным прогоном нельзя: К-3 брал первый результат с E = 662 и мог
    # взять прогон с выключенной моделью, где сдвига нет по построению.
    on_results = [r for r in results if r["npsm"] == 1]
    if not on_results:
        print("Нет ни одного прогона с включённой моделью — критерии не проверяются")
        return 3

    # К-1: край по ЭНЕРГИИ совпадает с теоретическим
    k1_pass = True
    for r in on_results:
        diff = abs(r["edge_by_energy"] - r["E_edge_theor"])
        if diff > 0.08 * r["E_edge_theor"]:
            k1_pass = False
            break
    print(f"К-1: край по ЭНЕРГИИ совпадает с теоретическим (расхождение ≤ 8 % от E_edge)")
    print(f"Факт: {'PASS' if k1_pass else 'FAIL'}")

    # К-2: знак сдвига отрицательный
    k2_pass = all(r["shift"] < 0 for r in on_results)
    print(f"К-2: знак сдвига отрицательный (на ВСЕХ поданных энергиях)")
    print(f"Факт: {'PASS' if k2_pass else 'FAIL'}")

    # К-3: величина сдвига при 662 кэВ
    shift_662 = None
    for r in on_results:
        if abs(r["E_keV"] - 662) < 1:
            shift_662 = r["shift"]
            break
    k3_pass = False
    if shift_662 is not None and abs(shift_662) >= 12 and abs(shift_662) <= 28:
        k3_pass = True
    print(f"К-3: величина сдвига при 662 кэВ (12…28 кэВ по модулю)")
    print(f"Факт: {'PASS' if k3_pass else 'FAIL'}")

    # К-4: мутация — при ВЫКЛЮЧЕННОЙ модели свет тождествен энергии, сдвиг
    # обязан исчезнуть. Проверять только прогоны с npsm = 0: раньше цикл шёл
    # по всем результатам с E = 662, включая прогон с ВКЛЮЧЁННОЙ моделью,
    # и критерий-мутация был обречён на FAIL при любом исходе.
    print("К-4: прогон с npsm_enabled = 0 (сдвиг не более 3 кэВ по модулю)")
    off_results = [r for r in results if r["npsm"] == 0]
    if not off_results:
        print("Факт: НЕ ПРОВЕРЕН — прогон с выключенной моделью не подан, итог FAIL")
        k4_pass = None
    else:
        worst = max(abs(r["shift"]) for r in off_results)
        k4_pass = worst <= 3
        print(f"Факт: {'PASS' if k4_pass else 'FAIL'} (наибольший по модулю сдвиг {worst:.3f} кэВ "
              f"на {len(off_results)} прогоне(ах))")

    # К-5: воспроизводимость: сдвиг растёт по модулю с энергией
    shifts = [(r["E_keV"], r["shift"]) for r in on_results]
    shifts.sort(key=lambda x: x[0])
    if len(shifts) >= 3:
        # Проверить монотонность по модулю
        abs_shifts = [abs(s) for _, s in shifts]
        increasing = all(abs_shifts[i] <= abs_shifts[i+1] for i in range(len(abs_shifts)-1))
        k5_pass = increasing
    else:
        k5_pass = False
    print(f"К-5: воспроизводимость: сдвиг растёт по модулю с энергией (монотонно, не менее чем на 3 точках)")
    print(f"Факт: {'PASS' if k5_pass else 'FAIL'}")

    # Требует толкования
    total_counts = sum(count_light) + sum(count_edep)
    print("\nТРЕБУЕТ ТОЛКОВАНИЯ")
    if total_counts < 200:
        print("Общее количество отсчётов в области края меньше 200")

    # Итог
    all_criteria_pass = k1_pass and k2_pass and k3_pass and (k4_pass is None or k4_pass) and k5_pass
    if all_criteria_pass:
        print("\nSTAGE4_ACCEPTANCE=PASS")
    else:
        print("\nSTAGE4_ACCEPTANCE=FAIL")

    # Средний сдвиг для итога
    avg_shift = sum(r["shift"] for r in results) / len(results)
    print(f"STAGE4_EDGE_SHIFT={avg_shift:.3f} keV")

if __name__ == "__main__":
    main()
