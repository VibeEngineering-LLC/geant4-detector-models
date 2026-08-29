# -*- coding: utf-8 -*-
import sys
import os
import csv
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

"""
Сложение независимых мюонных прогонов в один шаблон.

Назначение:
    Объединить несколько коротких прогонов (с разными зёрнами) в один длинный,
    так как длинные прогоны на данной машине часто прерываются.

Правила сложения:
    - counts: суммируются побинно.
    - n_events, n_hits_in_crystal, n_overflow: суммируются.
    - per_muon: пересчитывается как counts_сумма / n_events_сумма.
    - max_edep_keV: берётся максимум из всех прогонов.
    - Геометрические параметры (r_disk_mm, z_disk_mm, disk_area_cm2, pdg_expected_per_s, e_lo_gev, e_hi_gev):
      должны совпадать у всех файлов. При расхождении — ошибка.
    - Зёрна: если counts полностью совпадают у двух файлов, это признак общего зерна — ошибка.

Интерфейс:
    python merge_muon_runs.py <выходной.csv> <вход1.csv> <вход2.csv> [...]
"""

def parse_csv(filepath):
    """
    Парсит CSV файл с метаданными и данными.
    Возвращает словарь метаданных, массив bin_keV, counts, per_muon.
    """
    metadata = {}
    data_rows = []
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Файл не найден: {filepath}")
        sys.exit(1)
    except Exception as e:
        print(f"Ошибка чтения файла {filepath}: {e}")
        sys.exit(1)

    # Разделение на метаданные и данные
    # Ищем пустую строку, после которой идёт заголовок данных
    header_found = False
    data_start_index = -1
    
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not header_found:
            if stripped == "":
                # Следующая строка должна быть заголовком
                if i + 1 < len(lines) and lines[i+1].strip().startswith("bin_keV"):
                    data_start_index = i + 1
                    header_found = True
                    break
            else:
                # Это метаданные
                if ',' in stripped:
                    key, value = stripped.split(',', 1)
                    metadata[key.strip()] = value.strip()
    
    if not header_found or data_start_index == -1:
        print(f"Неверный формат файла {filepath}: не найден заголовок данных.")
        sys.exit(1)

    # Парсинг данных
    reader = csv.reader(lines[data_start_index:])
    try:
        header = next(reader)
        if header != ['bin_keV', 'counts', 'per_muon']:
            print(f"Неверный заголовок в {filepath}: ожидалось ['bin_keV', 'counts', 'per_muon'], получено {header}")
            sys.exit(1)
        
        bin_kev = []
        counts = []
        per_muon = []
        
        for row in reader:
            if not row or all(cell.strip() == '' for cell in row):
                continue # Пропуск пустых строк в данных
            try:
                b = float(row[0])
                c = int(float(row[1])) # counts могут быть записаны как float, но должны быть целыми
                p = float(row[2])
                bin_kev.append(b)
                counts.append(c)
                per_muon.append(p)
            except ValueError:
                print(f"Ошибка парсинга строки данных в {filepath}: {row}")
                sys.exit(1)
                
    except StopIteration:
        pass # Нет данных после заголовка

    if not bin_kev:
        print(f"Нет данных в файле {filepath}")
        sys.exit(1)

    return metadata, np.array(bin_kev), np.array(counts), np.array(per_muon)

def check_geometry_consistency(metadata_list):
    """
    Проверяет, что геометрические параметры совпадают у всех файлов.
    Возвращает первый словарь метаданных как эталон.
    """
    geo_keys = ['r_disk_mm', 'z_disk_mm', 'disk_area_cm2', 'pdg_expected_per_s', 'e_lo_gev', 'e_hi_gev']
    
    if not metadata_list:
        return {}

    reference_meta = metadata_list[0]
    
    for i, meta in enumerate(metadata_list[1:], 1):
        for key in geo_keys:
            val_ref = reference_meta.get(key)
            val_curr = meta.get(key)
            
            if val_ref is None or val_curr is None:
                # Если ключ отсутствует в одном из файлов, это может быть проблемой, 
                # но по условию задачи мы проверяем совпадение. 
                # Если ключа нет, считаем его несовпадающим или игнорируем?
                # Условие: "обязаны СОВПАДАТЬ". Если одного нет, а другое есть - расхождение.
                if val_ref != val_curr:
                    print(f"Расхождение в ключе '{key}': файл 0 имеет '{val_ref}', файл {i} имеет '{val_curr}'")
                    sys.exit(1)
            else:
                # Сравниваем как строки, чтобы избежать проблем с точностью float при чтении из CSV
                if val_ref != val_curr:
                    print(f"Расхождение в ключе '{key}': файл 0 имеет '{val_ref}', файл {i} имеет '{val_curr}'")
                    sys.exit(1)

    return reference_meta

def check_seed_uniqueness(counts_list, filenames):
    """
    Проверяет, что counts не совпадают полностью у разных файлов.
    Если совпадают - это признак общего зерна.
    """
    for i in range(len(counts_list)):
        for j in range(i + 1, len(counts_list)):
            if np.array_equal(counts_list[i], counts_list[j]):
                print(f"Обнаружено полное совпадение counts между файлами: {filenames[i]} и {filenames[j]}. Вероятно, использовано одно зерно.")
                sys.exit(1)

def check_shield_consistency(filenames):
    """Мюонный CSV НЕ содержит признака домика — в отличие от run_field, где в
    метаданные пишется shield=0/1 именно ради того, чтобы пару «с защитой / без»
    нельзя было перепутать постфактум. Здесь единственная доступная защита —
    имя файла: сложить прогон с защитой и без значит получить бессмыслицу,
    внешне неотличимую от правильного результата."""
    kinds = set()
    for f in filenames:
        base = os.path.basename(f)
        kinds.add("on" if "_on" in base else ("off" if "_off" in base else "?"))
    if len(kinds) > 1:
        print("Входные файлы разного типа по домику (%s):" % ", ".join(sorted(kinds)))
        for f in filenames:
            print("   ", os.path.basename(f))
        sys.exit(1)
    if "?" in kinds:
        print("ПРЕДУПРЕЖДЕНИЕ: по именам не видно, с домиком прогоны или без.")


def main():
    if len(sys.argv) < 4:
        print("Использование: python merge_muon_runs.py <выходной.csv> <вход1.csv> <вход2.csv> [...]")
        sys.exit(2)

    output_file = sys.argv[1]
    input_files = sys.argv[2:]

    metadata_list = []
    bin_kev_list = []
    counts_list = []
    per_muon_list = []
    
    # Парсинг всех входных файлов
    for filepath in input_files:
        meta, b, c, p = parse_csv(filepath)
        metadata_list.append(meta)
        bin_kev_list.append(b)
        counts_list.append(c)
        per_muon_list.append(p)

    # Проверка геометрии
    check_shield_consistency(input_files)
    reference_meta = check_geometry_consistency(metadata_list)

    # Проверка уникальности зёрен (по совпадению counts)
    check_seed_uniqueness(counts_list, input_files)

    # Проверка согласованности bin_kev (должны быть одинаковыми у всех файлов для корректного сложения побинно)
    reference_bin = bin_kev_list[0]
    for i, b in enumerate(bin_kev_list[1:], 1):
        if not np.allclose(reference_bin, b):
            print(f"Расхождение в bin_keV между файлом 0 и файлом {i}. Сложение невозможно.")
            sys.exit(1)

    # Суммирование
    total_counts = np.sum(counts_list, axis=0)
    
    # Суммирование метрик
    total_n_events = sum(int(m.get('n_events', 0)) for m in metadata_list)
    total_n_hits = sum(int(m.get('n_hits_in_crystal', 0)) for m in metadata_list)
    total_n_overflow = sum(int(m.get('n_overflow', 0)) for m in metadata_list)
    
    # Максимум max_edep_keV
    max_edep_values = [float(m.get('max_edep_keV', 0)) for m in metadata_list]
    final_max_edep = max(max_edep_values)

    # Пересчёт per_muon
    if total_n_events > 0:
        final_per_muon = total_counts / total_n_events
    else:
        final_per_muon = np.zeros_like(total_counts, dtype=float)

    # Подготовка метаданных для вывода
    output_metadata = reference_meta.copy()
    output_metadata['n_events'] = str(total_n_events)
    output_metadata['n_hits_in_crystal'] = str(total_n_hits)
    output_metadata['n_overflow'] = str(total_n_overflow)
    output_metadata['max_edep_keV'] = str(final_max_edep)
    output_metadata['n_merged'] = str(len(input_files))
    
    # Формирование списка исходных файлов (базовые имена)
    base_names = [os.path.basename(f) for f in input_files]
    output_metadata['merged_from'] = ';'.join(base_names)

    # Вывод в stdout
    print("Результаты по входным файлам:")
    for i, filepath in enumerate(input_files):
        n_ev = metadata_list[i].get('n_events', 'N/A')
        n_hi = metadata_list[i].get('n_hits_in_crystal', 'N/A')
        print(f"  {os.path.basename(filepath)}: n_events={n_ev}, n_hits_in_crystal={n_hi}")
    
    print(f"\nИтог:")
    print(f"  Суммарные n_events: {total_n_events}")
    print(f"  Суммарные n_hits_in_crystal: {total_n_hits}")
    print(f"  Выходной файл: {output_file}")

    # Запись выходного файла
    try:
        with open(output_file, 'w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f)
            
            # Метаданные
            for key, value in output_metadata.items():
                writer.writerow([key, value])
            
            # Пустая строка
            writer.writerow([])
            
            # Заголовок данных
            writer.writerow(['bin_keV', 'counts', 'per_muon'])
            
            # Данные
            for b, c, p in zip(reference_bin, total_counts, final_per_muon):
                # bin_kev может быть float, counts int, per_muon float
                writer.writerow([f"{b:.6f}", str(int(c)), f"{p:.10e}"])
                
    except Exception as e:
        print(f"Ошибка записи файла {output_file}: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
