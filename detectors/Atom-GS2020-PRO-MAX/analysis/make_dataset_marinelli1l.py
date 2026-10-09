#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Скрипт для создания набора данных детектора GS2020 PRO MAX + Маринелли 1 л
из сетки откликов.
"""

import sys
import os
import json
import glob
import hashlib
import shutil
import csv
import re
import datetime
import numpy as np

# Первая строка после импортов
sys.stdout.reconfigure(encoding="utf-8")

# Импорт функции разбора сетки
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_eff_becqmoni import parse_grid_file
from fep_def import fep_net   # #GS-66

# Константы
GRID_DIR = os.environ.get("GS_GRID_DIR") or r"C:\g4work\gs2020\run_marinelli\out_v5_oisn10"   # v1.1: сборная папка с узлами #GS-65
OUT_DIR = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "dataset"))   # в репозитории: analysis/ → dataset/
GDML_SRC = r"C:\g4work\gs2020\gs20_matrix\mx_oisn10.gdml"
TAG = "marinelli1l_oisn10"
_canon = lambda p: re.sub(r'(\d\.\d*?)0+(?=")', r'\1', open(p, encoding="utf-8").read())   # «42.50» == «42.5»
GDML_EQUIV = [p for p in (os.environ.get("GS_GDML_EQUIV") or "").split(";") if p]   # #GS-64: тот же GDML, иная запись чисел
for _p in GDML_EQUIV:
    if _canon(_p) != _canon(GDML_SRC):
        raise SystemExit("ОТКАЗ: GDML " + _p + " не эквивалентен " + GDML_SRC)


def sha256_file(filepath):
    """Вычисляет SHA256 хеш файла."""
    sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def parse_header_and_light(filepath):
    """
    Читает шапку файла (key,value до строки bin_keV) и столбец count_light.
    Возвращает словарь с ключами: seed, em_cut_mm, em_deex, em_option, particle, npsm_enabled,
    и список count_light (список целых чисел).
    """
    header_keys = {}
    count_light = []
    bins = []
    
    with open(filepath, "r", encoding="utf-8") as f:
        lines = f.readlines()
    
    # Пропускаем первую строку (комментарий)
    # Ищем строку "bin_keV"
    bin_kev_idx = -1
    for i, line in enumerate(lines):
        if line.strip().startswith("bin_keV"):
            bin_kev_idx = i
            break
    
    if bin_kev_idx == -1:
        raise ValueError(f"Не найдена строка 'bin_keV' в файле {filepath}")
    
    # Читаем шапку (строки до bin_kev_idx)
    for i in range(bin_kev_idx):
        line = lines[i].strip()
        if not line or line.startswith("#"):
            continue
        if "," in line:
            parts = line.split(",")
            if len(parts) >= 2:
                key = parts[0].strip()
                value = parts[1].strip()
                header_keys[key] = value
    
    # Читаем данные (строки после bin_kev_idx)
    # Формат: bin_keV, edep, count_light (предполагаем, что count_light — третий столбец)
    for i in range(bin_kev_idx + 1, len(lines)):
        line = lines[i].strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(",")
        if len(parts) >= 3:
            try:
                light_val = int(parts[2])
                bins.append((float(parts[0]), float(parts[1])))   # (центр бина, отсчёты edep) — parse_grid_file бинов не отдаёт
                count_light.append(light_val)
            except ValueError:
                # Если не удалось преобразовать, пропускаем
                pass

    return header_keys, count_light, bins


def main():
    # Создаем выходную директорию
    os.makedirs(OUT_DIR, exist_ok=True)
    
    # Находим все файлы сетки
    grid_files = sorted(glob.glob(os.path.join(GRID_DIR, "grid_mar_E*.csv")))
    
    if not grid_files:
        raise SystemExit(f"ОТКАЗ: не найдено файлов сетки в {GRID_DIR}")
    
    print(f"Найдено файлов сетки: {len(grid_files)}")
    
    # Хранилище данных для каждого узла
    nodes_data = []
    
    # Общие параметры, которые должны совпадать
    sample_matrix = None
    sample_rho_g_cm3 = None
    physics_params = None
    bin_centers = None
    
    for filepath in grid_files:
        # Читаем первую строку (комментарий с GDML)
        with open(filepath, "r", encoding="utf-8") as f:
            first_line = f.readline().strip()
        
        # Извлекаем путь GDML из комментария
        gdml_path_in_file = None
        if "GDML " in first_line:
            gdml_path_in_file = first_line.split("GDML ")[1].strip()
        
        # Проверяем, что GDML совпадает с GDML_SRC
        if gdml_path_in_file is not None and os.path.normcase(os.path.normpath(gdml_path_in_file)) not in [os.path.normcase(os.path.normpath(x)) for x in [GDML_SRC] + GDML_EQUIV]:
            raise SystemExit(f"ОТКАЗ: GDML в файле {filepath} ({gdml_path_in_file}) отличается от GDML_SRC ({GDML_SRC})")
        
        # Разбираем файл сетки
        grid_data = parse_grid_file(filepath)
        
        if grid_data is None:
            print(f"Предупреждение: пропускаю файл {filepath} (parse_grid_file вернул None)")
            continue
        
        E = grid_data["E"]  # энергия узла в кэВ
        N = grid_data["N"]  # число первичных фотонов
        k_fep_ref = grid_data["counts"]  # отсчёты пика по определению compare_eff_becqmoni (бин с E + сосед при |E−round(E)|<0,02) — для сверки
        sample_matrix_node = grid_data["sample_matrix"]
        sample_rho_node = grid_data["sample_rho_g_cm3"]
        
        # Читаем шапку и count_light
        header_keys, count_light, counts = parse_header_and_light(filepath)   # counts — список (центр, edep); grid_data["counts"] — число отсчётов пика
        
        # Проверяем, что sample_matrix и sample_rho_g_cm3 совпадают
        if sample_matrix is None:
            sample_matrix = sample_matrix_node
            sample_rho_g_cm3 = sample_rho_node
        else:
            if sample_matrix != sample_matrix_node:
                raise SystemExit(f"ОТКАЗ: sample_matrix отличается в файле {filepath}")
            if sample_rho_g_cm3 != sample_rho_node:
                raise SystemExit(f"ОТКАЗ: sample_rho_g_cm3 отличается в файле {filepath}")
        
        # Проверяем, что physics параметры совпадают
        physics_node = {
            "em_cut_mm": header_keys.get("em_cut_mm"),
            "em_deex": header_keys.get("em_deex"),
            "em_option": header_keys.get("em_option"),
            "particle": header_keys.get("particle"),
            "npsm_enabled": header_keys.get("npsm_enabled")
        }
        
        if physics_params is None:
            physics_params = physics_node
        else:
            for key in physics_params:
                if physics_params[key] != physics_node[key]:
                    raise SystemExit(f"ОТКАЗ: параметр {key} отличается в файле {filepath}")
        
        # Строим словарь центр -> отсчёты edep
        h = {}
        for center, edep in counts:
            h[center] = edep
        
        # Определяем bin_centers (должны совпадать для всех узлов)
        bin_centers_node = sorted(h.keys())
        if bin_centers is None:
            bin_centers = bin_centers_node
        else:
            if bin_centers_node != bin_centers:
                raise SystemExit(f"ОТКАЗ: bin_centers отличаются в файле {filepath}")
        
        # #GS-66 (v1.1): пик без континуума под ним — общее определение fep_def.fep_net (v1: окно «бин с E + бин ниже»)
        k_fep, var_fep = fep_net(h, E)
        k_fep = max(k_fep, 0.0)   # узлы 12–16 кэВ: пика нет, вычет не уводит ε ниже нуля
        eps_fep = k_fep / N
        sig_fep = np.sqrt(var_fep) / N
        
        # k_tot = Σ h[c] для c ≥ 1.0
        k_tot = sum(h[c] for c in h if c >= 1.0)
        
        # eps_tot = k_tot / N, sig_tot = sqrt(k_tot)/N
        eps_tot = k_tot / N
        sig_tot = np.sqrt(k_tot) / N
        
        # Отказ, если eps_tot > 0.6 или k_fep > k_tot
        if eps_tot > 0.6:
            raise SystemExit(f"ОТКАЗ: eps_tot > 0.6 для узла E={E} (бин нулевого отложения попал в сумму)")
        if k_fep > k_tot:
            raise SystemExit(f"ОТКАЗ: k_fep > k_tot для узла E={E}")
        
        # Сохраняем данные узла
        nodes_data.append({
            "E": E,
            "N": N,
            "h": h,
            "k_fep": k_fep,
            "eps_fep": eps_fep,
            "sig_fep": sig_fep,
            "k_tot": k_tot,
            "eps_tot": eps_tot,
            "sig_tot": sig_tot,
            "count_light": count_light,
            "seed": header_keys.get("seed", "")
        })
    
    if not nodes_data:
        raise SystemExit("ОТКАЗ: нет данных узлов")
    
    # Сортируем узлы по возрастанию E
    nodes_data.sort(key=lambda x: x["E"])
    
    # Проверяем, что bin_centers не None
    if bin_centers is None:
        raise SystemExit("ОТКАЗ: не удалось определить bin_centers")
    
    # Проверяем, что ширина бина = 1.0
    if len(bin_centers) > 1:
        bin_width = bin_centers[1] - bin_centers[0]
        if abs(bin_width - 1.0) > 1e-9:
            raise SystemExit(f"ОТКАЗ: ширина бина {bin_width} != 1.0")
    
    # Создаем массивы для npz
    M = len(nodes_data)
    B = len(bin_centers)
    
    E_nodes = np.array([node["E"] for node in nodes_data], dtype=np.float64)
    n_primary = np.array([node["N"] for node in nodes_data], dtype=np.int64)
    seed = np.array([node["seed"] for node in nodes_data], dtype=str)
    
    # Создаем матрицы counts_edep и counts_light
    counts_edep = np.zeros((M, B), dtype=np.int64)
    counts_light = np.zeros((M, B), dtype=np.int64)
    
    # Создаем словарь центр -> индекс
    bin_center_to_idx = {center: idx for idx, center in enumerate(bin_centers)}
    
    for i, node in enumerate(nodes_data):
        # Заполняем counts_edep
        for center, edep in node["h"].items():
            if center in bin_center_to_idx:
                counts_edep[i, bin_center_to_idx[center]] = edep
        
        # Заполняем counts_light
        # count_light — это список, длина которого должна быть равна B
        if len(node["count_light"]) == B:
            counts_light[i, :] = node["count_light"]
        else:
            # Если длины не совпадают, заполняем нулями (предупреждение)
            print(f"Предупреждение: длина count_light ({len(node['count_light'])}) != B ({B}) для узла E={node['E']}")
            min_len = min(len(node["count_light"]), B)
            counts_light[i, :min_len] = node["count_light"][:min_len]
    
    # Выход 1: efficiency_<TAG>.csv
    csv_path = os.path.join(OUT_DIR, f"efficiency_{TAG}.csv")
    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["E_keV", "eps_fep", "sig_fep", "eps_total", "sig_total", "n_primary", "k_fep", "k_total"])
        for node in nodes_data:
            writer.writerow([
                f"{node['E']:.3f}",
                f"{node['eps_fep']:.6e}",
                f"{node['sig_fep']:.6e}",
                f"{node['eps_tot']:.6e}",
                f"{node['sig_tot']:.6e}",
                int(node["N"]),
                "%.1f" % node["k_fep"],   # #GS-66: за вычетом φ·континуума — нецелое
                int(node["k_tot"])
            ])
    
    # Выход 2: response_matrix_<TAG>.npz
    npz_path = os.path.join(OUT_DIR, f"response_matrix_{TAG}.npz")
    np.savez_compressed(
        npz_path,
        E_nodes=E_nodes,
        bin_centers_keV=np.array(bin_centers, dtype=np.float64),
        counts_edep=counts_edep,
        counts_light=counts_light,
        n_primary=n_primary,
        seed=seed
    )
    
    # Выход 4: Копия GDML
    gdml_dest = os.path.join(OUT_DIR, f"geometry_{TAG}.gdml")
    shutil.copy2(GDML_SRC, gdml_dest)
    gdml_sha256 = sha256_file(gdml_dest)
    
    # Выход 3: meta_<TAG>.json
    meta = {
        "detector": "Atom GS2020 PRO MAX",
        "vessel": "Маринелли 1 л (РАДЭК)",
        "tag": TAG,
        "created": datetime.datetime.now().isoformat(),
        "geometry_gdml": f"geometry_{TAG}.gdml",
        "geometry_sha256": gdml_sha256,
        "grid_dir": os.path.basename(GRID_DIR),   # локальная папка прогонов, не в репозитории
        "n_nodes": M,
        "E_min_keV": float(E_nodes[0]),
        "E_max_keV": float(E_nodes[-1]),
        "bins_keV": 1.0,
        "n_bins": B,
        "sample_matrix": sample_matrix,
        "sample_rho_g_cm3": sample_rho_g_cm3,
        "physics": physics_params,
        "n_primary_min": int(n_primary.min()),
        "n_primary_max": int(n_primary.max()),
        "fep_definition": "v1.1 (#GS-66): k_fep = пик без континуума — целое E: [E-1,E+1) минус бин [E-2,E-1); дробное E: бин с E минус φ·(бин ниже), φ = дробная часть E; eps_fep = k_fep/N. v1: бин с E + бин ниже (включал 1–1,5 кэВ континуума, ниже 100 кэВ +3–8 %)",
        "total_definition": "сумма бинов ≥ 1 кэВ / N",
        "files_sha256": {}
    }
    
    # Вычисляем SHA256 для файлов
    meta["files_sha256"]["efficiency"] = sha256_file(csv_path)
    meta["files_sha256"]["response_matrix"] = sha256_file(npz_path)
    meta["files_sha256"]["geometry"] = gdml_sha256
    
    meta_path = os.path.join(OUT_DIR, f"meta_{TAG}.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=1)
    
    # Самопроверка
    print("\n=== Самопроверка ===")
    
    # Читаем npz обратно
    npz_data = np.load(npz_path, allow_pickle=True)
    counts_edep_loaded = npz_data["counts_edep"]
    E_nodes_loaded = npz_data["E_nodes"]
    
    # Находим узел с E ближайшей к 1460
    idx_1460 = np.argmin(np.abs(E_nodes_loaded - 1460))
    E_1460 = E_nodes_loaded[idx_1460]
    
    # Находим узел с E ближайшей к 238
    idx_238 = np.argmin(np.abs(E_nodes_loaded - 238))
    E_238 = E_nodes_loaded[idx_238]
    
    # Для узла 1460: проверяем сумму counts_edep
    sum_counts_edep_1460 = counts_edep_loaded[idx_1460].sum()
    
    # Находим исходный файл для узла 1460
    # Ищем файл с E близким к 1460
    target_file_1460 = None
    for filepath in grid_files:
        with open(filepath, "r", encoding="utf-8") as f:
            first_line = f.readline().strip()
        # Пропускаем проверку GDML для самопроверки
        grid_data = parse_grid_file(filepath)
        if grid_data is not None and abs(grid_data["E"] - E_1460) < 1e-6:
            target_file_1460 = filepath
            break
    
    if target_file_1460 is not None:
        grid_data_1460 = parse_grid_file(target_file_1460)
        if grid_data_1460 is not None:
            sum_original_1460 = sum(c for _, c in parse_header_and_light(target_file_1460)[2])   # бины — из собственного разбора
            print(f"Узел E={E_1460:.3f} кэВ (близко к 1460):")
            print(f"  Сумма counts_edep из npz: {sum_counts_edep_1460}")
            print(f"  Сумма из исходного файла: {sum_original_1460}")
            if sum_counts_edep_1460 != sum_original_1460:
                print("  ПРЕДУПРЕЖДЕНИЕ: суммы не совпадают!")
    
    # Для узла 238: проверяем сумму counts_edep
    sum_counts_edep_238 = counts_edep_loaded[idx_238].sum()
    
    # Находим исходный файл для узла 238
    target_file_238 = None
    for filepath in grid_files:
        grid_data = parse_grid_file(filepath)
        if grid_data is not None and abs(grid_data["E"] - E_238) < 1e-6:
            target_file_238 = filepath
            break
    
    if target_file_238 is not None:
        grid_data_238 = parse_grid_file(target_file_238)
        if grid_data_238 is not None:
            sum_original_238 = sum(c for _, c in parse_header_and_light(target_file_238)[2])   # бины — из собственного разбора
            print(f"Узел E={E_238:.3f} кэВ (близко к 238):")
            print(f"  Сумма counts_edep из npz: {sum_counts_edep_238}")
            print(f"  Сумма из исходного файла: {sum_original_238}")
            if sum_counts_edep_238 != sum_original_238:
                print("  ПРЕДУПРЕЖДЕНИЕ: суммы не совпадают!")
    
    # eps_fep узла ближайшего к 1460
    node_1460 = nodes_data[idx_1460]
    print(f"\nУзел E={node_1460['E']:.3f} кэВ (близко к 1460):")
    print(f"  E = {node_1460['E']:.3f} кэВ")
    print(f"  k_fep = {node_1460['k_fep']}")
    print(f"  N = {node_1460['N']}")
    print(f"  eps_fep = {node_1460['eps_fep']:.6e}")
    
    # Печать: число узлов, E_min…E_max, путь каждого выходного файла и его размер в КБ
    print(f"\nЧисло узлов: {M}")
    print(f"E_min = {E_nodes[0]:.3f} кэВ")
    print(f"E_max = {E_nodes[-1]:.3f} кэВ")
    
    print("\nВыходные файлы:")
    for filepath in [csv_path, npz_path, gdml_dest, meta_path]:
        size_kb = os.path.getsize(filepath) / 1024.0
        print(f"  {filepath} ({size_kb:.1f} КБ)")
    
    print("\nГотово.")


if __name__ == "__main__":
    main()
