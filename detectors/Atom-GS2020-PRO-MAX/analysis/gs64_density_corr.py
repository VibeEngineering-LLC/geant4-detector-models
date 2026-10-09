# -*- coding: utf-8 -*-
#GS-64: f(E, rho, V) = eps_fep(вариант) / eps_fep(ОИСН-10 rho 0.9061, 1161 мл); контроль против набора.

import sys
import os
import re
import csv
import math
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

from compare_eff_becqmoni import parse_grid_file
from make_dataset_marinelli1l import parse_header_and_light   # бины узла
from fep_def import fep_net   # #GS-66: пик без континуума — то же определение, что в наборе v1.1


def main():
    if len(sys.argv) != 4:
        raise SystemExit("Использование: python gs64_density_corr.py <out_dir> <efficiency.csv> <dst.csv>")

    out_dir = sys.argv[1]
    efficiency_csv = sys.argv[2]
    dst_csv = sys.argv[3]

    # --- 1. Чтение файлов из out_dir ---
    pattern = re.compile(r'^(oisn|kcl)_r([\d.]+)_v(\d+)_E([\d.]+)_s\d+\.csv$')
    
    data = {}  # ключ: (mat, rho, V, E) -> (eps, sig)
    
    if not os.path.isdir(out_dir):
        raise SystemExit("ОТКАЗ: нет файлов в " + out_dir)
    
    files = [x for x in os.listdir(out_dir) if x.endswith(".csv")]   # рядом лежат .log прогонов
    if not files:
        raise SystemExit("ОТКАЗ: нет файлов в " + out_dir)
    
    for fname in files:
        m = pattern.match(fname)
        if not m:
            raise SystemExit("ОТКАЗ: не разобран " + os.path.join(out_dir, fname))
        
        mat = m.group(1)
        rho_str = m.group(2)
        V_str = m.group(3)
        E_str = m.group(4)
        
        rho = float(rho_str)
        V = int(V_str)
        E_name = float(E_str)
        
        path = os.path.join(out_dir, fname)
        g = parse_grid_file(path)
        
        if g is None:
            raise SystemExit("ОТКАЗ: не разобран " + path)
        
        E_val = g["E"]
        N = g["N"]
        counts, var = fep_net(dict(parse_header_and_light(path)[2]), g["E"])
        
        if abs(E_val - E_name) > 1e-9:
            raise SystemExit("ОТКАЗ: не разобран " + path)
        
        if counts <= 0:
            raise SystemExit("ОТКАЗ: не разобран " + path)
        
        eps = counts / N
        sig = math.sqrt(var) / N
        
        key = (mat, rho, V, E_val)
        
        if key in data:
            raise SystemExit("ОТКАЗ: дубль " + str(key))
        
        data[key] = (eps, sig)
    
    if not data:
        raise SystemExit("ОТКАЗ: нет файлов в " + out_dir)
    
    # --- 2. Чтение efficiency.csv ---
    E_keV_list = []
    eps_fep_list = []
    rel_sig_list = []
    
    with open(efficiency_csv, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if float(row['eps_fep']) <= 0:   # 6 узлов 12–15 кэВ без пика — в log-интерполяцию не берутся
                continue
            E_keV_list.append(float(row['E_keV']))
            eps_fep_list.append(float(row['eps_fep']))
            rel_sig_list.append(float(row['sig_fep']) / float(row['eps_fep']))

    if not E_keV_list:
        raise SystemExit("ОТКАЗ: нет данных в efficiency.csv")
    
    # Сортировка по E_keV
    sorted_indices = sorted(range(len(E_keV_list)), key=lambda i: E_keV_list[i])
    E_keV_sorted = [E_keV_list[i] for i in sorted_indices]
    eps_fep_sorted = [eps_fep_list[i] for i in sorted_indices]
    
    log_E_keV = np.log(np.array(E_keV_sorted))
    log_eps_fep = np.log(np.array(eps_fep_sorted))
    rel_sig_sorted = np.array([rel_sig_list[i] for i in sorted_indices])
    
    # --- 3. Определение опорного варианта REF ---
    REF_KEY = ("oisn", 0.9061, 1161)
    
    # Собираем все E, встречающиеся у любого варианта
    all_E = set()
    for key in data:
        all_E.add(key[3])
    
    # Проверяем, что для каждого E существует ref
    ref_data = {}
    for key, (eps, sig) in data.items():
        if key[0] == REF_KEY[0] and key[1] == REF_KEY[1] and key[2] == REF_KEY[2]:
            ref_data[key[3]] = (eps, sig)
    
    for E in all_E:
        if E not in ref_data:
            raise SystemExit("ОТКАЗ: нет опорного узла E=" + str(E))
    
    # --- 4. Расчёт f и sig_f ---
    results = []
    
    for key in sorted(data.keys()):
        mat, rho, V, E = key
        eps, sig = data[key]
        
        eps0, sig0 = ref_data[E]
        
        f = eps / eps0
        sig_f = f * math.sqrt((sig / eps) ** 2 + (sig0 / eps0) ** 2)
        
        results.append((mat, rho, V, E, eps, sig, f, sig_f))
    
    # --- 5. Контроль против набора ---
    print("КОНТРОЛЬ:")
    max_dev = 0.0
    out_of_range = 0
    n_nodes = 0
    
    for E in sorted(ref_data.keys()):
        eps, sig = ref_data[E]
        
        # Интерполяция в логарифмическом масштабе
        log_E = math.log(E)
        j = int(np.searchsorted(log_E_keV, log_E))
        if j == 0 or j == len(log_E_keV) or log_E_keV[j] - log_E_keV[j - 1] > math.log(1.5):   # #GS-65: дыра набора 15,8–70,8 кэВ
            print("КОНТРОЛЬ E=%7.1f eps=%.5e ±%.2f%% набор: НЕТ УЗЛОВ (соседи %.1f–%.1f кэВ) — не сравнивается" % (
                E, eps, sig / eps * 100, math.exp(log_E_keV[max(j - 1, 0)]), math.exp(log_E_keV[min(j, len(log_E_keV) - 1)])))
            continue
        log_eps_set = np.interp(log_E, log_E_keV, log_eps_fep)
        eps_set = math.exp(log_eps_set)

        ratio = eps / eps_set
        dev_sig = (ratio - 1) / math.hypot(sig / eps, float(np.interp(log_E, log_E_keV, rel_sig_sorted)))

        print("КОНТРОЛЬ E=%7.1f eps=%.5e ±%.2f%% набор=%.5e отношение=%.4f откл_сигм=%.1f" % (
            E, eps, sig / eps * 100, eps_set, ratio, dev_sig
        ))
        
        n_nodes += 1
        abs_dev = abs(ratio - 1)
        if abs_dev > max_dev:
            max_dev = abs_dev
        if abs_dev > 0.005:
            out_of_range += 1
    
    print("КОНТРОЛЬ ИТОГ: узлов %d, max|отношение−1| = %.2f %%, узлов вне ±0,5 %%: %d" % (
        n_nodes, max_dev * 100, out_of_range
    ))
    
    # --- 6. Запись dst.csv ---
    with open(dst_csv, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['material', 'rho_g_cm3', 'V_ml', 'E_keV', 'eps_fep', 'sig_eps', 'f', 'sig_f'])
        
        for (mat, rho, V, E, eps, sig, f_val, sig_f) in results:
            writer.writerow([
                mat,
                '%.4f' % rho,
                str(V),
                '%g' % E,
                '%.6e' % eps,
                '%.3e' % sig,
                '%.5f' % f_val,
                '%.5f' % sig_f
            ])
    
    print("записан %s строк %d" % (dst_csv, len(results)))


if __name__ == '__main__':
    main()
