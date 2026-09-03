# -*- coding: utf-8 -*-
"""
Сравнение прогонов программы run_field (Geant4-модель RadiaCode-103),
различающихся ТОЛЬКО физикой (порог продукции `em_cut_mm` и/или режим атомной деэкситации `em_deex`),
по энергетическим полосам. Дата: 03.09.2026.
Задача — проверка утверждения «Geant4 упрощает физику на краях спектра ради скорости счёта».

CLI:
    python compare_physics_runs.py <base.csv> <test.csv> [<test2.csv> ...]
Без аргументов или с одним — напечатать docstring модуля и вернуть код 2.

Причины отказа:
1. Один из файлов не содержит ключей em_cut_mm или em_deex.
2. Постановки различаются по ключам shield, stand_mm, screen_up, n_events, flux_total_cm2_s.
"""
import os
import sys
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
import fit_two_criteria as ftc


def load(path):
    meta, cps, counts = ftc.read_template(path)
    if 'em_cut_mm' not in meta or 'em_deex' not in meta:
        print(f"ОТКАЗ: {path} не содержит em_cut_mm/em_deex — старый формат CSV, физика прогона не записана")
        sys.exit(2)
    return meta, cps, counts


def check_same_setup(meta_a, meta_b, name_b):
    keys = ['shield', 'stand_mm', 'screen_up', 'n_events', 'flux_total_cm2_s']
    for key in keys:
        a_val = meta_a.get(key)
        b_val = meta_b.get(key)
        if a_val != b_val:
            print(f"ОТКАЗ: постановки различаются по '{key}': base={a_val} vs {name_b}={b_val}")
            sys.exit(2)


def band_table(base_meta, test_meta, base_counts, test_counts):
    print(f"=== {os.path.basename(test_meta['file'])}: em_cut_mm {base_meta['em_cut_mm']}→{test_meta['em_cut_mm']}, em_deex {base_meta['em_deex']}→{test_meta['em_deex']} ===")
    
    base_hits = np.sum(base_counts)
    test_hits = np.sum(test_counts)
    print(f"hits base={int(base_hits)} test={int(test_hits)}")
    
    bands = list(ftc.BANDS) + [(20, 2999)]
    print("полоса, кэВ\tc_base\tc_test\ttest/base\t±")
    
    for e1, e2 in bands:
        c_base = np.sum(base_counts[e1:e2])
        c_test = np.sum(test_counts[e1:e2])
        
        if c_base == 0:
            ratio = np.inf
            error = np.inf
        else:
            ratio = c_test / c_base
            error = ratio * np.sqrt(1/c_base + 1/c_test)
        
        stat_note = ""
        if (c_base < 100) or (c_test < 100):
            stat_note = "  (мало статистики)"
        
        print(f"{e1}–{e2}\t{c_base:.0f}\t{c_test:.0f}\t{ratio:.4f}\t{error:.4f}{stat_note}")


def main(argv):
    sys.stdout.reconfigure(encoding="utf-8")
    
    if len(argv) < 3:
        print(__doc__)
        return 2
    
    base_path = argv[1]
    test_paths = argv[2:]
    
    base_meta, base_cps, base_counts = load(base_path)
    base_meta['file'] = base_path
    
    for test_path in test_paths:
        test_meta, test_cps, test_counts = load(test_path)
        test_meta['file'] = test_path
        
        check_same_setup(base_meta, test_meta, os.path.basename(test_path))
        band_table(base_meta, test_meta, base_counts, test_counts)
    
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
