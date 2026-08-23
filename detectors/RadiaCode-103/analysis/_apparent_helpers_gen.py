"""
#FIT-1 Stage-3: apparent_method1_compare.py
Проверяет гипотезу о влиянии LCE-tailing на 62%-дефицит континуума.
Читает боевые шаблоны метода 1 и Stage-2 Y-распределения, строит
"apparent"-версии с учётом карты LCE(Y), применяет те же NNLS-амплитуды,
сравнивает метрики качества подгонки и сохраняет результаты в CSV.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")

import io, contextlib, math, os
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import numpy as np
import fit_nuclides as fn
import fit_coverage as fc
import rcspec
import read_rcxml

LCE_MAP = {-4.5: 25.0, -3.0: 25.7, -1.5: 23.6, 0.0: 20.0,
           1.5: 18.8, 3.0: 18.2, 4.5: 15.8}
STAGE2_DIR = "D:/Claude_files/repos/geant4-detector-models/build/RadiaCode-103/_stage2_ypos"

def load_ypos(path):
    """Читает _ypos.csv. Возвращает (y_cols, dist):
    y_cols: list[float] — Y-координаты колонок в порядке файла.
    dist: dict[int, np.ndarray] — ключ: индекс 1-кэВ канала (i = round(E_keV - 0.5)),
    значение: np.array counts по Y-колонкам (той же длины, что y_cols).
    Если файл не существует — возвращает (None, None)."""
    if not os.path.exists(path):
        return None, None
    y_cols = None
    dist = {}
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(',')
            if y_cols is None:
                y_cols = [float(p[1:]) for p in parts[1:]]
                continue
            e_kev = float(parts[0])
            i = round(e_kev - 0.5)
            counts = np.array([int(x) for x in parts[1:]])
