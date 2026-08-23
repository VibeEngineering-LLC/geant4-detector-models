Задача: написать ТОЛЬКО вспомогательные функции для Python-модуля (без main(),
без блока `if __name__`). Верни ПОЛНЫЙ код от первой строки (docstring) до
последней функции — без markdown fences, без пояснений до/после, чистый Python.

## Шапка модуля

```python
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
```

Скопируй эту шапку ДОСЛОВНО в начало файла, затем добавь 5 функций ниже.

## Функция 1: `load_ypos(path)`

Формат входного файла (`_ypos.csv`):
```
# ...шапка с #-комментариями (строки начинаются с "#")...
E_keV,y-4.50,y-3.00,y-1.50,y0.00,y1.50,y3.00,y4.50
351.5,3,5,4,2,6,1,0
...
```
Первая непустая, не-`#`-строка — заголовок колонок (`E_keV,y-4.50,...`).
Каждая колонка после первой называется `y<Y>` (Y может быть отрицательным,
формат `%.2f`, например `y-4.50` или `y1.50`). Дальше идут строки данных:
первое число — `E_keV` (float, центр 1-кэВ канала, например `351.5`), дальше —
целочисленные counts по каждой Y-колонке.

```python
def load_ypos(path):
    """Читает _ypos.csv. Возвращает (y_cols, dist):
    y_cols: list[float] — Y-координаты колонок в порядке файла.
    dist: dict[int, np.ndarray] — ключ: индекс 1-кэВ канала (i = round(E_keV - 0.5)),
    значение: np.array counts по Y-колонкам (той же длины, что y_cols).
    Если файл не существует — возвращает (None, None)."""
```
Пропускать пустые строки и строки-комментарии (`#`). Заголовок — первая строка
с данными (`E_keV,...`), она НЕ идёт в `dist`. `y_cols[k] = float(имя_колонки[1:])`
(убрать первую букву `y`, распарсить остаток как float — работает и для
отрицательных, `float("-4.50")` даёт -4.5).

## Функция 2: `global_lce_ref(all_dists)`

```python
def global_lce_ref(all_dists):
    """all_dists: list[tuple[list[float], dict[int, np.ndarray]]] — пары
    (y_cols, dist) по всем нуклидам с непустыми Stage-2 данными.
    Суммирует ВСЕ counts по ВСЕМ каналам ВСЕХ нуклидов в один вектор длины
    len(y_cols) (Y-колонки одинаковы у всех нуклидов — взять y_cols первого
    элемента как канон). Возвращает (lce_ref, n_total):
    lce_ref = взвешенное среднее LCE_MAP[y] по суммарным counts (float, %).
    n_total = сумма всех counts (для диагностической печати).
    Если all_dists пуст — вернуть (0.0, 0)."""
```

## Функция 3: `build_apparent_hist(raw_hist, y_cols, dist, lce_ref)`

```python
def build_apparent_hist(raw_hist, y_cols, dist, lce_ref):
    """raw_hist: np.ndarray длины rcspec.NBINS — RAW (до rcspec.fold)
    энерговыделение, канал i соответствует энергии i+0.5 кэВ.
    Возвращает np.ndarray той же длины — "apparent" версию (LCE-размазанную).

    Алгоритм:
    1. global_p — резервное относительное Y-распределение по ВСЕМ каналам
       этого нуклида: сумма всех dist.values() в один вектор, нормировать на
       сумму (если сумма 0 -> global_p = None).
    2. out = np.zeros_like(raw_hist).
    3. Для каждого канала i (0..len(raw_hist)-1): если raw_hist[i] <= 0 -
       пропустить.
    4. counts = dist.get(i). Если counts is None или counts.sum() < 20 -
       p = global_p. Иначе p = counts / counts.sum().
    5. Если p is None: out[i] += raw_hist[i] (не размазывать), продолжить
       со следующим i.
    6. Для каждого k, y в enumerate(y_cols): если p[k] <= 0 - пропустить.
       e_center = i + 0.5; e_app = e_center * LCE_MAP[y] / lce_ref;
       j = int(e_app). Если 0 <= j < len(out): out[j] += raw_hist[i] * p[k].
       Иначе (вышло за границу) - ничего не делать, событие теряется молча.
    7. Вернуть out."""
```

## Функция 4: `load_column(nuc)`

```python
def load_column(nuc):
    """Аналог тела fn.load_templates() для ОДНОГО нуклида, но возвращает RAW
    (до fold) массив, а не cps. Возвращает (raw_norm, meta) или (None, None)
    если файлов нет."""
    wf = os.path.join(fn.BUILD, "%s_%s.csv" % (fn.WF_PREFIX, nuc))
    bg = os.path.join(fn.BG_DIR, "%s_%s.csv" % (fn.BG_PREFIX, nuc))
    if not (os.path.exists(wf) and os.path.exists(bg)):
        return None, None
    flu = fn.read_wallfield_total(wf)
    r, hz = fn.CYL["r"] / 10.0, 0.5 * (fn.CYL["z1"] - fn.CYL["z0"]) / 10.0
    area = 2 * math.pi * r * (r + 2 * hz)
    rate = flu * area / 4.0
    meta, hist = rcspec.read_spec(bg)
    t_run = float(meta["N_primaries"]) / rate
    raw_norm = hist / t_run
    return raw_norm, meta
```
Скопируй эту функцию ДОСЛОВНО как есть (уже готовый код, ничего не менять).

## Функция 5: `predict(names, A, amp_dict)`

```python
def predict(names, A, amp_dict):
    pred = np.zeros(A.shape[0])
    for name, a in amp_dict.items():
        if name in names:
            pred += a * A[:, names.index(name)]
    return pred
```
Скопируй эту функцию ДОСЛОВНО как есть.

Верни ПОЛНЫЙ файл: шапка модуля + все 5 функций по порядку, ничего больше
(НЕ добавляй main(), НЕ добавляй `if __name__`).
