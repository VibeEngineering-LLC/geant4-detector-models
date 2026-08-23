# Спека: apparent_method1_compare.py (#FIT-1, Stage-3)

## Контекст и цель

Проверяем гипотезу: light-collection tailing (D-007/D-008, разброс LCE(Y) по
объёму кристалла 15,8%…25,7%) частично объясняет 62%-дефицит континуума в
методе 1 (`fit_nuclides.py`).

Для 8 нуклидов метода 1 (`K40, Ra226, Pb214, Bi214, Pb212, Ac228, Bi212, Tl208`)
уже посчитаны Stage-2 прогоны: файлы
`D:/Claude_files/repos/geant4-detector-models/build/RadiaCode-103/_stage2_ypos/s2_<нуклид>_ypos.csv`
дают распределение по 7 Y-бинам (energy-weighted координата взаимодействия) для
каждого 1-кэВ канала энерговыделения (формат ниже).

Задача: взять КАЖДЫЙ боевой (30М статистика) шаблон метода 1
(`bg_bare_field_m1_<нуклид>.csv`), «размазать» его по энергии согласно карте
LCE(Y) и относительному Y-распределению из Stage-2 (P(y|E)), получить
"apparent"-версию шаблона, прогнать её через ТЕ ЖЕ NNLS-амплитуды, что уже
подобраны боевым `fit_nuclides.fit_by_lines`, и сравнить итоговый
предсказанный спектр (baseline vs apparent) с измерением — обеими метриками
(`fit_coverage.form_residual_pct`, `fit_coverage.chi2_of`) и по полосам энергии.

## Формат входных файлов

`_ypos.csv` (Stage-2, `main.cc`/`rc_curves.exe` вывод):
```
# ...шапка с #-комментариями...
# y_bin_mm = 1.50, y_min_mm = -5.25, y_bins = 7
E_keV,y-4.50,y-3.00,y-1.50,y0.00,y1.50,y3.00,y4.50
351.5,3,5,4,2,6,1,0
...
```
Первая колонка — центр 1-кэВ канала (`E_keV`, float, но фактически целое+0.5).
Остальные 7 колонок — integer counts, заголовок `y<Y>` (Y может быть
отрицательным, формат `%.2f`). Строки со всеми нулями не пишутся (разрежённый
формат) — отсутствие строки для канала E значит 0 событий во всех Y-бинах
этого канала.

`bg_bare_field_m1_<нуклид>.csv` (боевой шаблон, читается через
`rcspec.read_spec(path)` -> `(meta: dict, hist: np.ndarray[NBINS])`,
`rcspec.NBINS = 3201`, канал i соответствует энергии `i + 0.5` кэВ.
`meta["N_primaries"]` — строка, привести к float.

## Файл `apparent_method1_compare.py` (новый, в той же папке, что fit_coverage.py)

Шапка модуля (docstring, 5-8 строк): что делает, зачем (#FIT-1 Stage-3), какие
файлы читает, что печатает. Кодировка: `sys.stdout.reconfigure(encoding="utf-8")`
в начале файла (до любых import, как в continuum_deficit.py/tail_check.py —
переиспользовать тот же паттерн, эти файлы лежат рядом, посмотри их структуру
как образец стиля: `sys.path` настройка, suppress stdout через
`io.StringIO()`/`contextlib.redirect_stdout`, финальный `if __name__ == "__main__": main()`).

### Импорты и путь

```python
import io, contextlib, math, os, sys
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import numpy as np
import fit_nuclides as fn
import fit_coverage as fc
import rcspec
import read_rcxml
```

### Константа LCE(Y), мм -> %, из D-007/D-008 (opticalcheck.cc, N=200000)

```python
LCE_MAP = {-4.5: 25.0, -3.0: 25.7, -1.5: 23.6, 0.0: 20.0,
           1.5: 18.8, 3.0: 18.2, 4.5: 15.8}
STAGE2_DIR = "D:/Claude_files/repos/geant4-detector-models/build/RadiaCode-103/_stage2_ypos"
```

### Функция `load_ypos(path)`

Читает `_ypos.csv`. Возвращает `(y_cols: list[float], dist: dict[int, np.ndarray])`
где `y_cols` — список Y-координат колонок в порядке файла (float, из заголовка
`y<Y>` — распарсить как `float(col_name[1:])`), `dist` — словарь: ключ —
целочисленный номер канала (`int(round(float(E_keV_str) - 0.5))`, т.е. индекс
канала для `rcspec` — E_keV в файле это ЦЕНТР канала `i+0.5`, поэтому
`i = round(E_keV - 0.5)`), значение — `np.array` длины `len(y_cols)` с counts
(int -> float). Пропускать пустые строки и строку заголовка. Если файл не
существует — вернуть `(None, None)` (вызывающий код обрабатывает как «нет
Stage-2 данных для этого нуклида», не падает).

### Функция `global_lce_ref(all_dists)`

`all_dists` — список `(y_cols, dist)` по всем нуклидам, у которых Stage-2 есть
(не None). Суммирует ВСЕ counts по ВСЕМ каналам ВСЕХ нуклидов в один вектор
длины 7 (по Y-бинам — все `y_cols` одинаковы между нуклидами, 7 точек
`-4.5..4.5`, можно взять `y_cols` первого элемента как канон). Возвращает
взвешенное среднее `sum(counts[k] * LCE_MAP[y_cols[k]] for k) / sum(counts)`
— это LCE_ref (float, %). Также вернуть общее число просуммированных событий
(для печати в диагностику).

### Функция `build_apparent_hist(raw_hist, y_cols, dist, lce_ref)`

`raw_hist` — `np.ndarray` длины `rcspec.NBINS`, RAW (до `rcspec.fold`)
энерговыделение канал->counts (нормированное на cps, как в `fn.load_templates`
до вызова `rcspec.fold`, см. ниже как получить). Возвращает `np.ndarray` той
же длины — "apparent" версию.

Алгоритм:
1. Посчитать `global_p` — резервное (fallback) относительное Y-распределение
   по ВСЕМ каналам этого нуклида: просуммировать все `dist.values()` в один
   вектор длины 7, нормировать на сумму (если сумма 0 — `global_p = None`).
2. `out = np.zeros_like(raw_hist)`.
3. Для каждого индекса канала `i` от 0 до `len(raw_hist)-1`, если
   `raw_hist[i] <= 0` — пропустить (`continue`).
4. `counts = dist.get(i)`. Если `counts is None` или `counts.sum() < 20`
   (мало статистики в этом конкретном канале Stage-2 для надёжного P(y|E)) —
   использовать `p = global_p`. Иначе `p = counts / counts.sum()`.
5. Если `p is None` (совсем нет Stage-2 данных) — `out[i] += raw_hist[i]`
   (не размазывать, оставить как есть) и `continue`.
6. Для каждого `k, y в enumerate(y_cols)`: если `p[k] <= 0` — пропустить.
   `e_center = i + 0.5`; `e_app = e_center * LCE_MAP[y] / lce_ref`;
   `j = int(e_app)` (индекс канала apparent, floor через int() от
   неотрицательного числа). Если `0 <= j < len(out)`:
   `out[j] += raw_hist[i] * p[k]`. Иначе (ушло за границу диапазона) —
   ничего не делать, событие теряется молча ЭТО ОЖИДАЕМО и ДОПУСТИМО (редкий
   краевой эффект, не варнинг).

### Функция `load_column(nuc)`

Аналог тела `fn.load_templates()` для ОДНОГО нуклида, но возвращает RAW
(до fold) массив, а не cps. Скопировать логику нормировки один-в-один:
```python
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
(`fn.CYL`, `fn.WF_PREFIX`, `fn.BG_PREFIX`, `fn.BUILD`, `fn.read_wallfield_total`
уже существуют в `fit_nuclides.py` — переиспользовать as-is, НЕ пересчитывать
формулу площади самостоятельно.)

### Функция `predict(names, A, amp_dict)`

```python
def predict(names, A, amp_dict):
    pred = np.zeros(A.shape[0])
    for name, a in amp_dict.items():
        if name in names:
            pred += a * A[:, names.index(name)]
    return pred
```

### `main()`

1. Печать заголовка: `"=== #FIT-1 Stage-3: apparent-спектр метода 1 (LCE-tailing) ==="`.
2. Для каждого `nuc in fn.NUCS` (порядок как в `fit_nuclides.NUCS`):
   вызвать `load_ypos(os.path.join(STAGE2_DIR, "s2_%s_ypos.csv" % nuc))`.
   Собрать список `(nuc, y_cols, dist)` для тех, где `dist is not None`.
3. Если список пуст — напечатать `"НЕТ STAGE-2 ДАННЫХ, СТОП"` и `return`.
4. `lce_ref, n_total = global_lce_ref([(yc, d) for _, yc, d in список])`.
   Напечатать `"LCE_ref (глобальный, %d событий всех нуклидов) = %.2f %%" % (n_total, lce_ref)`.
5. Для каждого `nuc in fn.NUCS`:
   - `raw, meta = load_column(nuc)`. Если `raw is None` — `print("[--] %s: нет шаблона" % nuc)`, `continue`.
   - Найти в списке шага 2 запись для этого `nuc` (по имени). Если есть —
     `apparent_raw = build_apparent_hist(raw, y_cols, dist, lce_ref)`; иначе
     (Stage-2 для этого нуклида не считался) — `apparent_raw = raw.copy()`
     и `print("[--] %s: нет Stage-2 Y-данных, apparent = baseline" % nuc)`.
   - `cps_baseline = rcspec.fold(raw, "103")`
   - `cps_apparent = rcspec.fold(apparent_raw, "103")`
   - Накопить в списки `names_raw`, `cols_baseline`, `cols_apparent`.
6. `names_b, cols_b = fn.merge_by_chain(names_raw, cols_baseline)`
   `names_a, cols_a = fn.merge_by_chain(names_raw, cols_apparent)`
   (`assert names_b == names_a`, иначе `raise SystemExit` с понятным текстом —
   не должно случиться, но проверить явно, не молчать).
7. `mu, pdg = fn.load_muons()`; если `mu is not None`: `names_b.append("mu")`,
   `cols_b.append(mu)`; `names_a.append("mu")`, `cols_a.append(mu)`.
8. Измерение (как в `fit_nuclides.main()`, скопировать 5 строк один-в-один):
```python
smp = read_rcxml.read(fn.MEASURED)[0]
cnt = smp.counts[:-1].astype(float)
ch = np.arange(len(cnt))
e_meas = np.asarray(sum(c * ch ** i for i, c in enumerate(fn.CAL_ROOM)))[:len(cnt)]
cps_meas = cnt / smp.live
```
9. Построить `A_b`, `A_a` через `fn.fl.rebin_model_to_meas` (тот же паттерн,
   что в `fit_nuclides.main()` — `A[:, k] = fn.fl.rebin_model_to_meas(np.arange(len(c)) + 0.5, c, e_meas)`
   для каждого `k, c in enumerate(cols_b)` / `enumerate(cols_a)`).
10. Получить амплитуды из БОЕВОГО (baseline) fit — подавляя stdout:
```python
mu_idx = names_b.index("mu") if "mu" in names_b else None
with contextlib.redirect_stdout(io.StringIO()):
    amp_dict, a_mu = fn.fit_by_lines(names_b, A_b, e_meas, cps_meas, smp.live, mu_idx, pdg)
```
11. `pred_baseline = predict(names_b, A_b, amp_dict)`
    `pred_apparent = predict(names_a, A_a, amp_dict)`
    (используются ОДНИ И ТЕ ЖЕ amp_dict — сравниваем эффект замены ФОРМЫ
    шаблонов при фиксированных амплитудах, это первый порядок оценки, явно
    так и написано в печати ниже).
12. Печать блока сравнения:
```
print("")
print("=== СРАВНЕНИЕ (те же NNLS-амплитуды, разные шаблоны baseline/apparent) ===")
print("метрика                    baseline      apparent")
print("form_residual_pct %%      %10.2f    %10.2f" % (
    fc.form_residual_pct(pred_baseline, cps_meas),
    fc.form_residual_pct(pred_apparent, cps_meas)))
print("fraction_covered          %10.4f    %10.4f" % (
    fc.fraction_covered(pred_baseline, cps_meas),
    fc.fraction_covered(pred_apparent, cps_meas)))
chi2_b, ndf_b = fc.chi2_of(pred_baseline, cps_meas, smp.live, e_meas)
chi2_a, ndf_a = fc.chi2_of(pred_apparent, cps_meas, smp.live, e_meas)
print("chi2/ndf                  %10.2f    %10.2f" % (chi2_b / ndf_b, chi2_a / ndf_a))
```
13. По полосам (тот же набор границ, что в `fit_nuclides.fit_by_lines`):
```python
print("")
print("=== ПО ПОЛОСАМ ===")
print("%-12s %10s %10s %10s %10s" % ("polosa,keV", "izmereno", "baseline", "apparent", "apparent/izm"))
for lo, hi in ((20, 100), (100, 300), (300, 700), (700, 1500),
               (1500, 2000), (2000, 2400), (2400, 2830)):
    m = (e_meas >= lo) & (e_meas < hi)
    ym = cps_meas[m].sum()
    bm = pred_baseline[m].sum()
    am = pred_apparent[m].sum()
    print("%5d-%-6d %10.5f %10.5f %10.5f %10.3f" % (lo, hi, ym, bm, am, am / ym if ym else float("nan")))
```
14. Блок «требует толкования» (§31.A #SA-4), печатать ВСЕГДА, буквально:
```python
print("")
print("=== ТРЕБУЕТ ТОЛКОВАНИЯ (§31.A #SA-4) ===")
print("1. LCE_ref = глобальное взвешенное среднее по ВСЕМ событиям всех")
print("   нуклидов — не единственно возможный выбор (альтернатива: мода,")
print("   максимум LCE как 'истинная' калибровочная точка). Смена референса")
print("   сдвигает АБСОЛЮТНУЮ шкалу apparent-спектра, форма размазывания")
print("   между Y-точками не меняется.")
print("2. Одни и те же NNLS-амплитуды применены к разным шаблонам (baseline")
print("   и apparent) — это НЕ полный re-fit с apparent-шаблонами. Если")
print("   разница метрик существенна, следующий шаг — пересчитать NNLS на")
print("   apparent-шаблонах напрямую (может дать другие амплитуды).")
print("3. rcspec.fold() применяет ПАСПОРТНОЕ FWHM(E), откалиброванное на")
print("   реальном приборе — оно МОЖЕТ уже частично включать вклад LCE-")
print("   дисперсии в ширину пиков. Двойного учёта тут стараемся избежать")
print("   тем, что LCE-kernel применяется к RAW (до fold) энерговыделению,")
print("   а не поверх уже свёрнутого спектра — но это не строгое")
print("   доказательство отсутствия пересечения эффектов.")
print("4. Stage-2 статистика (N=1e6 на нуклид) даёт шумное P(y|E) на редких")
print("   каналах — fallback на глобальное P(y) нуклида (порог 20 событий)")
print("   сглаживает это, но остаётся источником неопределённости.")
```
15. Сохранить `pred_baseline`, `pred_apparent`, `e_meas` в CSV для дальнейшего
    разбора: `D:/Claude_files/repos/geant4-detector-models/detectors/RadiaCode-103/analysis/results/apparent_vs_baseline_20260823.csv`
    (создать директорию `results/` если не существует — `os.makedirs(..., exist_ok=True)`),
    колонки `E_keV,measured,baseline,apparent`, кодировка utf-8.

## Требования к приёмке

1. Скрипт запускается `python apparent_method1_compare.py` без аргументов,
   печатает всё описанное выше в stdout (кроме подавленного вывода
   `fit_by_lines`), завершается кодом 0.
2. **Мутационная проверка (#SA-3):** после первого успешного прогона —
   намеренно испортить `build_apparent_hist` (например, зафиксировать
   `lce_ref` внутри функции на константу, равную LCE каждой точки — то есть
   заставить `e_app == e_center` всегда) и прогнать снова: результат ДОЛЖЕН
   дать `pred_apparent` ПОБАЙТНО (или численно, `np.allclose`) равным
   `pred_baseline` (никакого размазывания). Если совпадение НЕ идеальное —
   значит основной прогон тоже не размазывал по-настоящему, разбирать ДО
   доклада результата. Откатить порчу после проверки.
3. Числа `form_residual_pct`/`fraction_covered`/`chi2/ndf` для baseline
   ДОЛЖНЫ совпадать (в пределах округления, `%.2f`) с тем, что печатает
   штатный `fit_nuclides.py` при обычном запуске (`SVERKA PO POLOSAM`,
   секция после `AMPLITUDY PO LINIYAM`) — если баланс не сходится, значит
   где-то разошлась нормировка/матрица A_b, искать до доклада числа.

## Не делать

- Не менять `fit_nuclides.py`/`fit_coverage.py`/`rcspec.py` — только читать
  (импортировать) их функции.
- Не писать в `bg_bare_field_m1_*.csv` (боевые файлы) — только читать.
- Не пытаться заново считать NNLS на apparent-шаблонах в этой версии скрипта
  (это явный «следующий шаг», не делать сейчас, см. п.14.2 печати).
