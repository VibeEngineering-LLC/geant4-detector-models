Задача: написать ТОЛЬКО функцию `main()` и финальный блок `if __name__` для
Python-модуля `apparent_method1_compare.py`. Верни ТОЛЬКО код функции `main()`
целиком (с сигнатурой `def main():`) плюс две строки в конце:
```python
if __name__ == "__main__":
    main()
```
Без markdown fences, без пояснений до/после, чистый Python. НЕ пиши импорты,
НЕ пиши другие функции — они уже существуют в том же файле выше:

Уже существуют (сигнатуры, для справки, НЕ переопределять):
- `fn` = импортированный модуль `fit_nuclides` (алиас `import fit_nuclides as fn`)
- `fc` = импортированный модуль `fit_coverage`
- `rcspec`, `read_rcxml`, `np` (numpy), `io`, `contextlib`, `math`, `os`, `sys` — уже импортированы
- `LCE_MAP: dict[float, float]` — карта Y(мм) -> LCE(%)
- `STAGE2_DIR: str` — путь к папке со Stage-2 файлами `s2_<нуклид>_ypos.csv`
- `load_ypos(path) -> (y_cols: list[float] | None, dist: dict[int, np.ndarray] | None)`
- `global_lce_ref(all_dists: list[tuple]) -> (lce_ref: float, n_total: int)`
- `build_apparent_hist(raw_hist: np.ndarray, y_cols: list, dist: dict, lce_ref: float) -> np.ndarray`
- `load_column(nuc: str) -> (raw_norm: np.ndarray | None, meta: dict | None)`
- `predict(names: list, A: np.ndarray, amp_dict: dict) -> np.ndarray`

`fit_nuclides` (модуль `fn`) уже содержит: `fn.NUCS` (list из 8 имён нуклидов),
`fn.merge_by_chain(names, cols)`, `fn.load_muons()`, `fn.MEASURED` (путь к XML
измерения), `fn.CAL_ROOM` (list из 3 коэфф. калибровки), `fn.fl` (модуль
`fit_lines`, содержит `fn.fl.rebin_model_to_meas(x, y, x_new)`),
`fn.fit_by_lines(names, A, e_meas, cps_meas, live, mu_col_idx, pdg) ->
(amp_dict: dict[str, float], a_mu: float)`.

`fc` (модуль `fit_coverage`) содержит: `fc.form_residual_pct(model, meas) ->
float`, `fc.fraction_covered(model, meas) -> float`,
`fc.chi2_of(model, meas, live, e_meas) -> (chi2: float, ndf: int)`.

## Тело `main()`, по шагам

1. `print("=== #FIT-1 Stage-3: apparent-спектр метода 1 (LCE-tailing) ===")`.
2. Для каждого `nuc in fn.NUCS`: вызвать
   `load_ypos(os.path.join(STAGE2_DIR, "s2_%s_ypos.csv" % nuc))`. Собрать
   список троек `(nuc, y_cols, dist)` для тех, у кого `dist is not None`
   (назови список `stage2_have`).
3. Если `stage2_have` пуст: `print("НЕТ STAGE-2 ДАННЫХ, СТОП"); return`.
4. `lce_ref, n_total = global_lce_ref([(yc, d) for _, yc, d in stage2_have])`.
   `print("LCE_ref (глобальный, %d событий всех нуклидов) = %.2f %%" % (n_total, lce_ref))`.
5. Списки `names_raw = []`, `cols_baseline = []`, `cols_apparent = []`.
   Для каждого `nuc in fn.NUCS`:
   - `raw, meta = load_column(nuc)`. Если `raw is None`:
     `print("[--] %s: нет шаблона" % nuc); continue`.
   - Найти в `stage2_have` запись с таким `nuc` (например через генератор
     `next((t for t in stage2_have if t[0] == nuc), None)`).
   - Если найдена `(_, y_cols, dist)`: `apparent_raw = build_apparent_hist(raw, y_cols, dist, lce_ref)`.
     Иначе: `apparent_raw = raw.copy()`,
     `print("[--] %s: нет Stage-2 Y-данных, apparent = baseline" % nuc)`.
   - `cps_baseline = rcspec.fold(raw, "103")`
   - `cps_apparent = rcspec.fold(apparent_raw, "103")`
   - `names_raw.append(nuc)`, `cols_baseline.append(cps_baseline)`, `cols_apparent.append(cps_apparent)`.
6. `names_b, cols_b = fn.merge_by_chain(names_raw, cols_baseline)`
   `names_a, cols_a = fn.merge_by_chain(names_raw, cols_apparent)`
   Если `names_b != names_a`: `raise SystemExit("names_b != names_a - разошлись списки нуклидов, разбирать")`.
7. `mu, pdg = fn.load_muons()`. Если `mu is not None`:
   `names_b.append("mu"); cols_b.append(mu)`
   `names_a.append("mu"); cols_a.append(mu)`.
8. Измерение (скопировать ДОСЛОВНО этот код, ничего не менять):
```python
    smp = read_rcxml.read(fn.MEASURED)[0]
    cnt = smp.counts[:-1].astype(float)
    ch = np.arange(len(cnt))
    e_meas = np.asarray(sum(c * ch ** i for i, c in enumerate(fn.CAL_ROOM)))[:len(cnt)]
    cps_meas = cnt / smp.live
```
9. Построить матрицы:
```python
    A_b = np.zeros((len(e_meas), len(cols_b)))
    for k, c in enumerate(cols_b):
        A_b[:, k] = fn.fl.rebin_model_to_meas(np.arange(len(c)) + 0.5, c, e_meas)
    A_a = np.zeros((len(e_meas), len(cols_a)))
    for k, c in enumerate(cols_a):
        A_a[:, k] = fn.fl.rebin_model_to_meas(np.arange(len(c)) + 0.5, c, e_meas)
```
10. Получить амплитуды из baseline-фита, подавляя stdout:
```python
    mu_idx = names_b.index("mu") if "mu" in names_b else None
    with contextlib.redirect_stdout(io.StringIO()):
        amp_dict, a_mu = fn.fit_by_lines(names_b, A_b, e_meas, cps_meas, smp.live, mu_idx, pdg)
```
11. `pred_baseline = predict(names_b, A_b, amp_dict)`
    `pred_apparent = predict(names_a, A_a, amp_dict)`.
12. Печать сравнения метрик:
```python
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
13. Печать по полосам:
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
14. Блок «требует толкования» — печатать ДОСЛОВНО этот текст (скопировать
    как есть, это не место для творчества):
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
15. Сохранить результат в CSV:
```python
    out_dir = os.path.join(_HERE, "results")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "apparent_vs_baseline_20260823.csv")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("E_keV,measured,baseline,apparent\n")
        for i in range(len(e_meas)):
            f.write("%.2f,%.6e,%.6e,%.6e\n" % (e_meas[i], cps_meas[i], pred_baseline[i], pred_apparent[i]))
    print("")
    print("Сохранено: %s" % out_path)
```

Верни ТОЛЬКО `def main(): ...` (все 15 пунктов внутри, по порядку) и в самом
конце `if __name__ == "__main__": main()`.
