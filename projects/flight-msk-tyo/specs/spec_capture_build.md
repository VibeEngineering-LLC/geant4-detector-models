Write ONE complete Python 3 script `capture_build.py`. Output ONLY the code (no markdown fences, no prose). Standard library only, plus the two project modules imported as below. Comments in Russian, short. Keep every expression simple and on one line (earlier attempts failed on unbalanced parentheses, on quantities reused between loop iterations, and on confusing a line index with a level index).

## Purpose
Build a database of neutron-capture gamma cascades. For a target isotope T = (Z, A) the capture product is P = (Z, A+1) with excitation energy `Sn` (neutron separation energy). PRIMARY transitions (from the capture state down to a low level) are taken from the measured IAEA PGAA line list; everything BELOW the primary level is taken from the ENSDF level scheme (branching, internal conversion). The fraction of captures not covered by measured primary lines is left to the standard Geant4 model.

## Modules you MUST import (already written and tested; do not rewrite them)
```python
from capture_levels import read_levels, normalise
#  read_levels(path) -> {idx: {"E": float, "trans": [(daughter_idx, Eg, weight, alpha), ...]}}
#  normalise(levels) -> same structure, with weight replaced by prob (sums to 1 within a level, or all 0.0)
from capture_data import read_masses, sn_keV, read_targets, read_prompt, merge_lines, product_label
#  read_masses(path) -> {(Z, A): mass excess keV};  sn_keV(masses, Z, A_target) -> keV or None
#  read_targets(path) -> {symbol: [(Z, A, abundance_percent, sigma0_barn), ...]}
#  read_prompt(path) -> {product_label: [(E_keV, sigma_gamma_barn), ...]}   e.g. label '28-Al'
#  merge_lines(lines, tol=0.3) -> [(E, w)] (sorted, neighbours within tol merged; w summed)
#  product_label(symbol, A_target) -> f"{A_target+1}-{symbol}"
```
Add `sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))` before the imports so that they resolve.

## Command line
    python capture_build.py <data_dir> <levels_dir> <elements> <out_db> <out_report> [tol_a=1.5] [tol_b=3.0e-4]
`data_dir` holds `promptgammas.xls`, `isotope.xls`, `mass.mas20.txt`; `levels_dir` holds files `z<Z>.a<A>`; `elements` is a comma-separated list of symbols. `tol_a`, `tol_b` (keV, relative) define the matching tolerance `tol(E) = tol_a + tol_b*E`. Wrong argument count -> usage to stderr, exit 2.

## Algorithm — for EVERY target isotope of every requested element (build ALL structures inside the per-isotope loop; nothing may leak from one isotope to the next)
1. Skip (record a reason) if: `sn_keV(...)` is None or <= 0; the level file `z<Z>.a<A+1>` does not exist or `read_levels` raises; `product_label(symbol, A)` has no PGAA lines. Note the level file is named by the PRODUCT mass number `A+1`.
2. `levels = normalise(read_levels(path))`.
3. Yields per capture: for every PGAA line `(E, sigma_gamma)` of the product label compute `y = sigma_gamma / sigma0` where `sigma0` is that target's thermal cross section from `read_targets` (skip the isotope if `sigma0 <= 0`). Then `lines = merge_lines([(E, y), ...])` (yields add when merged).
4. Recoil mass `M = (A+1) * 931494.0` keV. For every merged line `(E, y)`: the level it would feed is `Lt = Sn - E - E*E/(2*M)`; find the level index `j` (over ALL levels of the scheme) with the smallest `|levels[j]["E"] - Lt|`; `dist` = that distance; `tol = tol_a + tol_b*E`; the line is a PRIMARY CANDIDATE iff `dist <= tol`. It is a SECONDARY CANDIDATE iff some transition of ANY level has `|Eg - E| <= tol`.
5. Decision: PRIMARY if `E >= 0.5*Sn` and it is a primary candidate; or `E < 0.5*Sn` and it is a primary candidate and NOT a secondary candidate. Otherwise NOT primary. Keep, for every primary line, the LEVEL index `j` it feeds (NOT the index of the line in the list).
6. `primary[j]` = sum of `y` of the primary lines feeding level `j`. `Wp = sum(primary.values())`. If `Wp > 1.0` scale every weight by `1/Wp` and set `Wp = 1.0` (write a warning into the report).
7. Reachable levels: breadth-first from the keys of `primary` following `daughter_idx` of transitions; include level 0 if reached.
8. MODEL YIELDS (analytic propagation, independent of the PGAA secondary lines): visit probability `v[j]` starts as `primary[j]` for primary levels (0 otherwise); process levels in DESCENDING order of level energy `levels[j]["E"]` (a level is always fed only by higher ones): for every transition `(d, Eg, prob, alpha)` of level `j`: `flow = v[j]*prob`; the gamma yield of that transition is `flow/(1+alpha)` (record `(Eg, yield)`); `v[d] += flow`. Also record the primary gammas: for every primary level `j`, energy `Sn - levels[j]["E"]` minus recoil, yield `primary[j]`. Merge model yields within 0.3 keV (`merge_lines`).
9. Self-check: for every PGAA line NOT decided primary with `y >= 0.005`, find the model line whose energy is within `tol(E)`; count `n_check` (lines for which any model line exists OR the model yield is 0) and `n_within2` (model yield / y within `[0.5, 2.0]`). Lines with no model line count as checked and NOT within.
10. Print one stdout line per isotope: `ISO <Z>-<A> Sn=<%.3f> sigma0=<%.4g> Wp=<%.4f> primaries=<len(primary)> levels=<len(reachable)> check=<n_check> within2x=<n_within2>`.

## Output database `out_db` (UTF-8 text, `\n`)
First line: `# capture_db v1 targets=<count written> elements=<elements>` (write the file at the END, after all isotopes are processed, so that the count is real). Then for every written isotope in ascending `(Z, A)`:
```
ISO <Z> <A_target> <Sn %.4f> <sigma0 %.6g> <Wp %.6f> <n_reachable> <n_primary>
PRI <level_idx> <weight %.6f>
LEV <idx> <E_level %.4f> <n_transitions>
TR <daughter_idx> <Eg %.4f> <prob %.6f> <alpha %.6g>
END
```
`PRI` lines: one per key of `primary`, ascending level index. `LEV` blocks: reachable levels only, ascending `idx`, each followed by its `n_transitions` `TR` lines (the normalised `prob`).

## Output report `out_report` (Markdown)
A title; a table `isotope | Sn | sigma0 | PGAA lines | primaries | Wp | fallback 1-Wp | check | within2x`; then for every written isotope a list of its PGAA lines with `y >= 0.02` that are NOT primary and whose model yield is below half of `y` (energy, y, model yield) — store these per isotope inside the loop; then the list of skipped isotopes with reasons. Last stdout line: `DONE isotopes=<written> skipped=<skipped>` with the REAL numbers (use an f-string).

## Pitfalls of the first attempt (ALL must be avoided)
1. Read every input file ONCE in `main` (`masses`, `targets`, `prompt`), pass the dicts to the functions. The first attempt re-read the Excel files for every isotope and called `process_isotope` three times per isotope. Call it exactly ONCE per isotope and keep the result dicts in a list `results`; `write_db` and `write_report` only format that list.
2. The database format is LITERAL: `ISO 13 27 7725.1730 0.231 0.912345 12 3` — bare numbers separated by single spaces, NO `key=value` tokens, NO `Sn=`, `weight=`, `E_level=`, `prob=`, `alpha=`, `Eg=`. The C++ reader parses `tag int int double ...`.
3. `sigma0` and `Z` come from a dict built ONCE from `read_targets`: `sigma0_of[(Z, A)] = sigma0` and `symbol_of[Z] = symbol`. Never index the list of a symbol by `A_target - 1` or by position.
4. Implement the primary decision EXACTLY as in step 5 (two cases: `E >= 0.5*Sn` needs only the primary candidate; `E < 0.5*Sn` needs primary candidate AND NOT secondary candidate). `classify_lines` must return, for every merged line, a record `(E, y, is_primary, level_index_or_None)`; the caller builds `primary[j]` from the records.
5. Primary gamma energy in the model list: `Sn - levels[j]["E"]` minus recoil `Eg*Eg/(2*M)` (as in the emitter).
6. Self-check (step 9) uses ONLY the records with `is_primary == False` and `y >= 0.005`; `n_check` = their count; `n_within2` = those with a model line within `tol(E)` whose yield ratio model/y is in `[0.5, 2.0]` (choose the CLOSEST model line in energy, not the first found).
7. The unexplained-lines list (report) uses the same records: non-primary, `y >= 0.02`, model yield (closest model line within tol, else 0) below `y/2`.
8. Print one stdout line per written isotope exactly as in step 10 (`ISO ...`), then `DONE ...`.
9. `write_db` writes the header line LAST-KNOWN count: build all text first, then write once with the header `# capture_db v1 targets=<written> elements=<elements>`.
10. Levels of an isotope are read from `levels_dir/z<Z>.a<A+1>` with `read_levels` (wrap in try/except and record the skip reason `no level file`).

## Requirements
* `sys.stdout.reconfigure(encoding="utf-8")`, `sys.stderr.reconfigure(encoding="utf-8")` at the top; write files with `encoding="utf-8", newline="\n"`.
* Functions: `classify_lines`, `propagate_yields`, `reachable_levels`, `process_isotope`, `write_db`, `write_report`, `main`; guard with `if __name__ == "__main__":`. `process_isotope` returns a dict with everything needed for the outputs, so that `main` only loops and writes.
* Never use the loop variable of one isotope after the loop; never write a list of element symbols.
