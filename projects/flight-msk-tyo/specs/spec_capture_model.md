Write ONE complete Python 3 module `capture_model.py`. Output ONLY the code (no markdown fences, no prose). Standard library only, plus `from capture_data import merge_lines`. Comments in Russian, short. Keep every expression simple and on one line.

## Purpose
Given a level scheme and the weights of the levels fed by the primary transitions, compute (a) the set of reachable levels and (b) the ANALYTIC gamma yield per capture of every transition of the cascade (used both as an independent self-check against measured secondary lines and for the report).

## Input structures
* `levels`: dict `idx -> {"E": float_keV, "trans": [(daughter_idx, Eg_keV, prob, alpha), ...]}`; `prob` sums to 1 within a level (or all are 0.0 for a level without transitions). Level 0 is the ground state.
* `primary`: dict `level_idx -> weight` (probability per capture that the primary transition feeds that level).
* `Sn`, `M`: neutron separation energy of the product (keV) and its mass in keV.

## Interface (exactly)
```python
def reachable_levels(levels, primary):
    """Множество индексов уровней, достижимых из ключей primary по переходам (включая сами ключи). Индексы дочерних уровней, которых нет в levels, пропускаются."""
def propagate_yields(levels, primary, Sn, M):
    """Список (E_keV, yield_per_capture) всех гамма модели: первичных и каскадных, склеенный merge_lines (допуск 0.3 кэВ)."""
```

## Rules
* `reachable_levels`: breadth-first search (use `collections.deque`), never raises for a missing daughter index.
* `propagate_yields`:
  * `v[j]` = probability that level `j` is visited: initialise `v = {j: 0.0 for j in levels}`; then `v[j] += w` for every `(j, w)` in `primary` (ignore keys not in `levels`).
  * Process ONLY the reachable levels, in DESCENDING order of `levels[j]["E"]` (a level is fed only by higher ones). For every transition `(d, Eg, prob, alpha)` of level `j`: `flow = v[j] * prob`; the gamma yield of this transition is `flow / (1.0 + alpha)`, appended as `(Eg, flow/(1+alpha))`; if `d` is in `v` then `v[d] += flow`.
  * Primary gammas: for every `(j, w)` in `primary` with `j` in `levels`: `E0 = Sn - levels[j]["E"]`; subtract recoil once: `E0 = E0 - E0*E0/(2.0*M)`; append `(E0, w)` (skip if `E0 <= 0`).
  * Return `merge_lines(list_of_pairs)`.
  * Do not modify the inputs.

## Command-line self-test
`python capture_model.py` runs a built-in check and prints `SELFTEST PASS` (exit 0) or `SELFTEST FAIL <reason>` (exit 1). Scheme: `levels = {0: {"E": 0.0, "trans": []}, 1: {"E": 300.0, "trans": [(0, 300.0, 1.0, 1.0)]}, 2: {"E": 1000.0, "trans": [(1, 700.0, 0.75, 0.0), (0, 1000.0, 0.25, 0.0)]}}`, `primary = {2: 0.6, 1: 0.2}`, `Sn = 8000.0`, `M = 11*931494.0`. Expected yields per capture: line 700 keV: `0.45`; line 1000 keV: `0.15`; line 300 keV: `(0.2 + 0.45)/(1+1.0) = 0.325`; primary gammas about 6997.9 keV with `0.6` and about 7697.5 keV with `0.2`. Check each with tolerance `1e-9` on yields and `2.0` keV on the two primary energies (find the model line closest in energy). Also check that `reachable_levels(levels, {2: 0.6})` equals `{0, 1, 2}`.

## Requirements
* `sys.stdout.reconfigure(encoding="utf-8")` at the top; guard the self-test with `if __name__ == "__main__":`. Add `sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))` before the import of `capture_data`.
