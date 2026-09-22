Write ONE complete Python 3 module `capture_classify.py`. Output ONLY the code (no markdown fences, no prose). Standard library only. Comments in Russian, short. Keep every expression simple and on one line.

## Purpose
Given a level scheme of the capture PRODUCT nucleus and the measured gamma lines of thermal neutron capture, decide which lines are PRIMARY transitions (from the capture state at energy `Sn` down to a low level) and which level each primary line feeds.

## Input structures
* `levels`: dict `idx -> {"E": float_keV, "trans": [(daughter_idx, Eg_keV, prob, alpha), ...]}` (level 0 is the ground state, `E == 0`).
* `lines`: list of `(E_keV, y)` — gamma energy and yield per capture, already merged.
* `Sn`: neutron separation energy of the product (keV). `M = (A_product) * 931494.0` keV (mass of the product, used for recoil). `tol_a`, `tol_b`: tolerance `tol(E) = tol_a + tol_b*E` keV.

## Interface (exactly)
```python
def tol_of(E, tol_a, tol_b):
    """Допуск совпадения в кэВ."""
def nearest_level(levels, Lt):
    """(idx, dist): индекс уровня, ближайшего по энергии к Lt, и расстояние |E_level - Lt|. levels не пуст."""
def is_secondary_candidate(levels, E, tol):
    """True, если у КАКОГО-ЛИБО уровня есть переход с |Eg - E| <= tol."""
def classify_lines(levels, Sn, M, lines, tol_a, tol_b):
    """Список записей (E, y, is_primary, level_idx_or_None) в порядке lines."""
```

## Rules for `classify_lines` (for every line `(E, y)`)
1. `Lt = Sn - E - E*E/(2*M)` — energy of the level the line would feed if it were a primary transition.
2. `j, dist = nearest_level(levels, Lt)`; `tol = tol_of(E, tol_a, tol_b)`; `primary_candidate = dist <= tol`.
3. `secondary = is_secondary_candidate(levels, E, tol)`.
4. Decision: if `E >= 0.5*Sn`: `is_primary = primary_candidate`. If `E < 0.5*Sn`: `is_primary = primary_candidate and not secondary`.
5. The record is `(E, y, is_primary, j if is_primary else None)`. Lines with `y <= 0` get `(E, y, False, None)`.

## Command-line self-test
`python capture_classify.py` (no arguments) runs a built-in check on a tiny synthetic scheme and prints one line `SELFTEST PASS` (exit 0) or `SELFTEST FAIL <reason>` (exit 1):
scheme `levels = {0: {"E": 0.0, "trans": []}, 1: {"E": 300.0, "trans": [(0, 300.0, 1.0, 0.0)]}, 2: {"E": 1000.0, "trans": [(1, 700.0, 0.75, 0.0), (0, 1000.0, 0.25, 0.0)]}}`, `Sn = 8000.0`, `M = 11*931494.0`, `tol_a = 1.5`, `tol_b = 3e-4`, lines: `(6997.9, 0.5)` must be PRIMARY feeding level 2 (`Lt` is about 1000), `(7697.8, 0.2)` PRIMARY feeding level 1, `(700.0, 0.3)` NOT primary (it is a secondary candidate and below `0.5*Sn`), `(50.0, 0.1)` NOT primary. Compute the recoil-consistent energies yourself in the check by using `E = Sn - Elevel - recoil` with `recoil = E*E/(2*M)` solved by two fixed-point iterations, rather than the rounded literals above, so that the test is exact.

## Requirements
* `sys.stdout.reconfigure(encoding="utf-8")` at the top; guard the self-test with `if __name__ == "__main__":`.
* Do NOT import any other project module.
