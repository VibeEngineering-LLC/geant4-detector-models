Write ONE complete Python 3 script `capture_build.py`. Output ONLY the code (no markdown fences, no prose). Standard library only, plus the project modules below. Comments in Russian, short. Keep every expression simple and on one line. This script is ONLY GLUE: all physics lives in the imported, already tested modules. Do not re-implement anything they provide.

## Modules you MUST import (already written and tested)
```python
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # BEFORE the imports below
from capture_levels import read_levels, normalise
from capture_data import read_masses, sn_keV, read_targets, read_prompt, merge_lines, product_label
from capture_classify import classify_lines
from capture_model import reachable_levels, propagate_yields
```
* `read_levels(path)` -> `{idx: {"E": float, "trans": [(daughter_idx, Eg, weight, alpha), ...]}}`; `normalise(levels)` -> same with `prob` instead of weight. May raise `OSError`/`ValueError`.
* `read_masses(path)` -> `{(Z, A): mass_excess_keV}`; `sn_keV(masses, Z, A_target)` -> keV or `None`.
* `read_targets(path)` -> `{symbol: [(Z, A, abundance_percent, sigma0_barn), ...]}`.
* `read_prompt(path)` -> `{product_label: [(E_keV, sigma_gamma_barn), ...]}`; `product_label(symbol, A_target)` -> `f"{A_target+1}-{symbol}"`.
* `merge_lines(lines, tol=0.3)` -> merged `[(E, w)]`.
* `classify_lines(levels, Sn, M, lines, tol_a, tol_b)` -> list of records `(E, y, is_primary, level_idx_or_None)` in the order of `lines` (`lines` = list of `(E, y)`, `y` = yield per capture).
* `reachable_levels(levels, primary)` -> set of level indices; `propagate_yields(levels, primary, Sn, M)` -> merged `[(E, yield)]` model gamma list; `primary` = dict `level_idx -> weight`.

## Command line
    python capture_build.py <data_dir> <levels_dir> <elements> <out_db> <out_report> [tol_a=1.5] [tol_b=3.0e-4]
`data_dir` holds `promptgammas.xls`, `isotope.xls`, `mass.mas20.txt`; `levels_dir` holds `z<Z>.a<A>` files; `elements` comma-separated symbols. Wrong argument count -> usage to stderr, exit 2.

## Algorithm (in `main`, all inputs read ONCE)
1. `masses = read_masses(...)`, `targets = read_targets(...)`, `prompt = read_prompt(...)`.
2. For every requested symbol (skip and warn on stderr if not in `targets`), for every `(Z, A, abundance, sigma0)` of `targets[symbol]` call `process_isotope(symbol, Z, A, sigma0, masses, prompt, levels_dir, tol_a, tol_b)` EXACTLY ONCE. It returns either `{"skip": reason_string, "label": f"{Z}-{A}"}` or a result dict. Append EVERY return value (skipped ones too) to one list `results`.
3. `process_isotope` (all its state local; every parameter passed explicitly, no globals):
   * `Sn = sn_keV(masses, Z, A)`; if `None` or `<= 0` -> skip reason `"no Sn"`.
   * level file path `os.path.join(levels_dir, f"z{Z}.a{A+1}")`; `levels = normalise(read_levels(path))` inside `try/except (OSError, ValueError)` -> skip reason `"no level scheme"`. Do NOT use a bare `except Exception`.
   * `label = product_label(symbol, A)`; if `label not in prompt` -> skip reason `"no PGAA lines"`. If `sigma0` is `None` or `<= 0` -> `"no sigma0"`.
   * `lines = merge_lines([(E, sg / sigma0) for E, sg in prompt[label]])`; `M = (A + 1) * 931494.0`.
   * `records = classify_lines(levels, Sn, M, lines, tol_a, tol_b)`.
   * `primary = {}`: for `(E, y, is_p, j)` in records with `is_p` true: `primary[j] = primary.get(j, 0.0) + y`. `Wp = sum(primary.values())`. If `Wp > 1.0`: divide every weight by `Wp`, set `Wp = 1.0`, remember `scaled = True`.
   * If `primary` is empty -> skip reason `"no primary lines"`.
   * `reachable = reachable_levels(levels, primary)`; `model = propagate_yields(levels, primary, Sn, M)`.
   * Self-check: for every record with `is_p` false and `y >= 0.005`: find the model line closest in energy; `tol = tol_a + tol_b*E`; `n_check += 1`; if a model line exists within `tol` and `0.5 <= y_model / y <= 2.0` then `n_within2 += 1`. Also, for records with `y >= 0.02` that are not primary and whose closest model yield within `tol` is missing or below `y/2`, append `(E, y, y_model_or_0.0)` to `unexplained`.
   * Return a dict with keys `label, Z, A, symbol, Sn, sigma0, Wp, scaled, n_lines, primary, reachable, levels, n_check, n_within2, unexplained`.
4. After the loop: `written = [r for r in results if "skip" not in r]`, `skipped = [r for r in results if "skip" in r]`. Print one stdout line per written isotope: `ISO <Z>-<A> Sn=<%.3f> sigma0=<%.4g> Wp=<%.4f> primaries=<len(primary)> levels=<len(reachable)> check=<n_check> within2x=<n_within2>`. Then write the DB and the report and print the LAST stdout line `DONE isotopes=<len(written)> skipped=<len(skipped)>` (f-string with the real numbers). If `len(written) == 0` print `ERROR no isotope processed` to stderr and exit with code 3.

## Output database (UTF-8, `newline="\n"`, written ONCE at the end)
First line `# capture_db v1 targets=<len(written)> elements=<elements string>`. For every written isotope sorted by `(Z, A)`:
```
ISO <Z> <A> <Sn %.4f> <sigma0 %.6g> <Wp %.6f> <len(reachable)> <len(primary)>
PRI <level_idx> <weight %.6f>          (ascending level index)
LEV <idx> <E %.4f> <n_transitions>      (reachable levels only, ascending idx)
TR <daughter_idx> <Eg %.4f> <prob %.6f> <alpha %.6g>
END
```
LITERAL format: bare numbers separated by single spaces, NO `key=value` tokens.

## Output report (Markdown)
Title; a table `isotope | Sn (keV) | sigma0 (b) | PGAA lines | primaries | Wp | fallback 1-Wp | check | within2x` (`PGAA lines` = `n_lines`); a line `scaled` warning for isotopes with `scaled` true; a section `Unexplained lines` with, per isotope, a list `E, y, model yield`; a section `Skipped isotopes` listing every skipped `label` with its reason.

## Requirements
* `sys.stdout.reconfigure(encoding="utf-8")`, `sys.stderr.reconfigure(encoding="utf-8")` at the top; guard `main` with `if __name__ == "__main__":`.
* Functions: `process_isotope`, `write_db`, `write_report`, `main`. `write_db` and `write_report` take the `results` list (and the elements string) and never call `process_isotope`.
* No element list, no re-reading of files inside loops, no `except Exception`.
