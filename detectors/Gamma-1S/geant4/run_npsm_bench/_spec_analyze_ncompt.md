Write a single self-contained Python 3 script. Output ONLY the code, no prose, no
markdown fences. Standard library only (csv, math, glob, argparse, pathlib, sys).

# Purpose

Collect the results of a series of Geant4 benchmark runs (one CSV per photon
energy), tabulate the mean number of Compton scatterings before absorption as a
function of energy, and check it against a published reference with
pre-declared tolerances. Tier-0 acceptance: the script prints PASS/FAIL per
criterion and returns a nonzero exit code on any FAIL.

# Input files

Each CSV was written by the benchmark and has this layout (key,value lines,
then a table):

```
# npsm_bench geometry/transport benchmark
energy_keV,500
n_events_requested,2000000
n_events_processed,2000000
seed,11
em_cut_mm,0.05
em_deex,max
lowest_electron_energy_keV,...
crystal_mm,102x102x406
crystal_volume_cm3,4224.7
r_src_mm,245.4
n_absorbed_phot,...
n_full_edep_1keV,...
n_conv,...
n_escaped,...
n_other,...
mean_ncompt_absorbed,1.4123
sem_ncompt_absorbed,0.0012
n_compt,count_absorbed,count_all
0,...,...
1,...,...
```

Parse the key,value section into a dict (strings; convert numerically where
needed). Lines starting with `#` are comments. The table after the
`n_compt,count_absorbed,count_all` header may be ignored except for one use
below.

# CLI

```
--glob PATTERN     (default "out/ncompt_*.csv")  files to collect
--ref-lo E MEAN TOL   (default 500 1.4 0.15)   reference at the low energy
--ref-hi E MEAN TOL   (default 2000 2.3 0.20)  reference at the high energy
--json PATH        optional machine-readable report
```

# Behaviour

1. Load every matching file. Skip a file with a warning if it lacks
   `energy_keV` or `mean_ncompt_absorbed`. Sort by energy.
2. Sanity per file, printed as a table row:
   energy, n_events_processed, n_absorbed_phot, n_full_edep_1keV,
   n_escaped, mean ± sem, and the ratio `n_absorbed_phot / n_full_edep_1keV`
   (print `nan` if the divisor is 0). Also compute, from the histogram table,
   the mean of `n_compt` weighted by `count_absorbed` and print it next to the
   file's own `mean_ncompt_absorbed` — the two must agree to 1e-3; if not,
   print a warning line `ПРОТИВОРЕЧИЕ В ФАЙЛЕ` for that file.
3. Criteria (each printed as `PASS` or `FAIL` with the numbers):
   - **K1 low energy:** the file whose energy is closest to `--ref-lo` E must
     have |mean − MEAN| ≤ TOL. Print measured, reference, tolerance, difference.
   - **K2 high energy:** same for `--ref-hi`.
   - **K3 monotonic:** means sorted by energy must be non-decreasing (allow a
     decrease no larger than 2·sem of the neighbouring points, since statistics
     is finite). Print the first violating pair if any.
   - **K4 process/energy consistency:** for every file,
     `n_absorbed_phot / n_full_edep_1keV` must lie in [0.90, 1.10]; print the
     worst file. Explain in a comment: the two are independent definitions of
     "fully absorbed" (photoelectric end-of-track vs. deposited energy equal to
     the primary within 1 keV); a large mismatch means the counter is keyed on
     the wrong thing.
4. Final line: `ИТОГ: PASS` if all four pass, else `ИТОГ: FAIL (<list of failed>)`.
   Exit code 0 on PASS, 1 on FAIL, 2 if fewer than 2 files were loaded.
5. `--json`: write {"files":[...per-file dicts...], "criteria":{...}, "pass":bool}
   with `ensure_ascii=False, indent=1`.

# Style

`sys.stdout.reconfigure(encoding="utf-8")` right after imports; all file I/O
with `encoding="utf-8"`. Module docstring and comments in Russian. Table printed
with fixed-width columns, numbers with `format`, no external libraries.
