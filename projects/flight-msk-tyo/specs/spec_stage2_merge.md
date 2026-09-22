Write ONE Python 3 script `stage2_merge.py` (numpy allowed, no other third-party packages). Output ONLY the code (no markdown fences, no prose). Comments in Russian, short. Start with `import sys, os, glob, argparse, tempfile`, `import numpy as np`, `sys.stdout.reconfigure(encoding="utf-8")`.

## Purpose
Merge the outputs of many stage-II processes (one per stage-I file) into per-component and total spectra.

## Input files in a directory `<dir>` (all written by a C++ program)
For every stage-I file `<base>` (e.g. `neutron_0`, `gamma_2`, `mum_1`; the component is `base.rsplit("_", 1)[0]`) there are:
* `<base>_meta.txt`: `key=value` lines; needed keys: `T_sim_s` (float) and, if present, `triggers`, `used`, `draws`.
* `<base>_total.csv`: first line `# ...` (comment), second line header `bin_keV,sumw,sumw2,light_sumw,light_sumw2`, then N+1 rows (bin index 0..N-1 and a last overflow row). Values are RATES per second (`sumw`) and squared uncertainties of rates (`sumw2`, s^-2) for THIS file's own simulated time T_f.
* `<base>_cat.csv`: header line `pg,proc,vol,pgname,procname,volname,rate_total,bin_0,...,bin_N`, then one row per category with rates (s^-1) per bin.

## Merging rules (exact)
For one component with files f = 1..F and times T_f (T = sum T_f): `rate = sum_f rate_f * T_f / T`; `var = sum_f var_f * T_f^2 / T^2`. The same rule applies to the light columns and to every bin of the category rows (categories are matched by the triple (pg, proc, vol); a category missing in a file counts as zero there; `rate_total` merged by the same rule). The TOTAL over components = sum of the component rates, sum of the component variances, and (for categories) the union of categories with rates summed over components that have them. Statistical independence of files and components is assumed.

## Command line
`stage2_merge.py <dir> <out_prefix>` and `stage2_merge.py --selftest`.
Outputs: `<out_prefix>_total_<comp>.csv` per component, `<out_prefix>_total_all.csv` (same format as the input `_total.csv`: first line `# rates per second; ...`, header, N+1 rows, `%.10g`), and `<out_prefix>_cat_all.csv` (same format as the input cat files, categories sorted by descending `rate_total`).
Print a table to stdout: component, number of files, T_sim total (s), total rate in bins 0..N-1 (s^-1), its statistical uncertainty (sqrt of the sum of the bin variances), and the same for the TOTAL. Files with a missing meta or total file -> message on stderr, skip them; if no file was merged exit code 2.

## Self-test (`--selftest`)
Create a temporary directory with two components, with synthetic files written by your own helper in the exact format above with N = 5 bins (so 6 rows): component `aa` with files `aa_0` (T=1, rates 2,0,0,0,0 and overflow 0, sumw2 = 4,0,...) and `aa_1` (T=3, rates 6,0,0,0,0, sumw2 = 12 in bin 0); component `bb` with one file `bb_0` (T=2, rates 0,1,0,0,0, sumw2 0.5 in bin 1). Expected merged `aa` bin 0: rate = (2*1 + 6*3)/4 = 5, var = (4*1 + 12*9)/16 = 7; `bb` bin 1: rate 1, var 0.5; total: bin 0 = 5 and bin 1 = 1. Also one category row `0,1,2,gamma,nCapture,Cargo,<rate_total>,<bins>` in `aa_0_cat.csv` with bin 0 rate 2 and in `aa_1_cat.csv` with bin 0 rate 6 (same triple) -> merged category bin 0 = 5. Check all of these with `abs(diff) < 1e-9`, print `SELFTEST PASS` / `SELFTEST FAIL` and return exit code 0/1. The self-test must run the same merge function that the normal mode uses (use a function `merge_dir(directory, out_prefix, nbin_expected=None)`).

## Structure (use exactly these small functions; every array has the FULL length of the file, i.e. N+1 entries including the overflow row / overflow bin)
* `read_meta(path) -> dict[str,str]`.
* `read_total(path) -> np.ndarray` of shape (N+1, 5) (`np.loadtxt(..., delimiter=",", skiprows=2)`; columns: bin index, sumw, sumw2, light_sumw, light_sumw2).
* `read_cat(path) -> dict[(int,int,int) -> (list_of_3_names, np.ndarray of length N+1)]`: the FIRST LINE of a cat file is the header (there is NO comment line in cat files); every further non-empty line is `pg,proc,vol,pgname,procname,volname,rate_total,bin_0..bin_N`: key = ints of the first three fields, names = fields 3..5, array = floats of fields 7.. (the rate_total field 6 is ignored on reading and recomputed as the sum of the bins 0..N-1 when writing).
* `merge_dir(directory, out_prefix, nbin_expected=None)`: group bases by component; for each component accumulate `sum_f rate_f*T_f` (columns sumw and light_sumw), `sum_f var_f*T_f**2` (sumw2, light_sumw2) and for every category `sum_f rate*T_f` (keep the three names from the first occurrence), then divide by `T` (rates) and `T**2` (variances); the category rates of a component are divided by that component's T. The total over components: sum of the component arrays. NEVER shadow a dictionary with a scalar variable of the same name inside a loop. If `nbin_expected` is given, a total file with a different number of rows than `nbin_expected + 1` must be skipped with a message on stderr.
* Output rows of `_total_*.csv`: bin index i = 0..N-1 as `bin_keV`, the overflow row gets the index N; `_cat_all.csv` header exactly `pg,proc,vol,pgname,procname,volname,rate_total,bin_0,...,bin_N` built from the actual length, names taken from the stored names (never hard-coded).
* The self-test writes its synthetic cat files WITHOUT a comment line, with 6 bins (bin_0..bin_5; N = 5), and additionally checks the merged category row: bin 0 must be 5.0 (read `test_cat_all.csv`, first data line, 8th field) and `rate_total` = 5.0. It must not create any subdirectories.

## Requirements
* Include every import you use. Read files with `np.loadtxt` for total files (`delimiter=","`, `skiprows=2`) and `open` for cat files. Write with `encoding="utf-8"`, `newline="\n"`. Do not use pandas.
