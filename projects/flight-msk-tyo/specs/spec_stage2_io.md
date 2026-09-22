Write ONE Python 3 module `stage2_io.py` (numpy allowed). Output ONLY the code (no markdown fences, no prose). Comments in Russian, short. Imports: `import numpy as np`.

## Purpose
Readers and writers for the CSV files of a spectrum-simulation pipeline. No command line, no main.

## Functions (exactly these names and signatures)
1. `read_meta(path) -> dict`: `key=value` per line -> dict of strings (skip lines without `=`). Raises the normal OSError if the file is missing.
2. `read_total(path) -> np.ndarray`: file format: line 1 = comment starting with `#`, line 2 = header `bin_keV,sumw,sumw2,light_sumw,light_sumw2`, then rows `index,sumw,sumw2,light_sumw,light_sumw2`. Return an array of shape (rows, 4) with the columns sumw, sumw2, light_sumw, light_sumw2 (i.e. WITHOUT the index column). Use `np.loadtxt(path, delimiter=",", skiprows=2, ndmin=2)` and take `[:, 1:5]`.
3. `write_total(path, arr, note)`: writes line 1 = `# ` + note, line 2 = header above, then one row per array row: index `i`, then the four columns, all with `"%.10g"`; `encoding="utf-8", newline="\n"`.
4. `read_cat(path) -> dict`: category file format: FIRST line is a header `pg,proc,vol,pgname,procname,volname,rate_total,bin_0,...,bin_N` (NO comment line); every other non-empty line: three integers, three names, `rate_total`, then the bin values. Return `{(pg, proc, vol): {"names": (pgname, procname, volname), "bins": np.ndarray of the floats from field 7 to the end}}`. Duplicate keys inside one file are added.
5. `write_cat(path, cats)`: `cats` has the same structure as the return of `read_cat`; write the header `pg,proc,vol,pgname,procname,volname,rate_total,bin_0,...,bin_{n-1}` where n = length of the bin arrays, then the categories sorted by descending `sum(bins[:-1])` (the last bin is the overflow and is NOT part of `rate_total`); `rate_total` = `sum(bins[:-1])`; numbers `"%.10g"`; `encoding="utf-8", newline="\n"`.
