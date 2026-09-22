Write ONE Python 3 script `spec_smear.py` (numpy allowed, no other third-party packages). Output ONLY the code (no markdown fences, no prose). Comments in Russian, short. Start with `import sys, argparse`, `import numpy as np`, `sys.stdout.reconfigure(encoding="utf-8")`.

## Purpose
Convert a stage-II energy-deposit spectrum (rates per second in 1 keV bins) into an instrument spectrum with Gaussian energy resolution, in counts per hour.

## Input file (CSV, produced by a C++ program)
First line starts with `#` (comment, skip it). Second line is the header `bin_keV,sumw,sumw2,light_sumw,light_sumw2`. Then N+1 data rows: `bin_keV` = 0..N-1 = lower edge of a 1 keV bin (bin i covers [i, i+1) keV), the LAST row is an overflow bin (energy >= N keV). Columns: `sumw` = rate in s^-1 in that bin, `sumw2` = squared statistical uncertainty of the rate in s^-2 (sum of squared weights divided by T_sim^2), `light_*` = the same for the quenched light output (may be ignored unless `--light` is given).

## Command line
`spec_smear.py <total.csv> <out.csv> [--fwhm662 41.6] [--emin 20] [--emax 10000] [--light]`
`--fwhm662`: FWHM in keV at 661.657 keV; FWHM(E) = fwhm662 * sqrt(E / 661.657) (energy of a bin = its centre i + 0.5). `--light` selects the light columns instead of the deposit columns.

## Algorithm
1. Read the columns with `np.loadtxt(..., delimiter=",", skiprows=2)` (the comment line and the header are the first two lines); drop the overflow row from the smearing but remember its rate (`overflow_rate_per_h`).
2. Energy centres `E_i = i + 0.5` (keV). For each source bin i with rate r_i (only bins with r_i > 0) spread it onto the output bins j (centres E_j) with a Gaussian of sigma_i = FWHM(E_i) / 2.35482 truncated at +-5 sigma: weight_ij = `erf`-based integral over the output bin [j, j+1) (use `math.erf` vectorised through `np.vectorize` or `scipy`-free formula `0.5*(1+erf(x/sqrt2))` computed with `np.frompyfunc`; correctness matters more than speed, but the whole run must finish in under 60 s for N = 10000 — use a vectorised approach over j for each i, skipping bins with zero rate).
   The uncertainty: the variance of the smeared bin j is `sum_i weight_ij^2 * var_i` (var_i = sumw2_i); this treats source bins as independent, which is the intended approximation.
3. Output CSV header: `E_keV,counts_per_h_per_keV,sigma_counts_per_h_per_keV` for output bins with `emin <= E_j <= emax`; values = smeared rate * 3600 and sqrt(variance) * 3600.
4. Print to stdout: the total counts per hour inside [emin, emax] of the SMEARED spectrum, the same for the UNSMEARED input restricted to the same range (they must agree within 1 % unless the range cuts a wide part; print the ratio), the statistical uncertainty of the total (sqrt of the sum of variances, in counts per hour) and `overflow_rate_per_h`.
5. Exit code 0; on a malformed file print an error to stderr and exit 2.

## Requirements
* Include every import you use (`math`, `argparse`, `numpy`). Do not use pandas or scipy.
