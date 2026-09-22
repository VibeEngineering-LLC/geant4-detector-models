Write ONE Python 3 script `lines_report.py` (numpy allowed). It imports `parse_lines` from the module `bline_table` in the same directory (`sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))`; `parse_lines(path)` returns dicts `{"E_keV", "text", "reaction", "yield", "place", "cls"}`). Output ONLY the code (no markdown fences, no prose). Comments in Russian, short. Imports: `import sys, os, math`, `import numpy as np`, `sys.stdout.reconfigure(encoding="utf-8")`.

## Purpose
Find every expected gamma line in a simulated instrument spectrum, measure its net area, and attribute it to physical origins using a per-category breakdown.

## Inputs
`lines_report.py <smeared.csv> <cat.csv> <table.md> <out_prefix>` and `lines_report.py --selftest`.
* `smeared.csv`: header `E_keV,counts_per_h_per_keV,sigma_counts_per_h_per_keV`, rows on a 1 keV grid, `E_keV` = bin centre (e.g. 20.5, 21.5, ...): a spectrum in counts per hour per keV with its statistical uncertainty. Load with `np.loadtxt(delimiter=",", skiprows=1)`.
* `cat.csv`: header `pg,proc,vol,pgname,procname,volname,rate_total,bin_0,bin_1,...,bin_N`, one row per category with UNSMEARED rates in s^-1 per 1 keV bin (bin i = energies [i, i+1) keV; the last bin is the overflow). Parse with plain `open` (the names are strings; the numeric bins are floats).
* `table.md`: the report with the line table, given to `parse_lines`.

## Functions (exactly these names)
* `fwhm_keV(E) = 41.6 * math.sqrt(E / 661.657)`.
* `net_area(E, cph, sig, E0, fwhm) -> (net, sigma, signif)`: arrays `E`, `cph` (counts/h per keV), `sig` (their uncertainties). Signal window = bins with `|E - E0| <= 1.2*fwhm`; two baseline windows `E0 - 4*fwhm <= E <= E0 - 2.2*fwhm` and `E0 + 2.2*fwhm <= E <= E0 + 4*fwhm` (a window must contain at least 2 bins, else return `(0.0, inf, 0.0)`); `base` = mean of `cph` over both baseline windows together (per keV), `nwin` = number of bins in the signal window (each bin is 1 keV wide); `net = cph[window].sum() - base * nwin`; `sigma = sqrt( (sig[window]**2).sum() + nwin**2 * (sig[baseline]**2).sum() / nb**2 )` with `nb` = number of baseline bins; `signif = net / sigma` (0 if sigma is 0 or inf).
* `attribute(cat, E0, fwhm, top=3) -> str`: `cat` = list of `{"names": (pgname, procname, volname), "bins": np.ndarray (rates per s per keV)}`; for every category compute `w = sum of bins j with |j + 0.5 - E0| <= 1.2*fwhm` and `b` = mean of the bins in the two baseline windows (same definition as above, by bin centre `j + 0.5`) times the number of bins in the signal window; `excess = max(w - b, 0)`; return a string with the `top` largest excesses as `"pg/proc/vol NN%"` (share of the total excess of all categories, one decimal), separated by `"; "`; empty string if the total excess is 0.
* `find_unlisted(E, cph, sig, listed_E, minsig=4.0)`: local maxima of `cph` (bin larger than all bins within `+-fwhm(E0)/2 keV`), at least `20` keV above the start of the array and `20` keV below its end, for which `net_area` gives `signif >= minsig` and `net > 0` and whose distance to the nearest energy in `listed_E` is more than `1.0 * fwhm` — return a list of `(E0, net, sigma, signif)`, merging maxima closer than `fwhm` (keep the most significant).
* `main()`.

## main behaviour
Load the spectrum; for every listed line with `30 <= E_keV <= 9900` compute `net_area`; write `<out_prefix>_lines.csv` with header `E_keV,reaction,place_class,net_cph,sigma_cph,signif,found,origin` (`found = 1` if `signif >= 3 and net > 0`; text fields must not contain commas: replace `,` by `;`); write `<out_prefix>_unlisted.csv` with `E_keV,net_cph,sigma_cph,signif,origin` from `find_unlisted` (origin from `attribute`). Print: number of listed lines in range, number found, number found with class `A` or `B` (the last field of `text` after the final `|`), number of unlisted peaks. Exit code 0.

## Self-test (`--selftest`), no files needed
1. Synthetic spectrum on a 1 keV grid from 20.5 to 3000.5: flat `100` counts/h/keV plus a Gaussian at `E0 = 1000` with area `5000` counts/h and `sigma = fwhm_keV(1000)/2.35482` (use the exact Gaussian bin integral or the density at the bin centre — the density is fine); uncertainties `sig = sqrt(cph)`. `net_area` must return `net` within 3 % of 5000 and `signif > 10`. The same call at `E0 = 2000` (no line) must give `abs(signif) < 3`.
2. `attribute` with two categories: A `("gamma","nCapture","Skin")` with bins array of length 3001 with value `0.5` at bins 998..1002 and zero elsewhere; B `("neutron","hadElastic","Pax")` with a constant `0.01` in all bins: the string must start with `gamma/nCapture/Skin` and A's share must be above 90 %.
3. `find_unlisted` on the spectrum of item 1 with `listed_E = [500.0]` must return exactly one entry near 1000 (|E0 - 1000| < 5) and with `listed_E = [1000.0]` an empty list.
Print `SELFTEST PASS` and return 0, else `SELFTEST FAIL: <what>` and return 1.
