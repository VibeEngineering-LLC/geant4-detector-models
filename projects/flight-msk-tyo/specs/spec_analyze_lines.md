Write ONE complete Python 3 script `analyze_lines.py`. Output ONLY the code (no markdown fences, no prose). Allowed imports: `sys`, `math`, `pandas` (with `xlrd` installed for `.xls`), no other third-party modules.

## Purpose
Compare the capture-gamma histogram produced by the Geant4 application `line_bench` with the IAEA PGAA prompt-gamma database, line by line, for a thermal-neutron run on a target made of given chemical element(s).

## Command line
    python analyze_lines.py <bench_csv> <elements> <data_dir> [top_n]
* `elements`: comma-separated element symbols of the target that produce the gammas to check, e.g. `Al` or `H,C`.
* `data_dir`: directory holding `promptgammas.xls` and `isotope.xls`.
* `top_n` optional integer, default 12: number of strongest expected lines to check per element.
* Wrong arguments -> usage to stderr, exit code 2.

## Inputs
**bench_csv** (text, UTF-8): line 1 `# material=... E_MeV=... N=... hp_mode=... seed=... radius_cm=... n_capture_gamma=... n_inelastic_gamma=... out_of_range=...`; line 2 `capture,<12000 integers>`; line 3 `inelastic,<12000 integers>`. Column j (0-based) of the histogram is the 1-keV bin `[j, j+1)` keV, centre `j+0.5` keV. Use only the `capture` line. Parse the header `key=value` tokens after the `#`.

**promptgammas.xls**: read with `pd.read_excel(path, header=None, skiprows=2)`. Columns by POSITION: `1` = product isotope label like `2-H` or `28-Al` (mass number of the PRODUCT nucleus, i.e. target mass number + 1; format `<A>-<Symbol>`), `4`... ignore; `5` = gamma energy in keV (numeric, "Value"), `7` = partial gamma cross section `sigma_gamma` in barn (numeric, "Value"). Drop rows where the energy or the cross section is NaN or non-numeric (use `pd.to_numeric(..., errors="coerce")`).

**isotope.xls**: `pd.read_excel(path, header=None, skiprows=2)`. Columns by POSITION: `0` element symbol (may be NaN on continuation rows - forward fill), `1` isotope label like `27-Al` (target isotope), `5` = natural abundance in percent (numeric, "Value").

## Expected lines (per element symbol `S`)
For every prompt-gamma row of an element with symbol `S` (parse the label: product mass number `Ap`, symbol; target label = `f"{Ap-1}-{S}"`), find the abundance `ab` (percent) of that target isotope in `isotope.xls` (match the isotope label string exactly; if not found, skip the row). Element-weighted line strength `w = ab/100 * sigma_gamma`. Merge rows of the same element whose energies differ by less than 0.3 keV by summing `w` and using the strength-weighted mean energy. Keep lines with `30 <= E <= 11500` keV. Sort by `w` descending and take the `top_n` strongest. `ref_w` = the largest `w` of the kept lines; `exp_rel = w / ref_w`.

## Observed line areas
For each expected line at energy `E`:
* half window `hw = 3 + 0.0006*E` keV (so ~3 keV at low energy, ~9 keV at 10 MeV);
* peak bins: all bins with centre in `[E-hw, E+hw]`; `peak_sum` = sum of counts there;
* baseline from two side bands, each 8 keV wide: `[E-hw-8, E-hw)` and `(E+hw, E+hw+8]` (bins whose centres are inside); `base_per_bin` = mean count per bin over both side bands (0 if empty); `area = peak_sum - base_per_bin * (number of peak bins)`;
* `sigma_area = sqrt(peak_sum + (base_per_bin**2 * nbins**2) / max(1, number of side-band bins))` (simple estimate);
* centroid = counts-weighted mean of bin centres in the peak window after subtracting `base_per_bin` from every bin (clip negative to 0); if the total is 0 print `nan`.
Observed reference: the observed area of the expected line with the largest `exp_rel` (the strongest expected line); `obs_rel = area / area_ref`.

## Output (stdout)
First line: `# file=<bench_csv> elements=<elements> hp_mode=<from header> n_capture_gamma=<from header>`.
Then a table, one line per expected line, separated by ` | `, exactly these columns in this order:
`element | E_exp_keV | exp_rel | area | sigma_area | obs_rel | obs_over_exp | centroid_keV | dE_keV`
where `obs_over_exp = obs_rel/exp_rel` and `dE_keV = centroid - E_exp`. Print floats with 4 significant digits (`%.4g`). After the table print one summary line per element: `SUMMARY <element> lines=<n> detected=<count with area > 3*sigma_area> median_obs_over_exp=<median over detected lines, nan if none>`.
Exit code 0 always when parsing succeeded.

## Pitfalls of the first attempt (all MUST be avoided)
1. The reference for `obs_rel` is the OBSERVED area (counts) of the strongest expected line — NOT its expected weight `w` (barn). Dividing counts by barn is meaningless. Compute all `area` values first, then `area_ref = area of the line with the largest exp_rel`, then `obs_rel = area / area_ref` (if `area_ref <= 0` print `nan` for `obs_rel` and `obs_over_exp`).
2. Merge lines closer than 0.3 keV (function `merge_close_lines`) BEFORE sorting and truncating to `top_n`, and actually call it.
3. Compute the peak quantities ONLY ONCE, inside `analyse`, which returns a list of dicts with the keys `element, E_exp, exp_rel, area, sigma_area, obs_rel, obs_over_exp, centroid, dE, detected` (`detected = area > 3*sigma_area`). `main` only prints; it must NOT recompute areas. The SUMMARY line uses the same dicts (`median_obs_over_exp` over lines with `detected == True` and finite `obs_over_exp`).
4. Bin index arithmetic: bin `i` has centre `i+0.5` keV; a bin belongs to a window `[lo, hi]` when its centre is inside, i.e. `lo <= i+0.5 <= hi`. Use `range(int(math.ceil(lo-0.5)), int(math.floor(hi-0.5))+1)` clipped to `0..len(hist)-1` for the bin indices of a window.
5. Use `df.ffill()` (not `fillna(method=...)`, removed in pandas 3) for forward filling.
6. The strongest expected line of an element is the reference for THAT element; the observed reference is per element.

## Requirements
* Functions: `read_bench(path)`, `read_expected(data_dir, symbols, top_n)` returning a dict symbol -> list of (E, w), `analyse(hist, expected)`, `main()`; guard with `if __name__ == "__main__":`.
* At the top `sys.stdout.reconfigure(encoding="utf-8")`, `sys.stderr.reconfigure(encoding="utf-8")`. Read files with `encoding="utf-8"`.
* Any parse failure: message to stderr, exit 2. Comments in Russian, short. Keep every expression on one line and simple: previous scripts of this project failed `ast.parse` on unbalanced parentheses in long expressions.
