Write ONE complete Python 3 script `yield_lines.py`. Output ONLY the code (no markdown fences, no prose). Allowed imports: `sys`, `math`, `pandas` (with `xlrd` for `.xls`).

## Purpose
Absolute yield check of neutron-capture gamma lines: for a `line_bench` run on a single-element target compare, for the strongest expected lines, the observed number of gammas PER CAPTURE with the tabulated yield of the IAEA PGAA database (gammas per 100 captures).

## Command line
    python yield_lines.py <bench_csv> <element> <Sn_keV> <data_dir> [top_n]
* `element`: symbol, e.g. `Al`. `Sn_keV`: neutron separation energy (keV) of the capture product, i.e. the total energy released per capture (Al-28: 7725.2; Fe-57: 7646.2; Ti-49: 8142.4; Cu-64: 7915.9; Si-29: 8473.6). `top_n` default 10 (number of strongest expected lines by yield).
* Wrong arguments -> usage to stderr, exit 2.

## Inputs
**bench_csv** (UTF-8 text): line 1 `# material=... hp_mode=... radius_cm=... n_capture_gamma=...` (`key=value` tokens after `#`); line 2 `capture,<12000 integers>` (1-keV bins, bin `j` = `[j, j+1)` keV, centre `j+0.5`); line 3 `inelastic,...` (ignore).
**promptgammas.xls** (`pd.read_excel(path, header=None, skiprows=2)`, columns by POSITION): `1` = product isotope label like `'  28-Al'` (STRIP whitespace!; format `<A>-<Symbol>`, mass number of the PRODUCT), `5` = energy keV (numeric), `7` = partial cross section `sigma_gamma` (barn). Convert with `pd.to_numeric(errors="coerce")`, drop NaN.
**isotope.xls** (`header=None, skiprows=2`): `0` element symbol (forward fill with `.ffill()`), `1` target isotope label like `'  27-Al '` (STRIP whitespace!), `5` abundance percent, `8` total thermal cross section `sigma0` (barn, numeric).

## Method
1. Number of captures `Ncap = sum over bins j of (j+0.5)*count[j] / Sn_keV` (energy conservation; conversion electrons are not counted, so this is a lower bound; say so in the output).
2. Expected yield per 100 captures for a line of product isotope `P` (target isotope `T = f"{Ap-1}-{symbol}"`): `Y = 100 * sigma_gamma / sigma0(T)`. Element-weighted expected yield for the natural element: `Y_el = sum over isotopes T of (ab_T/100 * sigma0_T * Y_T) / sum over isotopes T (ab_T/100 * sigma0_T)` — implement as: `num(line) = ab_T/100 * sigma_gamma`, `den = sum over ALL target isotopes of the element (ab_T/100 * sigma0_T)`, `Y_el = 100 * num / den`. Merge rows of the element whose energies differ by less than 0.3 keV (sum `num`, strength-weighted energy). Keep lines with `50 <= E <= 11500` keV, sort by `Y_el` descending, keep `top_n`.
3. Observed yield: half window `hw = 3 + 0.0006*E` keV; peak bins = bins whose centre lies in `[E-hw, E+hw]`; sidebands: two windows 8 keV wide adjacent to the peak window (bins by centre); baseline per bin = mean over the sideband bins; `area = peak_sum - baseline*n_peak_bins`; `sigma = sqrt(max(peak_sum,1))`; `obs_yield = 100 * area / Ncap`.
4. Output table (stdout): first line `# file=<bench_csv> element=<element> hp_mode=<..> n_capture_gamma=<..> Ncap_est=<%.6g>`; then one line per expected line, columns separated by ` | `: `E_keV | Y_exp_per100 | area | obs_per100 | sigma_per100 | obs_over_exp | dE_keV(centroid-E)` with `%.4g` formatting (`obs_over_exp = obs_per100/Y_exp`, `nan` if `Y_exp <= 0`). Centroid = counts-weighted mean of bin centres in the peak window after baseline subtraction (clip at 0; `nan` if total 0).
5. Last line: `VERDICT lines=<n> ok=<count with 0.5 <= obs_over_exp <= 2.0 and area > 3*sqrt(max(peak_sum,1))> median_ratio=<median of obs_over_exp over lines with area > 3*sigma, nan if none>`.
Exit code 0 if parsing succeeded; if `Ncap_est <= 0` or no expected lines: message to stderr and exit 2.

## Pitfalls of the first attempt (ALL must be avoided)
1. The prompt-gamma label is the PRODUCT nucleus (`28-Al`), the isotope table label is the TARGET (`27-Al`). Do NOT merge on equal labels: build the target label `f"{Ap-1}-{symbol}"` for every prompt row (parse `Ap` and `symbol` from the stripped product label) and look it up in a dict `target_label -> (abundance_percent, sigma0)`.
2. Filter prompt rows to the requested element symbol BEFORE any computation.
3. The expected yield is PER LINE: `Y_el = 100 * num_line / den_total` where `num_line = ab/100 * sigma_gamma` (after merging rows closer than 0.3 keV) and `den_total = sum over all target isotopes of the element of (ab/100 * sigma0)` is ONE fixed number. Never divide the sum of all `num` by `den_total` and reuse it for all lines.
4. Compute each line's peak quantities ONCE, in one function `line_stats(counts, E, Ncap)` returning a dict with `area, sigma, peak_sum, obs_per100, centroid, detected`; `main` only prints and computes the VERDICT from those dicts. No variable from a previous line may leak into another (the first attempt reused a stale `baseline`).
5. Do NOT put a conditional expression inside an f-string format spec. Format numbers with a helper `def fmt(x): return "nan" if x != x else f"{x:.4g}"` and call it.
6. The prefix `capture,` has 8 characters: use `line2[len("capture,"):]`.
7. Merge close lines by sorting the rows by energy first, then merging neighbours within 0.3 keV (strength-weighted energy), and only THEN sort by `Y_el` and cut to `top_n`.

## Requirements
* Functions: `read_bench`, `read_expected`, `main`; guard with `if __name__ == "__main__":`.
* At the top `sys.stdout.reconfigure(encoding="utf-8")`, `sys.stderr.reconfigure(encoding="utf-8")`; read files with `encoding="utf-8"`.
* Comments in Russian, short. Keep every expression on one line and simple (previous scripts of this project failed on unbalanced parentheses in long expressions and on unstripped isotope labels).
