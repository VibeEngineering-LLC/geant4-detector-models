Write ONE complete Python 3 script `compare_smeared.py`. Output ONLY the code (no markdown fences, no prose). Allowed imports: `sys`, `math`, `pandas` (with `xlrd`), no other third-party modules.

## Purpose
Compare, per capture, the gamma spectrum produced by the Geant4 application `line_bench` (thermal-neutron capture on a single-element target) with the spectrum reconstructed from the IAEA PGAA line database, AFTER folding both with the energy resolution of the CsI scintillation detector. Line-by-line comparison is not meaningful when transitions differ by tens of keV (below the detector resolution); comparison in wide windows is.

## Command line
    python compare_smeared.py <bench_csv> <element> <Sn_keV> <data_dir> [fwhm662_keV=41.6] [window_keV=500]
Wrong arguments -> usage to stderr, exit 2.

## Inputs
**bench_csv** (UTF-8): line 1 `# key=value ...` tokens after `#` (need `hp_mode`, `n_capture_gamma`); line 2 `capture,<12000 integers>` (prefix has 8 characters: use `line2[len("capture,"):]`); bin `j` = `[j, j+1)` keV, centre `j+0.5`.
**promptgammas.xls** (`pd.read_excel(path, header=None, skiprows=2)`, columns by POSITION): `1` product label like `'  28-Al'` (STRIP whitespace; format `<A>-<Symbol>`, mass number of the PRODUCT nucleus), `5` gamma energy keV, `7` partial cross section `sigma_gamma` barn. `pd.to_numeric(errors="coerce")`, drop NaN.
**isotope.xls** (`header=None, skiprows=2`, columns by position): `0` element symbol (`.ffill()`), `1` TARGET isotope label like `'  27-Al '` (STRIP), `5` abundance percent, `8` thermal cross section `sigma0` barn.

## Method
1. `Ncap = sum_j (j+0.5)*count[j] / Sn_keV` (energy conservation). If `Ncap <= 0` -> stderr, exit 2. Geant4 per-capture spectrum `g[j] = count[j]/Ncap`.
2. PGAA lines of the element: keep prompt rows whose stripped product label ends with `f"-{element}"`; the TARGET label of a row is `f"{Ap-1}-{element}"` (`Ap` = product mass number); look up `(abundance, sigma0)` in a dict keyed by the isotope-table label (already a target label — do NOT subtract 1 again). Element-weighted strength `num = ab/100 * sigma_gamma`; `den_total = sum over all isotopes of the element of ab/100*sigma0` (one number). Expected yield per capture of a line: `y = num/den_total` (per capture, NOT per 100). Merge lines within 0.3 keV (sort by energy first, strength-weighted energy, sum of `y`). Keep lines with `20 <= E <= 12000`.
3. Energy-balance diagnostics: `E_model = sum over lines of E*y` and `mult_model = sum over lines of y`; for Geant4 `E_g4 = sum_j (j+0.5)*g[j]` (equals `Sn_keV` by construction) and `mult_g4 = sum_j g[j]`.
4. Folding. Resolution `FWHM(E) = fwhm662 * sqrt(E/661.657)` keV, `sigma(E) = FWHM/2.35482`. Fold on the 1-keV grid of 12000 bins: for every source with energy `E0` and weight `w` add to every bin `k` with `|k+0.5-E0| <= 4*sigma(E0)` the amount `w * exp(-0.5*((k+0.5-E0)/sigma)**2) / (sigma*sqrt(2*pi))` (bin width 1 keV). For the model, sources are the merged lines (`E0` = line energy, `w = y`). For Geant4, sources are all bins with `g[j] > 0` (`E0 = j+0.5`, `w = g[j]`); to keep it fast skip bins with `g[j] == 0`. Result arrays `S_model[k]`, `S_g4[k]` in units of gammas per capture per keV.
5. Windows: consecutive windows of `window_keV` from 0 to 12000 keV (`floor(12000/window)` windows). Window sums `W_model = sum of S_model[k]` and `W_g4` over bins with centre in the window.

## Output (stdout)
Line 1: `# file=<bench_csv> element=<element> hp_mode=<..> Ncap_est=<%.6g> lines=<n merged lines>`.
Line 2: `BALANCE E_model_keV=<%.6g> Sn_keV=<%.6g> mult_model=<%.4g> mult_g4=<%.4g>`.
Then one line per window: `WIN lo-hi | model_per_capture | g4_per_capture | ratio_g4_over_model` with `%.4g` (`ratio = nan` if model < 1e-6).
Last line: `WINDOWS n=<count of windows with model >= 0.005> within30=<count of those with 0.7 <= ratio <= 1.3> median_ratio=<median of their ratios>`.
Exit code 0 on success.

## Requirements
* Functions: `read_bench`, `read_lines`, `fold`, `main`; guard with `if __name__ == "__main__":`.
* Top of file: `sys.stdout.reconfigure(encoding="utf-8")`, `sys.stderr.reconfigure(encoding="utf-8")`; read files with `encoding="utf-8"`.
* Format numbers with a helper `def fmt(x): return "nan" if x != x else f"{x:.4g}"`; never put a conditional expression inside an f-string format spec.
* No variable may leak between iterations; compute each quantity once. Comments in Russian, short. Keep every expression simple and on one line (earlier scripts failed on unbalanced parentheses, unstripped labels, wrong merge keys and double subtraction of the mass number).
