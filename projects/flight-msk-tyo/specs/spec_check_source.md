Write ONE complete Python 3 script `check_source.py`. Output ONLY the code (no markdown fences, no prose). Standard library only (no numpy).

## Purpose
Compare the histograms written by the Geant4 application `src_check` (particles that entered a probe sphere in an empty world) with the PARMA source table the particles were drawn from.

## Command line
    python check_source.py <table_file> <src_check_csv>
Wrong argument count -> usage to stderr, exit 2.

## Table file (text, UTF-8)
Line 1: comment starting with `#`. Line 2: `nebin nabin ie511 flux511_cont flux511_line total_flux`. Line 3: `nebin+1` energy bin edges (MeV). Then `nebin` lines with `nabin` numbers `D[k][ia]` (k = 1..nebin, ia = 1..nabin), units /cm2/s/MeV per angular bin.
Model of the source (must be reproduced exactly):
* `mass_k = (sum over ia of D[k][ia]) * (edges[k]-edges[k-1])`; for k == ie511 (only if ie511 > 0) add `flux511_line`; `p_k = mass_k / sum(mass)`.
* conditional angular probability `q[k][ia] = D[k][ia] / sum over ia' of D[k][ia']` (0 if the row sum is 0); cosine marginal `c_ia = sum over k of p_k * q[k][ia]`.

## CSV file written by src_check (text)
Line 1: `# N=<N> hits=<hits> underflow=<u> overflow=<o> TotalFlux=<..> R_cm=<..> probe_r_cm=<..> seed=<..> particle=<name>` — parse `key=value` tokens after the `#`.
Line 2: `E_hist,<nebin integers>`; column j (0-based) is energy bin j+1.
Line 3: `C_hist,<nabin integers>`; column j is cosine bin j+1.

## Checks (print one line each, final exit code 0 only if all PASS)
* **H (hit fraction)**: expected fraction `f = (probe_r_cm/R_cm)^2`; `obs = hits`, `exp = N*f`; `z = (obs-exp)/sqrt(N*f*(1-f))`. PASS if `|z| < 4.5`. Print `CHECK H hits=<obs> exp=<exp> z=<..> PASS|FAIL`.
* **S (sums)**: the sum of `E_hist` plus `underflow` plus `overflow` must equal `hits`, and the sum of `C_hist` must equal `hits`. PASS if both hold. Print `CHECK S sum_E=<..> sum_C=<..> hits=<..> PASS|FAIL`.
* **E (energy groups)**: group the fine energy bins into groups of 20 consecutive bins (`nebin/20` groups). expected count of a group = `hits * sum of p_k` over its bins; observed = sum of `E_hist` over its bins. Use only groups with expected >= 25. `z = (obs-exp)/sqrt(exp)`. PASS if `max|z| < 4.5` and `chi2/ndf < 1.6`. Print `CHECK E groups=<ndf> max_abs_z=<..> chi2_ndf=<..> PASS|FAIL`.
* **C (cosine bins)**: expected count of bin `ia` = `hits * c_ia`, observed = `C_hist[ia-1]`; same rules as E (use bins with expected >= 25). Print `CHECK C bins=<ndf> max_abs_z=<..> chi2_ndf=<..> PASS|FAIL`.
Finally print `SOURCE_CHECK PASS` or `SOURCE_CHECK FAIL`.

## Requirements
* At the top: `import sys`, `import math`, `sys.stdout.reconfigure(encoding="utf-8")`, `sys.stderr.reconfigure(encoding="utf-8")`.
* Functions: `read_table(path)`, `read_csv(path)`, `main()`; guard with `if __name__ == "__main__":`.
* Read files with `encoding="utf-8"`. Any parse failure: message to stderr and exit 2.
* Comments in Russian, short. Beware of parenthesis balance in long expressions (a previous attempt at another script failed `ast.parse` on an unmatched `)`); keep expressions simple, one per line.
