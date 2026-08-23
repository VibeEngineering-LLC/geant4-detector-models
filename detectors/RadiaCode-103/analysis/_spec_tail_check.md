OLLAMA TASK SPEC — generate ONE Python module, output ONLY raw Python code, no markdown fences.

# Target file: tail_check.py (module, same folder as fit_coverage.py/fit_nuclides.py)

## Purpose
Diagnostic script (#FIT-1 follow-up, 2026-08-22): cheap empirical check of the
"incomplete light collection -> low-energy tailing" hypothesis for the CsI(Tl)
crystal, WITHOUT running full Geant4 optical photon transport. The idea: each
line-fit model channel with counts m at energy E0 gets an additional exponential
tail added at all energies E < E0 (standard low-energy-tail parametrization used
in gamma spectrometry, e.g. Genie2000-style, simplified to pure exponential decay,
no erfc term). Fit the two free parameters (relative tail amplitude, tail decay
length in keV) to best match the measured spectrum, and report whether this closes
the continuum deficit found earlier (~62% form residual).

## Behavior — script with a main() function, run under `if __name__ == "__main__":`

1. Imports: `import sys, os` then insert this dir into sys.path (same pattern as
   fit_coverage.py), then `import numpy as np`, `import fit_coverage as fc`,
   `import fit_nuclides as fn`, `from scipy.optimize import minimize_scalar,
   minimize`.

2. Helper function `add_tail(model, e_meas, a_frac, beta_kev)`:
   - `model` is a 1D numpy array of cps per channel (already gaussian-folded line-fit
     model, e.g. from A @ amp).
   - `e_meas` is the energy grid (keV) matching model, assumed roughly uniform
     spacing (use `de = np.median(np.diff(e_meas))` for the channel width).
   - For each channel i with model[i] > 0, add to ALL channels j with e_meas[j] <
     e_meas[i]: `tail_contrib = model[i] * a_frac * np.exp(-(e_meas[i]-e_meas[j])/beta_kev) * de / beta_kev`
     (the `de/beta_kev` factor keeps the total added tail area proportional to
     a_frac * model[i], independent of the energy grid step). Do this VECTORIZED
     with numpy broadcasting, not a nested Python double loop (this is called
     repeatedly inside an optimizer, must be reasonably fast for ~1400 channels).
     A clean way: build a lower-triangular kernel matrix once outside the optimizer
     loop is NOT required for this script — a broadcasted per-call computation is
     fine as long as it uses numpy operations (outer subtraction + np.where + einsum
     or matrix multiply), not Python-level loops over both i and j.
   - Return `model + tail_added` (the tail-augmented model), same shape as model.

3. In main():
   - Call `names, A, e_meas, cps_meas, live = fc.assemble()`.
   - Suppress stdout the same way fit_coverage.py does
     (`import io, contextlib; with contextlib.redirect_stdout(io.StringIO()):`)
     and call `amps_dict, a_mu = fn.main()`.
   - `amp = np.array([amps_dict.get(n, 0.0) for n in names])`
   - `model_base = A @ amp` (this is the line-fit model WITHOUT tail — same as
     model_chi2 used elsewhere).
   - Restrict to the analysis window used elsewhere in the project:
     `mask = (e_meas >= 20) & (e_meas < 2830)`, work with `e_meas[mask]`,
     `cps_meas[mask]`, `model_base[mask]` from here on (call them e, y, m0).

4. Define an objective function `objective(params)` where `params = [a_frac,
   beta_kev]`: compute `model_tail = add_tail(m0, e, a_frac, beta_kev)`, return the
   sum of squared relative residuals `np.sum(((model_tail - y) / np.maximum(y,
   1e-6)) ** 2)`. Guard: if a_frac < 0 or beta_kev <= 0, return a large number
   (e.g. 1e12) instead of calling add_tail (invalid physical parameters).

5. Run `scipy.optimize.minimize(objective, x0=[0.5, 200.0], method="Nelder-Mead",
   options={"xatol": 1e-3, "fatol": 1e-6, "maxiter": 500})`. Print the optimizer
   result: success flag, final a_frac, final beta_kev, final objective value.

6. Compute BOTH the baseline (no tail) and best-fit (with tail) versions of two
   metrics, reusing existing project functions — DO NOT reimplement them:
   `fc.fraction_covered(model, meas)` and `fc.form_residual_pct(model, meas)`.
   Print a small before/after table:
   ```
   === TAIL-FIT (#FIT-1 follow-up) ===
   podobrannye parametry: a_frac=<value>, beta_kev=<value>
   %% zapolneniya formoj:  baseline=<X>%%  s tail=<Y>%%
   nevyazka formy:         baseline=<X>%%  s tail=<Y>%%
   ```
   using `print("%-24s baseline=%6.1f%%  s tail=%6.1f%%" % (...))` for each metric
   row, computed on the (e, y, m0) already restricted to the 20-2830 keV window.

7. Band-wise before/after (reuse the SAME band tuple as elsewhere in the project):
   `bands = ((20, 100), (100, 300), (300, 700), (700, 1500), (1500, 2000),
   (2000, 2400), (2400, 2830))`. For each band print energy range, measured sum,
   baseline model sum, tail-added model sum, baseline deficit percent, tail-added
   deficit percent. Deficit percent = `100 * (meas_sum - model_sum) / meas_sum`.
   Print a header row first, then one row per band with
   `print("%5d-%-6d %10.5f %10.5f %10.5f %8.1f%% %8.1f%%" % (...))`.

8. Return 0 from main(), call `sys.exit(main())` under the `if __name__` guard.

## Constraints
- Python 3.10, numpy + scipy only besides project imports (fit_coverage, fit_nuclides).
- No classes, flat functions only (add_tail, objective as a closure or module-level
  with a nonlocal/default-arg trick — either is fine as long as it works standalone).
- All print() string literals ASCII only (transliterate Russian per project
  convention — e.g. "podobrannye" not "подобранные" inside code, but the print
  OUTPUT text itself follows the same transliteration convention seen in
  fit_nuclides.py / fit_coverage.py, e.g. "kontinuum", "nevyazka formy").
- Do not reimplement fraction_covered or form_residual_pct — import and reuse from
  fit_coverage.py exactly as they are.
