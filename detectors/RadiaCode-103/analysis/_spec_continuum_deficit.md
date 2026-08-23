OLLAMA TASK SPEC — generate ONE Python module, output ONLY raw Python code, no markdown fences.

# Target file: continuum_deficit.py (module, same folder as fit_coverage.py/fit_nuclides.py)

## Purpose
Diagnostic script (#FIT-1 follow-up, 2026-08-22): decompose the fit residual
(measured spectrum minus the line-fit model) spectrally, to find out WHAT physical
component is missing from the template basis (continuum deficit was found to be ~62%
of the total counts). No new Geant4 runs — pure post-processing of existing arrays.

## Behavior — script with a main() function, run under `if __name__ == "__main__":`

1. Imports: `import sys, os` then insert this dir into sys.path (same pattern as
   fit_coverage.py: `_HERE = os.path.dirname(os.path.abspath(__file__))`,
   `sys.path.insert(0, _HERE)`), then `import numpy as np`, `import fit_coverage as fc`,
   `import fit_nuclides as fn`.

2. In main():
   - Call `names, A, e_meas, cps_meas, live = fc.assemble()`.
   - Call `import io, contextlib` then, suppressing stdout the same way fit_coverage.py
     does (`with contextlib.redirect_stdout(io.StringIO()):`), call
     `amps_dict, a_mu = fn.main()`.
   - Build `amp = np.array([amps_dict.get(n, 0.0) for n in names])`.
   - Build `model = A @ amp`.
   - Build `diff = cps_meas - model` (this is the deficit spectrum, same units as
     cps_meas: counts per second per channel).

3. Print a header: `"=== DEFICIT SPEKTRA (izmerenie minus model iz linij) ==="`.

4. Band-wise summary — reuse the SAME band tuple already used elsewhere in the
   project: `bands = ((20, 100), (100, 300), (300, 700), (700, 1500), (1500, 2000),
   (2000, 2400), (2400, 2830))`. For each band print: band range, sum of cps_meas in
   band, sum of model in band, sum of diff in band, and diff as percent of measured
   sum in that band (`100 * diff_sum / meas_sum` if meas_sum > 0 else nan). Use
   `print("%5d-%-6d %10.5f %10.5f %10.5f %8.1f%%" % (...))` with a matching header
   line printed first.

5. Peak search in the deficit spectrum — reuse the DIAGNOSTIC LINE TABLE and the
   EXISTING window-net function from fit_nuclides.py, do NOT reimplement peak finding:
   `fn.DIAG` is a list of (nuclide_name, energy_keV) tuples; `fn.net_window(spec,
   e_grid, e0, live=None)` returns `(net_area, sigma)` for one line window, using the
   project's established window convention. For each `(nuc, e0)` in `fn.DIAG`: call
   `net, sd = fn.net_window(diff, e_meas, e0, live)`. If `sd > 0` and `net > 3.0 * sd`
   (statistically significant peak found IN THE DEFICIT), record it. Print a section
   header `"=== ZNACHIMYE PIKI V DEFICITE (net > 3 sigma) ==="` then for each
   significant line: `print("%-8s %8.1f keV   net=%10.3e cps   sigma=%9.1f" %
   (nuc, e0, net, net/sd))`. If none found, print
   `"net znachimyh pikov - deficit gladkij (kontinuum/rasseyanie/fon), ne otdelnye linii"`.

6. Save the full diff array to CSV for later inspection: write
   `results/deficit_spectrum_20260822.csv` (path relative to the project's existing
   `results/` folder — build the path as
   `os.path.join(_HERE, "..", "results", "deficit_spectrum_20260822.csv")`)
   with two columns `energy_keV,diff_cps`, one row per channel, header line first,
   using plain Python file I/O (open/write), NOT pandas.

7. Return 0 from main(), call `sys.exit(main())` under the `if __name__` guard.

## Constraints
- Python 3.10, numpy only besides project imports.
- No classes, one flat main() function.
- All print() string literals ASCII only (transliterate Russian per project convention
  visible in fit_nuclides.py / fit_coverage.py — e.g. "kontinuum" not "континуум").
- Do not reimplement net_window or the diagnostic line table — import and reuse them.
