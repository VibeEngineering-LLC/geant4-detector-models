OLLAMA TASK SPEC — generate ONE Python module, output ONLY raw Python code, no markdown fences, no explanations.

# Target file: coverage_windows.py (module in the same folder as fit_coverage.py)

## Purpose
Diagnostic printer: signed form-residual (%) of one or more fitted spectrum models
versus the measured spectrum, inside the windows of the diagnostic gamma lines
plus the muon window. Used by fit_coverage.py (#FIT-1, 2026-08-22).

## Module contents — exactly one function, module-level imports only

```python
import numpy as np

def print_window_residuals(e_meas, cps_meas, models, fwhm_fn):
    ...
```

Parameters:
- e_meas: 1D numpy array, energy (keV) of each measured channel.
- cps_meas: 1D numpy array, measured cps per channel, same length as e_meas.
- models: list of tuples (name: str, model: 1D numpy array same length as cps_meas).
  Typically 2 models, but code must work for any count >= 1.
- fwhm_fn: callable, fwhm_fn(E0_keV) -> FWHM in keV at that energy.

## Behavior

1. Hard-coded list of diagnostic line windows:
   WINS = [("Bi214 609", 609.3), ("Ac228 911", 911.2),
           ("K40 1460", 1460.8), ("Tl208 2614", 2614.5)]
2. For each (label, e0) in WINS:
   - sigma = fwhm_fn(e0) / 2.35482
   - window mask: (e_meas >= e0 - 2.5*sigma) & (e_meas <= e0 + 2.5*sigma)
     (SAME window convention nsig=2.5 as fit_lines.line_net_area — do not change 2.5)
   - ym = sum of cps_meas inside mask
   - for each model: rm = 100.0 * (sum(model inside mask) - ym) / ym  (signed, percent);
     if ym <= 0 use float("nan")
3. After the line windows, one extra row for the muon window:
   label "mu 2700-2828", mask (e_meas >= 2700) & (e_meas < 2828), same math.
4. Printing (plain print, ASCII only in code literals):
   - header line: "=== NEVYAZKA PO OKNAM LINIJ (gross, znak: +model vyshe izmereniya) ==="
   - column header: window label column width 18, then "izm,cps" width 11,
     then for each model two columns: "<name> cps" width 12 and "<name> %" width 9.
   - each row: label, ym as %11.4f, then per model: window sum as %12.4f and
     signed percent as %+9.1f (or "nan").
5. No side effects besides printing. Return a dict {label: {name: rm_percent}} with
   the signed residuals, so callers can reuse the numbers.

## Constraints
- Python 3.10, numpy only (no scipy).
- No if __name__ block, no argparse, no file IO.
- Keep it short and flat: one function, no classes, no helper functions.
- All string literals in code ASCII only (transliterate Russian).
