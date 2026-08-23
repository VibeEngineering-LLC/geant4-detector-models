Write a Python 3 script `make_m1_macros_v3.py` for a Geant4 detector-modeling project.
It generates Geant4 GPS macro pairs for METHOD 1 (single-nuclide-chain-link full decay).

CONTEXT: the previous version (v2) used `/gps/hist/type arb` + `/gps/hist/inter Lin`
(continuous-density interpolation between points). That was found to be statistically
UNRELIABLE in this Geant4 build (11.02-patch-01): identical input histograms randomly
produced either a silent fallback to the default 1 MeV mono-energetic source, or NaN
primary energies, depending on RNG state between runs. The fix: switch to
`/gps/ene/type User` + `/gps/hist/type energy` (a plain STEPPED histogram of bin
CONTENTS, no interpolation between points at all — the `/gps/hist/inter` command must
NOT be used for this type).

CRITICAL SEMANTIC DIFFERENCE for the "energy" histogram type (verified against the
official Geant4 GPS manual, section 2.7.2.5 / Table 2.7):
- Each `/gps/hist/point Ehi Weight` pair means: Ehi is the UPPER edge of a bin, Weight
  is that bin's CONTENT (not a density-curve node).
- EXCEPTION: the very FIRST data pair's first value is the LOWER edge of the FIRST bin,
  and its second value (weight) is IGNORED by Geant4.
- Because there is no interpolation, EVERY bin along the full energy range must be
  emitted explicitly, INCLUDING ZERO-CONTENT bins. Skipping a zero bin between two
  distant points would make Geant4 treat the gap as ONE WIDE bin whose content is the
  weight of the LAST point — silently redistributing flux into the wrong place. This is
  the opposite of the old `arb` method, where skipping zeros between sparse points was
  safe (continuous interpolation naturally went to zero there).
- The `add_walls()` trick from v2 is NOT needed any more and must be REMOVED entirely —
  it was a workaround specific to the `arb`+`Lin` interpolation bug.
- Geant4 histograms are limited to 1024 bins total (hard framework limit, verified by
  binary search on this build: 1024 -> OK, 1025 -> STATUS_STACK_BUFFER_OVERRUN).

REQUIRED SCRIPT STRUCTURE (mirror v2's proven layout closely, only the histogram-building
and macro-writing logic changes):

```python
# -*- coding: utf-8 -*-
"""(module docstring explaining the switch to `energy`/`User`+ND-003 rationale,
   pointing to DECISIONS.md D-003 for the full writeup)"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(_HERE, "..", "..", ".."))
BUILD = os.path.join(REPO, "build", "RadiaCode-103")
RESULTS = os.path.normpath(os.path.join(_HERE, "..", "results"))
WALLION = os.path.join(RESULTS, "wallion")

NUCS = ["K40", "Ra226", "Pb214", "Bi214", "Pb212", "Ac228", "Bi212", "Tl208"]
SRC_TAG = "m1"

MAX_HIST_POINTS = 950  # same hard-won safety margin as v2 (real Geant4 limit is 1024)
NATIVE_STEP_KEV = 2.0  # wallfield.cc's native bin width
```

Function 1: `read_wallfield_csv(csv_path)` -> returns two lists `(energies_kev, fluences)`
parsed from the CSV. CSV format: comment lines start with `#`, data lines are
`E_keV,fluence_cm2_s`. Skip blank lines and lines that don't parse as two floats. This is
IDENTICAL parsing logic to v2 — copy it verbatim (lines 133-149 in v2, reproduced here):
```python
        energies = []
        fluences = []
        with open(csv_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split(",")
                if len(parts) != 2:
                    continue
                try:
                    e_kev = float(parts[0])
                    fluence = float(parts[1])
                    energies.append(e_kev)
                    fluences.append(fluence)
                except ValueError:
                    continue
```
Only non-zero-fluence rows are present in the CSV (wallfield.cc only writes bins with
`fHist[i] != 0`), so the input is SPARSE — gaps between consecutive energies are real
gaps of zero-content bins that must be reconstructed.

Function 2: `build_full_grid(energies_kev, fluences, step_kev)` -> returns
`(grid_energies_kev, grid_fluences)`, two lists representing bin CENTER energies and
their contents on a COMPLETE, gap-free grid at the given `step_kev`, spanning from
`min(energies_kev)` to `max(energies_kev)` inclusive, with a zero fluence at every grid
point that wasn't present in the sparse input. Use a dict keyed by
`round(energy, 3)` to place input values onto the grid safely (float rounding). The
number of grid points is `round((max-min)/step_kev) + 1`.

Function 3: `coarsen_grid(grid_energies_kev, grid_fluences, native_step_kev, max_points)`
-> if `len(grid_energies_kev) <= max_points - 1` (reserve one slot for the leading edge
point), return the grid UNCHANGED. Otherwise, compute the smallest integer multiple of
`native_step_kev` (i.e. `native_step_kev * k` for k=2,3,4,...) such that regridding the
full span at that coarser step keeps total points <= `max_points - 1`, then SUM (not
sample) the fluences of all fine bins falling into each coarse bin — the total fluence
sum must be preserved exactly (assert `abs(sum(coarse_fluences) - sum(grid_fluences)) /
sum(grid_fluences) < 1e-9`, exit with code 3 and a clear Russian-language error message
printed to stdout if this check fails). Return `(coarse_energies_kev, coarse_fluences)`
as bin-center energies again.

Function 4: `write_energy_macro(spectrum_path, grid_energies_kev, grid_fluences,
step_kev, total_fluence, nuc)` -> writes the Geant4 macro file. Format:
```
# Единичный (1 Бк/кг) отклик звена <nuc>, wallfield.exe.
# МЕТОД 1: полный распад ТОЛЬКО этого звена, nucleusLimits
# отсекает дочерние (канон geant4-spectrum-pipeline, D-001).
# ВЕРСИЯ 3 (22.08): /gps/hist/type energy (ступенчатая гистограмма), НЕ arb+Lin —
# см. DECISIONS.md D-003 (нестабильность Arb+Lin, дефолт-откат/NaN на реальных прогонах).
# FLUENCE_TOTAL_CM2_S = <total_fluence with 6 decimals>
/gps/particle gamma
/gps/ene/type User
/gps/hist/type energy
/gps/hist/point <lower_edge_of_first_bin_in_MeV> 0.0
/gps/hist/point <upper_edge_of_bin_1_in_MeV> <content_of_bin_1>
/gps/hist/point <upper_edge_of_bin_2_in_MeV> <content_of_bin_2>
... (one line per grid point, in ascending energy order)
```
Energies are in MeV (divide keV values by 1000.0), formatted `%.6f`; fluence weights
formatted `%.6e`. The lower edge of the first bin is `grid_energies_kev[0] - step_kev/2`
converted to MeV (clip to 0.0 if it would go negative). Each subsequent point's energy is
`grid_energies_kev[i] + step_kev/2` converted to MeV (the upper edge of that bin). Do
NOT write a `/gps/hist/inter` line — it must be absent for the `energy` histogram type.

Function `main()`: for each nuclide in `NUCS`:
1. Build `csv_path = os.path.join(WALLION, "wf_%s_%s.csv" % (SRC_TAG, nuc))`; if missing,
   print a Russian warning (`u"ПРЕДУПРЕЖДЕНИЕ: Отсутствует файл %s" % csv_path`) and
   `continue`.
2. Parse it with `read_wallfield_csv`. Compute `total_fluence = sum(fluences)`; if
   `total_fluence <= 0`, print `u"ОШИБКА: Общая флюенция для %s равна нулю или меньше" %
   nuc` and `sys.exit(3)`.
3. Call `build_full_grid(energies, fluences, NATIVE_STEP_KEV)`.
4. Call `coarsen_grid(...)` with `MAX_HIST_POINTS`. Print a Russian info line reporting
   before/after point counts whenever coarsening actually happened (mirror v2's style:
   `u"  %s: точек %d -> %d (ребиннинг на шаг %.1f кэВ, сумма потока сохранена)" % (nuc,
   n_before, n_after, used_step_kev)`).
5. If after coarsening `len(grid) > MAX_HIST_POINTS - 1`, print a Russian error and
   `sys.exit(3)` (should not happen given function 3's guarantee, but keep as a hard
   safety check mirroring v2's defensive style).
6. Write the spectrum macro via `write_energy_macro` to
   `os.path.join(RESULTS, "field_spectrum_m1_%s.mac" % nuc)`.
7. Regenerate the run-macro EXACTLY like v2 does (copy lines 203-221 of v2 verbatim —
   reads a template from
   `os.path.join(BUILD, "_attic_table_method_20260821", "field_run_nucb_%s.mac" % nuc)`,
   rewrites its `/control/execute` and `/rc/outFile` lines, writes to
   `os.path.join(BUILD, "field_run_m1b_%s.mac" % nuc)`). Exit code 3 with a Russian error
   if the template is missing.
8. Track a summary list of `(nuc, total_fluence, hist_points, spectrum_path, run_path)`.

At the end, print `u"Сгенерировано пар макросов: %d из %d" % (written_pairs, len(NUCS))`,
a summary table (mirror v2's formatting), and `sys.exit(0)` if all 8 were written,
otherwise print how many are missing and `sys.exit(1)`.

Keep ALL user-facing print strings in Russian (this project's convention — see the
Russian text embedded above, copy that style exactly). Use plain stdlib only, no
external dependencies. Target Python 3.8+ syntax (no walrus-heavy tricks, this repo runs
on Windows with a standard CPython install). Output ONLY the complete Python file
contents, no explanation, no markdown code fences.
