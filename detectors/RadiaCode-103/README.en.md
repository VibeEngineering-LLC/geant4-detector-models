# RadiaCode 101/102/103 — CsI(Tl) 10×10×10 mm

Consumer scintillation spectrometer-dosimeter. CsI(Tl) crystal, 1 cm³ cube,
SiPM photodetector. The model was built for the sake of **efficiency curves
in Marinelli beakers** — a task the instrument was not designed for and for
which it has no factory curves.

Interactive article with an activity calculator — right here, in
[docs/radiacode-103](../../docs/radiacode-103/), published at
**https://vibeengineering-llc.github.io/geant4-detector-models/**. The former
separate repository `radiacode-curves` is no longer used; there is no need to
reference its address.

## Key results

- Detection efficiency curves for the author's Marinelli beakers, 200 and 500 mL,
  range 50–3000 keV, densities 0.2–1.8 g/cm³.
- Normalization to measurement `K_NORM = 0.833` (see caveat below).
- Effective thickness 19.6 mm (200 mL) and 28.8 mm (500 mL), χ²/dof ≈ 1.
- Cascade summing 0.06–0.33% — negligible, the crystal is small.
- The factory-rated sensitivity of 30 counts/s per 1 µSv/h for Cs-137 **is
  not reproduced by the model**: the calculation gives 22–27. Discussed in
  the [report](REPORT.en.md#on-the-datasheets-30-countss); in short, the conditions of
  the factory measurement are unknown, and the discrepancy has the opposite
  sign from the `K_NORM` systematic, so it cannot be explained by a single
  efficiency error. What should be trusted is the room-background
  cross-check — it was done in situ and agrees to within 11%.
- Beta contribution and empty-beaker background are calculated separately.

## Status of checks under the common protocol

See [common/docs/validation.en.md](../../common/docs/validation.en.md).

| item | status |
|---|---|
| 1. Masses against datasheet | NOT DONE: model mass ≈53 g, datasheet gives instrument weight — can be checked, simply not checked yet |
| 2. Convergence of drawing dimensions | done (drawings + teardown photos) |
| 3. Datasheet efficiency point | **NOT REPRODUCED**: datasheet 30 counts/s per 1 µSv/h, calculation 22–27 |
| 4. **Point geometry as detector reference** | **NOT DONE** |
| 5. Density trend, d_eff | done |
| 6. Summing check on a nuclide without cascade | not required: effect 0.06–0.33% |
| 7. Line yields from own calculation | done |
| 8. Cross-check of reference curve bypassing processing | not applicable: the reference here is an independent LSRM CALCULATION, not a measurement |
| 9. Whether a correction was introduced into the reference curve | not applicable: same reason |
| 10. MDA against datasheet | not applicable |
| 11. Cross-check of data formats | not applicable: measurements in a single format |
| 12. Ready-made numbers as files | done: tables and `curves.json` in `results/` |

## Caveat on K_NORM = 0.833 — must read

This coefficient was compared with 0.858 for [Gamma-1S](../Gamma-1S/) as
confirmation of a "common systematic of the method." **That reading is
wrong.** (The 0.858 value itself was retracted on 2026-07-31: it rested on a
stale export of the calculation grids; the current value is 0.795
(results/compare_lsrm_summary.csv), and the discrepancy with the RadiaCode
is 4.8%, not 3%. The retraction does not
affect the argument below — the case against a "common systematic" never
relied on how close the two numbers were.)

For the Gamma-1S, point geometries separated the detector model from the
vessel model: on the current geometry and the unified extraction convention,
calc/passport is 1.0365 ± 0.0112 and 1.0381 ± 0.0214 at the two distances
(the 0.971 and 0.928 of the earlier revision were retracted by the
entrance-face fix of 2026-07-28), whereas the excess in the Marinelli
beaker — 1.257 — turned out to be a property of the beaker geometry, not
of the method. The coincidence of two coefficients obtained in the SAME
Marinelli-beaker setup proves nothing about the method.

What follows from this for the RadiaCode: the value 0.833 remains valid as
an empirical normalization **for this particular setup**, but its physical
interpretation is open until item 4 of the protocol is completed — a
point-geometry calculation separating the detector model from the beaker
model. This is the main unfinished business for this instrument.

## Reliability: what is known as of 28.07.2026

The model is connected as a submodule to a third-party project (see
[common/docs/consumers.en.md](../../common/docs/consumers.en.md)), and the
consumer needs the **shape** of the curve, not the level: the level is
absorbed by a free multiplier during activity fitting. Therefore defects are
separated by exactly what they corrupt.

The list is maintained based on the results of an external audit. Work on
the instrument resumed on 28.07.2026; items are closed in order, and the
section is edited in the same commit as the code — otherwise it turns into
an outdated warning.

### Fixed

**Matrix `ash`, sum of mass fractions 1.020 instead of 1.000**
(`RCDetector.cc`). The neighboring `soil` has a sum of exactly 1.000,
meaning this is a typo, not a convention. Geant4 responds to this with a
warning in the general output stream, which nobody reads during a batch
run, and the composition remains set as-is. Cost is not equal to magnitude:
attenuation is energy-dependent, so the extra percent of mass produces a
SLOPE change in the curve, not a shift, and cannot be fixed by an overall
normalization multiplier.

Closed in `MakeMatrix`: the composition is normalized to its own sum
(proportions are preserved); if the deviation exceeds 1e−6, a `G4Exception`
is issued with the number; if the sum is non-positive, it fails. The guard
was verified by a run: on `ash` it prints "sum of mass fractions 1.02
instead of 1; composition normalized," on `soil` it stays silent.

**The numbers in `results/` did not change, and here is why.** The `ash`
matrix was not used in ANY run: grids are calculated on `water` and `soil`
(`drivers/run_grid.py`), background — on `air`, `organic`, `water`, `soil`
(`run_bg.py`), nuclides — on `water`, `soil`, `organic` (`run_nuc.py`),
attenuation coefficients — on `air`, `organic`, `soil`, `water`
(`geometry/mucalc.cc`). The defect sat in a branch reachable by the user but
not exercised by any published calculation. The fix protects future runs.

**Center of the sampling region in the manual macros** (`test.mac`,
`bench.mac`, `nuclides.mac`). It was `-0.44 mm`, correct is `-0.56`.
Analysis of the two gaps that the audit requested gave a simple answer:
there is no conflict, `seatGap = 0.32` refers to the m200 beaker, and 1.14
from `RCDetector.cc` refers to m500. The number in the macro lagged behind
the gap fix 0.20 → 0.32: the center of the cavity is 0.5·(z_top + z_bottom),
and at 0.20 it equals −0.44, at 0.32 it equals −0.56.

**Published results are not affected by this.** None of the drivers read
the manual macros: `run_grid.py`, `run_nuc.py`, and `run_bg.py` BUILD the
macro themselves from the `VESSELS` table, where the numbers have long been
correct — m200 (33.24; 33.25; −0.56) and m500 (43.34; 47.10; +2.06). The
manual macros serve for geometry checks and single runs.

Incidentally found and fixed along the way: the same three numbers were
outdated in the comments of `RCDetector.cc` (`35.61 / −33.69 / 32.81`
instead of `35.49 / −33.81 / 32.69`), and the macros did not state that
they are valid **only for m200** — for m500 the sampling region differs
substantially, and silently applying them to it would lose most of the
sample. The caveat has been added.

**Variance of the net count without `k²`** (`normalization.py`,
`fit_peak.py`, `validate_bgsub.py`). It was `N_sample + bg_scaled`, correct
is `N_sample + k²·bg_raw`, i.e. via the already-scaled background —
`N_sample + k·bg_scaled`. At `k < 1` the uncertainty was underestimated, and
the "model vs. certified measurement" cross-check showed better agreement
than actually exists. The multiplier is now carried as an explicit variable
in all three places.

**The 0.80 multiplier in `compare_lsrm.py`.** It printed a value declared
erroneous in `curves.py` itself (a multiplier applied to the grid point
where self-absorption at 662 keV came out to 1.0012 — clearly noise). One
number in two forms within a single report. Now it is taken via import
`from curves import K_NORM`; no literal remains.

### Checked and NOT confirmed

The remark that `analysis/read_rcxml.py` reads `MeasurementTime` as live
time in violation of the rule from `common/docs/pitfalls.md`. The rule is
correct, but it applies to the **BecqMoni** format. In the native RadiaCode
XML there is no `LiveTime` tag at all: there is `MeasurementTime`,
`ValidPulseCount`, and `TotalPulseCount`. In the parsed record `Valid =
Total = 255 743` at a rate of 2.2 counts/s — there are no losses, no
correction is needed, and there is nothing else to read. If dead time is
ever needed, it should be derived from the Valid/Total pair.

### Recalculated from primary sources

**`h*(10)` at 662 keV is calculated, not taken from memory**
(`analysis/analyze_sens.py`). The value 4.13 pSv·cm² that had been used was
labeled "ICRP 74 with interpolation," there was no way to verify the label,
and the conclusion about the factory sensitivity depends directly on this
number.

The quantity decomposes into two, known independently:

    h*(10)/Φ = [E · (μ_en/ρ)_air] · [h*(10)/K_air]

The first factor is air kerma per unit fluence — pure physics with a
tabulated NIST attenuation coefficient. The second is the ratio published
by ICRP 74 itself; for a Cs-137 field it equals 1.20 Sv/Gy and changes
slowly with energy. **Only this ratio** is taken from the reference.

The method was verified at a point where the tabulated value is well
known: at 1 MeV the calculation gives 5.23 pSv·cm² against the published
ICRP 74 value of **5.24** — agreement to within two tenths of a percent. At
662 keV it comes out to **3.73**, the fluence rate for 1 µSv/h equals
74.5 cm⁻²·s⁻¹.

The previous 4.13 overstated the coefficient by a factor of 1.107, i.e.
understated the calculated count rate. The calculated 20–24 counts/s
becomes **22–27 against the datasheet's 30**: protocol item 3 remains NOT
REPRODUCED. The dose coefficient explains part of the gap but does not
close it, and there is no more point waiting for a scan of the table — at
any reasonable value in the range 3.3…4.13 the count rate falls short of
thirty.

## What else is not done

- Protocol items 4, 8, 9 (see above).
- The beta grid and the Ra-226/Th-232/U-238 runs exist only for the 200 mL
  beaker; for the 500 mL there is only Cs-137, K-40, and Sr-90 in organic
  matter at 0.49.
- The background of the 200 mL beaker has not been recalculated with
  organic matter at 0.49.
- **RadiaCode 110 is not modeled:** crystal dimensions are needed, or
  paired background spectra of the 103 and 110 to reconstruct them. The 103
  unit is no longer in the fleet.

## What is here

```
geometry/   RCDetector.{hh,cc} — instrument and author's beakers, measured from STL
            main.cc, mucalc.cc (attenuation coefficients),
            wallfield.cc (natural background field spectrum in the room — a separate task)
macros/     grids, nuclides, background, reference run
drivers/    run_grid.py, run_bg.py, run_nuc.py
analysis/   curves.py (curve reduction), compare_lsrm.py, selfabs_fit.py,
            peaks.py, normalization.py, validate_*.py, plots.py etc.
reference/  LSRM curves for comparison (SpectraLine, .in and text exports)
results/    final tables, curves.json and figures; no raw spectra
            mu.csv          mass attenuation coefficients (mucalc)
            wallfield.csv   natural background field spectrum of the room (wallfield)
            m*/background/  background of empty beaker and with sample
            field_*.mac     field macros
```

Raw calculated spectra (about 400 files) are not included in the
repository — they are reproduced by the drivers in `drivers/`.

Some scripts work with the operator's PERSONAL measurements
(`G4MODELS_MEASURED`), which are not in the repository and never will be.
Such a script will say so directly and point out that its result is already
available ready-made in `results/`.

## Difference of the beakers from the Gamma-1S geometries

The beakers here are **custom, 3D-printed**, measured by laser scanning into
STL, not taken from the LSRM beaker table. The well follows the shape of the
instrument: a rounded rectangle with a 1.25 mm wall, not a cylinder. These
are different geometries from the 1 L Marinelli beaker used with the
Gamma-1S, and they must not be mixed — see
[geometries](../../geometries/README.en.md), which also contains the model
sources, the license, and a breakdown of what the measurement confirmed and
what remains an assumption.
