# Gamma-1S model in NuclideMaster format (LSRM)

Export of our Geant4 model to EffCalcMC/NuclideMaster format — for
independent cross-validation against someone else's Monte Carlo with the same layers.
Files were made following the operator's samples (`NaI63x63.din`, `Marinelli_0cm.sin`
from the NuclideMaster kit).

| file | what it describes |
|---|---|
| `G1S_NaI63x63.din` | УДС-Г-63×63 detector per `geometry/G1SDetector.hh` |
| `G1S_Marinelli_1L_OISN16.sin` | 1 L Marinelli with ОИСН-16 filler |

## Load check: Geometry Master, 28.07.2026

The detector file was loaded by the operator into Geometry Master (LSRM) — parsed and
rendered without errors, crystal volume 196 cm³ = π·3.15²·6.3. Layer cross-check:

| Geometry Master layer | thickness, cm | material | our layer (G1SDetector.hh) |
|---|---|---|---|
| Crystal Ø6.3×6.3 | — | NaI 3.667 | crystal Ø63×63 |
| Front/side reflector | 0.6 / 0.365 | MgO **2.0** | MgO 6.0 / 3.65 mm, bulk (poured) |
| Crystal packaging | 0.05 | Al 2.7 | sealed can 0.5 mm |
| Seal/gasket | 0.2 / 0.2 | ρ 0.92 | rubber shock absorber 2 mm |
| Cover | 0.15 | Al 2.7 | outer housing 1.5 mm |
| Bottom cup | 0.2 | — | seating (from LSRM sample) |

Display quirk: the program labeled the "Seal/gasket" as "Polypropy…" — it
picks the name from the composition itself. In the file it's natural
rubber C₅H₈ (H 0.118 / C 0.882, ρ 0.92); polypropylene would be H 0.144. Does not
affect the calculation.

The Marinelli file was loaded there as well: all dimensions were parsed correctly (Ø15 × 11,
hole Ø8 × 7.2, walls 0.2, source height 8.4, distance 0.15 cm),
**source volume 1.01·10³ cm³** vs. the target 1000 — an extra 1 % from
rounding the fill height 8.36 → 8.4 cm; the 0.4 mm layer above the sample affects
detection efficiency by essentially nothing. The program labeled the vessel material
"Polyethylene" — PP and PE have identical mass fractions of H/C (0.144), so by
composition they're indistinguishable; the density 0.90 in the file is that of polypropylene.

## What assumptions of the model the files carry inside them

- **MgO ρ = 2.0 g/cm³** — bulk (poured) powder, not sintered ceramic (3.58).
  A density sweep (28.07.2026) showed: across the whole plausible range
  1.3–2.0, detection efficiency at 59.5 keV shifts by only +9 % — less than a third
  of the observed soft-edge deficit.
- **Layers outside the cover** (0.5 mm air + 1.0 mm rubber protector) have no place
  in the format — replaced by the source distance 0.15 cm in `.sin`. The difference
  "rubber instead of air" costs ~1 % transmission at 60 keV.
- Marinelli fill level — from the target 1000 mL: 11.6 mm above the well bottom,
  full source height 8.4 cm.

## EffCalcMC calculation, 28.07.2026 — RESULT

A run on the pairing of these two files: distance 0.15 cm, grid of 50 points
(logarithmic, 30–3000 keV), 10⁸ histories. The output is right here:
`EffReg_G1S_Marinelli.efa` (raw EffCalcMC output) and
`EffReg_G1S_Marinelli.efr` (the same in format 1.7, as the operator saved it
for application to spectra). Cross-check — `analysis/compare_effcalcmc.py`.

**The our/EffCalcMC ratio is NOT flat, it cannot be reduced to a single number**
(the first draft did exactly that — "mean 1.000, slope −1.2 %" — and it
was a masking average: the soft-edge offset canceled the hard-edge rise; caught
by the auditor, checking at OUR nodes by interpolating the dense ECM curve). Structure:

    45.3–59.5 keV   ours is higher by 20–22 %      (three nodes, stat. 1.5–2.3 %)
    88–166 keV      ours is lower by 3–7 %
    above 90 keV    monotonic rise +7.5 % per decade, agreement ±3 % (RMS)

The main conclusion is not undermined by this structure, but reinforced by it: the
code-to-code discrepancy is ~10 % over the range and is OPPOSITE in sign to the
measured trend (measurement — the ratio falls with energy, codes — it rises). The
measured 18–52 % cannot be reproduced by choice of code: **the energy dependence
belongs to the DESCRIPTION OF THE LAYERS.**

A caveat on exactly what is confirmed: two codes on the SAME layers verify the
CODE, not the level. The level is verified by measurement, and for the Marinelli it
gives 0.780 ± 0.025 — the model overestimates the detection efficiency, and code
agreement does not refute this, it points to the layers.

Candidates for the structure below 100 keV (our/ECM 1.2 at the edge): filler
packing arrangement, cross-sections in ОИСН-16 (71 % iron), peak area
determination. The point geometry has none of these factors — which is why it was
the decisive step.

## Point-source EffCalcMC calculation, 28.07.2026 — FORK RESOLVED

**NOTE (28.07.2026): the cross-check below is stale.** `EffReg_G1S_Point5cm.efa`
was computed by EffCalcMC from the `.din` exported BEFORE the end-face and
reflector-density fix (the same `.din` as in the Marinelli comparison above).
The numbers below remain as a historical trace of the "code vs. code"
investigation; the up-to-date model-vs-measurement comparison is in
`REPORT.md` §2.4а and §5а, which also gives an honest account of the density
scan on the new geometry (the direction is right, the shape is not).

The operator assembled a POINT source (distance 5 cm, radial 0) and ran the same
grid. File — `EffReg_G1S_Point5cm.efa`, cross-check with `results/eff_p5cm.csv`
(also stale at the time of this note):

    45.3–59.5 keV   our/ECM = 1.017…1.024   (soft edge CONVERGED)
    88–166 keV      0.92…0.97
    all 22 nodes    mean 0.968, slope −0.4 % per decade, RMS 2.9 %
    E ≥ 90          +5.6 % per decade

Three conclusions:

1. **On a geometry where the energy dependence IS MEASURED, the two codes agree
   within an RMS of 3 %**, while the measured curve sits away from both by
   0.78 → 1.14 across the range. The dependence belongs to the description of the
   layers definitively; the "Marinelli vs. point" caveat is lifted.
2. **The soft structure of the Marinelli comparison was explained**: the +20 % at
   45–60 keV disappeared on the point geometry (1.02) — it belonged to the filler
   (self-absorption in ОИСН-16, fill packing arrangement), not the detector.
3. Inter-code difference of +5.6 %/decade above 90 keV (in the Marinelli, +7.5) —
   a Geant4-vs-EffCalcMC systematic, likely cross-sections; 6–8 times smaller than
   the sought effect.

**The explanation of the 0.968 level via a distance difference is WRONG and is withdrawn**
(caught by the auditor via the device passport/specification, §2.10: the distance is
measured "from the surface of the detector cover", while `Point5cm.sin` has
`Distance,cm=5` with no compensation). The argument rested on the assumption that
the end face includes the rubber protector (equivalent to 5.15 cm) — after the
geometry fix of 28.07.2026 (no protector on the end face, outer plane = Al cover)
this difference disappeared, and both models, ours and ECM's, measure from the same
physical plane. Where the 0.968 level then comes from remains an open question SPECIFICALLY for
the code-vs-code comparison (our Geant4 vs. EffCalcMC); it requires a fresh ECM
run against a re-exported `.din`, which does not exist yet. Separately from this
— the comparison of the model against MEASUREMENT on the new geometry has already
been done and does not depend on the `.din`: `REPORT.md` §2.4а gives p5cm
MC/exp = 1.094 (was 0.971 before the end-face fix), p25cm = 1.010 (was 0.928).
These are two different questions: "do the two codes agree" and "does the model
agree with the instrument" — the first awaits re-export, the second is already
answered.

Incidental confirmation: fitting this curve, Efficiency itself converged on
TWO zones (20–245 / 70–3010, degrees 4/4, χ² = 8.8) — matching our own scan via
sliding node exclusion (two zones 5/2).

## The 5.15 cm distance: why the dispute over it does not spoil the cross-check (the auditor's argument)

A distance error shifts the LEVEL, not the shape: the inverse-square factor does
not depend on energy. The slope check does not depend at all on the choice of
"5.00 or 5.15", and the free scale factor in the shape comparison absorbs it by
itself. The choice of 5.15 is physically justified (protector and air outside the
cover), but even if it were disputed — the conclusion on the slope would not
change. This also closes the objection "the distance was tuned to force agreement":
the distance can only be used to tune the level.

Converse corollary: since 5 and 25 cm agreed after the weight fix (1.179 vs.
1.184), there is no gross distance error at either distance — 1.5 mm gives
6 % at 5 cm and 0.6 % at 25 cm, a mismatch would be visible.

## Why this is needed (task 93)

Our model has an **energy dependence of the calculation/experiment ratio**
(the accepted term; colloquially "slope") against the measured LSRM curve:
underestimation at the soft edge (МК/эксп 0.78 at 59.5 keV) and overestimation
at the hard edge (1.14 at 2614.5), confirmed by two independent references
(source passports/specifications and the .efa curve). An EffCalcMC calculation
with THESE SAME layers discriminates between two hypotheses:

- their curve reproduces this dependence → the defect is in the DESCRIPTION of the
  layers (thicknesses/densities of the input end face, crystal dimensions) — to be
  sought in the drawing and the assumptions;
- they have no such dependence → the defect is in OUR code (tracking/sampling,
  area readout, solid angle) — to be sought in geant4-detector-models.

The decisive geometry is the point one (5 and 25 cm from the end face): there is
no vessel and no self-absorption, leaving only the detector and its layers. Compare
against `results/eff_p5cm.csv` / `eff_p25cm.csv`, column `eps_net`.
A point-source `.sin` is not available here yet — the POINT-block format is not
shown in the samples.
