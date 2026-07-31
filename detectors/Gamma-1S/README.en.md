# GAMMA-1S — NaI(Tl) 63×63 in lead shield

Gamma-1S scintillation gamma-ray energy spectrometer
(ZAO NPTs "ASPEKT", Dubna, DTsKI.412131.001): detection unit
UDS-GTs-63×63-USB in the "Ekran-1SG" shield-screen.

**Goal — spectrum analysis algorithms and their verification against a
certified reference kit.** Calibration verification, peak search, area
extraction, line selection by separation, coupled deconvolution of blends,
summing and pile-up corrections, and activity averaging by the LSRM rule are
all worked through here. The Geant4 model is a tool within the chain: it
supplies the detection efficiency. The measure is the reference kit: the
instrument has passport curves for each geometry, the sources have passports
with activities — that is, forty records with answers known in advance.

Detailed report with all the numbers and caveats — [REPORT.md](REPORT.en.md).
Spectra with found peaks, calibration verification and deconvolution insets —
on the [summary page](https://vibeengineering-llc.github.io/geant4-detector-models/gamma-1s/).

**Ready-made calculated curves — [results/](results/)**: nine full-energy-peak
efficiency curves (Marinelli, "Denta", Petri dish at two densities, water,
point sources at 5 and 25 cm), a run manifest, and cascade-summing
corrections. Building Geant4 to obtain the efficiency at a given energy is
not required; column descriptions and caveats are in
[results/README.md](results/README.en.md).

LSRM reference data are located in [reference/lsrm/](reference/lsrm/) in two
formats: the instrument's binary `.spe` and BecqMoni XML. They are
interchangeable — for 39 matching pairs the counts, times, and channel
numbers coincide exactly, verified by the script
[analysis/xml_vs_spe.py](analysis/xml_vs_spe.py).

## Status

The table below was synchronized with `results/` on 2026-07-31, after all
grids were recomputed on the corrected geometry; the previous revision of
this table rested on retracted values (the entrance-face fix of 28.07 and
the stale export, task 137).

| check | result | reference value | source |
|---|---|---|---|
| shield metal masses | Pb 167.1 / Cu 1.60 / Cd 1.58 kg | passport ≥165 / 1.6 / 1.2 | `geometry/G1SDetector.cc`, ReportMasses |
| point-source FEP, 25 cm, 662 keV | 0.121 ± 0.001 % | passport ≥0.1 % | `results/eff_p25cm.csv` |
| resolution at the 662 keV peak | 7.5 % | passport ≤8 % | `analysis/detector_params.py` |
| **point source 5 cm, calc/passport, 21 lines** | **1.0365 ± 0.0112** | χ²/ν = 1.79 | `results/kit_activity_point.csv` |
| **point source 25 cm, calc/passport, 7 lines** | **1.0381 ± 0.0214** | χ²/ν = 2.10; lid open | `results/kit_activity_point.csv` |
| Marinelli against `.efr`, 15 lines | 1.2526 ± 0.0120 | calculation overestimates | `results/compare_lsrm_summary.csv` |
| "Denta" against `.efr`, 13 points | 1.011 | agreement | `analysis/compare_cups.py` |
| Petri dish against `.efr`, 14 points | 0.941 | calculation underestimates | `analysis/compare_cups.py` |

Both point distances are now extracted with ONE peak-area convention (unified
on 2026-07-31) and agree with each other; the certified curve itself
reproduces the same passports with deviations of 3.9% at 5 cm and 4.4% at
25 cm — this sets the scale against which the calculation is meaningfully
compared. The jump of the discrepancy at the hard edge (2614.5 keV) belongs
to the peak-area extraction convention, not to the physics of the model;
after both sides are brought to one convention the residual is flat, about
7% (a lower bound) — details in [docs/report.md](docs/report.md) §5.3.

Independently of the `.efr`, **the entire kit was recalculated against the
source passports** — tables in
[results/kit_recalc_volume.csv](results/kit_recalc_volume.csv) and
[results/kit_recalc_point.csv](results/kit_recalc_point.csv). Poorly
separated lines are not fed into the activity calculation (`purity` column,
threshold 0.95); activity is computed by the LSRM rule, a weighted average
with weights 1/(ΔA)².

| geometry | Cs-137 | K-40 | Ra-226 | Th-232 |
|---|---|---|---|---|
| Marinelli 1 l | 0.782 ± 0.039 | 0.725 ± 0.074 | 0.788 ± 0.056 | 0.796 ± 0.048 |
| "Denta" 120 ml | 1.065 ± 0.055 | — | 1.149 ± 0.052 | 1.090 ± 0.068 |
| Petri dish 60 ml | 1.253 ± 0.065 | — | 1.352 ± 0.074 | 1.195 ± 0.072 |

In the Marinelli beaker all four nuclides agree with a single value of
≈0.78 within the uncertainties: the vessel has a SINGLE factor, not a
nuclide-dependent one. The sign is the same as in the comparison with
`.efr`, and differs between vessels — the discrepancy belongs to the vessel
models, not to a common normalization.

## Verification status per the common protocol

See [common/docs/validation.md](../../common/docs/validation.en.md).

| item | status |
|---|---|
| 1. Masses against passport | done: Pb 167.1 / Cu 1.60 / Cd 1.58 kg |
| 2. Convergence of drawing dimensions | done: both sums, 78.3 and 74.5 mm |
| 3. Passport efficiency point | done: 0.121 ± 0.001 % against the ≥0.1 % requirement |
| 4. Point geometry as detector reference | done twice: calc/passport 1.0365 ± 0.0112 and 1.0381 ± 0.0214 |
| 5. Density sweep, d_eff | done: two densities for each vessel |
| 6. Summing check on a nuclide without a cascade | done: Cs-137 0.983 ± 0.015, K-40 1.009 ± 0.025 |
| 7. Line yields from the model's own calculation | done: `*_emit.csv` of the same run |
| 8. Verification of the reference curve bypassing processing | done: ratio 0.999 |
| 9. Whether a correction was introduced into the reference curve | done: introduced, three independent tests |
| 10. MDA against passport | done, the passport falls on the Currie formula |
| 11. Cross-check of data formats | done: 39 XML/`.spe` pairs, zero difference |
| 12. Ready-made numbers as files | done: nine curves in [results/](results/) |

The only thing that remains open is not a protocol item, but the vessel-model
discrepancy, see below.

## What is known and what is open

Status as of 2026-07-31, after the entrance-face fix, the recomputation of
all seven grids, and the unification of peak-area conventions; details and
caveats in [docs/report.md](docs/report.md).

**The detector model agrees with the passports.** Both point distances,
extracted with one convention, give calc/passport 1.0365 ± 0.0112 and
1.0381 ± 0.0214 and agree with each other; the certified curve itself
reproduces the same passports with deviations of 3.9% and 4.4%.

**The jump of the discrepancy at the hard edge belongs to the peak-area
extraction method.** The convention correction `B(E)`, measured on clean
monoenergetic runs, absorbs the break at 2614.5 keV; what remains
unexplained is a flat plateau of about 7% (a lower bound — the calculation
uses a linear background against the second-degree polynomial of the
processing software). Decomposing the plateau is tasks 109/111/122/142.

**The response to a change of geometry was verified across seven grids**:
the agreement with the certification improves the farther the sample
geometry is from a point source (0.5σ on the Denta → Marinelli step, 4.1σ
on Petri → Denta, 13.5σ on point → Petri); the sign of the trend is set by
self-absorption in the sample (see the nuclide table above).

**The LSRM reference curve has been checked and is correct.** Reconstructed
from the raw spectrum of the calibration source and the passport activity:
1.8686·10⁻² against the recorded 1.8713·10⁻², ratio 0.999.

**The LSRM curve has already been corrected for cascade summing** — three
independent tests. Applying a correction to the comparison would be an error
of 10–20 %.

## Unfinished work

- Decomposition of the ~7% plateau into signed contributions and the
  acceptance criterion as a closed balance (tasks 109, 111); the
  peak-to-total ratio awaits the long Th-228 decay run (task 122); the
  depression of the Th-228 nodes in the certification is not yet separated
  between pile-up and cascade summing (task 142).
- Two estimates of the FWHM at 662 keV differ by 10% (the caesium record
  against the thorium session); the two values are deliberately kept
  separate (task 138).
- The alternative physics list (Livermore against option4) has not been
  cross-checked — it costs a full grid recomputation (tasks 101/141).

## Composition of the mixed source

The kit's source is named "Am-Ti-Eu-Cs"; its composition is not recorded in
the files. It was established from the spectra and confirmed by `.efr`
sections: **Am-241 + Ti-44/Sc-44 + Eu-152 + Cs-137**. "Ti" is indeed
titanium-44: it decays into Sc-44, which emits a positron and a 1157 keV
quantum, hence the strong 511 line (annihilation, not Na-22) and the sum
peak 511 + 1157 = 1668 keV. The matrix is RISN-379, composition taken from
the original `.spe` files, see [materials](../../materials/README.en.md).
