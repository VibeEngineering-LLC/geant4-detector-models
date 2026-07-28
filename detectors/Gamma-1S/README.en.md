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

Detailed report with all the numbers and caveats — [REPORT.md](REPORT.md).
Spectra with found peaks, calibration verification and deconvolution insets —
on the [summary page](https://vibeengineering-llc.github.io/geant4-detector-models/gamma-1s/).

**Ready-made calculated curves — [results/](results/)**: nine full-energy-peak
efficiency curves (Marinelli, "Denta", Petri dish at two densities, water,
point sources at 5 and 25 cm), a run manifest, and cascade-summing
corrections. Building Geant4 to obtain the efficiency at a given energy is
not required; column descriptions and caveats are in
[results/README.md](results/README.md).

LSRM reference data are located in [reference/lsrm/](reference/lsrm/) in two
formats: the instrument's binary `.spe` and BecqMoni XML. They are
interchangeable — for 39 matching pairs the counts, times, and channel
numbers coincide exactly, verified by the script
[analysis/xml_vs_spe.py](analysis/xml_vs_spe.py).

## Status

The model has been built and validated against five geometries of the
verification kit. **Two vessel-geometry defects have not been resolved** —
see below.

| check | result | reference value |
|---|---|---|
| shield metal masses | Pb 167.1 / Cu 1.60 / Cd 1.58 kg | passport ≥165 / 1.6 / 1.2 |
| point-source FEP, 25 cm, 662 keV | 0.116 ± 0.008 % | passport ≥0.1 % |
| resolution at the 662 keV peak | 7.5 % | passport ≤8 % |
| **point source 5 cm, 24 lines** | **0.971** | detector confirmed correct |
| **point source 25 cm, 20 lines** | **0.931** | confirmation, lid open |
| Marinelli, 15 lines | 1.165 | calculation OVERESTIMATES |
| "Denta", 13 lines | 0.886 | calculation UNDERESTIMATES |
| Petri dish, 14 lines | 0.811 | calculation UNDERESTIMATES |
| MDA (Currie) | Cs 1.55 / K 19.1 / Ra 2.84 / Th 3.88 Bq/kg | passport 1.5 / 25 / 3 / 3 |

The "Marinelli / Denta / Petri" rows were obtained by comparison with the
processed `.efr` curve. Independently of that, **the entire kit was
recalculated against the source passports** — 53 lines, tables in
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

See [common/docs/validation.md](../../common/docs/validation.md).

| item | status |
|---|---|
| 1. Masses against passport | done: Pb 167.1 / Cu 1.60 / Cd 1.58 kg |
| 2. Convergence of drawing dimensions | done: both sums, 78.3 and 74.5 mm |
| 3. Passport efficiency point | done: 0.116 ± 0.008 % against the ≥0.1 % requirement |
| 4. Point geometry as detector reference | done twice: 0.971 and 0.931 |
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

**The detector is confirmed correct.** Two point-source geometries in
different shield configurations give 0.97 and 0.93. This separates the
detector model from the vessel models.

**The vessel models diverge, and in opposite directions** — the Marinelli
beaker overestimates, the cuvettes underestimate. The decisive argument:
K-40, a clean isolated line with no chains or blends, gives 1.22 / 0.82 /
0.78 across the three vessels. The geometry does not know which nuclide is
inside, so the cause lies in the vessels themselves. The culprits are the
nodes where assumptions stand in for drawings: the depth of the Marinelli
well and the width of the cuvettes. No drawings exist for the cuvettes.

**The LSRM reference curve has been checked and is correct.** Reconstructed
from the raw spectrum of the calibration source and the passport activity:
1.8686·10⁻² against the recorded 1.8713·10⁻², ratio 0.999.

**The LSRM curve has already been corrected for cascade summing** — three
independent tests. Applying a correction to the comparison would be an error
of 10–20 %.

## Unfinished work

- Two targeted runs are prepared but not yet executed: the depth of the
  Marinelli well (74/65/55/45 mm) and the bulk density of the MgO reflector
  (1.3/1.5/2.0). Driver `drivers/run_probe.py`, analysis
  `analysis/probe_analyze.py`. The rejection criterion is specified in
  advance: the soft edge must be corrected without breaking the middle or
  the hard edge.
- The soft edge of the curve is underestimated (0.78 at 59.5 keV against
  ~0.95 in the middle). Three independent pieces of evidence point to
  excess absorber in front of the crystal; the accepted MgO density of
  2.0 g/cm³ is an assumption from the 1.5–2.4 range, and the Am-241 point
  favors 1.3–1.5.
- The activities of Ti-44 and Eu-152 in the mixed source were obtained by
  decomposition, but **are not reported**: their lines lie in the affected
  soft region, and until the reflector density is corrected, publishing
  them would be incorrect.
- Recalculation of the 25 point-source records of the kit has been written
  but not fully executed.

## Composition of the mixed source

The kit's source is named "Am-Ti-Eu-Cs"; its composition is not recorded in
the files. It was established from the spectra and confirmed by `.efr`
sections: **Am-241 + Ti-44/Sc-44 + Eu-152 + Cs-137**. "Ti" is indeed
titanium-44: it decays into Sc-44, which emits a positron and a 1157 keV
quantum, hence the strong 511 line (annihilation, not Na-22) and the sum
peak 511 + 1157 = 1668 keV. The matrix is RISN-379, composition taken from
the original `.spe` files, see [materials](../../materials/README.md).
