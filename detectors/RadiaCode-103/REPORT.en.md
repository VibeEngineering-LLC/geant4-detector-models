# radiacode-curves

Detection efficiency curves for **RadiaCode** detectors in Marinelli
beakers — a Geant4 11.2.1 calculation for determining isotope activity in a
sample, accounting for sample density, self-absorption, the beta-radiation
contribution, and background attenuation by the sample itself.

All original drawings, teardown photographs, and measurements are in
[`drawings/`](drawings/) and are open for use.

**Article with interactive graphs and an activity calculator:
[vibeengineering-llc.github.io/radiacode-curves](https://vibeengineering-llc.github.io/radiacode-curves/)**
— it also contains the full breakdown of checks, the uncertainty budget, and
the applicability limits. Page source —
[`docs/radiacode-103/article.tmpl.html`](../../docs/radiacode-103/article.tmpl.html),
built with `python analysis/build_article.py` (data substituted from
`results/curves.json`, figures embedded as data-URIs).

## What has been calculated

| | Instrument | Beaker | Status |
|---|---|---|---|
| Curves ε_p(E), 25 energies × 6 matrices | 101/102/103 | 200 and 500 mL | done |
| Reduction to formula: smooth limit × f(μρd) | 101/102/103 | 200 and 500 mL | done |
| Background and per-channel multiplier k, 5 matrices | 101/102/103 | 200 and 500 mL | done |
| Beta penetration, 16 energies × 5 matrices | 101/102/103 | 200 mL | done |
| Nuclides by full decay (K-40, Cs-137, Ra-226, Th-232, U-238) | 101/102/103 | 200 mL | done |
| Nuclides: Cs-137, K-40, Sr-90 in organic matter at 0.49 | 101/102/103 | 500 mL | done |
| Same as above | 110 | — | crystal dimensions needed |

## How to use this

```
eps_p(E, матрица, rho) = eps_lim(E) · f(mu(E)·rho·d) · 0.833
f = (1 − e^(−x)) / x,   x = mu(E)·rho·d,   d = 1.96 cm (200 mL) / 2.89 cm (500 mL)
A = S / (eps_p · p_gamma · t)
```

The form of the self-absorption multiplier is the accepted industry model:
it is exactly the one used in SpectraLine since version 1.4 in place of the
former form `exp(-mu*rho*d_eff)` (LSRM, "Algorithmic Foundations," §8.5.2),
and it is given for Marinelli beakers in Gilmore (ch. 10). The difference
here is only that the effective thickness is not taken from a table but is
fitted to the calculated points for the specific geometry: the tabulated
LSRM d_eff values are given for a 63×63 mm head (0.5 L Marinelli beaker —
15±2 mm), while the 1 cm³ crystal collects photons from a greater distance,
hence 2.89 cm.

The `eps_lim` coefficients (a polynomial in log E), the μ/ρ table, and
everything else are in [`results/curves.json`](results/curves.json),
computed by [`analysis/curves.py`](analysis/curves.py). Working range
**50–3000 keV**: below that lie the K-edges of iodine at 33.2 and cesium at
36.0 keV, where the curve has a discontinuity, and where the self-absorption
formula also stops being valid.

## Model checks

The calculation is absolute: not a single parameter was fitted to
measurements.

| Check | Calculation | Measured | Discrepancy |
|---|---|---|---|
| Sample volume, 200 mL beaker | 199.5 cm³ | 200.15 (STL measurement) | 0.3% |
| Sample volume, 500 mL beaker | 499.8 cm³ | 498.9 (STL measurement) | 0.2% |
| Room background H*(10) | 0.134 µSv/h | 0.10–0.20 (typical) | within range |
| Spectral hardness of background: ratio of count drop to dose drop in a lead castle | 4.25 | 4.16 | 2% |
| Peak areas: line convolution vs. full decay (662, 609, 352 keV) | — | — | 2–10% |
| **Empty-beaker background count rate** | 7.2 counts/s | 6.33 (RC-103, 7 days) | **model higher by 14%** |
| **ε_p(662) in 500 mL Marinelli beaker** | 3.85·10⁻⁴ | 3.21·10⁻⁴ (from real sample) | **model higher by a factor of 1.20** |
| Peak area 662: one sample with and without the lead castle | should be 1 | 1.006 ± 0.012 | 0.5 σ |
| Curve shape, 240–1000 keV | — | independent LSRM calculation | 5% |
| **Background distortion by the sample, k(E) by bands** | calculation | 3 datasets, no model | **1–3%** |
| Reduction of grid to formula, 800–3000 keV | residual 7.9% | point statistics 8.1% | noise only |

The last two rows are a systematic of one sign, and this must be taken into
account (see below). The lead-castle check compares the instrument against
itself, so it does not depend on either the crystal size or the field
normalization, and it checks purely the **shape** of the calculated field
spectrum.

### Normalization of the absolute scale: ε_p × 0.833 ± 0.03

The multiplier was chosen so that the **published smoothed curve**
reproduces the certified reference tie-point on the reference sample:
measured ε_p(662) = 3.2095·10⁻⁴ against calculated 3.8541·10⁻⁴. Consistency
between the curve and its normalization is mandatory, otherwise the
calculator would not reproduce the reference point on which the scale
itself is built. The multiplier is applied to both beakers: the response
deficit is a property of the crystal, not of the beaker.

**On light collection: it is not in the data.** Monte Carlo calculates
energy deposition, and only events with fully collected light end up in the
peak, so light collection should show up as a discrepancy between the
multipliers for the total count and for the peak area. This must be
checked within a SINGLE run, otherwise quantities with different statistics
are being compared (`analysis/normalization.py`):

| Quantity | Measured | Calculated | Ratio |
|---|---|---|---|
| Total sample count, counts/s | 1.8830 | 2.2486 | **0.837 ± 0.013** |
| Peak area 662, counts/s | 0.22441 | 0.28100 | **0.799 ± 0.035** |
| Peak/total — this is light collection itself | — | — | **0.954 ± 0.042** |

This is **1.0 σ** from unity: the data show no light-collection deficit,
the upper bound from this cross-check is 14%. The earlier conclusion of
"light collection 0.932" was a comparison error: the total count was taken
from the full-decay run, and the peak area from a mono-grid point where
self-absorption at 662 keV came out to 1.0012, i.e. a clearly noisy value
(absorption cannot increase efficiency).

So the entire systematic reduces to a **single** multiplier of 0.833,
identical for the peak and for the total count — and the total count gives
it independently of the curve, 0.837 ± 0.013. This corresponds to an
effective crystal volume of 0.84 cm³ instead of the nominal 1.0 — a
tolerance, a dead layer, or an unaccounted-for detail near the nose. The
cross-check points the same way: the LSRM calculation, whose crystal is
exactly 1.27 times smaller than mine, is higher than the measurement by
only 5.6% ± 6.7%, i.e. it practically coincides with it.

Reference point: dried blueberries, 246 g, in a 500 mL Marinelli beaker,
RC-103, 114121 s, in a lead castle. Reference activity 3340 Bq/kg from the
RadiaCode app, whose calibration for this geometry was cross-checked
against a certified MKS-AT1125A (4.39 vs. 4.37 kBq/kg on a different
sample, agreement 0.5%).

A direct cross-check on a second line did not succeed: K-40 at 1461 keV in
the same spectrum gave 37 ± 22 counts, i.e. at the detection limit.
Therefore the energy dependence was obtained by reasoning about the
mechanism, not measured.

**What the normalization does NOT affect:** the shape of the curve versus
energy, the self-absorption correction (d = 1.9 and 2.78 cm), the
per-channel multiplier k, the gamma/beta split. These are relative
quantities, checked by ratios.

Incidentally, regarding the specific RC-103 unit: the energy calibration is
shifted by 4.3 keV at 662 (the peak sits at 657.4), resolution 9.4% against
a datasheet value of 8.4%.

Separately, there is the lead-castle check: it compares the instrument
against itself, so it does not depend on either the crystal size or the
field normalization, and it checks purely the **shape** of the calculated
field spectrum.

The datasheet sensitivity of 30 counts/s per µSv/h (Cs-137) is not
reproduced by the model — the calculation gives 20–24. The discrepancy is
localized in the definition of the datasheet figure, not in the model: see
["On the datasheet's 30 counts/s"](#on-the-datasheets-30-countss).

### Cross-check with an independent LSRM calculation (BecqMoni)

The curves were cross-checked against the LSRM calculation from the
[BecqMoni](https://github.com/Am6er/BecqMoni) package — this is a second,
fully independent code that calculated the same instrument in the same two
author's Marinelli beakers. Copies of the models and exported curves with a
breakdown of the differences are in [`reference/`](reference/), the
cross-check — `analysis/compare_lsrm.py`. It must be compared against the
**water 1.00 g/cm³** configuration: that is exactly what the LSRM sample is.

![cross-check with LSRM](results/m500/figures/lsrm_compare.png)

**The curve shape is confirmed independently.** In the 240–1000 keV range,
where almost all the analytical lines lie (Pb-212 239, Pb-214 352, Bi-214
609, Cs-137 662, Ac-228 911), the ratio of the two curves, normalized to 1
at 662 keV, lies within **5%** — even though the curves were computed by
different codes for different geometries. This is exactly the quantity on
which activities computed from different lines depend.

**Absolute scale: the discrepancy is explained, not fitted.** My
calculation is higher than LSRM by a factor of 1.15 (0.2 L) and 1.12
(0.5 L) in the 100–1500 keV region. The expectation from the crystal
geometry alone is 1.27: for a convex body the direction-averaged projection
equals S/4, which for a cube is 150 mm², for their cylinder 118 mm², and
volume gives exactly the same ratio. The remainder (1.27 → 1.15) is due to
their well, Ø20 mm instead of rectangular, which brings the sample closer
to the crystal. At 80 keV the ratio grows to 1.38: their heavy reflector is
added (1 mm TiO₂ — that is 0.43 g/cm² against my 0.18).

**The direction of the normalization is confirmed by a second code.**
Rescaling the LSRM curve to the real sample's density (×1.114, via
self-absorption — the thing the models have in common):

| ε_p(662), 0.5 L Marinelli beaker, 0.49 g/cm³ sample | Value | vs. measurement |
|---|---|---|
| LSRM, rescaled from water to 0.49 | 3.54·10⁻⁴ | higher by 1.10 |
| this calculation | 4.11·10⁻⁴ | higher by 1.28 |
| **measurement (certified reference)** | **3.21·10⁻⁴** | — |

The independent calculation is **also higher than the measurement**, even
though its crystal is 21% smaller. So the remaining systematic is not the
crystal volume, but something absent from both models: light collection
into the peak. The 0.80 multiplier remains valid.

**Above 1.5 MeV the discrepancy is not resolved by the LSRM data.** The
ratio drifts to 0.72 (0.2 L) and 0.65 (0.5 L) at 2614 keV, but the
uncertainty there is 18–22% **for them** and 6–12% for mine, meaning their
curve sets the precision limit, and pushing my own statistics further is
pointless. Plus their export is strictly monotonic across 150 points at
that stated uncertainty — that is a smoothed fit, not calculated points.
For the Tl-208 2614 keV line (thorium) an explicit calculation is more
reliable.

## Instruments

RadiaCode 101/102/103 — case **123.0 × 34.0 × 17.5 mm**, scintillator
**CsI(Tl) 10×10×10 mm** on a SiPM. The geometry is the same for all three,
differing only in energy resolution (9.5% for the 102, 8.4% for the 103 at
662 keV), which is introduced by convolving the spectrum in
post-processing.

From the cross-sections
([longitudinal](drawings/rc101-103_case_side_section.png),
[plan](drawings/rc101-103_case_top_section.png)): the **center** of the
crystal is 12.00 mm from the outer plane of the nose, 8.20 mm from one wide
face and 9.30 from the other (sum exactly 17.5), centered exactly across
the width.

From the [teardown photographs](drawings/rc_teardown_overview.jpg): the
crystal sits in a white reflective cup with a window for the
photodetector, the cup is in a separate module at the nose, connected to
the board by a ribbon cable, meaning **the main board does not extend
under the crystal**; the Li-Po battery is marked 602560, i.e.
6.0 × 25 × 60 mm.

## Marinelli beakers

Measured by laser-scanning STL of the column mesh.

| | 200 mL | 500 mL |
|---|---|---|
| Inner radius | 33.24 mm | 43.34 mm |
| Sample height | 66.50 mm | 94.20 mm |
| End wall | 2.80 mm | 3.30 mm |
| Cup wall | 2.81 mm | 2–12 mm (barrel-shaped) |
| Well cavity | 34.70 × 18.14 mm | 36.68 × 19.78 mm |
| Well sleeve wall | 1.25 mm | 1.91 mm |
| Well depth | 47.81 mm | 65.60 mm |
| Crystal offset from sample center | 0.44 mm | 2.06 mm |
| Sample volume | 200.15 cm³ | 498.9 cm³ |
| Plastic (cup + lid) | 83.0 cm³ | 301.3 cm³ |

In both, the well follows the shape of the instrument — a rounded
rectangle, not a cylinder. On the 500 mL beaker, the barrel wall is thick
and variable, so **the empty beaker itself attenuates external background
by 10–17%**.

## How to calculate activity

**A [Bq] = (S − B·k) / (0.80 · ε_p · p_γ · t)**,  **a [Bq/kg] = A / m**

where S is the peak area in the sample, B is the area of the same window in
the blank run with the empty beaker, k is the background-attenuation
multiplier due to the sample, p_γ is the line yield, t is the count time, m
is the sample mass.

### Efficiency and self-absorption

Instead of a table by matrix — a formula with a single geometric parameter:

**ε_p(E, matrix, ρ) = ε_p^limit(E) · (1 − e^(−μd)) / (μd)**,  μ = ρ·(μ/ρ)(E)

For the 200 mL beaker, **d = 1.9 cm**. A fit over 125 points gives an RMS of
6.5% — this is within the statistics of the runs themselves, i.e. the
formula describes the data as accurately as it is known at all. Verified on
organic/soil/water matrices at densities 0.5–1.6 g/cm³ and energies
30–3000 keV.

ε_p^limit(E) — in
[`results/m200/efficiency.csv`](results/m200/efficiency.csv), a table by
natural-background lines — in
[`efficiency_table.en.md`](results/m200/efficiency_table.en.md); same for the
500 mL vessel — [`results/m500/efficiency.csv`](results/m500/efficiency.csv)
and [`efficiency_table.en.md`](results/m500/efficiency_table.en.md).

Self-absorption separates the curves **only below 150 keV**, but there
radically: at 30 keV, 1.6 g/cm³ soil gives four times less than the
limiting case. Above 300 keV all matrices converge within 15%. That is, for
Pb-210 (46.5) and Am-241 (59.5) density is essential, for Cs-137 and K-40 a
rough estimate is sufficient.

### Beta-radiation contribution

The path from the sample to the crystal — well sleeve, gap, case wall, and
reflective cup — totals **0.49 g/cm²** for the 200 mL beaker (0.57 for the
500 mL). By practical range, this corresponds to a threshold of
**Emax ≈ 1.15 MeV** (1.29 for the 500 mL).

Below the threshold there is still a signal — this is **bremsstrahlung**
from beta particles stopped in the sample; in this region the response
grows with the matrix's Z. Above the threshold, the electrons themselves
arrive, and there the response falls with density.

Fraction of charged particles in the total count per decay (1.6 g/cm³
soil):

| Nuclide | total/decay | β fraction |
|---|---|---|
| Cs-137 | 4.35·10⁻³ | **3%** |
| Ra-226 (equilibrium) | 1.50·10⁻² | 32% |
| U-238 (equilibrium) | 1.96·10⁻³ | 35% |
| K-40 | 5.65·10⁻⁴ | 42% |
| **Th-232 (equilibrium)** | 2.51·10⁻² | **50%** |

For the thorium chain, more than half the count comes from beta particles,
not gamma: it has three hard beta emitters per decay (Ac-228 2.08 MeV at
100% yield, Bi-212 2.25 MeV at 64%, Tl-208 1.80 MeV at 36%). This continuum
underlies all the peaks and degrades the detection limit. Cs-137 is clean
against this background.

### Correction for background distortion by the sample

The blank run is taken with an **empty** beaker, while under the sample the
background is already different. **The sign of the correction is not
universal**, and this is the main thing to understand here:

- in the **line window**, the sample attenuates the background line itself,
  k < 1, and subtracting the blank run overcorrects by 8–15%;
- in the **total count**, the opposite happens: scattering in the sample
  adds a soft continuum, and for light samples k > 1. For berries at
  0.49 g/cm³ in the 500 mL beaker, the total background count grows by 2%.

Line windows, 500 mL Marinelli beaker (`analysis/analyze_bg.py`):

| Line | organic 0.49 | water 1.0 | soil 1.2 | soil 1.6 |
|---|---|---|---|---|
| 351.9 | 0.882 | 0.851 | 0.872 | 0.842 |
| 609.3 | 0.894 | 0.870 | 0.831 | 0.783 |
| 1460.8 | 0.890 | 0.799 | 0.929 | 0.898 |

Per-channel — [`background_k.csv`](results/m500/background_k.csv).

**Verified by measurements, without a model**
(`analysis/validate_bgsub.py`). The operator has three datasets: sample
without a lead castle, empty beaker without a castle, and the same sample
in the castle, where background is suppressed by a factor of 19. The
sample's own emission is taken from the third one, and then k is computed
from measurements alone:

| Band, keV | k measured | k calculated | Discrepancy | Fraction of own emission |
|---|---|---|---|---|
| 20–100 | 1.069 ± 0.002 | 1.067 | +0.2% | 17% |
| 100–300 | 0.995 ± 0.002 | 0.989 | +0.6% | 18% |
| 300–600 | 0.968 ± 0.006 | 0.942 | +2.8% | 51% |
| 600–750 | 0.973 ± 0.022 | 0.894 | +8.8% | 81% |
| 750–1400 | 0.927 ± 0.008 | 0.916 | +1.2% | 5% |
| 1400–1550 | 0.901 ± 0.022 | 0.880 | +2.3% | 3% |
| 1550–2700 | 0.941 ± 0.017 | 0.953 | −1.3% | 0.3% |

Agreement 1–3% across the whole spectrum. Only the 600–750 keV band stands
out, and the reason is in the last column: there 81% of the count is the
cesium peak itself, and after subtracting it only crumbs remain of the
background. This is the most substantive check in the whole work: it
checks what otherwise cannot be checked — the calculation of the
natural-background field of the room plus photon transport through the
sample.

The trouble is that for natural-background samples, the background lines
are exactly the same lines being measured in the sample. For natural
potassium in soil in an open room, the correction comes out **larger than
the signal itself**, and without it K-40 cannot be measured at all. In the
lead castle, the background in the 1461 window drops by a factor of three,
and potassium becomes measurable.

## Background field

The field spectrum was not assigned from reference tables but calculated:
fluence in the air cavity inside concrete with natural background
radionuclides (K-40 400, Ra-226 40, Th-232 30 Bq/kg per UNSCEAR). This is
the standard approximation for the middle of a room surrounded by
enclosing structures over the full solid angle.

**Key result: the lines account for only 33% of the fluence, the
remaining 67% is scattered continuum**, concentrated below 400 keV, where
CsI has its efficiency maximum. A field assembled by hand from lines alone
would understate the count rate by a large factor.

## On the datasheet's 30 counts/s

The model gives 20–24 counts/s per µSv/h for Cs-137 against the
datasheet's 30, i.e. here, conversely, it **understates**. The direction is
opposite to the systematic seen for background and for the Marinelli
beaker, so it cannot be explained by a single efficiency error — the
response is smooth in energy. Likely causes of the discrepancy lie in the
definition of the datasheet figure:

- if it is referenced to **air kerma** rather than to H*(10), the
  calculation gives 26.7;
- calibration with a point source **in an ordinary room** adds scattering
  from the walls and floor to the count, whereas the dose in the
  denominator is calculated from the direct beam only.

The datasheet figure should not be used as a reference: the conditions of
its measurement are unknown. The normalization above is built on a
cross-check with a certified instrument.

## Build and run

```bash
cmake -S sim -B sim/build -G Ninja && cmake --build sim/build
```

```bash
python drivers/run_grid.py gamma m200
```

```bash
python analysis/analyze_gamma.py m200 && python analysis/plots.py m200
```

`rc_curves <macro> [full|bare|empty] [matrix] [density] [m200|m500]`.
Matrices: water, soil, ash, organic, air. Individual targets: `wallfield` —
room field spectrum, `mucalc` — matrix attenuation coefficients.

## Assumptions

Marked in the code with the word ДОПУЩЕНИЕ (Russian for "assumption"). Main
ones: instrument case material (ABS), composition of the white reflective
cup, battery surrogate, nose wall thickness, board and display dimensions.
Not modeled: beaker threads and ribs, the black detector module around the
cup, the well's entrance chamfer. The outer surface of the 500 mL beaker is
not axially symmetric (ribs) and is replaced by a median circumferential
contour — agrees in plastic volume to within 2.7%.

## Structure

```
drawings/          drawings, teardown photos, beaker photos and renders
sim/               Geant4 model, run drivers, post-processing
reference/         independent LSRM calculation (BecqMoni) + breakdown of model differences
docs/              interactive article (template + built page)
results/curves.json  working-curve coefficients: limit, d, mu/rho
results/m200/      200 mL beaker: spectra, curves, plots
results/m500/      500 mL beaker
```

Post-processing and checks:

| Script | What it does |
|---|---|
| `run_grid.py`, `run_bg.py`, `run_nuc.py` | runs: energy grid, background field, full decay |
| `analyze_gamma.py`, `analyze_beta.py`, `analyze_nuc.py` | curves, beta contribution, gamma/beta separation |
| `analyze_bg.py` | per-channel k(E) and recount in line windows |
| `validate_bgsub.py` | **cross-check of k(E) against measurements** (three datasets, no model) |
| `compare_lsrm.py` | **cross-check of curves against independent LSRM calculation** |
| `curves.py` | reduction of grid to formula; export of `results/curves.json` |
| `build_article.py` | build of the interactive article in `docs/` |
| `normalization.py`, `fit_peak.py`, `peaks.py` | scale normalization, peak areas, activities |
| `read_rcxml.py` | reading of RadiaCode/AtomSpectra spectra |

All drivers and analyzers accept the beaker as an argument: `python
run_bg.py m500`.
