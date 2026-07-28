# GAMMA-1S: developing analysis algorithms against a certified verification kit

## What this does and why

**The goal of this work is gamma-spectrum analysis algorithms and their verification, not a model of the instrument.** What is being developed and tested: checking the energy calibration against found centroids, peak search, area extraction with continuum and background subtraction, line selection by separation, coupled deconvolution of blends, corrections for cascade summing and for pulse-pileup losses, activity averaging over lines by the LSRM rule.

The Geant4 calculation is **a tool inside this chain**, not its goal: it supplies the detection efficiency, without which activity cannot be obtained from a peak area. The advantage of the model is that any intermediate quantity can be taken from it — the line yield, a neighbor's contribution to a blend, the coincidence fraction — without taking anything from a reference table.

**The certified verification kit serves as the measure, and that is the main point here.** The spectrometer has efficiency curves for each geometry, recorded during verification, and each source in the kit has a passport with an activity and its uncertainty. So the algorithm is checked not by eye but by running it against forty records with answers known in advance. There are three checks, each with a different weak point:

1. **curve against curve** — computed efficiency against the LSRM passport fit (§2.2, §2.3). Weak point: our errors add to someone else's, and LSRM's method of area extraction is unknown;
2. **kit recalculation** — activity is obtained from a measured spectrum and checked against the source's passport ("Recalculating the whole kit" in §2). Nothing to fit, but the whole chain is needed at once;
3. **algorithm against algorithm** — blend deconvolution against windowed extraction on the same lines ("Coupled deconvolution" in §2). It checks nothing external, but it separates algorithm error from model error.

Kit spectra with found peaks, calibration checks, and deconvolution insets for each group are on the summary page
<https://vibeengineering-llc.github.io/geant4-detector-models/gamma-1s/>.

Below is everything that agreed and everything that did not, with an indication of exactly where.

## The model everything is computed on

The GAMMA-1S unit (ZAO NPTs "ASPEKT", Dubna, DTsKI.412131.001): detection assembly UDS-GTs-63x63-USB in the lead shield "Ekran-1SG" DTsKI.305179.038, 1 L Marinelli beaker.

The calculation was carried out in Geant4 11.2.1, physics list `EmStandardPhysics_option4` (with fluorescence — needed for the X-ray emission of lead, cadmium, and copper) plus `RadioactiveDecay`. Secondary-particle production threshold 0.05 mm.

**The computed energy range is set by the range of the passport curves.** The energy grid must cover the domain of the LSRM fit in every geometry — otherwise there is nothing to compare at the edges. The zones are declared as `Zone_i` blocks in the `.efa` files themselves; their union is 45.3…3552.5 keV, and the grid edges are set exactly to that. The remaining nodes are placed on the actual lines of the kit's reference sources, so that the comparison happens at nodes rather than through interpolation. The list is in `drivers/grid_energies.py`, one file shared by all drivers.

---

## 1. Where every dimension comes from

The rule the calculation adhered to: **no number is entered "by common sense" — every one has a source, and wherever there is no source the word ASSUMPTION is written.** Numerical data — line yields, branching ratios, attenuation coefficients — are not taken from reference tables from memory, but computed by the same Geant4 that does the transport.

### 1.1 The detector — from the technical drawing

The drawing "BDS-G, UDS-G-63x63" labels the layers and three overall dimensions: Ø78, head length 74, and total length 315, plus Ø71 and the MgO thicknesses — 3.65 radially and 6 at the entrance face.

The labels allowed two readings: whether the 1 mm rubber and the 0.5 mm air gap are included in the radial stack. This was resolved arithmetically — both sums only agree for one variant:

| Direction | Stack from crystal outward | Sum | On drawing |
|---|---|---|---|
| radial | MgO 3.65 → Al 0.5 → rubber 2 → Al 1.5 | Ø78.3 | Ø78 |
| axial (end face) | MgO 6 → Al 0.5 → rubber 2 → Al 1.5 → air 0.5 → rubber 1 | 74.5 | 74 |

So the 1 mm rubber and the air gap are not present radially — they occur only at the entrance face, and the outer housing is 1.5 mm aluminum. This also explains Ø71: it is both the PMT envelope and the outer diameter of the crystal can (63 + 2·(3.65+0.5) = 71.3) — the can and the photomultiplier are butted together.

The drawing was checked for scale accuracy by measuring the raster at 1200 dpi from the fraction of dark pixels per row: crystal diameter 90.42 pt with the centerline in the middle, i.e., 1.435 pt/mm, and the layer boundaries fall on the labeled thicknesses.

Detector ASSUMPTIONS: MgO reflector bulk/packing density 2.0 g/cm³ (not sintered ceramic 3.58; the usual range for such assemblies is 1.5–2.4), 5 mm glass light guide, 1.5 mm PMT envelope wall, 120 mm envelope length, the voltage-divider block and PCB behind the envelope replaced with a homogeneous mixture of fiberglass-laminate and copper at a density of 0.5 g/cm³.

### 1.2 The shield — from the metal masses in the device passport

Composition specified by the operator: steel 3 mm, lead 50 mm, cadmium, copper. Layer order from the cavity outward — Cu → Cd → Pb → steel: this is a standard graded shield, the 72–88 keV lead X-rays are absorbed by cadmium, the 23 keV cadmium X-rays by copper, and the 8 keV copper X-rays are already below the working range (passport: from 50 keV).

**The cavity size is not given in any available document** — the manual refers to a separate operating manual (RE) for the "Ekran-1SG" shield, which is not available. The size was reconstructed from Table 2.2 of the passport, "Non-ferrous metal content, kg, not less than":

| metal | model | passport |
|---|---|---|
| lead | 167.1 kg | not less than 165 |
| copper | 1.60 kg | not less than 1.6 |
| cadmium | 1.58 kg | not less than 1.2 |

A cavity of Ø200 × 190 mm with 1 mm liners gives all three masses at once, with copper landing exactly on the passport value. The cadmium and copper thicknesses were not specified by the operator and were taken as 1 mm each (Budyka: Cd 0.5–1, Cu 1–2) — and the passport copper mass confirms exactly this choice.

### 1.3 The Marinelli beaker — from the LSRM cuvette table

External dimensions Ø150 × 110 — from the table of measurement cuvettes, LSRM "Precision Measurements," p. 11. Well Ø80 for the Ø78.3 head, 2 mm polypropylene wall. The vessel is larger than the nominal volume, so a 1000 mL sample does not fill it to the top: the fill level is computed from the target volume, with air above the sample up to the lid. The constructed solid yields exactly 1000.0 cm³.

The OISN-16 matrix — mass composition from the LSRM `.efa`: H 0.022, C 0.206, N 0.009, O 0.049, **Fe 0.714**, ρ = 1.6. This is a bulk-sample simulant on an iron base, neither soil nor organic matter. For the MDA calculation a second matrix is used — distilled water, as in the passport — and it is computed directly, not by rescaling from OISN-16: transferring between a 71% iron matrix and water using a self-absorption formula would be an extrapolation.

---

## 2. Model checks

None of the checks is a fit: in every case the calculation was compared against a number that was not built into it.

### 2.1 Passport efficiency for a point source

Passport, item 2.10: FEP efficiency at the 662 keV line for a point source at 25 cm from the surface of the detector lid — not less than 0.1%.

Calculated: **0.116 ± 0.008%**. The intrinsic efficiency that follows from this, 0.32, is the textbook value for a 63×63 NaI crystal at this energy.

Separately: with the shield lid CLOSED, the count is zero. This is not a glitch but physics — 50 mm of lead attenuates 662 keV by roughly a factor of five hundred. The passport 25 cm geometry is measured with the lid open, which is also confirmed by the name of the LSRM kit's background file (`open_lid_point25cm`).

**Passport limits (items 2.2, 2.11), checked against a scan dated 2026-07-28.** The registered energy range is **50…3000 keV**; the measured point-source curve at 5 cm starts at **59.5 keV (Am-241)** and does not go lower — the passport floor is respected. The model-grid nodes at 45.3 and 56.1 keV were therefore only compared against another code (EffCalcMC, §5a), not against measurement: there is nothing to compare against there — these points are a model extrapolation beyond the passport range, not a verified result.

The maximum permissible error of the FEP efficiency (item 2.11) is **±10%**, for a source certified to no worse than ±5%. This is a passport CEILING, not the actual accuracy of a given calibration: the real `dp` values in the SN-01 instrument's curve are 4.63% at 59.5 keV, 3.88% at 88.0, 1.58% at 661.7 (maximum over the whole curve 8.28% at 238.6). The observed model deficit at the soft edge (22–27% at 59.5–88 keV, §5a) is several times larger than both the actual measurement uncertainty and the passport ceiling — it cannot be attributed to calibration uncertainty.

### 2.2 Comparison against the measured efficiency curve

Experiment — file `УДС-ГЦ-63х63-USB__SN-01_-_Маринелли.efr`, the FEP curve built by LSRM from certified sources, 15 points, 238–2615 keV.

Calculation — a grid of 24 mono-energies, 45.3…3552.5 keV, 400 thousand events each, sampled over the whole sample volume. The FEP efficiency is extracted with a ±6 keV window, subtracting the continuum from the left shelf E₀−30…E₀−10 (to the right of E₀, for a monoenergetic source with no pileup, there are no counts).

**Weighted-average MC/experiment = 1.165 ± 0.011, equivalent to K_NORM = 0.858.** Without continuum subtraction the result is 1.196, i.e., the result depends only weakly on the peak-extraction method.

The main point here is not the number itself but the comparison: in the previous project the same coefficient for RadiaCode 103 in the author's own Marinelli beaker came out to **0.833**. Different instruments — a 1 cm³ CsI crystal versus a 196 cm³ NaI crystal — different vessels, different matrices, different laboratories, different years. The excess of calculation over measurement agreed to within 3%. So a systematic factor of about 1.2 belongs not to the model of the specific instrument, but to the method itself: Monte Carlo scores energy deposition, while the experiment measures peak area.

The shape of the curve is reproduced worse than the normalization: after normalization the scatter is RMS 11.3%, χ²/dof = 9.8. The outliers are examined in section 3.

### 2.3 Independent check bypassing .efr

The `.efr` curve is an already processed LSRM result. The same answer was obtained bypassing it: the peak area taken directly from the kit spectrum (BecqMoni XML), minus the measured background of the same geometry, divided by the source activity from the passport, recalculated to the measurement date and to the line yield.

Cs-137, Marinelli, activity on the measurement date 1025 Bq:

| window | recovered ε | vs. .efr 1.871·10⁻² |
|---|---|---|
| ±1.00 FWHM | 1.901·10⁻² | +1.6% |
| ±1.25 FWHM | 1.988·10⁻² | +6.3% |
| ±1.50 FWHM | 2.024·10⁻² | +8.1% |

Spectrum reading, times, and calibrations are correct. Along the way, the resolution was measured from the peak itself: **FWHM 7.5% at 662 keV** against the passport's "not more than 8%."

This also reveals the main systematic of this procedure — window width, up to 8%. A Gaussian cannot account for that much (±1.0 FWHM contains 98.2% of the area, ±1.5 — 99.98%), so the remainder comes from under-subtraction of the baseline at the wide window.

### 2.4 Effective sample thickness

The dependence of efficiency on density is reduced to the standard industry form f(x) = (1 − e⁻ˣ)/x, x = μ(E)·ρ·d, where d is the only parameter. The ratio of two densities does not depend on the absolute efficiency, so d is determined independently at each point, and the agreement between points is a check, not a fit.

**d = 31.5 mm, Δχ² = 1 band: 31.0…32.5, χ²/dof = 0.91 over 20 points.**

LSRM recorded **Thick = 31 ± 2 mm** in the `.efa` for this specific detector unit — an exact match. The tabulated value for a typical assembly (26 ± 2 mm, "Precision Measurements," p. 11) fits worse; this makes sense, as the passport of the specific instrument is more accurate than the table.

A χ²/dof near unity means that a single thickness describes the whole 45.3…3552.5 keV range, i.e., the shape of the correction itself is correct.

The matrix attenuation coefficient is computed by a separate program (`mucalc`) using the same physics and the same mixture as the transport. A manually entered XCOM table is kept only for cross-checking: the discrepancy is −2.1…+3.4%.

### 2.4a Point geometry: the DETECTOR model separated from the sample geometry

A point source at 5 cm is the cleanest reference: there is no vessel, no self-absorption in a sample, leaving only the crystal, reflector, housing, and shield. LSRM's range here is the widest: 24 lines from 59.5 to 2614.5 keV. The calculation uses a grid with a cone around the direction toward the detector (the efficiency is recovered by dividing by the solid-angle fraction, 0.25); direct at 8 nodes, interpolated elsewhere with a 5th-degree polynomial over its own 20 nodes (fit χ²/dof 2.68).

**Weighted-average MC/experiment = 0.971 over 24 points** at 5 cm and **0.931 over 20 points** at 25 cm.

| geometry | MC/exp | points | shape RMS | χ²/dof |
|---|---|---|---|---|
| point 5 cm, lid CLOSED | 0.971 | 24 | 8.9% | 9.0 |
| point 25 cm, lid OPEN | 0.928 | 20 | 6.1% | 4.0 |
| Marinelli 1 L | 1.165 | 15 | 11.3% | 9.8 |

The two point geometries agree with each other within 4% and both sit near unity, in DIFFERENT shield configurations — with the lid closed and with the lid open. So the detector model is verified twice, independently, while the Marinelli result is off by 20–26%.

This forces us to RETRACT the conclusion drawn in section 2.2 from a single geometry. There, the 1.165-fold excess of calculation over measurement was attributed to the approach itself — "Monte Carlo scores energy deposition, the experiment scores peak area." If that were so, the point geometry would give the same excess. It gives unity instead. So the discrepancy belongs to the **volume-source geometry**, not the detector, and the agreement with RadiaCode (0.833 vs. 0.858) is likewise a property of the Marinelli setup, not a universal systematic of the method.

The prime suspect is the Marinelli beaker: its overall dimensions are taken from a table giving two numbers (Ø150, H = 110), while the 2 mm wall, the Ø80 well, and its 74 mm depth are ASSUMPTIONS. If the sample actually sits farther from the crystal (a shallower well or a thicker wall), the computed efficiency drops. This is to be checked with a direct run at several well depths — a cheap calculation for one or two energies, which is planned.

### Recalculating the whole kit: 53 lines against the passports

The comparison above relies on the `.efr` curve — an already processed LSRM result. Recalculating the kit relies on nothing but the sources' passports: the count rate in a peak is taken from the measured spectrum, divided by the computed efficiency per decay, yielding an activity, which is compared to the passport activity referred to the measurement date. There is nothing to fit here.

Tables — `results/kit_recalc_volume.csv` (22 rows) and `results/kit_recalc_point.csv` (31 rows). `ratio` = A_measured / A_passport; below unity means the model OVERSTATES the efficiency.

**Lines for the activity are selected, not all taken at once.** A poorly resolved line does not qualify: its area is the sum of several transitions and cannot be attributed to just one. The measure is the fraction of the yield within the window that belongs to the line itself, from the emission spectrum of the same run (the `purity` column); threshold 0.95. Rejected: Ra-226 1120.3 and 1764.5, Th-232 238.6, 583.2, and **911.2 (purity 0.48: nearby Ac-228 968 — 29%, 964 — 9%, Tl-208 860 — 8%)**. The 911 keV line is usable only as material for developing the deconvolution. The threshold was cross-checked against the Sparrow limit of 0.85·FWHM — both measures give the same set.

**The activity is computed by the LSRM rule** ["Algorithmic Foundations of SpectraLine" §5; the "Activity in Counting Samples" methodology]: a weighted average over lines with weights 1/(ΔA)², with the uncertainty taken as the maximum of the weighted-average uncertainty and the scatter. The median, which was used here previously, satisfied no rule: it gives equal weight to a precise line and a crude one.

| geometry | Cs-137 | K-40 | Ra-226 | Th-232 |
|---|---|---|---|---|
| Marinelli 1 L | 0.782 ± 0.039 | 0.725 ± 0.074 | 0.788 ± 0.056 | 0.796 ± 0.048 |
| "Denta" 120 mL | 1.065 ± 0.055 | — | 1.149 ± 0.052 | 1.090 ± 0.068 |
| Petri dish 60 mL | 1.253 ± 0.065 | — | 1.352 ± 0.074 | 1.195 ± 0.072 |

**This is the decisive result of the section.** In the Marinelli beaker all four nuclides agree with ONE value ≈0.78 within their uncertainties. So the vessel has a single multiplier, not a nuclide-dependent one — and a single cause per vessel is what should be sought. By the same rule applied to the ratios across nuclides: Marinelli **0.780 ± 0.025**, "Denta" **1.105 ± 0.033**, Petri dish **1.264 ± 0.044**.

The summary numbers are in `results/kit_activity_volume.csv` and `kit_activity_point.csv` and are taken by the report and the page FROM THERE. Previously the page combined the same table of lines with the median, while the report used the LSRM rule, so one and the same recalculation carried two different numbers; now there is one rule and one source.

The point geometries are combined with the SAME rule and the same purity selection: **1.179 ± 0.019** over 21 qualifying lines at 5 cm and **1.184 ± 0.048** over 7 lines at 25 cm. Three of twenty-four lines were rejected: Ba-133 356.0 (purity 0.87, nearby 384 keV — 13% of the window yield), Co-57 122.1 (0.89, nearby 136 — 11%), Eu-152 964.1 (0.91). At 25 cm the selection removed nothing beyond what was already excluded due to contaminated shelves.

**The combination proceeds in two steps, and this is not a formality.** The first version combined all qualifying lines of a geometry into a single weighted average with PURELY STATISTICAL weights: the passport uncertainty of the kit's activities (2–5%) did not enter at all. Within a single nuclide this is legitimate — the passport uncertainty enters all of its lines equally and does not affect their relative weight. Across nuclides it is independent for each source and is the DOMINANT term there: statistics give 0.1–0.6% per line, the passport gives percent-level uncertainty. The result got worse the stronger the line: at 25 cm the combined value came out at 1.093, i.e., below five of the seven values — effectively just Cs-137 661.7 (±0.002) with an add-on.

The correct order: each nuclide's series is combined by statistics first, its passport share is added in quadrature, and only then are the series combined with each other. After this, **the two distances agree with each other** (1.179 and 1.184, versus the earlier 1.183 and 1.093) — in itself this is an argument for the new order: the distance to the source cannot change the ratio to the passport, and previously it only changed because of the weights.

Next to the combined value, **χ²/ν — the measure of agreement across the series** — is printed, and it says honestly that the set cannot be reduced to a single number: 4.0 at 5 cm and 7.9 at 25 cm (for the vessels: 0.23 / 0.64 / 1.17, where combining is legitimate). The LSRM rule "maximum of the weighted estimate and the scatter" is exactly the Birge expansion — Δ_avg = Δ_weighted·√(χ²/ν) — but the expansion only rescues you when the weights are right, not when they are understated.

The reason for the disagreement has been identified: the **energy dependence of the calculation/experiment ratio** (hereafter "the slope"; the associated number is the slope of the linear regression of ln R against ln E). At 25 cm the ratio falls monotonically — 1.41 at 88 keV, 1.24 at 121.8, 1.14 at 238.6, 1.04 at 661.7, and 0.88 at 2614.5. That is, the model underestimates the efficiency at the soft edge and overestimates it at the hard edge: this is a distortion of the SHAPE of the curve, not a shift of its level, and it is precisely the shape that matters to anyone using the curve to recover activities — a level shift is absorbed into an overall multiplier, whereas the energy dependence spreads apart the lines of one nuclide. At 5 cm the same trend is twice as weak (1.29 at 59.5 → about 1.15 at the hard lines). The difference between the two distances points to the 5→25 cm calculation path — the mono-energy grid with the cone and the solid angle — rather than to the lines themselves.

The purity SELECTION itself barely shifted the result (at 5 cm, the third significant digit). This is expected and serves as a check: at 5 cm the efficiency is taken from the decay run with the SAME window, so a neighbor's contribution enters both sides of the ratio and cancels out. A line is excluded not because it gives a wrong number, but because a group's area cannot be attributed to a single transition — one rule applies to the whole kit.

Purity is computed from the emission spectrum of its OWN run, and the same line can be usable in one geometry and unusable in another: Th-228 238.6 and 583.2 are clean for the point source (0.99 and 1.00), but give 0.89 in the Marinelli beaker — there the source is Th-232, and above Th-228 sits Ac-228 with its own 209 and 562 keV lines, which the point source does not have at all. It is the sources that differ, not the thresholds.

### Coupled deconvolution: blends resolved by fitting

The purity selection left Th-232 with just one line — 2614.5 keV, with no margin. Lines 583.2 and 911.2 were dropped not because they are bad, but because windowed extraction cannot separate them. Fitting can: a group of lines has a single free parameter — activity — with the areas tied as S_k = A·I_k·ε_k·t [LSRM, "Algorithmic Foundations," formula 5.2-7 in its limiting form]. Implementation — `analysis/deconv.py`, table — `results/deconv_lines.csv`.

Three decisions everything rests on:

- **normalization — by a second, identical fit** to the model decay spectrum broadened to the instrument's resolution: A·t = N_decay·S_meas/S_model. This ratio cancels out the yield, the FEP efficiency, cascade summing, and the neighbor contribution. The first version took ε from the mono-energy grid — that grid knows nothing about coincidences and overstated ε by 12% for Tl-208 2614.5 in the Marinelli beaker; this produced a 20% discrepancy on a single line, where deconvolution has nothing to do;
- **FWHM — from an instrument calibration** FWHM² = −495 + 4.40·E (keV²), taken from three strong single lines in the Marinelli beaker (662, 1461, 2614.5). The square-root law from a single reference point gives a biased value: at 583 keV it overstates (46.8 vs. 40.4), at 2614 it understates (99 vs. 105). Control — the 238.6 line, which was not part of the calibration: 23.8 vs. the model's 23.6. The window is nearly indifferent to a width error, but it cost the fit 9%;
- **neighbors from just past the edge of the fit region** (within 3σ) are included in the group model: without them, the fit of 351.9 landed on the shoulder of the neighboring 295.2 line — χ²/dof 14–39 and an overstatement of up to 24%; with them, χ²/dof 0.8–1.4.

The numerical machinery — weights 1/√counts, non-negativity bounds (lsq_linear), covariance via SVD — is taken from `gamma/peaks/coupled_multiplet.py` in the gamma-spectrum-analysis package; the tie via ε·t is ours (in their version the areas are tied by bare intensities, which does not give becquerels and does not account for summing).

**The normalization check passes**: on single lines the deconvolution matches windowed extraction — 662: 0.773 vs. 0.775; 1461: 0.712 vs. 0.720; 2614.5 for the Petri dish 1.196 vs. 1.195. **What remains open** is the systematic excess of groups 583 and 911 over the single 2614.5 line of the same nuclide (+12…26% across vessels) — this is a mismatch between the model's and the measurement's continuum shape in the middle of the spectrum, not the weights and not the normalization (checked, including by masking 511 keV). Until this is resolved, activity is published from windowed extraction with clean-line selection, and the deconvolution numbers remain method development.

#### Exactly where the discrepancy goes: peak-to-continuum balance

Here, activity is a ratio of amplitudes of ONE shape, taken from the measurement and from the broadened model. The ratio is correct exactly to the extent that the shapes match. So the check is set up directly: in each window, we compute what fraction of the total area the fit assigns to the peaks and what fraction to the baseline and step — separately for the measurement and for the model (`analysis/deconv_balance.py`).

| geometry | line | peak fraction, measurement | peak fraction, model | meas/model |
|---|---|---|---|---|
| Marinelli 1 L | 583.2 | 47.0% | 46.0% | 1.021 |
| Marinelli 1 L | 911.2 | 68.9% | 68.2% | 1.011 |
| Marinelli 1 L | 2614.5 | 69.8% | 73.1% | 0.955 |
| "Denta" 120 mL | 583.2 | 52.3% | 46.4% | 1.126 |
| "Denta" 120 mL | 911.2 | 70.4% | 64.1% | 1.099 |
| "Denta" 120 mL | 2614.5 | 73.4% | 75.6% | 0.971 |
| Petri dish 60 mL | 583.2 | 53.1% | 45.9% | 1.156 |
| Petri dish 60 mL | 911.2 | 70.1% | 72.4% | 0.968 |
| Petri dish 60 mL | 2614.5 | 70.3% | 75.3% | 0.934 |

The last column carries over into the activity as a one-to-one multiplier, and it reproduces the observed excess: for the Petri dish, 1.156 vs. 0.934 — that is exactly the 24%; for "Denta," 1.126 vs. 0.971 — 16%; for the Marinelli beaker, 1.021 vs. 0.955 — 7%.

**The sign of the discrepancy is opposite at the two ends.** Under the 583.2 group the model gives MORE continuum than measured, under the single 2614.5 line — LESS. This is a rotation, not a shift: the computed spectrum's baseline has a different slope. An overall efficiency multiplier cannot fix this, and this also explains why the check on single lines showed nothing — it was performed on 662, 1461, and 2614.5, i.e., at one end only.

**The effect grows as the sample shrinks:** Marinelli 7%, "Denta" 16%, Petri dish 24%. The direction points to scattering outside the crystal — in the vessel, the holder, and the shield: the more compact the source and the closer it is to the crystal, the larger the fraction of quanta that arrive by an indirect path, and the more strongly an inaccuracy in describing them shows up. The tail of the measured peak cannot produce this: it would work in the opposite direction, taking area away from the peak in the measurement.

What needs checking next is not the deconvolution but the computed spectrum itself — by comparing the model's continuum with the measurement in a wide window, away from the peaks.

#### Continuum away from the peaks: exactly where the computed spectrum diverges

The check is set up without any deconvolution (`analysis/continuum.py`). The model gives counts per N sampled decays, the measurement gives counts per live time, so the ratio `N·(measurement count)/(model count)/t` is the activity that a spectral segment would give if the model were correct. Everything is normalized to the clean single 2614.5 line; in segments WITHOUT LINES the deviation from unity shows the baseline discrepancy directly. The segments are chosen automatically — everything within 2.5σ of any line with a yield brighter than 2% is discarded.

Marinelli, Th-232 record (ratio to the anchor):

| segment, keV | to anchor | segment, keV | to anchor |
|---|---|---|---|
| 108…172 | 1.04…1.08 | 1668…2092 | 1.05…1.27 |
| 388…652 | 0.97…0.98 | 2108…2252 | 1.17…1.38 |
| 1028…1172 | 1.00…1.14 | 2308…2452 | 0.83…0.89 |
| 1188…1292 | 0.91…0.97 | 2468…2492 | 1.03 |

**Below 700 keV the computed continuum is correct** — agreement within 3%, against a statistical uncertainty of 0.2…0.5%. This rules out the simple explanation of an error in the low-energy part of the model.

**The discrepancy sits above 1600 keV and has structure.** The model gives 5…38% too little baseline in the 1668…2252 band and 11…17% too much in the 2308…2452 band. The second band lies exactly on the Compton edge of the 2614.5 line (2381 keV): the calculation places the edge sharper than the instrument sees it.

**The order of the vessels here is REVERSED** relative to what the peak balance gives: the Marinelli beaker diverges the most, the Petri dish the least (there the same bands give 0.88…1.08). The Marinelli source is the strongest one in the kit, and the dependence on count rate points to pulse pileup, which is entirely absent from the calculation. Two mechanisms, not one: pileup explains the upper band and the order of the vessels, but does not explain the peak balance at 583.2, where the order is reversed.

Neither has been resolved down to a number, and both remain open.

### Independent cross-check: the same record processed by someone else's pipeline

The gamma-spectrum-analysis package processes the very same Th-232 record in the Marinelli beaker with its own pipeline and checks the result against the source's passport (`cert_zcheck`). Importantly, **it takes the efficiency from the LSRM passport curve** (`УДС-ГЦ-63х63-USB_-_Маринелли.efr`), not from our calculation. Agreement with our numbers is not built into this scheme.

| line | nuclide | our A/A_passport | their A/A_passport |
|---|---|---|---|
| 238.6 | Pb-212 | — (purity 0.89) | 0.948 |
| 911.2 | Ac-228 | 0.886 (purity 0.48) | 0.749 |
| 2614.5 | Tl-208 | **0.796** | **0.832** |

**Both pipelines give an activity below the passport value**, on one and the same record, with different sources of efficiency. On the 2614.5 keV line — the only one where both sides work with a clean single line — the discrepancy between us is 4.5%, while the shortfall relative to the passport is 17…20%.

This shifts the weight of the hypotheses on the section's main open question. The underestimate obtained when using the **instrument's own verification curve** cannot be explained by an error in our vessel model: that model does not enter that calculation. What remains are explanations common to both pipelines — for example, the actual mass or packing of the fill in this particular source specimen, or the referral of the passport activity to the measurement date.

**A caveat, without which the conclusion would be premature.** Their own certificate gate on this record is not passed (1 of 3 nuclides), the cause is not established, and the χ²/ν of the multiplet fit reaches three-digit values; their scatter across three lines of one series in secular equilibrium is 20% (0.75…0.95), against our 9% (0.80…0.89). So their numbers are accepted as an independent indication of the sign and order of magnitude of the discrepancy, but not as a benchmark. Checking the fill-mass hypothesis requires access to the source itself and has been split off into a separate task.

### BACKGROUND calibration is checked separately from the sample calibration

Operator's remark: "it looks like your background isn't calibrated." This was checked — and the check gave not quite what was expected, but something useful.

**There was a basis for the suspicion.** The kit's background records have a LINEAR calibration, with two coefficients, whereas the samples have four or five, and one and the same background file serves almost all records. The LSRM rule on this point is direct: the background calibration may differ from the sample calibration even on the same instrument (a different session, temperature, gain), and the calibration steps are applied to the background INDEPENDENTLY.

**What the check against the background's own natural lines showed** (K-40 1460.8; Tl-208 2614.5; Pb-214 351.9; Bi-214 609.3 — present in any lead-shielded background). For the main background file the residuals were −0.05 / −0.15 / −0.29 / −0.14 FWHM: the worst, 0.29, against the LSRM threshold of 0.30 — the calibration is right at the edge of acceptability, but does not need correcting.

**A separate case is the background of the thorium record**, where the residual came out to 0.49 FWHM. Investigation showed the culprit is not the scale but the ANCHOR: 1460.8 and 2614.5 sit exactly right (0.00 and −0.01 FWHM), while 609.3 is off — in a thorium background, 583.2 Tl-208 sits nearby, and on NaI the peak search returns a single peak for both. After discarding this anchor, the fit over the remaining three gives a residual of 0.00. This is exactly the case for which LSRM ranks anchors by RECOGNIZABILITY, not by intensity.

**The cost of the issue was measured, not eyeballed:** recalculating the background calibration changes the count rates in the kit's analytical peaks by less than 0.15%, and the activity table did not shift in a single digit. The kit's sources are strong, and the background within their windows is negligible. But for MDA and weak lines, where the background is what determines the result, this is no longer true, so the check is built into the pipeline: all analysis scripts read records through `becqmoni.read_checked()`, which checks the background calibration and, if the threshold is exceeded, recomputes it from the anchors after discarding the misbehaving ones. The verdict for each record is shown on the summary page.

### Acquisition time: the uncertainty is computed from the SHORTER arm of the pair

Operator's instruction: synchronize the acquisition time to the smaller of the sample–background pair. The rule is simple and correct in substance — **one cannot claim a precision better than the short arm allows**.

There is one subtlety, without which a straightforward implementation would change nothing. The count rate itself does not depend on the time base at all: R = S/t_sample − B/t_background, and "bringing both spectra to a common time" by dividing the counts changes nothing — neither the rate nor its variance. Reducing to a common time T is **thinning** of the longer acquisition, and for a Poisson quantity thinned with a fraction k, the variance drops by a factor of k, not k². Hence

    D(R) = σ_sample²/(T·t_sample) + σ_background²/(T·t_background),   T = min(t_sample, t_background)

At equal times this is exactly the usual formula; for any inequality it is more conservative than it.

**Times in the kit.** The background is almost everywhere the same, 54,000 s, with samples ranging from 300 s (Am-241, point source) to 62,000 s. In three records the background is SHORTER than the sample, by 5…15%: K-40 in "Denta" and the Petri dish, Th-232 in the Petri dish.

**The cost of the issue.** Line uncertainties grew by 2…12% — there, and only there, where the background makes a noticeable contribution: Ra-226 1120.3 in the Petri dish ×1.12, Cs-137 662 in the Petri dish ×1.11, Ra-226 1120.3 in "Denta" ×1.10. The activities themselves did not shift; in the summary table, individual uncertainty digits changed. The old formula is left available via the `sync=False` parameter — so the difference can be measured rather than argued about.

### Found during the check: a background shelf can sit on someone else's peak

Cross-checking the page's tables against `results/deconv_lines.csv` uncovered a defect not in the deconvolution but in the **windowed extraction**. The baseline is taken from the left shelf, E−2·FWHM … E−FWHM. For the Ra-226 351.9 keV line this is 279…316 keV, and sitting right there is **295.2 keV of the same Ra-226** with an intensity of 18%: part of another line is subtracted as "background."

The tell-tale sign was this: the result at 351.9 shifts by 5…12% across vessels when the FWHM law is changed (square root from a single point vs. the calibration FWHM² = a + b·E), while at all other lines the shift is under 1.5%. Capturing the peak itself cannot account for that much — ±1 FWHM contains 98.2% of a Gaussian. What changes is precisely the fraction of the 295.2 peak that falls into the shelf.

The defect is partly cancelled, because the area from the model spectrum is extracted with the same shelf, but not completely. **Consequence: the purity measure is insufficient.** It only looks inside the window and gives 0.99 for 351.9, i.e., "the line is clean," even though the baseline is corrupted. A second measure is needed — shelf contamination from the emission spectrum of the same run; this has been filed as a task. The deconvolution does not have this defect: its continuum is fitted, not taken from a shelf.

This also explains why **two** FWHM laws coexist in the calculation and must not be confused: windowed extraction works with the square-root law (that's how the published numbers were obtained), the deconvolution uses the calibration, because the fit is sensitive to the width's accuracy. The page's calibration table lists both, next to the measured width.

### Kit spectra: the algorithm's three steps, visible

All the numbers above are compiled from the kit's records, and the summary page shows exactly what they come from: each record's spectrum with its own background and difference, found peaks against the tabulated energies, measured FWHM against the calibrated law, and a deconvolution inset for each group — the data, the full model, the continuum with the step, and the contribution of each line separately.

The figures are built by `analysis/spectra_figs.py` **using the same run of the same functions** that compute the published numbers: the deconvolution curves are taken from the same design-matrix columns that solved the problem. The figure cannot diverge from the table — and if it does, that means one of the two was computed differently, and this is immediately visible.

What is worth seeing with your own eyes there:

- **the calibration shift is real**: in the cesium record the peak sits at 658.6 instead of 661.657 keV, for Tl-208 — at 2610.7 instead of 2614.5. An area window placed at the tabulated energy clips part of the peak and moves the background shelves onto the slope;
- **the blend's "centroid" is not a calibration shift**: for Ac-228's 911 keV it drifts by +18.7 keV, because it is the center of gravity of the 911+965+969 group. Shifting the window to it worsened the result (the scatter for "Denta" grew from 1.37 to 1.66) — hence the purity limiter;
- **the square-root-from-one-point FWHM law is unfit for fitting**: at 583 keV it gives 46.8 against a measured 40.4. For a ±1 FWHM window this is nearly immaterial; for the fit it cost 9% of the area at 2614.5.

Two records did not make it into the section and are called out explicitly on the page: K-40 in "Denta" and in the Petri dish — the spectra exist, but there is no K-40 decay run for these geometries, so there is nowhere to take the per-decay efficiency from. A silent omission reads as "everything has been calculated," and one such omission has already cost weeks.

**These numbers were obtained after fixing a defect that made them meaningless.** In the chain runs, Geant4 carried the whole series through to completion within a single event: the long-decay threshold was raised to 1e30 ns, since otherwise long-lived members would not decay at all. The energy depositions of Ac-228, Tl-208, and Bi-212 were then summed, producing coincidences between nuclei that in nature decay years apart. This depleted the peaks by 25–46% in the volume geometries and by 7% in the point geometry — there the probability of a double hit is much smaller. Now the energy depositions are cut into groups by global time with a 1 μs threshold; single nuclides still give byte-identical spectra to before.

What changed in the ratios (Marinelli, three thorium lines): it was 1.26 / 0.99 / 1.05, now it is 0.82 / 0.89 / 0.80. Cesium and potassium did not shift — their runs are single-nuclide.

**The main point here is not the level shift but the vanished contradiction inside the vessel.** Previously, in the Marinelli beaker, cesium said 0.78 and thorium said 1.26 — a 60% difference on one and the same geometry, which the geometry cannot explain, since it does not know which nuclide is inside. Now the whole vessel fits within 0.69–0.89, and each geometry is left with ONE multiplier instead of a scatter across nuclides. The discrepancy between vessels, meanwhile, persisted and kept opposite signs — so this is still not an overall normalization but the vessel models.

**Shape is a separate finding.** The trend of the ratio with energy:

| E, keV | 59.5 | 81 | 88 | 122 | 245–360 | 583–1115 | 1173–2615 |
|---|---|---|---|---|---|---|---|
| MC/exp | 0.78 | 0.82 | 0.90 | 0.88–0.92 | 0.90–1.00 | 0.87–0.97 | 1.01–1.14 |

The soft edge is understated by 10–22%: there is too much absorber in front of the crystal in the model. The first suspect was named during construction — **the MgO reflector's bulk/packing density was taken as 2.0 g/cm³ as an ASSUMPTION** (the usual range is 1.5–2.4).

This hypothesis received independent confirmation from a comparison of the two distances: at 59.5 keV the shortfall is 0.78 at 5 cm and only 0.86 at 25 cm. At five centimeters the source is seen at a wide angle and the quantum's path through the front stack is longer; at twenty-five centimeters the entry is nearly normal. An excess of absorber must give a larger deficit precisely at oblique entry, which is what is observed. A mechanism unrelated to the absorber (light collection, for example) would not produce this dependence on distance.

At 59.5 keV, six millimeters of MgO at ρ = 2.0 transmit 71%, while at ρ = 1.3 they already transmit 80%. In other words, the americium point effectively measures this density, and it favors 1.3–1.5 rather than the assumed 2.0. This should not be "fitted," but recomputed with the corrected density and checked that the soft edge is fixed WITHOUT breaking the hard edge.

The hard edge is overstated (2614.5 → 1.137). Here a known mechanism is at work, which the model fundamentally does not reproduce: in a real scintillator, part of the full-energy events are lost to incomplete light collection and move from the peak into the tail, and the effect grows with energy, because the interaction points spread over the whole volume of the crystal.

### 2.5 The same procedure on three vessels: where the parameter stops being a length

The fit was repeated for all geometries in the kit. The density pairs used were 0.60 and 1.60 for the cuvettes (the kit's sources have densities from 0.57 to 1.60 — mass at one nominal volume) and 1.00/1.60 for the Marinelli beaker.

| vessel | fill height h, mm | fit d, mm | χ²/dof | LSRM `.efa` | d/h | deviation |
|---|---|---|---|---|---|---|
| Marinelli 1 L | 31.0 | 31.5 (30.8–32.5) | 0.91 | 31 ± 2 | 1.02 | matches |
| "Denta" 120 mL | 29.5 | 27.0 (26.8–27.2) | 0.73 | 33 ± 3 | 0.92 | −2.0 σ |
| Petri dish 60 mL | 10.6 | 12.5 (12.2–12.5) | 1.08 | 10 ± 1 | 1.18 | +2.5 σ |

The main point here is χ²/dof near unity in ALL three cases: a single thickness describes the whole 45.3…3552.5 keV range, i.e., the shape of the correction f(μρd) is correct for all geometries. This is the model's operability.

But **the numerical value of d is a fitting parameter, not a geometric length**, and the exact match for the Marinelli beaker cannot be presented as confirming its physical meaning. Arguments:

1. The deviations from `.efa` for the two cuvettes have OPPOSITE SIGNS (−2.0 σ and +2.5 σ). An error in the model's physics would push them the same way.
2. The ratio d/h is likewise not constant (1.02 / 0.92 / 1.18), nor is LSRM's own `.efa`/h ratio (1.00 / 1.12 / 0.94).
3. For a flat cuvette on the detector's end face, the mean chord through the layer under isotropic emission is KNOWN to exceed the fill height — oblique paths. So the physical "effective thickness" must satisfy d > h, while the fit for "Denta" gives d < h. These are different quantities.

The discrepancies for the cuvettes point to the weakest spot in their geometry: the outer diameters are taken from the cuvette table, and the 1.5 mm wall is an ASSUMPTION. Notably, the table describes a Petri dish of 0.075 L at Ø88, while the kit's Petri dish is 60 mL. The fitted 12.5 mm would correspond to an inner diameter of about 78–80 mm instead of the assumed 85; for "Denta" the shift goes the other way and corresponds to an inner Ø closer to 75 than the assumed 72. In other words, the fit acted as a measuring instrument and showed exactly where the assumption is weakest. This can only be verified with cuvette drawings, which are not available.

---

## 3. Cascade summing

In a mono-energetic run, a single quantum arrives at the crystal, so there is nothing to sum with. In a full-decay run the cascade quanta fly out simultaneously, and if two of them interact in the crystal, the event moves out of its own peak. Hence

    eps_decay(E) = A_peak(E) / N_emitted(E),      C(E) = eps_mono(E) / eps_decay(E)

The denominator is the number of quanta of that energy ACTUALLY emitted during the run. The line yield is nowhere entered by hand: it comes from the same PhotonEvaporation database as the transport. Agreement with reference values — Cs-137 661.7 → 0.8513 (reference 0.851), Bi-214 609.3 → 0.4574 (0.4549), Tl-208 2614.5 → 0.9974 (0.9975) — serves as a check of the counter.

This calculation is done by a run, not by the formula C = 1/(1 − Σpε), DELIBERATELY: the formula only knows about losses and in principle cannot give C < 1, whereas in a real decay scheme the sum of two cascade quanta is often exactly equal to the energy of a cross-over transition and lands as a sum peak right on a real line.

### 3.1 Method control

Cs-137 and K-40 have no cascade, so their correction must come out to unity. For cesium, the 661.7 quantum and the 31–37 keV barium X-ray never coincide in time at all: they are alternatives — the transition proceeds either by gamma or by internal conversion.

| nuclide | C |
|---|---|
| Cs-137, 661.657 | 0.983 ± 0.015 |
| K-40, 1460.822 | 1.009 ± 0.025 |

Both agree with unity. The procedure is correct.

### 3.1a Chain-truncation windows: do they even work

The entire branching method rests on the `/process/had/rdm/nucleusLimits` windows: a run of a "single" nuclide must terminate the chain at it, otherwise daughter lines will contaminate the spectrum and the line yield will end up computed against the wrong parent. The long-decay threshold is raised to 10³⁰ ns, meaning **nothing but the window holds the chain in place**. At the same time, it is known that in analog mode — which is the default here and is not overridden anywhere — the windows can fail to trigger, in which case they are purely decorative.

Silence here is not proof, so the check produces a number (`analysis/check_nucleuslimits.py`). Am-241, with the window `241 241 93 95`, decays into Np-237 (2.14 million years), then further into Pa-233 (27 days) with lines at 300.1 and 311.9 keV, yields 6.6 and 38.6%. Np-237 is cut off by mass number: A = 237 is outside the window 241…241. If the window is not working, with the time threshold raised these lines would come through at full strength — about 180,000 counts for 400,000 primary nuclei.

| quantity | obtained |
|---|---|
| Am-241's own 59.5 keV line, yield per decay | 0.361 (reference 0.359) |
| Pa-233 lines 300.1 + 311.9 keV | **1 count** against an expected 180,800 |

**The windows work**: suppression by a factor of 1.8·10⁵, the single count being a random coincidence with the continuum. The first row serves as a control that the run is meaningful at all. The branching method legitimately relies on the windows.

### 3.2 Magnitude of the effect

| nuclide | E, keV | C |
|---|---|---|
| Tl-208 | 583.187 | 1.146 ± 0.020 |
| Tl-208 | 2614.511 | 1.185 ± 0.030 |
| Bi-214 | 609.320 | 1.115 ± 0.023 |
| Bi-214 | 768.360 | 1.149 ± 0.124 |
| Bi-214 | 1120.294 | 1.156 ± 0.049 |
| Bi-214 | 1764.491 | 0.954 ± 0.035 — unreliable, see below |

Losses of 11–16%. For comparison: for RadiaCode with a one-cubic-centimeter crystal, the same effect was 0.06–0.33%, i.e., here it is thirty to forty times larger. The reason is the total efficiency in the Marinelli geometry, about 9%, against a fraction of a percent for the small crystal.

Robustness was checked by varying the peak window and the continuum shelf. The controls stay put (0.983…0.984 and 1.007…1.009), lines 609.3, 768.4, and 1120.3 hold within 2%. **Line 1764.5 varies from 0.95 to 1.19** — there, the sum peak 609.3 + 1155.2 = 1764.5 lands exactly on the line itself, the "peak" becomes a superposition of the line and the sum, and they cannot be separated by fitting the continuum. This correction is excluded from further conclusions.

### 3.3 The LSRM curve is already corrected for summing

Before applying the correction to the comparison, it had to be established whether it was already built into the measured curve. This is not an idle question: an error here would introduce a 10–20% systematic where none exists. The check was performed on LSRM's own data, without reference to our model — in three independent ways.

**A cascade-free point against cascade-bearing neighbors.** Cs-137 sits above the interpolation between 609.3 and 768.4 (both Bi-214) by a factor of **1.011 ± 0.036**. This is 0.3 σ from "corrected" and 2.9 σ from "not corrected."

**Smoothness of the curve.** The efficiency curve is physically smooth. A polynomial fit in log-log over seven points with the known correction:

| degree | curve as-is | curve × C |
|---|---|---|
| 2 | χ²/dof = 2.50 | 12.91 |
| 3 | χ²/dof = 2.79 | 17.09 |

Multiplying by the correction worsens the smoothness fivefold.

**Cascade-free points against the curve built from cascade lines.** Cs-137 gives 0.961, K-40 — 0.959, against an expectation of 1.00 for a corrected curve and 1.13 and 1.16 for an uncorrected one.

All three tests agree: **the LSRM curve is corrected for cascade summing**, and applying the correction C to the comparison would be an error. Section 2.2 remains valid unchanged.

A small detail for the future: both cascade-free points landed 4% BELOW the curve built from cascade lines. Individually each is a bit more than one sigma, but both lean the same way — a slight over-correction on LSRM's part cannot be excluded. This cannot be asserted from just two points.

### 3.4 Radon: where equilibrium is assumed, and where it is checked

Everything above that relies on Bi-214 lines assumes secular equilibrium Ra-226 → Rn-222 → Bi-214, and this holds only to the extent that the sample RETAINS radon (operator's remark). Emanation is not modeled in the calculation — equilibrium there is exact by construction, and the "chain vs. direct run" comparison checks the code, not the physics of the sample.

For the real M_ra source, leakage is checked by the comparison itself: with noticeable emanation, the count in the Bi-214 lines is understated, the efficiency recovered from it is understated, and the MC/exp ratios at Ra-226 points would stand systematically higher than for the cascade-free Cs-137 and K-40. In fact: 609.3 → 1.21; 768.4 → 1.21; 1120.3 → 1.08, against 1.19 (Cs) and 1.22 (K) — there is no excess, radon is mostly retained. This agrees with the source's history: sealed in 1999, equilibrium with radon is established over ~a month.

For FUTURE samples the caveat remains in force: the Ra-226 activity from Bi-214 lines (including the MDA from 609.3) is correct only to the extent of the particular sample's emanation; loose and moist samples lose radon noticeably.

### 3.5 Practical implication

The summing calculation does not explain the shape scatter (the curve is already corrected), but it is needed by anyone who will compute activities from real spectra of this unit: neglecting it will understate the result for thorium and radium by 11–16%.

---

## 4. Minimum detectable activity

An end-to-end check of the whole chain: the unit's **measured** background (Marinelli beaker with distilled water, 54,000 s, embedded in the kit's XML) + the **computed** efficiency in water (a separate grid, water in the model, not rescaled from OISN-16) + line yields per parent decay from the chain runs. Passport, item 2.2: 1 L Marinelli beaker with distilled water, 2 h.

| nuclide | line | n_bg, cps | MDA (4√2 rule) | MDA (Currie) | passport |
|---|---|---|---|---|---|
| Cs-137 | 661.7 | 0.324 | 1.86 | **1.55** | 1.5 |
| K-40 | 1460.8 | 0.199 | 22.9 | **19.1** | 25 |
| Ra-226 | 609.3 | 0.365 | 3.42 | **2.84** | 3 |
| Th-232 | 2614.5 | 0.033 | 4.55 | **3.88** | 3 |

Bq/kg. An observation that was not built in: the passport values fall on the **Currie** formula (2.71 + 4.65·√B)/(t·η), not on the 4√2·√n_bg/(√t·η) from the LSRM methodology — for cesium, 1.55 against a passport value of 1.5, versus 1.86 from 4√2. That is, judging by the numbers, ASPEKT normalized its MDA by Currie/GOST. Both columns are given.

Caveats. (1) The efficiency here is on the MC scale; on the experimental scale both columns are multiplied by 1.165 — then cesium becomes 1.8 against a passport value of 1.5, and the remainder is explained by a difference in backgrounds: our background is from 2016, the passport verification was performed in 2016 at the manufacturer's site with its own background. (2) K-40: the background window contains potassium's own line from the surroundings — the potassium MDA is determined by the background potassium, and this is visible in the numbers. (3) The Ra-226 MDA is from the Bi-214 line, correct only to the extent the sample retains radon; water retains radon (operator's remark). (4) MDA is a characteristic of the method, not an upper bound for a specific sample.

---

## 5. What was not modeled

- The PMT voltage divider, the USB-ADC board, and connectors — replaced with a homogeneous block behind the envelope. Far from the crystal, they affect the result only through backscattering.
- The permalloy magnetic shield — absent from the drawing.
- The carriage, hinges, and mechanism of the shield lid; the cable entry in the lead.
- The intrinsic activity of the shield materials (Pb-210 in the lead, K-40 in the steel) — the background is taken as measured, not computed.
- The non-proportionality of NaI(Tl) light yield: scoring is done by energy deposition, not by collected light. Above 300 keV the effect is ≤3–4%, below that it grows.
- Angular correlations of cascade quanta: Geant4 reports at startup `Enable correlated gamma emission 0`, and the quanta are sampled isotropically and independently. For the Marinelli geometry, where the sample surrounds the detector on almost all sides, the contribution is small; for the point geometry it would be more significant.

---

## 5a. Independent check: the model in NuclideMaster format

The model was exported to the EffCalcMC/NuclideMaster (LSRM) format — layer by layer, following the operator's samples: `nuclidemaster/G1S_NaI63x63.din` (detector) and `nuclidemaster/G1S_Marinelli_1L_OISN16.sin` (Marinelli 1 L with OISN-16). The point: a different Monte Carlo code with the SAME layers discriminates between two hypotheses about the model's energy dependence (§2) — if it reproduces it, the defect is in the layer description; if it does not, the defect is in our code.

Both files were checked by the operator in Geometry Master on 2026-07-28: they parsed and rendered without errors, crystal volume 196 cm³ (= π·3.15²·6.3), source volume 1.01·10³ cm³ against a target of 1000 (1% — rounding of the fill height 8.36 → 8.4 cm). The layer-correspondence table, the assumptions the files carry (MgO ρ = 2.0 bulk density; layers outside the lid replaced with a 0.15 cm distance), and two minor Geometry Master display quirks are in `nuclidemaster/README.md`.

**Result of the first calculation (Marinelli, 10⁸ trials).** The ratio ours/EffCalcMC is not flat, and it cannot be reduced to a single number: ours is higher by 20–22% at three nodes below 70 keV, lower by 3–7% in 88–166, and above 90 keV there is a monotonic rise of +7.5% per decade with agreement of ±3% (RMS). A straight line through all points comes out near zero only because the soft-end offset cancels the hard-end rise — a masking average (caught by the auditor; the comparison was redone: the dense ECM curve is interpolated onto our nodes, not the other way around).

This does not weaken the conclusion but strengthens it: the code-to-code discrepancy is ~10% across the range and is opposite in sign to the measured trend — the measured 18–52% cannot be reproduced by choice of code — **the energy dependence belongs to the layer description, not the code**. At the same time, two codes on the same layers check the CODE, not the level: the level is checked by measurement, and it gives 0.780 ± 0.025 for the Marinelli beaker — agreement between the codes does not refute this, it points to the layers. The EffCalcMC curves and the comparison are in `nuclidemaster/EffReg_G1S_Marinelli.{efa,efr}`, `analysis/compare_effcalcmc.py`.

**The EffCalcMC point-source calculation (5 cm, 10⁸ trials) closed the fork for good.** On the geometry where the energy dependence is MEASURED, the two codes agree within an RMS of 2.9% (all 22 nodes: mean 0.968, slope −0.4% per decade; soft nodes 45–60 keV — 1.017…1.024), whereas the measured curve departs from both by 0.78 → 1.14 across the range. Along the way, the soft structure of the Marinelli comparison was explained: the +20% at 45–60 keV in the point-source case disappeared — it belonged to the fill (self-absorption in OISN-16), not to the detector. The inter-code residual of +5.6%/decade above 90 keV (+7.5 in the Marinelli case) is probably cross-sections; it is 6–8 times smaller than the effect being sought. Files — `nuclidemaster/EffReg_G1S_Point5cm.{efa,efr}`. Next to look for is a discrepancy between the drawing assumptions and the real instrument: the entrance-face stack (the MgO density has already been ruled out by the scan) and the crystal dimensions.

**The explanation of the 0.968 level via a difference in distances is WRONG and has been retracted** (caught by the auditor: passport item 2.10 sets the distance "from the surface of the detector lid," while `.sin` has `Distance,cm=5` with no compensation; the argument rested on a rubber protector which, after the geometry fix (see below), does not exist). The level requires a new explanation on the recomputed grid.

**Dimensional uncertainties of the sources** (operator's remark) enter the same budget, but differently. For point OSGI sources they shift the LEVEL, not the shape: the active-spot position ±1 mm — ±4% of the level at 5 cm by the inverse square, energy-independent; the encapsulation film ≤0.3% even at 59.5 keV — they cannot account for the energy dependence measured in the point geometry. For the Marinelli beaker and the cuvettes, the dimensions bend the SHAPE at the soft edge through self-absorption: the fill level, the actual mass versus nominal, the packing density and its non-uniformity — exactly what the +20% discrepancy between the two Monte Carlo codes at 45–60 keV belonged to. The kit's sources are UNAVAILABLE for physical measurement, so the budget is estimated only by calculation — sensitivity scans via runs (level, density, mass per the inventory) — and enters the final numbers as an interval of systematic uncertainty, not as a correction: there is nothing to justify a correction without a physical measurement.

---

## 6. Reproducing the results

Paths are taken from environment variables; the code contains not a single machine-specific path (see `common/py/paths.py`).

```bash
export GEANT4_ROOT=<Geant4 root>
cmake -S detectors/Gamma-1S/geometry -B build/Gamma-1S -G Ninja
cmake --build build/Gamma-1S
```

| step | command |
|---|---|
| reference data | `python reference/fetch_efr.py`, `fetch_kit_xml.py`, `fetch_bg.py` |
| masses and geometry | `g1s masses.mac shield` |
| passport point 25 cm | `g1s test.mac open` |
| grids for all geometries | `python drivers/run_all_grids.py` |
| attenuation coefficients | `mucalc` |
| decay runs | `g1s macros/decay_all.mac vessel` |
| decays in cuvettes and point source | `python drivers/run_cup_chains.py`, `run_point_decay.py` |
| comparison against LSRM | `python analysis/compare_lsrm.py`, `compare_point.py`, `compare_cups.py` |
| effective thickness | `python analysis/selfabs_fit.py` |
| summing and control | `python analysis/summing.py` |
| checking the reference curve | `python analysis/tcc_evidence.py`, `resolve_tension.py` |
| recalculating kit records | `python analysis/kit_recalc.py`, `point_recalc.py` |
| mixture unfolding | `python analysis/mix_unfold.py` |
| MDA | `python analysis/mda.py` |

`resolve_tension.py` needs `SPECTRAVIBE_ROOT`: it reads raw LSRM `.spe` files with the standard `gamma.io.lsrm_spe` reader — the counts in them are packed in binary and cannot be parsed with a custom reader.

The reference data are stored here in full, under `reference/lsrm/`: the verification kit in two formats (the instrument's binary `.spe` and BecqMoni XML), the LSRM working tree with the 2016 and 2024 verifications, `.efr`/`.efa` curves, averaged backgrounds, reference-source passports, and nuclide libraries. Calculating from the XML and from the `.spe` gives the same result — checked by `analysis/xml_vs_spe.py` over 39 matching pairs.

The data can be refreshed from the original source with the loaders in `tools/`; the original source itself is
[spectravibe-toolkit](https://github.com/VibeEngineering-LLC/spectravibe-toolkit).

The finished computed curves are in [results/](results/), with their description and caveats in [results/README.md](results/README.md).

---

## Statement on the use of artificial intelligence

The text of the report, the code of the computation and analysis scripts, the model geometry, and the layout of the summary page were prepared with the participation of AI (Claude, model Opus 5) in the roles of analyst and verifier: building the geometry from drawings and passports, implementing the analysis algorithms, finding and investigating discrepancies, self-checks, and cross-checks against independent implementations.

The operator set the tasks, provided the verification-kit measurements and passport data, accepted the results, and bears responsibility for them. The review of the gamma-spectrum-analysis package's source code, on which the methodological borrowings rely, was carried out by a separate auditor agent.

All numbers were obtained by calculation and are presented as they are, including those that did not agree.
