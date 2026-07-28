# Model validation protocol

A model without validation is a picture. Below is the minimum that makes
models of different instruments comparable with one another. Each check
compares the calculation against a number that was **not built into it**.

## 1. Geometry: masses against the nameplate

An instrument's certificate usually states the content of non-ferrous
metals. This is the only way to recover the dimensions of components that
are not in the drawings.

*Example.* The size of the Gamma-1S shield cavity is not specified in any
document. A cavity of Ø200 × 190 mm with 1 mm liners gives 167.1 kg of lead,
1.60 of copper, 1.58 of cadmium against the nameplate's "not less than
165 / 1.6 / 1.2" — three masses at once.

## 2. Drawing: do the stated overall dimensions add up

Layer labels usually admit more than one reading. The correct one is the
one for which ALL stated summary dimensions add up.

*Example.* Radially, `63 + 2·(3.65+0.5+2+1.5) = 78.3` against Ø78; axially,
from the crystal centre to the outer plane of the entrance face,
`31.5 + (6+0.5+1+2) = 41.0`.

*And an example of this check NOT working.* The first reading of the
entrance face put a 1 mm rubber protector and a 2 mm shock absorber there,
giving axially `63 + (6+0.5+2+1.5+0.5+1) = 74.5` against the stated 74 — it
agreed to within half a millimetre, and the reading was taken as confirmed.
The operator later pointed out that there is no rubber on the entrance face
at all (it is radial only), and that the face stack is
`Al 2 / air 1 / Al can 0.5 / MgO 6`. The sum of overall dimensions is
consistent with both variants: the check rejects grossly wrong readings, but
it cannot distinguish readings that produce the same sum from different
terms. The composition of the layers is confirmed only by someone who has
held the instrument in their hands.

## 3. Nameplate efficiency point

There is almost always at least one: "full-absorption-peak efficiency at
the 662 keV line for a point source at 25 cm — not less than X." Check it
with the shielding configuration in which it was measured (for the
Gamma-1S — with the lid OPEN; with the lid closed, 50 mm of lead gives
zero, and that is not a malfunction).

## 4. Point-source geometry as a reference for the DETECTOR

**The single most important check.** There is no vessel and no
self-absorption — only the crystal, reflector, housing and shielding
remain. Only this check separates a detector error from a sample-geometry
error.

*Lesson.* For the Gamma-1S, a discrepancy of 1.17 in the Marinelli geometry
was initially attributed to a systematic of the method. Point-source
geometries gave 0.971 and 0.931 — the conclusion had to be withdrawn: the
detector is correct, it is the vessel models that disagree. Without this
check, the wrong explanation would have stayed in the report.

Measure at TWO distances: the difference between them is sensitive to
excess absorber in front of the crystal (at oblique incidence the path is
longer).

## 5. Density trend: effective thickness

The ratio of efficiencies at two densities of the same matrix does not
depend on the absolute normalization, so it gives an independent
parameter:

`eps(E, rho) = eps_lim(E) · f(mu(E)·rho·d)`, `f(x) = (1-e^-x)/x`

The check is not the value of `d` itself, but the **chi2/dof of the fit**:
if it is around unity, one thickness describes the entire range, i.e. the
shape of the correction is correct.

Caution: `d` is a fit parameter, not a geometric length. For the three
Gamma-1S vessels, the ratio of `d` to the fill height came out as
1.02 / 0.92 / 1.18 — there is no constancy, and agreement with the
nameplate value must not be presented as confirmation of physical meaning.

## 6. Coincidence summing: a control check on a nuclide WITHOUT a cascade

Compute the correction with a full-decay run, not with the formula
`1/(1-Σpε)`: the formula only knows about losses and in principle cannot
give `C < 1`, whereas the sum of two cascade quanta often equals exactly
the energy of a cross-over transition and lands as a sum peak right on a
real line.

**Mandatory control:** Cs-137 and K-40 have no cascade, the correction for
them must come out to 1.00. If it does not — there is an error in the
procedure (window, shelves, emission counting), and the numbers for the
cascade lines cannot be trusted.

## 7. Line yields — from your own calculation

The counter of quanta produced by radioactive decay gives `p_gamma` from
the same database as the transport. Cross-checking against a reference is
a check on the counter, not a data source: Cs-137 661.7 → 0.8513
(reference 0.851), Bi-214 609.3 → 0.4574 (0.4549), Tl-208 2614.5 → 0.9974
(0.9975).

For chains, run the whole series through — the branching ratios come out
on their own (Bi-212 → Tl-208 came out at 0.3594 against a reference of
0.3594, and it is not written in anywhere).

## 8. Checking the reference curve bypassing its processing

The manufacturer's curve is already a processed result. It is worth
reconstructing it at least once from the raw spectrum and the source's
nameplate activity:

`eps = R_peak / (A · p_gamma)`

*Example.* For the Gamma-1S Marinelli, this gave 1.8686·10⁻² against the
recorded 1.8713·10⁻² — a ratio of 0.999. This removed doubt about the
reference curve and unambiguously attributed the discrepancy to the
calculation.

## 9. Whether the correction has already been folded into the reference curve

Before applying the summing correction to a comparison, check whether it is
already in the curve. Three independent tests
(`analysis/tcc_evidence.py`):

- the point of a cascade-free nuclide against interpolation between
  neighboring points of a cascade nuclide;
- smoothness of the curve before and after multiplying by the correction;
- the position of cascade-free points relative to a curve drawn only
  through cascade points.

*Result for the Gamma-1S:* the LSRM curve **has already been corrected**,
and applying the correction would introduce a systematic of 10–20% where
none exists.

## 10. End-to-end check: MDA against the nameplate

Measured setup background + calculated efficiency + line yields → MDA,
compared against the nameplate. Checks the whole chain at once.

*Side finding for the Gamma-1S:* the nameplate values fall on the Currie
formula `(2.71 + 4.65√B)/(t·eta)`, not on the LSRM procedure's
`4√2·√n/(√t·eta)`. State both.

## 11. Data in two formats — cross-check, don't just trust

If measurements exist both in the instrument's native format and in a
converted one, the calculation must be run on one of them, with the other
serving as a check. Do not compare "by eye": check the channel count, the
total sum of counts, live and real time, and the channel-by-channel
difference.

*Example for the Gamma-1S:* the calibration verification set exists both as
LSRM binary `.spe` and as BecqMoni XML; the whole analysis was written
against XML. Over 39 matching pairs, the per-channel difference came out
to zero for every channel — the formats are interchangeable, and this is
recorded by the script `analysis/xml_vs_spe.py`, not asserted in the report.

The silent assumption that "the converter doesn't corrupt anything" is
dangerous because, if a discrepancy exists, ALL conclusions will diverge at
once, and there will be nowhere to notice it.

## 12. Finished numbers must exist as files

A report with numbers embedded in the text is not a result, it is a story
about one. Curves, corrections and tables are delivered as machine-readable
files, with a description of every column: what quantity it is, in what
units, with what window the area was taken, what enters the uncertainty,
and by which run it is reproduced.

*Sign of a violation:* to find out the efficiency at a single energy, one
has to build Geant4 and compute for a day.
