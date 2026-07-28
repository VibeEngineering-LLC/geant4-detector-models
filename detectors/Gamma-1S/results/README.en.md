# Computed curves for Gamma-1S

Ready-made numbers: to take the efficiency at the energy you need, there is no
need to build Geant4 and compute for a day. The computed spectra themselves
are not committed to the repository — they are reproduced by the drivers, and
what is stored here is the thing they were computed for.

Everything in this directory is produced by a single script:

```bash
python detectors/Gamma-1S/analysis/export_curves.py
```

It needs the directory of computed spectra — the variable
`G4MODELS_BUILD_GAMMA_1S` (or the shared `G4MODELS_BUILD`, if there is only
one instrument). If there are no spectra yet, first run:

```bash
python detectors/Gamma-1S/drivers/run_grid.py
python detectors/Gamma-1S/drivers/run_all_grids.py
```

## What eps means in these tables

The absolute detection efficiency in the **full-energy peak** (FEP): the
fraction of quanta of a given energy emitted in the sample that produced a
count in the peak. It is not divided by activity or by the line yield — it is
a characteristic of the "sample + vessel + detector + shielding" assembly, not
a property of the nuclide.

To get the count rate from it:

    N_peak [counts/s] = A [Bq] · p_gamma · eps(E)

where `p_gamma` is the line yield per decay. For cascading nuclides this also
includes the true coincidence summing correction, see `summing_C.csv`.

## Files

### Format: two safeguards, not one

All the tables here are **standard CSV**, written via `csv.writer`. On top of
that, **no value contains a single comma**, so the file can also be parsed by
a naive `split(",")` by position. The second safeguard is needed because the
reader is not under our control: quoting saves `csv.DictReader`, but does not
save whoever slices the line by hand.

Before 28.07.2026 there was neither of the two: `export_curves.py` joined
fields with commas without escaping, and the `geometry` of the point-source
grids was "point-source, 25 cm". The 48 rows of both point-source curves in
`efficiency_curves.csv` ended up with 14 fields against a 13-field header. The
parser did not crash, it just shifted everything one field to the right —
`E_keV` came out as the word "ОТКРЫТА" (OPEN), `eps_net` as the energy value.
The volume curves were read correctly, which is why the defect survived.
**The numbers themselves were correct**: only the layout across columns was
broken, and after the fix not a single value changed.

The geometry labels became `точечный 5 см` and `точечный 25 см` (point-source
5 cm / point-source 25 cm, without a comma), and the note in
`runs_manifest.csv` became "рабочая кривая; сверялась с .efr ЛСРМ" (working
curve; checked against the LSRM .efr).

The `python tools/check_csv.py` check is part of `recalc_all.py` and of the
writing step itself: `write_csv` re-reads its own file and fails if the row
length does not match the header.

| file | contents |
|---|---|
| `efficiency_curves.csv` | all nine curves in one long table — for processing |
| `eff_<label>.csv` | one file per curve per geometry — for reading and plotting |
| `runs_manifest.csv` | what and with which assumptions each grid was computed |
| `summing_C.csv` | true coincidence (cascade) summing corrections per line |
| `kit_recalc_volume.csv` | recalculation of the verification kit's volume entries (`kit_recalc.py`) |
| `kit_recalc_point.csv` | recalculation of the kit's point-source entries (`point_recalc.py`) |
| `kit_activity_volume.csv` | summary per the LSRM rule, volume: one row per nuclide plus `*` for the vessel |
| `kit_activity_point.csv` | the same for point-source geometries |

The last two tables are the **single source of the summary numbers**: both
the report and the `docs/gamma-1s/` page read them, rather than each
combining the line table with its own formula. While the page combined by
median and the report by the LSRM rule, the same recalculation carried two
different numbers.

### Kit recalculation: what it is and how to read it

The inverse problem. From the measured spectrum of a reference source, the
count rate in the peak is taken, divided by the computed efficiency per
decay — giving an activity. It is compared with the certified value, brought
to the measurement date. Column `ratio` = A_measured / A_certified:

* `ratio < 1` — the model **overestimates** the efficiency (the computed eps
  is higher than the true one, so the recovered activity comes out lower than
  the certified one);
* `ratio > 1` — the model **underestimates** it.

This is the strictest possible check of the model: it is not about the shape
of the curve but about the absolute value, and it cannot be fitted — the
certificate is set independently.

### Which lines go into the activity calculation and which do not

**A poorly resolved line is not usable for the activity calculation.** Its
area is the sum of several lines, and it cannot be attributed to a single
transition, no matter how many corrections are applied. The measure of
resolution is the `purity` column: the fraction of the yield within the
±1 FWHM window (±3 keV) attributable to the line itself. It is computed from
the emission spectrum of the same Geant4 run, so it accounts for the
**intensity** of neighboring lines, not just their distance. The `usable`
threshold is 0.95.

| line | purity | goes into activity | contaminated by |
|---|---|---|---|
| Cs-137 661.7 | 1.00 | yes | — |
| K-40 1460.8 | 1.00 | yes | — |
| Ra-226 295.2 / 351.9 / 609.3 | 0.96 / 0.99 / 0.98 | yes | — |
| Th-232 2614.5 | 1.00 | yes | — |
| Ra-226 1120.3 | 0.85 | no | Bi-214 1155 — 9 % |
| Ra-226 1764.5 | 0.82 | no | Bi-214 1730 — 15 % |
| Th-232 238.6 | 0.89 | no | Ac-228 209 — 8 % |
| Th-232 583.2 | 0.89 | no | Ac-228 562 — 3 % |
| Th-232 911.2 | **0.48** | no | Ac-228 968 — 29 %, 964 — 9 %, Tl-208 860 — 8 % |

The threshold has been checked against an independent measure — the Sparrow
limit of 0.85·FWHM, formal resolvability by distance: both measures give the
same set of usable lines. Unusable rows are not removed from the table
(`usable=0`) — they are used to test the deconvolution, and the 911 keV line
is well suited for exactly that.

**Consequence:** Th-232 is left with only ONE usable line, 2614.5 keV. There
is no margin of lines there.

### Activity per the LSRM rule

The activity is computed neither as a median nor from a single line, but as a
**weighted mean**, with weights `w = 1/(ΔA)²` [Algorithmic Foundations of
SpectraLine, §5; the "Activity in Counting Samples" procedure], with the
uncertainty taken as the **maximum** of the two estimates — the weighted-mean
uncertainty and the spread. With such weights, a line with a threefold smaller
error gets nine times the weight, i.e. the result is determined by the best
line — but this is recorded honestly, through averaging.

| geometry | Cs-137 | K-40 | Ra-226 | Th-232 |
|---|---|---|---|---|
| Marinelli 1 L | 0.782 ± 0.039 | 0.725 ± 0.074 | 0.788 ± 0.056 | 0.796 ± 0.048 |
| "Denta" 120 mL | 1.065 ± 0.055 | — | 1.149 ± 0.052 | 1.090 ± 0.068 |
| Petri 60 mL | 1.253 ± 0.065 | — | 1.352 ± 0.074 | 1.195 ± 0.072 |

In the Marinelli beaker all four nuclides agree with a single value of
≈0.78 within their uncertainties — i.e. the vessel really has ONE single
factor, not a nuclide-dependent one. For "Denta" it is ≈1.10, for the Petri
dish ≈1.25. K-40 is absent for the cuvettes: no potassium-decay run has been
computed for them.

The point-source geometries are combined by THE SAME rule and with the same
purity selection (`kit_recalc_point.csv`): 1.179 ± 0.019 from 21 usable lines
at 5 cm and 1.184 ± 0.048 from 7 lines at 25 cm. Excluded were Ba-133 356.0
(0.87), Co-57 122.1 (0.89), Eu-152 964.1 (0.91) — at 25 cm these same lines
had already been excluded earlier, due to contaminated shelves.

The combination is two-step: the nuclide's lines are combined by statistics,
then the certified-source fraction is added in quadrature, and only then are
the series combined with each other. The certified uncertainty (2–5 %) is
independent across nuclides and dominant — statistics gives 0.1–0.6 % per
line. The `chi2_dof` column on the `*` rows shows the agreement of the
series: 4.0 and 7.9 for the point-source geometries (the set does not
collapse to a single number, there is an energy trend), 0.23 / 0.64 / 1.17
for the vessels.

Purity is computed from the emission spectrum of its OWN run, so the same
line can be usable for a point source and unusable for a volume source:
Th-228 238.6 and 583.2 are clean in the point-source geometry (0.99 and
1.00), but in the Marinelli beaker — 0.89, because there the source is
Th-232 and above Th-228 stands Ac-228 with lines at 209 and 562 keV. This is
not different thresholds, but different sources.

Differing signs across vessels are not an overall normalization: a single
efficiency error cannot behave like that. The discrepancy belongs to the
**vessel models**, whose wall thickness and well are taken as an operator
assumption (there are no cuvette drawings, see `geometries/README.md`).

### Curve columns

| column | meaning |
|---|---|
| `E_keV` | quantum energy, keV |
| `eps_net` | FEP efficiency, peak area **net of** the left-hand continuum shelf |
| `d_eps` | statistical uncertainty of `eps_net` (Poisson on the peak plus shelf noise) |
| `eps_gross` | the same without subtracting the shelf |
| `N_primaries` | how many quanta were sampled in the run |
| `net_counts` | net peak area in counts |
| `solid_angle_fraction` | fraction of the solid angle (see below); 1.0 for volume geometries |

**Why there are two efficiency columns.** The FEP window is ±6 keV around the
line, the continuum shelf runs from E−30 to E−10 keV. `eps_net` subtracts the
shelf (this is how the `compare_lsrm.py` curve is built, and these numbers
were used for the comparison with LSRM), `eps_gross` does not (this is how
`compare_point.py` computes). The difference between the columns is the
systematics of the area-taking method, not a spread of the calculation.
Publishing one number while staying silent about the other would amount to
passing off a processing choice as a property of the detector. For soft
lines the difference reaches 10 %, for hard lines it is under 1 %.

The shelf is taken only on the left: for a monoenergetic source with no
pile-up there are no counts to the right of the line, and the continuum under
the peak is formed by events with almost full energy deposition.

**The uncertainty here is statistical only.** The uncertainty of the geometry
(MgO reflector thickness, well depth, vessel dimensions), of the cross
sections and of the matrix composition is not included in it and is
deliberately larger. Comparison with the measured LSRM curve gives an idea of
the full uncertainty — see `../REPORT.md`.

### Geometries

| label | geometry | matrix | density | shielding lid |
|---|---|---|---|---|
| `rho1.00`, `rho1.60` | Marinelli 1 L | OISN-16 | 1.00 and 1.60 g/cm³ | closed |
| `water1.00` | Marinelli 1 L | water | 1.00 g/cm³ | closed |
| `denta0.60`, `denta1.60` | Denta 120 mL | OISN-16 | 0.60 and 1.60 g/cm³ | closed |
| `petri0.60`, `petri1.60` | Petri 60 mL | OISN-16 | 0.60 and 1.60 g/cm³ | closed |
| `p5cm` | point source, 5 cm from the end face | — | — | closed |
| `p25cm` | point source, 25 cm from the end face | — | — | **open** |

Two densities are computed for each vessel not for completeness, but so that
fitting `f(x) = (1−e^−x)/x`, `x = μ(E)·ρ·d` gives an effective thickness `d`
and thereby lets any density be predicted: the verification-kit sources have
different densities at the same nominal cuvette volume.

`water1.00`, at the same density as `rho1.00`, shows the contribution of the
matrix composition separately from the contribution of density.

**The lid at 25 cm is open, not by oversight.** With the lid closed, 50 mm of
lead attenuate the flux by about e⁻⁶·², and the count is practically zero;
the certified "Point-source-25 cm" geometry is measured with the shielding
open, which is confirmed by the name of the LSRM background file
(`open_lid_point25cm`).

**The point-source geometries were sampled into a cone.** An isotropic source
at 25 cm wastes 99.4 % of events, so the quanta were sampled into a cone
around the direction toward the detector, and the efficiency is scaled to the
full solid angle by multiplying by its fraction (column
`solid_angle_fraction`). This is valid for the FEP: practically only direct
quanta go into the peak. The total (non-peak) efficiency under such sampling
would come out underestimated — quanta scattered off the shielding walls are
not sampled.

### True coincidence summing corrections

`summing_C.csv`. In the monoenergetic run a single quantum arrives at the
crystal, there is nothing to sum with — this is the "pure" efficiency
`eps_mono`. In the full-decay run the cascade quanta fly simultaneously, and
if two of them register, the event leaves its own peak:

    C(E) = eps_mono(E) / eps_decay(E)

The number of emitted quanta `N_emitted` is taken from the same Geant4 run
(`*_emit.csv`), not from a reference table: otherwise numbers from different
sources would be compared.

**Method control.** Cs-137 and K-40 have no cascade, so `C` must come out as
unity for them. The values obtained are 0.983 ± 0.015 and 1.009 ± 0.025 —
consistent. The uncertainty must be read: without it, 0.983 would look like a
1.7 % bias, which it is not.

`C` below unity (Bi-214, 1764 keV: 0.954 ± 0.035) means summing **into** the
peak: the event arrives in it from the addition of two other quanta.

**Caveat.** At startup Geant4 prints "Enable correlated gamma emission 0":
angular correlations between the cascade quanta are disabled by default, the
quanta are sampled isotropically and independently. In the Marinelli
geometry the sample surrounds the detector from almost all sides and the
contribution of correlations is small, but for the point-source geometry it
would be more significant.

## What is not here

Curves in the LSRM format (`.efr`/`.efa`) for loading into the instrument's
software. Producing them is a separate piece of work: the format contains
not only the points but also approximation parameters and a binding to the
geometry, and passing off a computed curve as a certified one is not
permitted. The measured LSRM curves themselves are in
`../reference/lsrm/efficiency/`.
