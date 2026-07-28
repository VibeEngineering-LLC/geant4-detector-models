# Sample geometries

Common to all detectors: the same 1 L Marinelli vessel is used both for
Gamma-1S and for any subsequent instrument. That's why the dimensions live
here, not inside the detector folder.

## Dimension sources

The primary one is the LSRM table of measurement vessels, "Precision
Measurements," p. 11. It gives the overall dimensions and effective
thickness, but does NOT give the wall thickness, well diameter, and well
depth — these have to be taken as assumptions, and this is the weakest point
of all volumetric geometries.

| vessel | nominal, L | dimensions | d_eff LSRM, mm |
|---|---|---|---|
| Marinelli | 0,5 | Ø125 H=100 | 15 (2) |
| Marinelli | 1,0 | Ø150 H=110 | 26 (2) |
| Marinelli | 3,0 | Ø180 H=200 | 60 (5) |
| "Denta" | 0,12 | Ø75 H=35 | 36 (2) |
| Petri | 0,075 | Ø88 H=14 | 15 (2) |

**Note:** the table values refer to a typical assembly. The passport of a
specific detector unit (`.efa`) gives its own: for Gamma-1S — Marinelli 31
(2), "Denta" 33 (3), Petri 10 (1). The discrepancy between the table and the
`.efa` is noticeable, and you should rely on the `.efa` of the instrument you
are actually working with.

## Model status

| geometry | where implemented | comparison with measurement |
|---|---|---|
| Marinelli 1 L | `detectors/Gamma-1S/geometry` | calculation OVERESTIMATES by 17 %; suspect the well depth (assumption) |
| "Denta" 120 mL | same | calculation UNDERESTIMATES by 11 % |
| Petri 60 mL | same | calculation UNDERESTIMATES by 19 % |
| point source, 5 and 25 cm | same | agreement 0,97 and 0,93 |

The deviations of the volumetric geometries are of **opposite sign**, which
means this is not a common normalization. The decisive argument: K-40 — a
clean isolated line with no chains or blends — gives 1,22 / 0,82 / 0,78
across the three vessels; the geometry doesn't know which nuclide is inside,
so it is specifically the vessel models that diverge.

**Confirmed independently, 28.07.2026.** Recalculating the entire kit
against the source passports (53 lines,
`detectors/Gamma-1S/results/kit_recalc_*.csv`) gives medians A_meas/A_pass:
Marinelli 0,80, "Denta" 1,19, Petri 1,30, point sources 1,17 and 1,18. The
sign is the same; the `.efr` curve is not involved in this calculation.
Significantly, after fixing the time-based event separation (see the
Gamma-1S report), the spread WITHIN each vessel collapsed: previously in the
Marinelli, cesium gave 0,78 and thorium 1,26 — now the whole vessel lies
within 0,69–0,89. Each geometry is now left with a single factor, and the
cause to look for is one per vessel, not a nuclide dependence.

**What's needed to close this:** drawings of the vessels and of the
Marinelli beaker. They don't exist. An intermediate step is a run at several
well depths (`detectors/Gamma-1S/drivers/run_probe.py`).

## Printed RadiaCode vessels: sources and license

The vessels for RadiaCode 101/102/103 are not laboratory vessels but printed
models by the author **dnpro**, published on Thingiverse under the **CC BY**
license (Creative Commons Attribution). The geometry in
`detectors/RadiaCode-103` is derived from these STL files, so it is a
derivative work and requires attribution.

| vessel | publication | date | files |
|---|---|---|---|
| 200 mL, v.2 | [thing:6562353](https://www.thingiverse.com/thing:6562353) | 03.04.2024 | `v.2_-_Can_200_ml_-_RC-Mr.stl`, `v.2_-_Cap_200_ml_-_RC-Mr.stl` |
| 500 mL, v.1 | [thing:6325102](https://www.thingiverse.com/thing:6325102) | 20.11.2023 | `Marinelli_Beaker_for_RadiaCode_101.stl`, `Marinelli_Beaker_cap_for_RadiaCode_101.stl` |

The STL files themselves are not included in the repository: they are
obtained from the author via the link. Dimensions are obtained with
`detectors/RadiaCode-103/analysis/measure_stl.py` — by ray scanning the
triangle mesh, without assumptions about the body's shape.

**How this differs from the LSRM vessels above.** There, the wall and well
dimensions are absent and have to be assigned as assumptions — hence the
unresolved discrepancy across the three vessels. Here the geometry is taken
from the very same file the vessel was printed from, so it contains almost
no assumptions; what remains an assumption is marked in `RCDetector.hh` with
the word `ДОПУЩЕНИЕ` (ASSUMPTION) (fit clearance, unmodeled thread and
ribs).

**Independent confirmation of the measurement.** The author, in the
thing:6562353 discussion (11.05.2024), states the design thicknesses: outer
wall 2,8 mm, inner wall 1,2 mm. The STL measurement gave 2,81 and 1,25 mm.
The agreement is significant: these are two different sources — the
author's statement and the measurement of the delivered geometry.

**Material.** The author recommends Pet-G, but the 500 mL vessel used for
the measurements was printed by the operator from PLA — the model's default
corresponds to this (for the 200 mL vessel, the specimen's material is not
separately confirmed). Anyone who printed following the author's advice
needs to set the sixth argument of `rc_curves` to `PETG`. The choice has
little effect: the linear attenuation coefficient of the two plastics
differs by 0,8 % at 100 keV, which over the path through the 2,8 mm end wall
gives 0,04 % transmission — below the statistics of the runs.

**Infill.** The author prints the well wall at **100 % infill**, and this is
the only wall between the sample and the crystal, so solid density there is
justified — previously this was an unstated assumption. The outer wall is
printed with less than 100 % infill; the model slightly overestimates its
attenuation, but the path through it is of secondary importance for the
response.

## Effective thickness: what it is not

The parameter `d` in the correction `f(mu·rho·d)` is a fitting parameter,
not a geometric length. For a flat vessel, the mean chord through the layer
under isotropic emission is known to be greater than the fill height, while
the fit for the "Denta" gives less. The ratio of `d` to fill height across
the three vessels: 1,02 / 0,92 / 1,18 — there is no constancy. The parameter
does its job (a single number describes the entire energy range at
chi2/dof ≈ 1), but it cannot be called a physical length.
