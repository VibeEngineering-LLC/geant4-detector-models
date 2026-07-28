# Independent LSRM calculation (from the BecqMoni package)

Files copied from the open repository
[Am6er/BecqMoni](https://github.com/Am6er/BecqMoni), directory `LSRM
Geometries`, and are kept here unmodified — so that the cross-check is
reproducible and so it is visible exactly what the curves were compared
against. Rights to them belong to the BecqMoni/LSRM authors.

| File | What it is |
|---|---|
| `RadiaCode_AuthorMarinelli0.2.in` | model: RadiaCode in the author's 0.2 L Marinelli beaker |
| `RadiaCode_AuthorMarinelli0.5.in` | model: same, 0.5 L |
| `RadiaCode_Marinelli0.5.in` | model: RadiaCode in a **classic** 0.5 L Marinelli beaker |
| `RadiaCode - author marinelli 0.2.txt` | exported ε(E) curve, 20–3000 keV |
| `RadiaCode - author marinelli 0.5.txt` | same |
| `RadiaCode - marinelli 0.5.txt` | same, classic Marinelli beaker |
| `RadiaCode - cilinder.txt` | same, cylindrical beaker |

## How the LSRM model differs from mine

The breakdown is needed so that curve discrepancies read as the consequence
of specific differences, rather than as an "unexplained residual."

| Parameter | LSRM | here |
|---|---|---|
| Crystal | **cylinder** Ø10 × 10 mm, 0.785 cm³ | cube 10 × 10 × 10 mm, 1.000 cm³ |
| Well | axially symmetric **Ø20 mm** | 36.7 × 19.8 mm (0.5 L), 37.2 × 20.7 (0.2 L) |
| Instrument case | absent | ABS 1.5–2.7 mm + board, display, battery |
| Crystal wrap | 1 mm TiO₂ (ρ 4.26) + 1 mm shell + 1 mm mount | 1.25 mm polymer with TiO₂ (ρ 1.45) |
| Sample | water **1.0 g/cm³** | grid 0.0012…1.60 g/cm³ |
| 0.5 L cup wall | 1.8 mm | 2.0–12.0 mm (STL measurement, variable) |
| 0.5 L inner diameter | 89.2 mm | 86.7 mm (STL measurement) |

The key to these differences: **the LSRM code is axially symmetric**. A
rectangular instrument in a rectangular well cannot be expressed in it, so
the crystal became a cylinder, and the well became a Ø20 mm pipe. In other
words, the cylinder in their file is not information about the instrument,
but a limitation of the model: in the teardown photo
(`drawings/rc_teardown_crystal_window.jpg`) the reflective cup is **square
in cross-section**, with a square 6 × 6 mm window for the SiPM.

Within LSRM itself, the description of the instrument is not consistent
between files: in the author's Marinelli beaker the wrap is specified as
polyethylene 0.93 + TiO₂ 4.26, while in the classic one it is aluminum 2.7 +
MgO 2.25. This limits their curves as an absolute reference.

The exported values are strictly monotonic across all 150 points at a
stated uncertainty of 15–22% at high energies. Raw Monte Carlo cannot look
like that, meaning what was exported is a **smoothed fit**, not calculated
points.

Cross-check: `python analysis/compare_lsrm.py m500`, result — in the
project README.
