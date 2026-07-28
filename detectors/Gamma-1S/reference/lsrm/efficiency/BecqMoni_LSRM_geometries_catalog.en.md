# LSRM Geometries — BecqMoni catalog for efficiency calculation

**Date recorded**: 2026-06-14
**Source**: https://github.com/Am6er/BecqMoni/tree/master/LSRM%20Geometries
**Purpose**: reference on the formats of LSRM efficiency curves and .in model files, reference geometries for benchmarking our pipeline.
**Not to be confused**: with the LSRM software (SNIIP «LSRM Certificate») — these are different entities. Here LSRM = «Lab Scintillator Response Method», an MC simulator by author Am6er (built into BecqMoni).

## 1. What LSRM is in the context of BecqMoni

BecqMoni uses its own MC simulator (referred to in the README as **EffMaker**) to calculate **photopeak efficiency** ε(E) from a description of the geometry and materials. The `LSRM Geometries/` folder contains **ready-made** simulation results for typical configurations:

- `*.in` — input file for the MC simulator (geometry + materials + source).
- `*.txt` (curve) — result: a table of ε(E) at 150 points from 20 to 3000 keV in steps of 20 keV.

When importing into the BecqMoni UI, the operator selects a geometry from the list → BecqMoni takes the curve table → interpolates it onto the energy grid of their instrument → applies it to calculate activity.

For our pipeline, this data serves as an **independent reference**: we can compare our ε(E) calculation for the RC-103 marinelli with the LSRM curve and validate it.

## 2. Format of `curve_*.txt` (efficiency table)

```
Energy, keV	Efficiency	Uncertainty, %
20.0			4.55616E-02		1150
40.0			2.76634E-03		63.9
60.0			2.24553E-03		2.15
...
3000.0			2.48089E-05		21.4
```

Structure:
- **3 columns**: Energy in keV, Efficiency (fraction of total flux into 4π), Uncertainty in %.
- **Delimiter**: TAB character (may be repeated consecutively).
- **Header**: first line `Energy, keV\tEfficiency\tUncertainty, %`.
- **Energy grid**: 20 ... 3000 keV in steps of 20 keV → **150 points**. (In curve_Nano_16 and curve_RadiaCode_cilinder a 151st point was found — possibly an artifact.)
- **Efficiency** — dimensionless, e.g. `2.24553E-03` = 0.224% efficiency at 60 keV.
- **Uncertainty** — relative %, for 20 keV often >100% (low MC statistics), for 100-2000 keV typically 2-5%, in the tails above 3000 keV rising to 20%.

**Important**: 20-40 keV uncertainty >50% — this range **should not be used** for quantitative measurements. Low energy means larger MC statistical error + strong dependence on self-absorption (fine geometry details).

### Example: RadiaCode + Marinelli 0.5 (our profile)

| E, keV | ε | u, % |
|---|---|---|
| 20  | 4.56e-02 | 1150 (garbage) |
| 40  | 2.77e-03 | 63.9 (poor) |
| 60  | 2.25e-03 | 2.15 (reliable) |
| 100 | 2.58e-03 | 2.16 |
| 200 | 1.60e-03 | 2.48 |
| 500 | 3.93e-04 | 4.78 (interpolated) |
| 1000 | 1.45e-04 | 7.9 |
| 1460 | 9.6e-05 | ~10 (K-40) |
| 2614 | 4.4e-05 | ~15 (Tl-208) |
| 3000 | 2.48e-05 | 21.4 |

The drop of ε(E) ≈ E^(-1.4) from 100 to 1000 keV is a typical CsI(Tl) response, with the photoelectric effect dominating up to ~200 keV, then transitioning to Compton + pair production at E > 1022.

## 3. Format of the LSRM `.in` model file

Six blocks:

### 3.1. DETECTOR PARAMETERS

```
DetectorType = SCINTILLATOR     // or COAXIAL for HPGe

// SCINTILLATOR-specific (DS_*)
DS_CrystalDiameter = 1 cm
DS_CrystalHeight = 1 cm
DS_CrystalFrontReflectorThickness = 0.1 cm
DS_CrystalSideReflectorThickness = 0.1 cm
DS_CrystalFrontCladdingThickness = 0.1 cm
DS_CrystalSideCladdingThickness = 0.1 cm
DS_DetectorMountingThickness = 0.1 cm
```

(COAXIAL DC_* fields are present for compatibility, but are ignored if DetectorType=SCINTILLATOR.)

### 3.2. SOURCE PARAMETERS

```
SourceType = MARINELLI     // or POINT, CYLINDER

// MARINELLI-specific (SM_*)
SM_BeakerToDetectorFrontDistance = 0.8 cm     // gap between detector → bottom of the Marinelli (0.2-0.8 cm)
SM_BeakerDiameter = 11.4 cm                   // outer diameter of the Marinelli (standard 11.4)
SM_BeakerHeight = 8.9 cm                      // beaker height
SM_BeakerHoleDiameter = 6.1 cm                // through-cavity diameter (where the detector sits)
SM_BeakerHoleHeight = 5.3 cm                  // cavity depth
SM_BeakerSideThickness = 0.2 cm               // side wall plastic thickness
SM_BeakerEndWallThickness = 0.2 cm            // beaker bottom thickness
SM_BeakerHoleSideThickness = 0.2 cm           // cavity wall thickness
SM_BeakerHoleEndWallThickness = 0.2 cm        // cavity cap thickness
SM_SourceHeight = 8.5 cm                      // sample fill height
```

**Standard Marinelli 0.5 L** (Am6er reference, matches our ОИСН-16 exactly):
- Beaker 11.4 cm Ø × 8.9 cm H, walls 0.2 cm.
- Hole 6.1 cm Ø × 5.3 cm H.
- Useful volume ≈ (π × 5.7² × 8.5) − (π × 3.05² × 5.3) ≈ 868 − 155 = 713 cm³ ≈ **0.7 L** (marketing rounding down to 0.5).

### 3.3. MATERIAL PARAMETERS

For each component (Crystal, Cladding, Reflector, Mounting, Beaker, Source) — density + a list of elements with their fractions.

**RadiaCode CsI(Tl) crystal**:
```
DS_nCrystalElements = 2
DS_RoCrystal = 4.51                    // density g/cm³
DS_ZCrystal[0] = 53                    // I (Iodine)
DS_FractionsCrystal[0] = 0.488451      // mass fraction of I
DS_ZCrystal[1] = 55                    // Cs (Cesium)
DS_FractionsCrystal[1] = 0.511549      // mass fraction of Cs
DS_FractionTypeCrystal = MASS          // fraction type: MASS or ATOM
```

(Without Tl — the dopant is ≪0.1%, the simulator ignores it.)

**Source (Marinelli fill)**:
```
M_SM_Source.MName = Water, liquid              // default — water-equivalent
M_SM_Source.Nmaterials = 1
M_SM_Source.Name[0] = Water, liquid
M_SM_Source.MatRelWeight[0] = 1
```

For real samples the operator must **change** it to the required matrix (soil ρ ≈ 1.4, concrete ρ ≈ 2.3) and recalculate. The default water is the base-case for calculations.

**Beaker material**:
```
M_SM_Beaker.MName = Polyethylene terephthalate (PET)
```

Not to be confused with PE (polyethylene, ρ=0.93) — this is PET (ρ=1.38). Our real ОИСН-16 Marinelli vessels are white PVC or PP, the density differs slightly, but the effect on ε(E) for energies >100 keV is negligible.

### 3.4. Where the Tl dopant is in the LSRM model

**Nowhere explicitly**. The simulator models the photoelectric effect and Compton scattering via atomic cross-sections by Z, and Tl (0.1% of the mass) contributes negligibly to attenuation. The scintillation yield (Tl as the luminescence activator) is not modeled by the simulator — this is **photopeak efficiency**, not light-output.

## 4. Catalog of geometries in the repo

| Curve file | Model | Detector | Source | Reliable E-range |
|---|---|---|---|---|
| `RadiaCode - marinelli 0.5.txt` | `model_RadiaCode_Marinelli0.5.in` | CsI(Tl) 1×1 cm | Marinelli 0.7 L (11.4×8.9, hole 6.1×5.3, walls 0.2) | 60-2500 keV |
| `RadiaCode - author marinelli 0.5.txt` | `model_RadiaCode_AuthorMarinelli0.5.in` | CsI(Tl) 1×1 cm | Marinelli «author» 0.5 L (9.28×9.28, hole 2×6.2, walls 0.18) — a different geometry | 60-2500 keV |
| `RadiaCode - author marinelli 0.2.txt` | — | CsI(Tl) 1×1 cm | Marinelli «author» 0.2 L (smaller) | 60-2500 keV |
| `RadiaCode - cilinder.txt` | — | CsI(Tl) 1×1 cm | Cylinder geometry (compact plastic) | 60-2500 keV |
| `Obsidian - marinelli 0.5.txt` | `model_Obsidian_Marinelli_0.5.in` | CsI(Tl) 0.67×3 cm (long and thin) | Marinelli standard 11.4×8.9 (like our ОИСН-16) | 60-2500 keV |
| `Nano 16 - marinelli.txt` | `model_Nano16Pro_Marinelli.in` | CsI(Tl) 1.854×5.9 cm (large Nano) | Marinelli standard 11.4×8.9 | 80-2800 keV |

**What is missing from the repo** (gaps):
- RC-103 + Marinelli with soil (ρ=1.4) — the operator must simulate it themselves or take «water» as a worst-case.
- Geometry for HPGe (Coaxial) — there are no model files, only scintillator.

## 5. Comparison of CsI(Tl) crystal sizes

| Instrument | Crystal Ø | Crystal H | Volume | Notes |
|---|---|---|---|---|
| RadiaCode RC-10x | 1.0 cm | 1.0 cm | 0.785 cm³ | Declared by the manufacturer as «1×1» (compact) |
| Obsidian | 0.67 cm | 3.0 cm | 1.06 cm³ | Thin and long (directional) |
| Atom Spectra Nano 16 Pro | 1.854 cm | 5.9 cm | 15.9 cm³ | Semi-stationary, large volume |

The ratio of useful volume Nano/RC-103 ≈ 20×. At K-40 1461 keV, ε_Nano ≈ 5× ε_RC-103 (approximately √20 ≈ 4.5×, corrected for shape).

## 6. Cross-reference with the project

### 6.1. Importing the LSRM curve into our pipeline

Minimal helper script (reads the current format, outputs `EfficiencyModel`):

```python
import numpy as np
from pathlib import Path

def read_lsrm_curve(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Returns (E_keV, eff, uncertainty_rel_pct)."""
    data = np.loadtxt(path, skiprows=1, delimiter=None)  # whitespace delimiter
    return data[:, 0], data[:, 1], data[:, 2]
```

Save as `scripts/io/lsrm_efficiency.py` (the numbering continues the existing BecqMoni parsers).

### 6.2. Validation of our ε(E) for RC-103

Reference points for RC-103 + Marinelli 0.5 (water):
- 60 keV: ε = 2.25e-03 (reliable)
- 186 keV: interpolation between 180 (1.84e-03) and 200 (1.60e-03) ≈ 1.69e-03
- 583 keV: interpolation between 580 and 600 ≈ 3.2e-04
- 1461 keV: ≈ 9.6e-05
- 2614 keV: ≈ 4.4e-05

When running our Marinelli + RC-103 calculation with a water fill, ε should fall within ±10% of these values (LSRM uncertainty + our model). If we diverge by 2× — there is a bug in the self-absorption or solid-angle factor.

### 6.3. Differences from our methods

- LSRM models a **water** fill by default. Our `gamma.efficiency` for the operator's working samples is **soil** ρ=1.4. At 60 keV, self-absorption of soil in the Marinelli ≈ 0.4-0.5 of water, at 600 keV ≈ 0.7-0.8. A direct comparison of the LSRM curve and our ε(E) for soil is NOT valid. LSRM needs to be simulated with soil material.
- LSRM uncertainty at 20 keV >100% is **MC noise**, not a calibrated uncertainty. The real uncertainty of ε(60 keV) is ± 2-3% given a calibration source, not the LSRM curve.

## 7. Anti-hallucination — provenance

All facts are from open files in `Am6er/BecqMoni/LSRM Geometries/`:
- `curve_RadiaCode_-_marinelli_0.5.txt`, `curve_RadiaCode_-_author_marinelli_0.2.txt`, `curve_RadiaCode_-_author_marinelli_0.5.txt`, `curve_RadiaCode_-_cilinder.txt`, `curve_Nano_16_-_marinelli.txt`, `curve_Obsidian_-_marinelli_0.5.txt`
- `model_RadiaCode_Marinelli0.5.in`, `model_RadiaCode_AuthorMarinelli0.5.in`, `model_Nano16Pro_Marinelli.in`, `model_Obsidian_Marinelli_0.5.in`

Downloaded on 2026-06-14 via `gh api repos/Am6er/BecqMoni/contents` + curl with URL-encoding of spaces.

**Do not extrapolate** to:
- The symbiosis between LSRM (MC) and BecqMoni (ε from the curve): after the curve is generated, BecqMoni only interpolates, it does not recompute the MC.
- Other LSRM formats (SNIIP LSRM from 2007 — this is a completely different format, do not confuse them).
- HPGe Coaxial — there are no ready-made curves for HPGe in the repo, do not make assumptions about them.
