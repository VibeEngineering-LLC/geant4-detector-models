# Sample matrices

Common to all detectors. Composition by mass, fractions.

## OISN-16

A bulk-sample simulant on an iron base — **not soil and not organic
material**: 71 % iron by mass. LSRM builds efficiency curves for volumetric
geometries on this matrix at ρ = 1,6 g/cm³.

| element | H | C | N | O | Fe |
|---|---|---|---|---|---|
| fraction | 0,022 | 0,206 | 0,009 | 0,049 | **0,714** |

Source: the `Material` field in the `.efa` files of the calibration kit.

Practical consequence: transferring the efficiency from OISN-16 to water or
organic material using the self-absorption formula is extrapolation, not
recalculation. Iron sharply changes the trend at the soft edge. Calculate
directly in the target matrix.

## RISN-379

The matrix of the mixed source of the calibration kit: a light
organo-mineral base with a noticeable amount of calcium, ρ = 1,0 g/cm³.

| element | H | C | N | O | Na | Mg | Ca |
|---|---|---|---|---|---|---|---|
| fraction | 0,043 | 0,330 | 0,012 | 0,348 | 0,041 | 0,022 | **0,203** |

Source: the `MATERIAL` field in the original `.spe` files ("Verification
2016"). The XML versions of the kit do not contain the composition — only
the original files do.

Calcium is significant at 59,5 keV (the Am-241 line): calculating the
mixture using water instead of RISN-379 is not permitted without
verification.

## Water

Distilled water — the matrix for which instrument passports specify the MDA
(minimum detectable activity). Implemented as H₂O with a settable density.

## How to add a matrix

`G1SDetector::MakeMatrix` takes the matrix name, density, and the Geant4
material name. Density is a parameter, composition is fixed. The material
name in the model is always `Sample` (a role), not the name of the
composition: the matrix is switchable.
