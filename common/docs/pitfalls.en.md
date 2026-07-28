# Pitfalls

Each of the items listed below has already cost an error of several-fold or
days of compute time. All of them are silent: the calculation runs, the
numbers look plausible, and the error is only found by cross-checking
against an independent source.

## Spectrum analysis

### Peak area from the model and from the measurement — with the SAME window

**Symptom.** The activity from one line of a nuclide is systematically
higher than from other lines of the same nuclide.

**Example.** Ac-228 911.2 and 968.97 keV are separated by 58 keV, and the
FWHM at 911 keV for NaI happens to be exactly 58 keV — in the spectrum this
is ONE peak. From the measurement, a ±1 FWHM window picks up both lines;
from the model (unbroadened) spectrum, a narrow ±6 keV window picks up only
one. The thorium activity from this line was overestimated by a factor of
**1.5**.

**Fix.** Broaden the model spectrum to the instrument's resolution and take
the area with the same window (`common/py/becqmoni.py`: `broaden`,
`area_broadened`). Even more reliable — take efficiency **per decay directly
from a full-decay run**: then line yield, cascade summing, blends and
continuum are accounted for automatically and identically on both sides.

### Trapezoidal background shelves contaminated by neighboring lines

**Symptom.** One line falls out of line with the rest, while the others are
consistent.

**Example.** Th-232 583.2 keV: the left continuum shelf picks up the Tl-208
510.77 keV line, the continuum is overestimated, and the area is cut by
**30%**.

**Fix.** `common/py/contam.py` checks the shelves against the nuclide's line
list; contaminated points are **excluded with a message**, not fitted
around. Important: if the area is taken with the same broadened window on
both sides, the contamination shrinks on its own — only schemes where the
model side is taken from a mono-energy grid or from the nameplate activity
are dangerous.

### Window width — the dominant systematic on area

For Cs-137 the recovered efficiency depends on the window: ±1.0 FWHM gives
+1.6% relative to the reference value, ±1.25 → +6.3%, ±1.5 → +8.1%. A
Gaussian does not account for this much (98.2% vs. 99.98% of the area) — the
remainder comes from under-subtracted continuum at a wide window. Fix the
window and state it in the report.

### Live time vs. real time

In BecqMoni XML, `MeasurementTime` is REAL time, `LiveTime` is live time.
On a recording with a 70 cps load they diverge by 2%.

### The decay correction can silently become unity

The `gamma.io.lsrm_spe` reader returns `None` in the `end_datetime` and
`file_created_datetime` fields. Without an explicit check, the decay
correction comes out equal to 1. For Cs-137, 9 years past its nameplate
date, this is **18%** — and this is exactly what created an apparent
contradiction between two reference points.

### Full-spectrum fitting: chi2/dof is not a quality measure

At a million counts, purely Poisson weights hand the fit over to a few of
the strongest peaks, and a three-percent shape inaccuracy inflates chi2/dof
into the hundreds. A systematic weight floor is needed:
`sigma^2 = y + background + (0.03*model)^2`. The quality measure is not
chi2, but the **stability of amplitudes under a change of weighting scheme**
(observed 9–13%) and consistency across geometries.

### A fit that has run to a parameter boundary is not a result

When part of the templates is missing, the optimizer runs to the boundary
and zeroes out components, producing amplitudes that look plausible on the
surface. Check the completeness of the input data BEFORE fitting and reject
output that has run to a boundary.

## Geant4

### A long-lived nucleus will not decay without raising the threshold

By default the "very long decay" threshold is 1 year, and the primary
nucleus is simply killed. This affects K-40 (1.25·10⁹ years), Ra-226 (1600
years), and the thorium-chain members Ra-228 (5.75 y) and Th-228 (1.9 y).

```
/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+30 ns
```

Chain boundaries must be set with
`/process/had/rdm/nucleusLimits Amin Amax Zmin Zmax`, otherwise, with the
threshold raised, things will decay that do not decay in the real source
(Np-237 in Am-241).

### Cascade angular correlations are disabled by default

Geant4 prints at startup `Enable correlated gamma emission 0`. For a
Marinelli geometry the contribution is small; for a point geometry it would
be more significant.

### The GPS sampling region must cover the ENTIRE source body

If the sampling cylinder is smaller than the sample volume, the source sits
only in its inner part, the efficiency comes out overestimated, and the
error is silent.

### A sampling cone is valid for the full-absorption peak but not for total efficiency

Restricting the solid angle saves compute time for distant point sources
(at 25 cm, 99.4% of events are wasted). Essentially only direct quanta reach
the full-absorption peak, so this is legitimate for that peak. Total
efficiency with such sampling is underestimated — quanta scattered off the
shielding are not sampled. For cascade runs a cone is FORBIDDEN: coincidences
depend on the directions of both quanta.

### Tessellated solids wreck performance

`G4ExtrudedSolid` is a descendant of the tessellated solid. Four such solids
gave 1.7 thousand events/s versus 12 thousand on primitives. Assemble a
rounded rectangle from boxes and cylinders instead.

## Tools

### `subprocess` with `text=True` and no `encoding`

Decodes output using the locale's codepage and fails on the first byte
outside the table. The secondary error `'NoneType' has no attribute
'splitlines'` masks the real cause. Always use `encoding="utf-8",
errors="replace"` and `(r.stdout or "")`.

### BOM in a macro

`Out-File -Encoding utf8` in PowerShell writes a BOM, on which Geant4
stumbles: the first command becomes unrecognized, and the message looks
like a syntax error. Write without a BOM.

### A one-off "in-place" edit not carried back into the script

Makes the result irreproducible: the numbers sitting on disk are correct,
but a rerun will silently produce different ones. Parameters pulled out for
individual runs (reflector density, well depth) must be **printed on every
run**.

### An import block inserted after its first use

**Symptom.** The script fails from a clean clone with `NameError: name
'paths' is not defined` or `ModuleNotFoundError: No module named
'becqmoni'`, even though the files are present.

**How it happened.** Converting paths to environment variables was done
with an automated replace, and the line `import paths` ended up BELOW the
first `paths.build(...)`; likewise with `common/py` in `sys.path` and
`import becqmoni`. On the author's machine this never surfaced, because
the packages were found via a different path there.

**Fix.** Place the block that adds `common/py` to `sys.path` immediately
after the standard-library imports, before any of the project's own
imports. Verify not by eye, but by running each script from a clean clone.

### Data is not where the script looks for it

**Symptom.** The script runs and prints "found 0 records," without failing.

**How it happened.** The layout of the downloaded set is flat, while the
committed one is organized by geometry and nuclide. A hard-coded path
`ref()/filename` finds nothing in a clean clone, and `glob("*/*.xml")`
does not see files one level down.

**Fix.** Look up a file by name (`paths.find_data`) and a curve by geometry
name (`paths.efficiency_curve`), rather than concatenating a path. Do not
treat an empty result as normal: print an explicit error stating where the
search looked.

### Instrument files are in CP-1251, and the script opens them as UTF-8

**Symptom.** `UnicodeDecodeError: 'utf-8' codec can't decode byte 0xd3 in
position 1` on an `.efr`/`.efa` curve.

**How it happened.** The downloader saved the downloaded data as UTF-8, and
on the working copy everything read fine. The original LSRM files are
CP-1251.

**Fix.** `paths.read_text()`: try UTF-8 first, then CP-1251. Order matters:
UTF-8 text "successfully" decodes as CP-1251 too, but produces mojibake.
