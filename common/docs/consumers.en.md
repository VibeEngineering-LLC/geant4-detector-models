# Repository consumers

Since 28.07.2026 the repository has been connected as a submodule to a
third-party project. A standalone repository can be rewritten freely; a
connected one cannot. Below are the four obligations that did not exist
before, and the editing procedure that follows from them.

## Who takes what

`Am6er/BecqMoni`, commit `cbc826a`, file `.gitmodules`:

```
[submodule "geant4-detector-models"]
	path = geant4-detector-models
	url = https://github.com/VibeEngineering-LLC/geant4-detector-models.git
```

Pinned at commit `b7867bc` (28.07, 11:58). Connected at the same level as
the certified curves from the instrument's factory package — that is, the
models and curves arrive together with the clone, `git submodule update
--init` and that's it.

Why the consumer needs this: for the RC101 and RC103 groups our CsI(Tl)
10×10×10 model is a **second independent source of shape** on top of the
certification; the Gamma-1S package is the very same one that fed the third
column of its comparison.

## Rule 1. History is frozen

No `push --force`, no `rebase` of published history, no history rewrite.

An externally pinned SHA disappears on a rewrite, and `git submodule update
--init` fails not for us but for every clone of the consumer. The 28.07
rewrite done for de-identification went unpunished only because it happened
before the connection was made.

From now on privacy is fixed with **new commits**: delete a file, extend the
template, add a rule. If a leak specifically requires history scrubbing —
that is an operator-level decision requiring agreement with the consumer,
not a routine edit.

## Rule 2. Fixes do not arrive on their own

The consumer is pinned to a commit, not a branch. Right now `2d74d79` and
`654f343` already sit on top of the pinned `b7867bc` — the consumer does not
see them.

Hence: once a block of fixes is closed, say so, so the pointer gets moved.
A silently fixed defect is not fixed at all as far as the consumer is
concerned.

## Rule 3. The order of fixes determines what the consumer reads

The consumer needs the curve's **shape**, not its level: the level is
absorbed by the free scaling factor during the activity fit. So a defect
that produces an offset is harmless to them, while a defect that produces a
slope corrupts exactly what the model is used for.

Hence the order for the remainder of the external audit on the RC side —
fix what bends the shape first:

**1. `.live` is read from `MeasurementTime`** — `detectors/RadiaCode-103/analysis/read_rcxml.py:70`.

The repository's own pitfall entry (`common/docs/pitfalls.md`, "Live time
vs. real time") states it plainly: `MeasurementTime` is REAL time,
`LiveTime` is live time; at a 70 cps load they diverge by 2%. The code reads
real time as live time, so all measured ε values are underestimated — and
underestimated **unevenly**: dead time grows with count rate, and count
rate differs across sources. A systematic effect turns into a slope. The
rule exists, the implementation is behind.

**2. `ash` matrix Σ = 1.020** — `detectors/RadiaCode-103/geometry/RCDetector.cc:232-235`.

Normalize the fractions and add a Σw check to `MakeMatrix` with
`G4Exception` — the failure mechanism is already used there. Attenuation in
the sample is energy-dependent, so a two-percent mass error produces a
slope, not an offset.

**3. `centre −0.44` in the macros** — all three: `macros/test.mac:14`,
`macros/bench.mac:13`, `macros/nuclides.mac:24`.

`seatGap` was changed 0.20 → 0.32, but the macros did not follow. The
correct value is −0.56. A 0.12 mm shift on a 10 mm crystal hits low
energies hardest — again a slope.

**4. Net-count variance `N_s + bgs` instead of `N_s + k²·bgs`** —
`analysis/normalization.py:58,66`, `analysis/fit_peak.py`,
`analysis/validate_bgsub.py`.

Does not bend the shape, but understates σ, so a "model vs. certification"
comparison will show better agreement than actually exists. The reference
for the correct derivation is `detectors/Gamma-1S/analysis/export_curves.py`,
where the same algebra has already been fixed.

**5. Normalization from two sources** — `analysis/compare_lsrm.py:260-261`
prints the multiplier 0.80, which is declared erroneous in `curves.py`.

Level, does not matter for shape, but two different numbers in one report
read as an error. Import from `curves`.

**6. `h*(10) = 4.13 pSv·cm²`** — `analysis/analyze_sens.py:7,41` — cross-check
against ICRP 74, Table A.21. The first-principles estimate ≈3.7 changes the
conclusion about the nameplate 30 cps figure (20–24 → 22–27). Concerns the
nameplate point in the README, not the curves, so it comes last.

## Rule 4. Known systematics are to be disclosed in advance

While items 1–3 remain open, a consumer comparing our RC model against
certification will attribute the discrepancy to the certification or to
their own algorithm — not to us. This is bad faith by default, not by
malice: they cannot know what we know.

Therefore: known unresolved systematics in the **pinned** commit are to be
communicated to the consumer before they run into them. Communicating
outward is the operator's job; the repository model's job is to keep this
list current and not let it go silently stale.

## What from the same audit is already closed

- **`nucleusLimits` windows in analog mode** — verified with a number, not
  by silence: `detectors/Gamma-1S/analysis/check_nucleuslimits.py`, discussed
  in `detectors/Gamma-1S/REPORT.md`. The Pa-233 lines at 300.1 + 311.9 keV
  gave 1 count against an expected 180,800, a suppression factor of
  1.8·10⁵. The control line (the Am-241 59.5 keV line itself: 0.361 against
  a reference value of 0.359) shows the run is sound. The branching method
  legitimately relies on the windows.
- **Privacy** — personal cloud paths and source serial numbers have been
  scrubbed, `check_paths` has been fixed in both places, history has been
  rewritten (the last time, see Rule 1).
- **Chain separation by time** — a global timestamp, sorting, a 1 µs cut;
  spurious inter-nuclide coincidences eliminated.
- **The sampling region follows the geometry**, `decay_control.mac` no
  longer overwrites good data with empty ones.

## What remains on privacy

`tools/anonymize.py`: `TOKEN` has been extended to cover "#" and wrapped
with `SKIP_EXT` rules, but `BARE` still requires the strict form `0NNN-NN`,
and there is no pattern at all for "nuclide plus bare number"
(`Th-232 420_3064`). The next import of factory certificates with such
numbers will pass straight through, and only a human will catch it.
Also: extend `--verify` with a second pass of `check_paths.scan()` over the
processed directory, with a nonzero exit code on findings — currently
verify only checks the integrity of `.spe` checksums.
