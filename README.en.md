# Geant4 Detector Models

Models of scintillation gamma spectrometers built from drawings and
datasheets, and **verified against measured data**, not just assembled.

The rule this repository holds to: every dimension has a source, and where
there is no source, the code contains the word `ДОПУЩЕНИЕ` (ASSUMPTION).
Numerical data — line yields, branching ratios, attenuation coefficients —
are not taken from reference books from memory, but are computed by the same
Geant4 that does the transport.

## Models

| detector | crystal | sample geometries | status |
|---|---|---|---|
| [Gamma-1S](detectors/Gamma-1S/) | NaI(Tl) 63×63 in 50 mm lead | 1 L Marinelli, "Denta", Petri dish, point sources at 5 and 25 cm | verified against 5 geometries; discrepancies localized in vessel models, cause not found — [report](detectors/Gamma-1S/REPORT.en.md) |
| [RadiaCode 101–103](detectors/RadiaCode-103/) | CsI(Tl) 10×10×10 | custom-authored 200 and 500 mL Marinelli vessels | curves built and cross-checked; **three protocol items not completed** — [details](detectors/RadiaCode-103/README.en.md) |

Both instruments were checked against LSRM curves in the Marinelli-vessel
setup and gave close normalization coefficients (0.858 and 0.833). **This is
not evidence of a systematic bias common to the method**: for the Gamma-1S,
point-source geometries showed that the detector model itself is correct,
and the excess belongs to the vessel geometry. Two numbers agreeing, obtained
in the same setup, proves nothing about the method — see the caveat in the
RadiaCode card.

## Repository layout

Sample geometries, matrices, and the parsing library **do not belong to any
one detector**: the 1 L Marinelli vessel and the OISN-16 matrix will be
needed by the next instrument, and spectrum reading and curve fitting are
needed by all of them. That is why they are placed at the top level rather
than inside the instrument's folder.

```
common/
  src/     model skeleton: physics list, spectrum accumulation, emitted-quanta counting
  cmake/   shared build template
  py/      parsing library: becqmoni.py, contam.py, paths.py
  docs/    methodology, pitfalls, verification protocol — shared across all models
detectors/<instrument>/
  geometry/ macros/ drivers/ analysis/ results/ reference/ drawings/
geometries/  shared sample geometries: Marinelli, "Denta", Petri dish, point source
materials/   shared sample matrices: OISN-16, RISN-379, water
```

```
tools/       de-identification, pre-publication checks, data downloaders
docs/        title page and the RadiaCode article
```

**Computed spectra are not committed.** There are hundreds of files and tens
of megabytes of them; instead, a run manifest and **ready-made curves** live
in `results/`. To read off the efficiency at a given energy, building Geant4
is not required.

**Reference data lives here in full.** The Gamma-1S verification kit, LSRM
curves, averaged backgrounds, and reference-source certificates are in
`detectors/Gamma-1S/reference/lsrm/`, in two formats: the instrument's native
binary `.spe` and BecqMoni XML. The formats are interchangeable, which is
checked by the `analysis/xml_vs_spe.py` script. The repository is
self-contained: external sources are needed only to update the data.

Names of measurement participants and serial numbers of instruments and
sources have been removed from the spectra and certificates. The mapping of
pseudonyms to real values is not published. Names in the reference lists are
retained — those are the authors of books and methodologies, not measurement
participants.

## Running

Paths are taken from environment variables; there is not a single
machine-bound path in the code. **All scripts run immediately after
cloning**, without a single variable set: if the data needed for a
calculation is missing, the script explains what is missing rather than
failing with a traceback.

| variable | what it sets | required |
|---|---|---|
| `GEANT4_ROOT` | Geant4 root for CMake | for building a model |
| `G4MODELS_BUILD_<DETECTOR>` | the calculation directory for one instrument, e.g. `G4MODELS_BUILD_GAMMA_1S` | no |
| `G4MODELS_BUILD` | the same, when there is only one instrument; also understood as a root with per-instrument subdirectories | no |
| `G4MODELS_REF` | directory with downloaded reference data; without it, the committed set is used | no |
| `G4MODELS_MEASURED` | the operator's personal measurements (not present in the repository) | only for RadiaCode |
| `SPECTRAVIBE_ROOT` | root of gamma-spectrum-analysis: the standard readers for LSRM formats | for reading raw `.spe` files |
| `ANON_MAP` | de-identification map, **outside the repository** | only when publishing |

```bash
cmake -S detectors/Gamma-1S/geometry -B build/Gamma-1S -G Ninja
cmake --build build/Gamma-1S
```

Next come the run drivers in `detectors/<instrument>/drivers/` and parsing in
`analysis/`. The order and meaning of each step is in the model's report.

## Before publishing

```bash
ANON_MAP=<map.json> python tools/check_paths.py
```

The check looks for local paths, names, serial numbers, email addresses, and
tokens — in file contents, including binary spectrum headers, and in file
names. Once such data lands in git history, it does not disappear from there
on its own. Details and de-identification rules are in
[tools/TOOLS_INVENTORY.md](tools/TOOLS_INVENTORY.en.md).

**History is frozen.** The repository is attached as a submodule to a
third-party project, and a history rewrite destroys the commit pinned
externally — it breaks not for us, but for everyone who clones the consumer.
The rules and editing order that follow from this are in
[common/docs/consumers.md](common/docs/consumers.en.md).

## What to read before building a new model

- [common/docs/validation.md](common/docs/validation.en.md) — the mandatory
  verification protocol. A model without it is not a model, just a picture.
- [common/docs/pitfalls.md](common/docs/pitfalls.en.md) — pitfalls, each of
  which has already cost a multi-fold error, with symptoms and remedies.
- [common/docs/consumers.md](common/docs/consumers.en.md) — who takes the
  repository as a submodule, what that forbids, and in what order to fix the
  rest of the audit findings.

## License and data provenance

Code — MIT. Instrument drawings and datasheets belong to their
manufacturers and are given to the extent needed for reproducibility of the
calculation. Reference spectra — LSRM and CJSC NPC "ASPECT", links in the
model reports.
