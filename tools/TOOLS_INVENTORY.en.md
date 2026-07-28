# Tool registry — geant4-detector-models

> The canonical list of applied tools in the repository. **Before writing a
> new script, check this file and the registry of the adjacent project**
> [spectravibe-toolkit/scripts/TOOLS_INVENTORY.md](https://github.com/VibeEngineering-LLC/spectravibe-toolkit/blob/main/scripts/TOOLS_INVENTORY.md):
> if a suitable tool already exists — use it, extending it if necessary,
> rather than rewriting it from scratch.
>
> Tools are not deleted. Obsolete ones are marked `LEGACY`; the files remain.

---

## Publication and security

### `tools/check_paths.py`
- **Purpose**: checking the repository before a commit to the public
  repository — absolute paths, surnames and initials, serial numbers,
  email addresses, tokens.
- **Why mandatory**: once a local path or surname gets into the git
  history, it stays there forever; it cannot be scrubbed after the fact.
- **Status**: ACTIVE.
- **Run**: `python tools/check_paths.py [root]`; returns 1 on findings.

### `tools/check_csv.py`
- **Purpose**: read every table in `detectors/*/results/` with the
  standard `csv.reader` and check that the number of fields in each row
  matches the header.
- **Why mandatory**: `results/` is what the external consumer reads, and
  it reads it with a plain `csv.DictReader`. An unquoted comma inside a
  value doesn't break the read — it SILENTLY shifts the fields: in
  `efficiency_curves.csv`, 48 rows of both point-source curves ("точечный,
  25 см" — "point, 25 cm") returned `E_keV` = "ОТКРЫТА" ("OPEN"), and
  `eps_net` = the energy value. The volumetric curves, meanwhile, read
  correctly, which is why the defect went unnoticed.
- **Where the check lives**: `common/py/csvio.py` — ONE implementation for
  the whole repository. It is also invoked on every write (`csvio.write`),
  and it is also what forbids a comma inside a value: the consumer parses
  the export by position, and a quoted field doesn't help it. A copy of
  the guard in `export_curves.py` had drifted from this one in how it
  handled comments — it has been removed.
- **Status**: ACTIVE.
- **Run**: `python tools/check_csv.py [file...]`; returns 1 on findings.
  Included in `recalc_all.py` ALWAYS, including with `--only`.

### `tools/anonymize.py`
- **Purpose**: anonymizing spectra and passports — removing the
  `OPERATOR` field, replacing instrument serial numbers with `SN-nn`,
  source certificate numbers with `SRC-nn`, both in the content AND in
  file names.
- **Feature**: numbers not listed in the map receive auto-generated
  aliases `SRC-Ann`; the mapping to the genuine numbers **is not
  published**.
- **Format protection**: before editing, it checks on one `.spe` that the
  sum of counts before and after matches — in an LSRM `.spe`, the counts
  are packed in binary.
- **Status**: ACTIVE.
- **Run**: `python tools/anonymize.py <directory>`; `SPECTRAVIBE_ROOT` is
  required for the format check.

---

## Obtaining reference data

### `tools/fetch_efr.py`
- **Purpose**: downloading efficiency curves `.efr`/`.efa` from
  spectravibe-toolkit, with CP1251 re-encoding and point parsing.
- **Status**: ACTIVE.

### `tools/fetch_kit_xml.py`
- **Purpose**: downloading the calibration kit in BecqMoni XML (all
  geometries).
- **Feature**: downloads via the git tree and blob SHAs, rather than via
  `contents/<path>` — the paths contain Cyrillic, and URL-encoding breaks
  it.
- **Status**: ACTIVE.

### `tools/fetch_bg.py`
- **Purpose**: downloading the averaged background spectra of the setup.
- **Status**: ACTIVE.

> The Gamma-1S reference data is stored in full in the repository
> (`detectors/Gamma-1S/reference/lsrm/`), so downloading is needed only for
> updates or for other detectors.

---

## Parsing library — `common/py/`

### `common/py/paths.py`
- **Purpose**: path roots from environment variables; there is not a
  single machine-bound path anywhere in the code.
- **Variables**: `GEANT4_ROOT`, `G4MODELS_BUILD`, `G4MODELS_REF`,
  `G4MODELS_MEASURED`, `SPECTRAVIBE_ROOT`.
- **Status**: ACTIVE.

### `common/py/becqmoni.py`
- **Purpose**: reading BecqMoni XML spectra, peak areas with a
  trapezoidal background, measuring the FWHM from the peak itself,
  **broadening the modeled spectrum to the instrument resolution**
  (`broaden`, `area_broadened`).
- **Key point**: the area from the model and from the measurement must be
  taken with the SAME window, otherwise blends are computed differently
  (an error of up to 1,5 times).
- **Adjacent counterpart**: `gamma/io/becqmoni_xml.py` in
  spectravibe-toolkit — cross-check against it when developing further.
- **Status**: ACTIVE.

### `common/py/contam.py`
- **Purpose**: checking the purity of background shelves — whether a
  strong neighboring line falls into them. Contaminated points are
  **excluded**, not fitted around.
- **No counterpart found in spectravibe-toolkit** — an addition specific
  to this repository.
- **Status**: ACTIVE.

---

## Model tools

Instrument-specific scripts live in `detectors/<instrument>/`: `drivers/`
— Geant4 runs, `analysis/` — analysis. The list and purpose are given in
the README and report of the corresponding model.

Worth noting separately are the ones that are portable between detectors
and should be promoted to `common/py/` for the next model:

| script | what it does | portability |
|---|---|---|
| `Gamma-1S/analysis/selfabs_fit.py` | effective thickness from a pair of densities | high |
| `Gamma-1S/analysis/summing.py` | cascade summing via a run + control | high |
| `Gamma-1S/analysis/tcc_evidence.py` | whether the correction was introduced into the reference curve | high |
| `Gamma-1S/analysis/loading.py` | dead time and pile-up from the inventory | high |
| `Gamma-1S/analysis/mda.py` | declared MDA and Currie's MDA | high |
| `Gamma-1S/analysis/mix_unfold.py` | NNLS unfolding of a mixture by templates | medium |
| `Gamma-1S/analysis/export_curves.py` | exporting full-energy-peak (FEP) efficiency curves to `results/` tables | high |
| `Gamma-1S/analysis/xml_vs_spe.py` | cross-checking BecqMoni XML against the original LSRM .spe | high |
| `RadiaCode-103/analysis/build_article.py` | rebuilding the article from the template and `curves.json` | medium |
| `RadiaCode-103/analysis/nucdata.py` | line yields and half-lives for nuclide analysis | high |
| `RadiaCode-103/analysis/make_table.py` | efficiency table in markdown | high |
| `RadiaCode-103/analysis/check_results.py` | cross-checking the final tables against each other | high |

> An adjacent counterpart for the unfolding is the quasi-template
> full-spectrum WLS (LSRM §13) in spectravibe-toolkit, `gamma/activity/`.
> Compare approaches when developing further.


## Anonymization: the map lives OUTSIDE the repository

`tools/anonymize.py` does not contain a single genuine surname or a single
genuine number — the correspondence map is the key to reverse
identification, and if published together with the data it defeats the
anonymization. The path to the map is set by the `ANON_MAP` variable; the
format is described at the top of the script.

```bash
ANON_MAP=<map.json> SPECTRAVIBE_ROOT=<root>     python tools/anonymize.py detectors/Gamma-1S/reference --verify
ANON_MAP=<map.json> python tools/check_paths.py
```

The directory must be specified explicitly and must be a DATA directory:
running it over the repository root corrupts the markup and the source
files (in markdown `#` is a heading, in CSS a color, in Python a comment)
and reaches even the tool itself. Both programs now refuse to run under
such conditions, but this should not be relied upon.

`--verify` requires `SPECTRAVIBE_ROOT`: the integrity of every edited
`.spe` is checked with the standard LSRM reader — the sum of counts before
and after must match.

Without `ANON_MAP`, the `check_paths.py` check works only from general
patterns and warns about this. The precise check against the map is the
only reliable one.
