# detectors/Gamma-1C/

Detector-specific assets for the **Gamma-1C** spectrometric complex:

- **Detector head:** УДС-ГЦ-63×63 (NaI(Tl) 63×63 mm crystal, БДЭГ-63×63-USB) by Aspect
- **DAQ software:** Lsrm SpectraLine
- **Aliases:** Колибри-1М, Гамма-1С, БДЭГ-63×63 (see `data/aliases.json` → `detector.Gamma-1C`)

---

## 1. Taxonomy lock — Crystal vs Station (2026-06-05, user)

Two distinct levels — never conflate:

| Level | What it is | Carries… |
|---|---|---|
| **Crystal-class** | scintillator physics (type + dimensions) | universal nuclide line patterns (relative positions, BR-normalized intensities after FWHM normalization). Templates are transferable between all instances of this class. |
| **Station-instance** | a specific detector instance + electronics + serial | individual FWHM polynomial, efficiency curve, drift over time, dead-time. Calibration is not transferable. |
| **Resolution (FWHM curve)** | a property of a specific instance | differs for every NaI 63×63 (different crystals — different light yield, PMT currents). Affects the width of peak-search windows, but not the intensity ratios. |

**Gamma-1C as an entity**: the only physical **station-instance** in this
project (LSRM, BДЭГ-63×63-USB №SN-01 + БДЭГ-63×63-USB — both
serial numbers refer to **one and the same** station across different epochs / verification cycles),
and its **crystal-class = NaI 63×63**. See `audit/_rag/visual_templates/SCHEMA.md:10`
and `scripts/gamma/io/lsrm_spe.py:29` («DETECTOR=… — e.g. "Гамма-1С" (NaI 63×63)»).

### Canonical name fixation (HARD-LOCK 2026-06-05, user)

> «Гамма-1с (кириллица) = Gamma-1C (латиница)»

**Canonical = `Gamma-1C`** (homoglyph mapping: Cyrillic «С» → Latin «C»).
`detectors/Gamma-1C/` — the **only** folder for this station.

**Do NOT** create `detectors/Gamma-1S/`, `detectors/Гамма-1С/`, or any other
variants — these are duplicates of the same sub-tree.

`data/aliases.json:detector` contains *two* canonical tokens (`Gamma-1C` and
`Gamma-1S`) — this is a legacy v1.11.1 registry decision (`_meta.version`), kept
for backwards compatibility with .spe headers, where the LSRM converter transliterates the Cyrillic «С»
by sound → Latin «S». **For our station the authoritative canonical is `Gamma-1C`**.
BUG-40 (`KNOWN_AND_FIXED_ISSUES.md:1401-1422`) — a defensive warning for when
the canonicalizer catches a cyrillic-latin homoglyph ambiguity in an .spe header.

**Methodological consequence**: visual-templates for the NaI-63×63 crystal-class
are transferable to any future NaI-63×63 station (not only ours); each
template record's provenance specifies the station-instance (for downstream
efficiency/FWHM attribution). Currently all 24 canonical templates are observed-on
= `Gamma-1C`.

---

## 2. Shielding configuration

- **Outer shield:** **Pb 50 mm** — primary attenuation of cosmic background and external γ
- **Inner liner (graded shield):** **Cd + Cu cup** — suppresses Pb K-X-rays (73–87 keV)
  generated in the Pb shield by cosmic-ray excitation, cascading through Cd
  (K-XR 23 keV) → Cu (K-XR 8 keV)
- **Observable consequence:** without the Cd+Cu liner, the Pb K-XR peak in BG
  dominates the 73–87 keV region (~1400 counts in maximum vs ~700 continuum in
  a side-by-side comparison). With the liner, residual Pb K-XR in BG is
  ~0.1–0.3 cps — at the level of natural background continuum.
- **Cd Kα 23 keV and Cu Kα 8 keV are not seen in the BG spectrum** — both sit
  below the ADC threshold of this NaI 63×63 + USB front-end configuration
  (effective cutoff ~10–15 keV). Absence is not evidence against the liner.
- **For sample analysis:** the Pb K-XR cluster (73–87 keV) seen in any sample
  spectrum is dominated by **internal conversion X-rays from the source decays**
  (Pb-212 238 keV IC, Bi-212 / Tl-208 transitions IC) plus direct chain γ
  (Th-228 84.37 keV) — *not* by Pb fluorescence of the shielding.

The graded liner is what makes Gamma-1C a **low-background** instrument suitable
for trace-level identification at LSRM/ISO 11929 sensitivities.

---

## 3. Vessel canonicalization (LSRM source-of-truth, p. 11)

**Vessel canonicalization rule** (LSRM source-of-truth «Precision
Measurements. Reference and Calibration Sources», p. 11 +
user locks 2026-06-05):

### Constructive vs operational — three separate fields, not one tag

Previously mistakenly conflated into a single `geometry` string. Correct separation:

| Field | What it is | Source | Per-detector? |
|---|---|---|---|
| `vessel_class` | Constructive specification of the vessel (capacity, dimensions, d_eff) | LSRM source-of-truth | NO — transferable between detectors of the same crystal-class |
| `effective_thickness_mm` | d_eff for self-attenuation correction | LSRM source-of-truth | NO |
| `useful_sample_volume_ml` | Actually loaded sample volume for a given measurement geometry | Operator / passport / .spe header | **YES** — may differ |
| `placement_distance_cm` | Sample-to-detector distance | Operator / passport | **YES** |

### Canonical vessel-classes (per LSRM-source p. 11)

| `vessel_class` | Capacity | Dimensions (mm) | `effective_thickness_mm` |
|---|---|---|---|
| `marinelli_0.5L` | 0.5 L | ⌀125, H=100 | `[15, 2]` |
| `marinelli_1L`   | 1.0 L | ⌀150, H=110 | `[26, 2]` |
| `marinelli_3L`   | 3.0 L | ⌀180, H=200 | `[60, 5]` |
| `denta_120ml`    | 0.12 L | ⌀75, H=35 | `[36, 2]` |
| `petri_75ml`     | 0.075 L | ⌀88, H=14 | `[15, 2]` |
| `point_source`   | — | non-table | — |
| `other_<spec>`   | non-LSRM | passport | passport |

### Operator-data canonicalization (this FS layout → vessel tags)

| Folder name | `vessel_class` | `useful_sample_volume_ml` | Detector scope |
|---|---|---|---|
| `Маринелли` (89 records) | `marinelli_1L` (lab convention: default = 1L) | default = vessel capacity 1000 ml (fill to ring mark) | universal |
| `Маринелли 1л` (37 records) | `marinelli_1L` | 1000 ml | universal |
| `MARINELLI` (11) / `Marinelli` (1) | `marinelli_1L` (presumed, verify per .spe) | 1000 ml (presumed) | universal |
| `marinelli_0cm` | `marinelli_1L` (suffix = placement, not volume; +`placement_distance_cm: 0`) | 1000 ml | universal |
| `Дента-120мл` | `denta_120ml` | 120 ml | universal |
| `Дента-120` (typo) | `denta_120ml` | 120 ml | universal |
| `Дента-100` (3 records, Поверка-2016) | **PENDING source-reconciliation** (see below) | — | — |
| `Петри-60мл` | `petri_75ml` (vessel = 75 ml LSRM) | **60 ml (useful)** — **Gamma-1C ONLY** | Gamma-1C **only** |

### `Дента-100` — pending reconciliation

User lock 2026-06-05 = «a genuinely separate 100 ml geometry», but
the LSRM source p. 11 does not mention a 100 ml vessel. Three candidate
interpretations:
- (a) non-LSRM vessel 100 ml → `vessel_class: "other_denta_100ml"`, flag `lsrm_standard: false`;
- (b) LSRM Denta 120 ml with partial-fill → `vessel_class: "denta_120ml"` + `useful_sample_volume_ml: 100`;
- (c) operator typo. Decision deferred pending explicit user instruction.

### Petri 60 ml = useful, not vessel — **Gamma-1C ONLY** (user lock 2026-06-05)

> «This is the total volume; the useful one is 60 ml, and that's what we use for Gamma-1C.
>  For other detector types it needs to be checked individually.» (user)

⇒ For other crystal-classes (HPGe-coaxial-20pct, LaBr3, CdZnTe, Si(Li),
Si-surface-barrier) with the same `vessel_class: "petri_75ml"` —
`useful_sample_volume_ml` is unpacked from the passport / .spe header
of the specific record, **not inherited** from the Gamma-1C default of 60 ml.

### Architectural rule

The visual-template shape is determined by the tuple (`crystal_class` ×
`vessel_class` × `useful_sample_volume_ml`). The constructive part
(vessel + d_eff) is transferable between stations of the same crystal-class.
The operational part (useful_volume + placement_distance) is per-station.

---

## 4. Runtime conflict-resolution rule (user lock 2026-06-05)

> «During specific measurements, geometry parameters can be changed
>  by the operator. Therefore, when analyzing a spectrum in case of discrepancies,
>  you must ask for clarification.» (user)

**Geometry parameters at measurement time may diverge from any nominal
spec** — the operator may load a partial fill, use a non-LSRM
vessel substitute, move the sample to a non-standard distance, and
so on. Folder name / .spe header / passport / operator metadata —
**four independent provenance layers**, and in the general case they can
disagree.

**Precedence for logging** (not for silent decision):

| Rank | Layer | Authority |
|---|---|---|
| 1 | `operator_explicit_metadata` | Operator explicitly specified in session metadata |
| 2 | `passport_pdf_block` | Source passport / certificate |
| 3 | `spe_header_sample_geometry` / `spe_header_sample_volume` | .spe header field |
| 4 | `folder_name_hint` | Folder name (easily renamed, not authoritative) |
| 5 | `lsrm_source_default` | LSRM-source constructive default (only for vessel_class fallback, **never** for operational useful_volume) |

**Conflict handling**:

- **Offline RAG-build pipeline (W4)**: when a discrepancy is detected in
  vessel/volume/distance between layers — log to `__geometry_conflicts`,
  set `__needs_operator_review: true`, the template goes to
  `_pending_review/` (NOT into the canonical pool, **not**
  eligible for the similarity API).
- **Online analyzer (production runtime)**: when analyzing a specific
  spectrum, if a discrepancy between provenance layers is detected —
  the analyzer **MUST prompt the operator**: «Folder name suggests X,
  but .spe header reports Y. What are the actual measurement parameters?»
  Efficiency calculation is **gated** until an explicit operator decision.
  **Do NOT resolve by precedence silently** — the precedence rank is only
  for structured logging.

This applies to all three operational fields (`vessel_class`
substitution / `useful_sample_volume_ml` / `placement_distance_cm`)
and extends the project's CLAUDE.md anti-hallucination rule «every
statement references a specific offset/line/table in the
source» — extending it to the multi-layer conflict case.

---

## 5. Layout

```
detectors/Gamma-1C/
├── README.md                          # this file
├── certificates/                      # passports of standard sources (.xls / .pdf / .src)
├── data/                              # secondary peaks catalogs + aliases overrides
│   ├── averaged_backgrounds/          # 5 averaged .spe background files per geometry
│   ├── secondary_peaks.json           # Cs-137 + K-40 secondary peak catalog
│   └── secondary_peaks_v2.json        # 9-isotope rich catalog (incl. chain proxies)
├── efficiency/                        # .efr efficiency curves per geometry
│   └── Gamma-1C_NaI_63x63_USB_SN-01/
├── lsrm-libraries/                    # LSRM SpectraLine nuclide libraries
├── reference_spectra/                 # .spe reference spectra (verification campaigns)
│   └── Gamma-1C_NaI_63x63_USB_SN-01/
├── references/
│   ├── 05_intrinsic_detector_activity.md   # NaI(Tl) 63×63-specific intrinsic signatures
│   └── 07_dead_time_correction.md          # A, B coefficients for the УДС-ГЦ
└── raw_lsrm/                          # LOCAL-ONLY operator LSRM-tree, gitignored (F-115)
    ├── Work/
    │   ├── BG/Gamma-1C/Spe/           ← working trunk of the LSRM station
    │   │   ├── Маринелли/             ┐ both subfolders → one canonical
    │   │   ├── Маринелли 1л/          ┘ tag `marinelli_1L` (default = 1 L,
    │   │   │                            user lock 2026-06-05)
    │   │   ├── Точечная-25см/
    │   │   ├── Background/
    │   │   └── Spe — поверки/Поверка YYYY/
    │   ├── Calibration/
    │   └── …
    └── passports/                     ← (optional) .pdf/.txt passports of sources
```

---

## 6. Path resolver

Python code accesses these assets via the resolver module:

```python
from gamma.detectors.gamma1c import (
    DETECTOR_ROOT,
    CERTIFICATES_DIR,
    EFFICIENCY_DIR,
    REFERENCE_SPECTRA_DIR,
    LSRM_LIBRARIES_DIR,
    AVERAGED_BACKGROUNDS_DIR,
    SECONDARY_PEAKS_PATH,
    SECONDARY_PEAKS_V2_PATH,
    DEFAULT_REFERENCE_DIR,
    DEFAULT_EFFICIENCY_DIR,
    DETECTOR_NAME,
)
```

Never hardcode `detectors/Gamma-1C/...` paths in calling code — always go through the
resolver so future detectors can swap to their own subtrees without invasive churn.

---

## 7. Local-only `raw_lsrm/` working copy (NOT committed)

The operator's complete LSRM tree is placed in the `raw_lsrm/` subfolder, and it
**never** enters git (pattern F-150/F-293 «books_library»).

### Exclusion guarantees (defence-in-depth)

`raw_lsrm/` is excluded from artifacts at **three** levels:

1. **`.gitignore`** — pattern `detectors/Gamma-1C/raw_lsrm/` →
   git tracking will never pick up the operator's .spe files / passports / certificates.
2. **`scripts/build_release_archive.py:EXCLUDE_DIRS`** — the basename
   `raw_lsrm` is excluded → even if a file accidentally ends up in the working tree
   between gitignore passes, the release archive will skip it.
3. **F-115 anonymizer** (`scripts/gamma/reporting/anonymize.py`) — any
   output artifact (JSON / Markdown / HTML / PDF) that references
   a path inside `raw_lsrm/` will have the absolute path scrubbed down to the basename
   before being written to disk.

---

## 8. Isolation policy (v1.12.0)

This folder is **isolated** from any other detector. Algorithms in `scripts/gamma/`
are shared, but data, certificates, .efr curves, intrinsic-activity references and
secondary-peak catalogues here are valid **only** for the Gamma-1C complex.

When the AtomSpectra / AtomNano / RadiaCode pipelines are added (deferred), each
gets its own `detectors/<canonical>/` folder. Scripts will be copied across only
after they are stabilized in the Gamma-1C branch (per user policy 2026-05-29).

---

## 9. Crystal-class map (locked on 2026-06-05)

What is precisely known from project source comments — for future stations:

| Station-instance | Crystal-class (exact) | Source-pin |
|---|---|---|
| `Gamma-1C` (the only one in this project) | NaI 63×63 | `audit/_rag/visual_templates/SCHEMA.md:10`, `scripts/gamma/io/lsrm_spe.py:29` |
| `GP_HPGe20` | HPGe coaxial 20% | scope-glob name `Work\GP\HPGe(20%)` (`build_spectra_index.py:60`) |
| `NM_HPGe20` | HPGe coaxial 20% | scope-glob name `Work\NM\HPGe(20%)` (`build_spectra_index.py:66`) |
| `Handy_NaI` | NaI (size TBD) | scope-glob name `Work\Handy\Handy(NaI)` |
| `Handy_HPGe` | HPGe (size TBD) | scope-glob name `Work\Handy\Handy(HPGe)` |
| `Handy_LaBr` | LaBr3 (size TBD) | scope-glob name `Work\Handy\Handy(LaBr)` |
| `Simple_NaI` | NaI (Demo, size TBD) | scope-glob name `Work\Simple\NaI(Demo)` |
| `Simple_HPGe` | HPGe (Demo, size TBD) | scope-glob name `Work\Simple\HPGe(Demo)` |
| `Simple_TeCd` | CdZnTe / CZT (size TBD) | scope-glob name `Work\Simple\TeCd(Demo)` |
| `Simple_SiLi` | Si(Li) (size TBD) | scope-glob name `Work\Simple\SiLi(Demo)` |
| `Simple_Alpha` | Si surface barrier (alpha) | scope-glob name `Work\Simple\Alpha(Demo)` |

«size TBD» — to be clarified with the operator (or from the passport DETECTOR=
field of specific .spe files) BEFORE the W5+ harness; add a
`station → crystal_class` lookup-table to the build script.

**Note**: «Gamma-1S» is NOT a separate station. If the operator brings a new
LSRM export with the DETECTOR header «Гамма-1С №NNNN-NN» — this is **the same**
Gamma-1C in a different verification epoch. It is placed in `raw_lsrm/Work/BG/Gamma-1C/`
(see §7); cross-epoch / cross-verification drift study is handled via
the `_drift_study/` mirror (see SCHEMA.md «drift-study isolation»).

---

## 10. F-070 W4 Use-case (extended Gamma-1C ingest)

See `audit/_plans/F-070_W4_gamma1c_visual_templates_TODO.md` (renamed from
`_gamma1s_` after 2026-06-05 user lock).

Brief contract:
1. The operator copies the LSRM tree into `detectors/Gamma-1C/raw_lsrm/Work/...`
2. The build script `scripts/rag/build_visual_templates_nai63x63.py` (TODO)
   reads `raw_lsrm/...` via **relative paths**.
3. Each emitted VT-*.json:
   - provenance basename only (F-115 anonymizer)
   - `detector_id` → `УДС-ГЦ-63×63-USB` (S/N `№NNNN-NN` stripped)
   - `sample_id` cert-S/N patterns (`420-7-XX`) → `None`
   - `crystal_class: "NaI-63x63"` field + `station_observed_on: "Gamma-1C"`
     field (new in schema 0.2 — added in the S0 retrofit).
4. Result — JSON templates in
   `audit/_rag/visual_templates/<class>/VT-<NUC>-<GEOM>-<EPOCH>.json`
   (with the station-instance noted in the provenance).

---

## 11. Cross-refs

- F-78 / F-78a (aliases): `data/aliases.json`, `scripts/gamma/data/aliases.py`
- F-83 (detector isolation): this README §8
- F-115 (anonymization): `scripts/gamma/reporting/anonymize.py`
- F-150 / F-293 (external-data working-copy pattern): `books_library/` precedent
- F-155 (root-folder allow-list includes `detectors/`)
- BUG-40 (cyrillic-latin homoglyph warning): `KNOWN_AND_FIXED_ISSUES.md:1401-1422`,
  `tests/step04_detector_type/test_bug40_cyrillic_latin_warning.py`
- Crystal-vs-station distinction (LOCK 2026-06-05): this file §1
- W1 / W2 / W3 visual templates harness for NaI-63×63 / Gamma-1C station:
  `audit/_rag/visual_templates/SCHEMA.md`, `SIMILARITY_POLICY.md`
- SPECTRA_INDEX Gamma-1S metadata records (~394, historical artifact
  of v1.11.1 aliases): `audit/_rag/SPECTRA_INDEX.json` (filter `by_detector.Gamma-1S`)
  — this is the same physical Gamma-1C, transliterated by the LSRM-header parser
  via cyrillic-latin homoglyph mapping; see §1 «Canonical name fixation»

---

## 12. Last release that touched this folder

- `v1.12.0` — initial isolation of Gamma-1C-specific assets (F-83).
- `v1.25.0` (in progress) — F-070 W4 schema 0.2 retrofit (crystal-class
  abstraction + vessel taxonomy + runtime conflict-resolution rule).
