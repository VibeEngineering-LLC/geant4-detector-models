Write ONE complete Python 3 script `build_capture_db.py`. Output ONLY the code (no markdown fences, no prose). Allowed imports: `sys`, `math`, `os`, `re`, `pandas` (with `xlrd` for `.xls`). Comments in Russian, short. Keep every expression simple and on one line (earlier scripts of this project failed on unbalanced parentheses, unstripped labels, wrong merge keys, double subtraction of the mass number, and quantities computed once and reused for all rows).

## Purpose
Build a data-driven database of neutron-capture gamma CASCADES for a list of target isotopes: PRIMARY transitions (from the capture state at energy `Sn` down to a low level) taken from the measured IAEA PGAA line list, and everything BELOW the primary level taken from the ENSDF level scheme (Geant4 PhotonEvaporation data), with branching and internal conversion. The Geant4 application will sample one cascade per capture from this database (energy conserving, correlated); the fraction of captures whose primary transitions are not covered by measured lines is left to the standard Geant4 model.

## Command line
    python build_capture_db.py <data_dir> <levels_dir> <elements> <out_db> <out_report>
* `data_dir` holds `promptgammas.xls`, `isotope.xls`, `mass.mas20.txt`. `levels_dir` holds files named `z<Z>.a<A>` (Geant4 PhotonEvaporation). `elements` is a comma-separated list of element symbols (e.g. `H,C,Al,Fe`). Every STABLE isotope (natural abundance > 0 in `isotope.xls`) of each element is processed as a target.
* Wrong argument count -> usage to stderr, exit 2.

## Input formats
**mass.mas20.txt** (AME2020, fixed columns, read with `encoding="utf-8", errors="replace"`): a data line has `len(line) >= 43` and `line[0]` in `" 01"`; then `N=int(line[4:9])`, `Z=int(line[9:14])`, `A=int(line[14:19])`, mass excess (keV) `= float(line[28:42].strip().replace("#", "."))`. Skip lines that fail to parse. Store `mass_excess[(Z, A)]`. The neutron is `(0, 1)`. Neutron separation energy of the PRODUCT of capture on target `(Z, A)`: `Sn = mass_excess[(Z,A)] + mass_excess[(0,1)] - mass_excess[(Z,A+1)]` keV.

**isotope.xls** (`pd.read_excel(path, header=None, skiprows=2)`, columns by POSITION): `0` element symbol (forward fill with `.ffill()`), `1` target isotope label like `'  27-Al '` (STRIP whitespace; format `<A>-<Symbol>`), `5` natural abundance percent, `8` thermal capture cross section `sigma0` in barn. Convert with `pd.to_numeric(errors="coerce")`.

**promptgammas.xls** (`header=None, skiprows=2`, columns by position): `1` PRODUCT isotope label like `'  28-Al'` (STRIP; mass number of the product = target A + 1), `5` gamma energy keV, `7` partial cross section `sigma_gamma` (barn). `to_numeric(errors="coerce")`, drop NaN. Yield per capture of a line: `y = sigma_gamma / sigma0(target)`.

**Level files `z<Z>.a<A>`** (product nucleus `(Z, A_product)`): text lines with whitespace-separated columns. A LEVEL line has a NON-numeric second token (a string such as `-`, `+X`): `idx  flag  E_keV  halflife_s  Jpi  n_gamma`. It is followed by `n_gamma` TRANSITION lines whose second token IS numeric: `daughter_idx  Eg_keV  Irel  multipolarity  mixing  alpha  [shell columns...]`. Level 0 is the ground state (E = 0). Total transition weight is `Irel*(1+alpha)`; gamma emission probability of a transition is `1/(1+alpha)`, internal conversion `alpha/(1+alpha)`.

## Algorithm for each target isotope T = (Z, A), product P = (Z, A+1)
1. Skip (record reason in the report) if: no `mass_excess` for T or P; no level file for P; `sigma0` missing or <= 0; no PGAA rows for the product label.
2. `Sn` as above. Recoil of the product after emitting a gamma of energy `E`: `E*E/(2*M)` with `M = (A+1)*931494.0` keV.
3. Parse the level scheme of P into arrays `levels[idx] = (E_level, [transitions])` where each transition is `(daughter, Eg, weight, alpha)` with `weight = Irel*(1+alpha)`. Normalise weights per level to sum to 1 (`prob`). A level with no transitions and idx > 0 is TERMINAL.
4. PGAA lines of P: merge lines closer than 0.3 keV (sort by energy; keep strength-weighted energy, sum `y`).
5. Classify each PGAA line `i` (energy `E`, yield `y`) as PRIMARY or NOT:
   * `Lt = Sn - E - E*E/(2*M)` (energy of the level the primary transition would feed);
   * `tol = 1.5 + 3.0e-4*E` keV; find the level `j` with the smallest `|E_level_j - Lt|`; it is a PRIMARY candidate only if that distance `<= tol`;
   * secondary test: the line is a SECONDARY candidate if some level scheme transition has `|Eg - E| <= tol`;
   * decision: if `E >= 0.5*Sn` and primary candidate -> PRIMARY; if `E < 0.5*Sn` and primary candidate and NOT secondary candidate -> PRIMARY; otherwise NOT primary.
6. Primary weight of level `j` = sum of `y` of the PRIMARY lines feeding it. `Wp = sum of all primary weights`. If `Wp > 1.0`: scale all primary weights by `1/Wp` and set `Wp = 1.0`, and write a warning in the report. The fallback fraction (capture left to Geant4's own model) is `1 - Wp`.
7. Reachable levels: breadth-first from the primary levels following transition daughters; only these are written.
8. Self-check (independent of how the DB was built): propagate the primary weights down the scheme analytically (probability of visiting a level = sum over feeders; probability of emitting the gamma of a transition = visit probability * `prob` * `1/(1+alpha)`); this gives the modelled yield per capture for every gamma energy (merge energies within 0.3 keV, ignore fallback fraction). Compare with the PGAA lines that were NOT used as primary and have `y >= 0.005`: count `n_check`, and `n_within2` = those whose modelled yield (at energy within tol) is within a factor 2 of `y` (`0.5 <= model/y <= 2.0`). Print to stdout one line per isotope: `ISO <Z>-<A> Sn=<%.3f> sigma0=<%.4g> Wp=<%.4f> primaries=<n> levels=<n> check=<n_check> within2x=<n_within2>`.

## Output database `out_db` (UTF-8 text, `\n` newlines, floats `%.6g` unless stated)
First line: `# capture_db v1 targets=<count written> elements=<elements>`. Then for every written isotope:
```
ISO <Z> <A_target> <Sn %.4f> <sigma0> <Wp %.6f> <n_levels> <n_primary>
PRI <level_idx> <weight>            (one line per primary level)
LEV <idx> <E_level %.4f> <n_transitions>
TR <daughter_idx> <Eg %.4f> <prob %.6f> <alpha %.6g>     (n_transitions lines, prob normalised to sum 1)
END
```
`LEV` blocks in ascending `idx`, only reachable levels; a terminal level has `n_transitions = 0`.

## Output report `out_report` (Markdown, UTF-8)
Title, then a table with one row per target isotope: `isotope | Sn | sigma0 | PGAA lines | primaries | Wp | fallback 1-Wp | check | within2x`. Below it, for every isotope, a list of PGAA lines with `y >= 0.02` that are NOT primary and NOT reproduced by the model (modelled yield below half of `y`), with energy and yield (this is the list of lines that remain unexplained). Then the list of skipped isotopes with reasons. Also print on the last stdout line `DONE isotopes=<written> skipped=<skipped>`.

## Pitfalls of the first attempt (ALL must be avoided)
1. NEVER write a list of element symbols inline (the first attempt built a 90-element list inside an f-string, used the wrong mass number, and produced unbalanced parentheses). The symbol and the atomic number `Z` of a target come from `isotope.xls`: column `0` (symbol, forward filled), column `2` = `Z`, column `3` = mass number `A` (numeric). Build a dict `targets[symbol] = [(Z, A, abundance, sigma0), ...]` once and iterate over the requested symbols; skip isotopes with abundance NaN or <= 0.
2. The PGAA label of the PRODUCT is `f"{A+1}-{symbol}"` (mass number of the product). Look up `prompt_lines[label]` (a dict built once from the stripped labels of `promptgammas.xls`) with exactly that label.
3. Put label construction in one small function `product_label(symbol, A)` returning `f"{A+1}-{symbol}"` and use it everywhere.
4. Compute `Sn`, `M`, the level scheme, the merged lines and all per-isotope arrays INSIDE the per-isotope loop only.
5. No expression longer than about 100 characters; split anything long into named intermediate variables.

## Requirements
* Functions: `read_masses`, `read_isotopes`, `read_prompt`, `read_levels`, `classify_primary`, `propagate`, `main`; guard with `if __name__ == "__main__":`.
* At the top `sys.stdout.reconfigure(encoding="utf-8")`, `sys.stderr.reconfigure(encoding="utf-8")`; read text files with `encoding="utf-8", errors="replace"`.
* Do NOT reuse a variable computed for one isotope in the next one; build every structure inside the per-isotope loop. Isotope labels always `.strip()`-ed before comparison. Do not subtract 1 from a mass number that is already a target mass number.
* Deterministic output (sort isotopes by `(Z, A)`).
