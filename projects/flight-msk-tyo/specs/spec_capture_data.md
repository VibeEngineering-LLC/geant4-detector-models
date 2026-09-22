Write ONE complete Python 3 module `capture_data.py`. Output ONLY the code (no markdown fences, no prose). Allowed imports: `sys`, `math`, `pandas` (with `xlrd` for `.xls`). Comments in Russian, short. Keep every expression simple and on one line.

## Purpose
Read the nuclear data needed to build neutron-capture gamma cascades: atomic mass excesses (AME2020), the isotope table (natural abundance, thermal capture cross section) and the IAEA PGAA prompt-gamma line list.

## Interface (exactly)
```python
def read_masses(path):
    """dict (Z, A) -> избыток массы в кэВ."""
def sn_keV(masses, Z, A_target):
    """Энергия связи нейтрона в продукте захвата на ядре (Z, A_target), кэВ, или None, если масс нет."""
def read_targets(path):
    """dict symbol -> list of (Z, A, abundance_percent, sigma0_barn), только изотопы с abundance > 0 и sigma0 > 0."""
def read_prompt(path):
    """dict product_label -> list of (E_keV, sigma_gamma_barn); product_label вида '28-Al' (без пробелов)."""
def merge_lines(lines, tol=0.3):
    """lines: list of (E, w). Сортирует по E и склеивает соседей ближе tol кэВ: энергия — средняя, взвешенная по w; w суммируются. Возвращает list of (E, w)."""
def product_label(symbol, A_target):
    """Метка ПРОДУКТА захвата: f"{A_target+1}-{symbol}"."""
```

## File formats
**mass.mas20.txt** (AME2020, fixed columns, read with `encoding="utf-8", errors="replace"`): a data line has `len(line) >= 43` and `line[0]` in `" 01"`; then `N=int(line[4:9])`, `Z=int(line[9:14])`, `A=int(line[14:19])`, mass excess `float(line[28:42].strip().replace("#", "."))` in keV. Lines that fail to parse are skipped. The neutron is `(0, 1)`.
`sn_keV(masses, Z, A_target) = masses[(Z,A_target)] + masses[(0,1)] - masses[(Z,A_target+1)]`; return `None` if any of the three is missing.

**isotope.xls**: `pd.read_excel(path, header=None, skiprows=2)`; columns by POSITION: `0` element symbol (forward fill with `.ffill()`), `2` = Z, `3` = mass number A, `5` = natural abundance in percent, `8` = thermal capture cross section sigma0 in barn. Convert numeric columns with `pd.to_numeric(errors="coerce")`. Drop rows where Z, A, abundance or sigma0 is NaN, and rows with abundance <= 0 or sigma0 <= 0. Symbols must be `str(...).strip()`. Return integers for Z and A.

**promptgammas.xls**: `pd.read_excel(path, header=None, skiprows=2)`; columns by POSITION: `1` product label like `'  28-Al'` (STRIP whitespace), `5` gamma energy in keV, `7` partial gamma cross section sigma_gamma in barn (both `to_numeric(errors="coerce")`, drop NaN rows). Keep the list of `(E, sigma_gamma)` per stripped label.

## Command-line self-test
`python capture_data.py <data_dir>` (directory holding the three files) prints, each on its own line:
* `masses=<count>`
* `Sn Al-27 = <sn_keV(13,27) %.3f>`, `Sn Fe-56 = <sn_keV(26,56) %.3f>`, `Sn H-1 = <sn_keV(1,1) %.3f>`
* `targets_H=<number of H isotopes with abundance>` and for `Al`: `Al27 sigma0=<sigma0>` (the isotope with `A == 27`)
* `prompt_labels=<count>` and `H-2 lines=<count of lines of label '2-H'>`, `first H-2 line E=<%.3f> sigma_gamma=<%.4g>`.
Wrong argument count -> usage to stderr, exit 2.

## Requirements
* `sys.stdout.reconfigure(encoding="utf-8")`, `sys.stderr.reconfigure(encoding="utf-8")` at the top.
* Guard the command-line part with `if __name__ == "__main__":`.
* Do NOT hard-code lists of element symbols anywhere.
