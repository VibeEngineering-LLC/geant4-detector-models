Write ONE complete Python 3 module `capture_levels.py`. Output ONLY the code (no markdown fences, no prose). Standard library only. Comments in Russian, short.

## Purpose
Parse one nuclear level-scheme file of the Geant4 `PhotonEvaporation` database (a file named like `z13.a28`, the nucleus Al-28) into a dictionary.

## File format
Each line is whitespace-separated. There are two kinds of lines:
* a LEVEL line: its SECOND token is NOT a number (it is a marker string such as `-`, `+X`, `+Y`). Columns: `idx  marker  E_keV  halflife_s  Jpi  n_gamma`. `idx` is an integer starting at 0 (the ground state, `E_keV = 0`). `n_gamma` is the number of transition lines that follow.
* a TRANSITION line: its SECOND token IS a number (float). Columns: `daughter_idx  Eg_keV  Irel  multipolarity  mixing  alpha  [more columns]`. It belongs to the most recent LEVEL line. `daughter_idx` is an integer, `Eg_keV` the gamma energy, `Irel` the relative intensity, `alpha` the total internal conversion coefficient (sixth column, index 5).

Decide the kind of a line ONLY by trying `float(tokens[1])`: if it raises `ValueError` the line is a LEVEL line, otherwise a TRANSITION line. Skip blank lines. Any other malformed line (too few tokens, non-numeric where a number is required) -> raise `ValueError` with the line number in the message. A TRANSITION line before the first LEVEL line -> `ValueError`.

## Interface (exactly)
```python
def read_levels(path):
    """Возвращает dict: idx -> {"E": float, "trans": [(daughter_idx, Eg, weight, alpha), ...]}."""
def normalise(levels):
    """Возвращает dict: idx -> {"E": float, "trans": [(daughter_idx, Eg, prob, alpha), ...]}, где prob = weight/сумма weight по уровню (0, если сумма 0)."""
```
* `weight = Irel * (1.0 + alpha)`; `alpha` is `0.0` if the transition line has fewer than 6 columns.
* The number of transitions actually read for a level may differ from `n_gamma`; do not treat that as an error.
* `normalise` must not modify its input.

## Command-line self-test
`python capture_levels.py <file>` prints exactly one line `levels=<n_levels> transitions=<n_transitions> ground_E=<E of idx 0> max_E=<largest level energy>` (floats with `%.4f`) and exits 0; on `ValueError` prints the message to stderr and exits 2; wrong argument count -> usage to stderr, exit 2.

## Requirements
* `sys.stdout.reconfigure(encoding="utf-8")` and `sys.stderr.reconfigure(encoding="utf-8")` at the top; read the file with `encoding="utf-8", errors="replace"`.
* Guard the command-line part with `if __name__ == "__main__":`.
* Keep every expression simple and on one line.
