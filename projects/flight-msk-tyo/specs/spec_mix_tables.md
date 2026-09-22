Write ONE complete Python 3 script `mix_tables.py`. Output ONLY the code (no markdown fences, no prose). Standard library only, no numpy.

## Purpose
Combine several PARMA source-table text files of the SAME particle into one weighted mixture file (used to build the route-average cosmic-ray source).

## Command line (exactly this)
    python mix_tables.py <out_file> <file1>:<w1> <file2>:<w2> ...
* At least one input. Each argument is `path:weight` (weight is a float; split on the LAST colon, because Windows paths contain `C:`). Weights are used as given (NOT renormalised); the caller passes weights that sum to 1.
* Any bad argument, unreadable file, or files that are inconsistent with each other -> print the reason to stderr and exit with code 2. Nothing is written in that case.

## Input file format (produced by `parma_tables.exe`, text, UTF-8)
Line 1: comment starting with `#` (e.g. `# ip=0 W=1.1e+01 Rc=... d=... g=...`).
Line 2: six fields: `nebin nabin ie511 flux511_cont flux511_line total_flux` (ints, ints, int, floats).
Line 3: `nebin+1` floats — the energy bin edges `ehigh[0..nebin]` in MeV.
Then `nebin` lines with `nabin` floats each: `D[k][ia]`, k = 1..nebin (line order), in /cm2/s/MeV per angular bin.

## Consistency checks (all inputs must agree, otherwise exit 2)
* `nebin`, `nabin`, `ie511` identical.
* the `ip=` token in line 1 identical.
* the energy edges of line 3 identical to relative tolerance 1e-8 element-wise.

## What to compute
Mixture (weights `w_j`, files `j`):
* `D_mix[k][ia] = sum_j w_j * D_j[k][ia]`
* `flux511_cont_mix = sum_j w_j * flux511_cont_j`, `flux511_line_mix = sum_j w_j * flux511_line_j`
* `total_flux_mix = sum_j w_j * total_flux_j`
* line 3 (energy edges) copied from the first file.

## Output file
Exactly the same format as the input (so that the same reader reads a single-point file and a mixture):
Line 1: `# mix ip=<ip> n_inputs=<n> weights=<w1,w2,...> sum_w=<sum of weights, %.9f>`
Line 2: `<nebin> <nabin> <ie511> <flux511_cont_mix> <flux511_line_mix> <total_flux_mix>` with floats as `%.9e`.
Line 3: energy edges as `%.9e` separated by single spaces.
Then `nebin` lines of `nabin` values `%.9e`, single-space separated.
Write with `encoding="utf-8", newline="\n"`.

## Also print to stdout (one line each)
* `total_flux_mix=<%.9e>`
* `continuum_check=<%.9e>` where `continuum_check = sum over k, ia of D_mix[k][ia]*(ehigh[k]-ehigh[k-1])` (k 1-based, edges as above).

## Requirements
* At the top `import sys` and `sys.stdout.reconfigure(encoding="utf-8")`, `sys.stderr.reconfigure(encoding="utf-8")`.
* Read files with `encoding="utf-8"`.
* Functions: `read_table(path)` returning a dict; `main()`; guard with `if __name__ == "__main__":`.
* Comments in Russian, short. No third-party imports.

## Syntax reminders (the previous attempt failed `ast.parse` on exactly this)
* The element-wise comparison of energy edges MUST be written with `any(...)`, for example
  `if any(abs(a - b) > 1e-8 * max(abs(a), abs(b), 1e-300) for a, b in zip(first["ehigh"], t["ehigh"])):`
  — a bare generator expression after `if` with an unmatched closing parenthesis is a syntax error.
* Output must be a complete, syntactically valid Python file.
