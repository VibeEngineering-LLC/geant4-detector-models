Write ONE Python 3 script `page_claims.py`. Output ONLY the code (no markdown fences, no prose). Comments in Russian, short. Start with `import sys, re` and `sys.stdout.reconfigure(encoding="utf-8")`.

## Purpose
List every numeric claim of a Russian markdown text so that a reviewer can check them one by one.

## Command line
`page_claims.py <text.md>` — read the file as UTF-8, print to stdout.

## Behaviour
* Split the text into lines; skip empty lines and the lines that contain only `@@NAME@@` placeholders (regex `^@@[A-Z]+@@$`).
* For each remaining line find numbers with the regex `\d+(?:[  ]\d{3})*(?:,\d+)?` (Russian decimal comma, space as thousands separator, e.g. `265 200`, `73,66`, `0,38`, `1460,5`).
* For every line that has at least one number print one row: `<line_number>|<numbers joined with ;>|<the first 160 characters of the line with markdown asterisks removed>`.
* At the end print `строк с числами: N, чисел всего: M`.
* Exit code 0; if the file is missing print an error message to stderr and exit with code 2.

## Self-test
`page_claims.py --selftest`: build an inline text with three lines: `Итог 73,66 ± 0,23 отсч./с`, an empty line, `@@CHART@@`, and `Пик 1460,5 кэВ, 265 200 отсч./ч`; run the same processing; require that exactly 2 lines are reported, the numbers of the first line are `73,66;0,23` and the numbers of the last line are `1460,5;265 200`. Print `SELFTEST PASS` and exit 0, else print `SELFTEST FAIL` and exit 1.
