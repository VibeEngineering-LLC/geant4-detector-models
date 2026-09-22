Write ONE Python 3 module `bline_table.py` (standard library only). Output ONLY the code (no markdown fences, no prose). Comments in Russian, short. No main, no command line.

## Purpose
Parse the big line table of a Markdown report into a list of lines.

## Input format
A Markdown file (UTF-8). The table starts after the heading line that begins with `## 9. БОЛЬШАЯ СВОДНАЯ ТАБЛИЦА ЛИНИЙ` and ends at the next line that begins with `## ` . Inside it the table rows look like `| 27,36 | ¹²⁷I(n,γ) | 6,94/100 захв. | сам кристалл NaI/CsI | C |` (five cells: energy in keV with a DECIMAL COMMA, reaction/nuclide, yield, place, visibility class A/B/C/D) or `| 57,608 / 58,11 | ¹²⁷I(n,n′) (...) и ¹²⁷I(n,γ) 58,11 | 4,52 (захват) | ... | B |` where the energy cell may contain SEVERAL energies separated by ` / `, or `| 511 | ... |`, or ranges like `1461` / `2 614,5` with a space as thousands separator. Header and separator rows (`| E, кэВ |`, `|---|`) and any row whose first cell has no number must be skipped.

## Function (exactly)
`parse_lines(path) -> list[dict]`, each dict: `{"E_keV": float, "text": str, "reaction": str, "yield": str, "place": str, "cls": str}` — one dict per ENERGY in the first cell (a cell with two energies gives two dicts with the same reaction text). Rules for the energy cell: split on `/`; for every part remove spaces (also non-breaking spaces ` `) and replace the decimal comma by a dot, then take the FIRST number matched by `[0-9]+(?:\.[0-9]+)?`; ignore parts without a number; discard energies outside `1 <= E <= 20000`. `text` = `"%s | %s | %s" % (reaction, place, cls)`. Rows with fewer than 5 cells: fill the missing ones with empty strings. Return the list sorted by energy. Raise `ValueError("table not found")` if the heading is absent.
