Write a single self-contained Python 3 script. Output ONLY the code, no prose,
no markdown fences. The script must run on Windows with `python script.py ...`.

# Purpose

Mechanically cross-check every NUMBER written in a Markdown document against the
text of one or more source PDF files. The script does not judge meaning. It sorts
numbers into three buckets and prints them.

# CLI

```
--doc  PATH   (required, exactly one) path to a Markdown file
--src  PATH   (required, repeatable via action="append") path to a PDF file
--json PATH   (optional) where to write a machine-readable report
```

Use `argparse`. Use `pathlib.Path` for all path arguments.

# Required behaviour

1. First line of executable code after imports:
   `sys.stdout.reconfigure(encoding="utf-8")`.
   Every file read/write must pass `encoding="utf-8"` explicitly.

2. Extract text from each PDF with PyMuPDF:
   `import fitz`, `doc = fitz.open(path)`, join `page.get_text()` over all pages
   with newlines. Concatenate the text of all `--src` files into one string.

3. Define a function that extracts all numbers from a string and returns them as
   a `set` of normalised strings:
   - first replace the Unicode minus `−` with ASCII `-`
   - regex: `\d+(?:[.,]\d+)?`
   - normalise each match by replacing `,` with `.`

4. Define a function `variants(num: str) -> set[str]` returning alternative
   spellings of the same number, because papers write numbers differently:
   - the number itself
   - if it contains `.`: the same with `,` instead of `.`; the value with
     trailing zeros and trailing dot stripped; and the integer part alone
   - if it is all digits and at least 2 chars long: the form with a decimal
     point/comma inserted after the first digit (so `322` also matches `3.22`
     and `3,22`, which is how `3.22 * 10^2` appears in text)
   - if it contains `.`: the same with the dot removed

5. A module-level tuple `DERIVED_MARKERS` of lowercase Russian/symbol substrings
   that mark a line as containing a DERIVED (computed, estimated, ours) number:
   "расчёт", "расчет", "вычисл", "оцен", "наш", "мы ", "получ", "размах",
   "отношение", "во сколько", "×", "раз", "%", "процент", "производн",
   "следств", "гипотез", "порядок", "приблиз", "≈", "~"

6. A module-level `set` named `NOISE` containing meaningless/noisy tokens to skip:
   the strings "0" through "14", "100", and the years
   "2009", "2011", "2012", "2014", "2015", "2022", "2023", "2026".

7. Iterate over the document lines with 1-based line numbers. Skip a line if it
   is blank after strip, or if the stripped line starts with any of
   "```", "|---", "---".
   For each line compute `is_derived_line` = True if any element of
   DERIVED_MARKERS occurs in the lowercased line.
   For every number found on the line, in sorted order:
   - skip it if it is in NOISE, or if the number with "." removed is shorter
     than 2 characters
   - build a record: dict with keys "line" (int), "num" (str), "text"
     (the stripped line truncated to 110 chars)
   - if `variants(num)` intersects the set of source numbers -> bucket `found`
   - elif `is_derived_line` -> bucket `derived`
   - else -> bucket `orphan`

8. Printing, in this exact order:
   - a line of 78 "=" characters
   - `ДОКУМЕНТ:` and the doc path
   - one `ИСТОЧНИК:` line per source, printing the file NAME only
   - a line of 78 "="
   - `чисел найдено в источнике дословно : {n}`
   - `чисел на строках с меткой производного: {n}`
   - `ЧИСЕЛ БЕЗ ОПОРЫ И БЕЗ МЕТКИ           : {n}`
   - blank line
   - if `orphan` is non-empty: a line of 78 "-", then the literal three lines
     `ТРЕБУЕТ ТОЛКОВАНИЯ — числа, которых нет в источнике и которые не`
     `помечены как производные. Каждое обязано получить письменный ответ:`
     `источник / вычислено из чего / ошибка.`
     then a line of 78 "-", then for each orphan record a line formatted as
     `  стр.{line:>4}  [{num}]  {text}`
   - else print `(чисел без опоры не найдено)`

9. If `--json` was given, write JSON with `ensure_ascii=False, indent=1` holding
   keys "found" (the integer count), "derived" (the list of records) and
   "orphan" (the list of records); then print `\nмашинный отчёт: {path}`.

10. `main()` returns 1 if the orphan bucket is non-empty, else 0. Guard with
    `if __name__ == "__main__": sys.exit(main())`.

# Style

Module docstring in Russian explaining that this is tier-0 mechanical auditing:
the tool splits numbers into three buckets and does not interpret them;
interpretation is the human's job. Keep comments in Russian, code identifiers in
English. No external dependencies beyond `fitz` (PyMuPDF) and the standard library.
