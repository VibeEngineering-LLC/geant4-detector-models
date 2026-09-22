You are a C++17 engineer (MSVC, Windows). Write ONE complete source file `replay_selftest.cpp` (contains `main`, NO Geant4). Output ONLY the code (no markdown fences, no prose). Comments in Russian, short.

## Purpose
Known-answer self-test of `replay.hh` (`#include "replay.hh"`, which includes `psp_record.hh`). API of `replay.hh`:
* `class PspEvents { bool AddFile(const std::string&); bool NextEvent(std::vector<PspRec>&); long long EventsRead() const; };`
* `template <class Rng> bool ChooseGroup(const std::vector<PspRec>& ev, double ell_cm, double L_cm, Rng&& rand01, std::vector<int>& idx, double& w, double& zd_cm);` — returns false when the drawn position is outside [-L/2, L/2] (then w = 0).
* `PspRec` (packed 52-byte POD, fields `evt`, `y_cm`, ...), `PspWriter{Open, Write, Close, Count}`, `PspReader` from `psp_record.hh`.

## Test 1: unbiasedness of the weights (known analytic answer)
For an event with records at axial positions y_1..y_n the expected value of the returned weight (0 when ChooseGroup returns false) over many draws equals `measure(union_j [y_j - ell/2, y_j + ell/2] ∩ [-L/2, L/2]) / L`. Use `ell = 30`, `L = 300` and these four events (only the `y_cm` field matters, set the others to 0): (a) y = {0} -> expected 30/300 = 0.1; (b) y = {0, 10, 100} -> union [-15,25] + [85,115] = 70 -> 0.233333...; (c) y = {140, 145} -> union clipped to [125, 150] = 25 -> 0.083333...; (d) y = {-149, -100, 0, 149} -> 16 + 30 + 30 + 16 = 92 -> 0.306667. Compute the expected values IN THE CODE by a small function `unionMeasure(ys, ell, L)` (sort the intervals, clip to [-L/2, L/2], merge, sum) — do not hard-code the numbers; ALSO print the hard-coded hand values above next to them for the reader. For each event draw `N = 4000000` samples with `std::mt19937_64 gen(12345 + eventIndex)` and `rand01 = [&]{ return std::uniform_real_distribution<double>(0.0, 1.0)(gen); }` (create the distribution object once outside the lambda), accumulate the sum and the sum of squares of w; PASS when `|mean - expected| <= 5 * sqrt(var / N)` AND additionally `|mean - expected| <= 0.01 * expected` (both must hold). Print one line per event: `T1 <label> expected=<..> mean=<..> se=<..> PASS|FAIL`.

## Test 2: event grouping across files
Write three files in the current directory with `PspWriter`: `rs_A.psp` with records whose `evt` values are 5,5,6,7,7,7 (in this order); `rs_B.psp` EMPTY (open and close the writer without writing); `rs_C.psp` with evt values 7,7,8. Then `PspEvents ev; ev.AddFile("rs_A.psp"); ev.AddFile("rs_missing.psp"); ev.AddFile("rs_B.psp"); ev.AddFile("rs_C.psp");` (the missing file must be skipped, `AddFile` still returns true). Read all events with `NextEvent` and require the sequence of event sizes to be exactly {2, 1, 3, 2, 1} (event 7 at the end of file A must NOT merge with event 7 at the start of file C) and `EventsRead() == 5`. Also check that the `evt` values in the first event are both 5. Print `T2 sizes=<list> expected=2,1,3,2,1 PASS|FAIL`. Delete the three temporary files with `std::remove` at the end.

## Exit code
Print a final line `SELFTEST PASS` and return 0 only if every check passed; otherwise print `SELFTEST FAIL` and return 1. Any exception -> message to stderr and return 2.

## Requirements
* Include every header you use (<random>, <vector>, <cmath>, <cstdio>, <algorithm>, <string>, <iostream>). Compile with MSVC `cl /EHsc /std:c++17 /utf-8`. Fill a `PspRec` with `PspRec r{};` then set fields.
