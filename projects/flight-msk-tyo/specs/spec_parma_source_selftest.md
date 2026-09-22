You are a C++ engineer. Write ONE complete C++17 source file `parma_source_selftest.cpp` and output ONLY the code (no markdown fences, no prose). It tests the header `parma_source.hh` (class `ParmaSource`, given below) by statistics. Include it with `#include "parma_source.hh"`.

## Command line
    parma_source_selftest <table_file> <N_events> <radius_cm> <seed>
Wrong argument count -> usage to stderr, exit 2. If `Load` fails -> print `err` to stderr, exit 2.

## Randomness
`std::mt19937_64 gen(seed); std::uniform_real_distribution<double> u01(0.0, 1.0);` and pass `[&]{ return u01(gen); }` as `rand01`.

## What to do
Draw `N` events with `src.Sample(radius_cm, rand01)`. Use `long long` counters. Every check prints one line and the program exits with code 0 only if all checks PASS, otherwise 1.

**Check A (geometry).** For each event with position `p=(x_cm,y_cm,z_cm)` and direction `d=(u,v,w)`:
* `| |d| - 1 | <= 1e-9`;
* `| p·d + radius_cm | <= 1e-6*radius_cm` (the start point lies at distance `radius_cm` behind the origin along the flight direction);
* `|p × d| <= radius_cm*(1+1e-9)` (the line passes within `radius_cm` of the origin).
Count events violating ANY of the three. PASS if the count is 0. Print `CHECK A geometry violations=<n> PASS|FAIL`.

**Check B (energy).** Let `nebin = src.NEBin()`, edges `src.EdgesMeV()`. Find the fine energy bin of each event (the smallest `k` in 1..nebin with `e <= edges[k]`, use `std::lower_bound` on the edges; events with `e` exactly equal to `edges[k]` belong to bin `k`). Group the fine bins into `nebin/20` consecutive groups of 20 fine bins (nebin is 1000). Expected count of a group = `N * sum of src.EnergyBinProb(k)` over its fine bins. Only groups with expected count >= 25 take part. For those compute `z = (obs-exp)/sqrt(exp)`. PASS if `max|z| < 4.5` and `chi2/ndf < 1.6` (chi2 = sum z², ndf = number of participating groups). Print `CHECK B energy groups=<ndf> max_abs_z=<..> chi2_ndf=<..> PASS|FAIL`.

**Check C (cosine).** `nabin = src.NABin()`, cosine bin edges `-1 + 2*ia/nabin`. Bin index of an event from `cos_down` (clamp to 1..nabin, `ia = min(nabin, max(1, (int)std::floor((cos_down+1)/2*nabin)+1))`). Expected = `N * src.CosBinProbMarginal(ia)`. Same participation rule (expected >= 25), same PASS criteria as B. Print `CHECK C cosine bins=<ndf> max_abs_z=<..> chi2_ndf=<..> PASS|FAIL`.

**Check D (511-keV line).** Only if `src.IE511() > 0`: count events with `e_MeV == 0.51099895` exactly (compare with `==` on the same double literal). Expected = `N * src.EnergyBinProb(src.IE511()) * src.LineFrac()`. `z = (obs-exp)/sqrt(exp)`; PASS if `|z| < 4.5` and (if `exp < 25`) print `SKIP` instead. Print `CHECK D line511 obs=<..> exp=<..> z=<..> PASS|FAIL|SKIP`. If `IE511() == 0` print `CHECK D line511 n/a`.

**Check E (means, information only, no PASS/FAIL).** Print `INFO mean_cos_down=<..> mean_energy_MeV=<..>`.

Finally print `SELFTEST PASS` or `SELFTEST FAIL`.

## Pitfalls of the first attempt (all four MUST be avoided)
1. The EXPECTED counts must NOT depend on the drawn events. Compute them BEFORE looking at the events: for energy group `g` (fine bins `k = 20*g+1 .. 20*g+20`, 1-based) `exp_g = N * sum_k src.EnergyBinProb(k)`; for cosine bin `ia` (1-based) `exp_ia = N * src.CosBinProbMarginal(ia)`. Only the OBSERVED counts are accumulated inside the event loop. Never add `expected[...]` inside the event loop.
2. Do NOT declare your own `Event` struct. Store events as `ParmaSource::Event` (or do not store them at all: fill the histograms and the geometry counter inside one loop over `i < N`, which also saves memory).
3. When `Load` fails print the `err` STRING to stderr (`std::cerr << err`), not the word "err".
4. Energy bin index: with `edges` (size nebin+1) and 1-based bin `k` covering `(edges[k-1], edges[k]]`, `k = std::lower_bound(edges.begin()+1, edges.end(), e) - edges.begin()`; if `e <= edges[0]` use `k = 1`; if `k > nebin` use `k = nebin`. The 0-based index into arrays is `k-1`. An event with `e` exactly equal to `edges[k]` therefore belongs to bin `k`.

## Requirements
* Standard library only. Comments in Russian, short. Must compile with MSVC `cl /EHsc /std:c++17 /W4 /utf-8` and g++ -std=c++17 -Wall.
* Use `double` for all statistics; avoid narrowing conversions (use `static_cast`).
* The class interface you can rely on (all `const` except `Load`):
  `bool Load(const std::string& path, std::string& err); double TotalFlux() const; int NEBin() const; int NABin() const; int IE511() const; double LineFrac() const; const std::vector<double>& EdgesMeV() const; double EnergyBinProb(int k) const /*1-based*/; double CosBinProbMarginal(int ia) const /*1-based*/; template<class Rng> Event Sample(double radius_cm, Rng&& rand01) const;`
  `struct Event { double e_MeV; double u, v, w; double x_cm, y_cm, z_cm; double cos_down; };`
