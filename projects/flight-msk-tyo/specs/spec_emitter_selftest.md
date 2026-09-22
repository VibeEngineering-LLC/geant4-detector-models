You are a C++ engineer. Write ONE complete C++17 source file `emitter_selftest.cpp` and output ONLY the code (no markdown fences, no prose). Comments in Russian, short. It tests the header `capture_emitter.hh` (class `CaptureEmitter`, struct `Cascade`; interface below) with `#include "capture_emitter.hh"`.

## Command line
    emitter_selftest <db_file> <Z> <A_target> <N> <seed> [dump_prefix]
Wrong argument count -> usage to stderr, exit 2. If `Load` fails print the `err` string to stderr and exit 2.

## Randomness
`std::mt19937_64 gen(seed); std::uniform_real_distribution<double> u01(0.0, 1.0);` and pass `[&]{ return u01(gen); }` as `rand01`.

## What to do
Draw `N` captures with `emitter.Sample(Z, A, rand01, c)`. Count `n_ok` (returned true) and `n_fallback` (returned false). For the `n_ok` cascades accumulate:
* a histogram of gamma energies with 1-keV bins over 0..12000 keV (bin = `floor(E)`), counts of ALL emitted gammas in all cascades, and `n_gamma_total`;
* `sum_E = sum over gammas of E`, and the per-cascade energy balance `E_total = sum(gamma_keV) + conv_keV + recoil` where `recoil = sum over gammas of E*E/(2*M)` with `M = (A_target+1)*931494.0`; track `max_abs_balance_err = max over cascades of |E_total - Sn|` where `Sn = emitter.Sn_keV(Z, A)`. (The primary gamma already includes recoil correction in the emitter; the emitter ignores the recoil of the cascade gammas, so allow a tolerance in the report but print the maximum.)
* the number of cascades whose primary level is each level index (map `int -> long long`).

## Output (stdout, exactly these lines)
`N=<N> n_ok=<n_ok> n_fallback=<n_fallback> Wp_expected=<emitter.Wp(Z,A)> Wp_observed=<n_ok/N>`
`gammas_per_capture=<n_gamma_total/(double)n_ok>`
`max_abs_balance_err_keV=<..>`
Then, sorted by descending count, one line per gamma-energy bin with more than 0.1% of `n_ok` counts: `LINE <bin_lo_keV> per_capture=<count/n_ok>` (per covered capture).
Then one line per primary level: `PRIMARY level=<idx> fraction=<count/n_ok>`.
If a third-to-sixth argument `dump_prefix` was given, additionally write `<dump_prefix>_hist.csv` with `energy_keV,count` for every non-empty bin.
Use `std::setprecision(8)`. Exit code 0.

## Interface you can rely on
```cpp
struct Cascade { std::vector<double> gamma_keV; double conv_keV; int primary_level; };
class CaptureEmitter { public:
  bool Load(const std::string& path, std::string& err);
  bool Has(int Z, int A_target) const; double Wp(int Z, int A_target) const; double Sn_keV(int Z, int A_target) const; int Count() const;
  template <class Rng> bool Sample(int Z, int A_target, Rng&& rand01, Cascade& out) const; };
```
Must compile with MSVC `cl /EHsc /std:c++17 /W4 /utf-8` without warnings and with g++ -std=c++17 -Wall. Standard library only. Use `static_cast` for narrowing conversions.
