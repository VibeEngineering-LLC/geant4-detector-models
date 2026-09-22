You are a C++ engineer. Write ONE complete C++17 header file `capture_emitter.hh` (header-only, NO Geant4 dependency, standard library only). Output ONLY the code (no markdown fences, no prose). Comments in Russian, short.

## Purpose
Sample neutron-capture gamma cascades from a prebuilt text database (`capture_db.txt`), energy conserving and correlated: a PRIMARY transition from the capture state down to a low level, then a cascade down the level scheme with branching and internal conversion.

## Database format (UTF-8 text)
Lines starting with `#` are comments (skip). For every isotope a block:
```
ISO <Z> <A_target> <Sn_keV> <sigma0_b> <Wp> <n_levels> <n_primary>
PRI <level_idx> <weight>            (n_primary lines; weights sum to Wp)
LEV <idx> <E_level_keV> <n_transitions>
TR <daughter_idx> <Eg_keV> <prob> <alpha>     (n_transitions lines; prob sums to 1 within the level)
END
```
`LEV` blocks appear for reachable levels only, in ascending `idx`; level 0 (ground state) may or may not be present; a level with `n_transitions == 0` and `idx > 0` is TERMINAL. Index levels by their `idx` through a `std::map<int, Level>` or an `unordered_map`.

## Interface (exactly)
```cpp
#pragma once
#include <string>
#include <vector>
#include <map>
struct Cascade {
  std::vector<double> gamma_keV;   // emitted gamma energies in emission order
  double conv_keV = 0.0;           // energy released as internal conversion (stays local)
  int primary_level = -1;          // level fed by the primary transition
};
class CaptureEmitter {
 public:
  bool Load(const std::string& path, std::string& err);   // false + err on any format problem
  bool Has(int Z, int A_target) const;
  double Wp(int Z, int A_target) const;                    // covered fraction, 0 if isotope unknown
  double Sn_keV(int Z, int A_target) const;                // 0 if unknown
  int Count() const;                                       // number of isotopes loaded
  // rand01: callable returning a uniform double in (0,1). Returns false if the capture must be left to the
  // standard model (isotope unknown, or the draw falls into the uncovered fraction 1-Wp). Otherwise fills `out`.
  template <class Rng> bool Sample(int Z, int A_target, Rng&& rand01, Cascade& out) const;
};
```
`Sample` must be a template defined inside the class.

## Sampling algorithm
1. Find the isotope; if absent return false. `u = rand01(); if (u >= Wp) return false;`
2. Choose the primary level: `r = rand01() * Wp`; accumulate the `PRI` weights in file order; the first level whose cumulative weight `>= r` is chosen (last one if none). `out.primary_level = idx`.
3. Primary gamma energy: `Eg0 = Sn - E_level`; correct for nuclear recoil `Eg0 = Eg0 - Eg0*Eg0/(2*M)` with `M = (A_target+1)*931494.0` keV (one iteration is enough); push it into `gamma_keV`. If `E_level <= 0` (ground state) the cascade ends here.
4. Cascade: current level `j = primary level`. Loop while `j != 0`: look up level `j`; if it has no transitions (terminal or missing) push the gamma `E_level_j` (keV, direct transition to the ground state) and stop; otherwise pick a transition by cumulative `prob` (`r = rand01()`, first transition with cumulative `>= r`, the last one if none); with probability `1/(1+alpha)` push its `Eg` into `gamma_keV`, otherwise add `Eg` to `conv_keV`; set `j = daughter`. Guard against infinite loops: stop after 200 steps.
5. Return true.

## Requirements
* `Load`: read with `std::ifstream`, parse with `std::istringstream`, check every read; reject with a meaningful `err` if a block is malformed, if `n_primary` lines are missing, if a `TR` line is outside a `LEV` block, if `Wp` is not in `[0, 1.000001]`, if any `prob` is negative or not finite. An empty database (no `ISO` block) is an error.
* Key the isotope map by `Z*1000 + A_target`.
* No global state; `Sample` must not allocate except for growing `out.gamma_keV`.
* Use `double` everywhere. `#include` everything you use (`<fstream>`, `<sstream>`, `<unordered_map>`, `<cmath>`, `<cstdlib>`...). Must compile with MSVC `cl /EHsc /std:c++17 /W4` without warnings and with g++ -std=c++17 -Wall.
