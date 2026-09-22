You are a C++ engineer. Write ONE complete C++17 header file `parma_source.hh` (header-only, NO Geant4 dependency, standard library only). Output ONLY the code (no markdown fences, no prose).

## Purpose
Sample cosmic-ray primary particles (energy, direction, start position) from a tabulated PARMA/EXPACS distribution, reproducing EXACTLY the sampling algorithm of the official PARMA `main-generator.cpp`, but reading the table from a text file.

## Table file format (text, UTF-8)
Line 1: comment starting with `#`.
Line 2: `nebin nabin ie511 flux511_cont flux511_line total_flux` (int int int double double double).
Line 3: `nebin+1` doubles: energy bin edges `ehigh[0..nebin]` in MeV.
Then `nebin` lines, each with `nabin` doubles: `D[k][ia]` for k = 1..nebin (file order), ia = 1..nabin, units /cm2/s/MeV per angular bin.
The angular mesh is fixed: `ahigh[ia] = -1 + (2/nabin)*ia`, ia = 0..nabin (cosine of zenith angle; +1 means the particle comes from above, i.e. moves downward).

## Interface (exactly)
```cpp
#pragma once
#include <string>
#include <vector>
class ParmaSource {
 public:
  struct Event { double e_MeV; double u, v, w; double x_cm, y_cm, z_cm; double cos_down; };
  // returns true on success; on failure returns false and sets `err`
  bool Load(const std::string& path, std::string& err);
  double TotalFlux() const;          // total_flux from line 2, /cm2/s (includes the 511-keV line)
  int NEBin() const; int NABin() const; int IE511() const;
  const std::vector<double>& EdgesMeV() const;        // ehigh[0..nebin]
  // probability mass of energy bin k (1-based k=1..nebin), normalised so that the sum over k is 1,
  // INCLUDING the line in bin ie511
  double EnergyBinProb(int k) const;
  // probability of cosine bin ia (1-based) integrated over all energies (weighted by EnergyBinProb), sums to 1
  double CosBinProbMarginal(int ia) const;
  // radius_cm: radius of the target sphere; rand01: callable returning a uniform double in (0,1)
  template <class Rng> Event Sample(double radius_cm, Rng&& rand01) const;
};
```
`Sample` must be a template defined inside the class (in the header).

## Algorithm (identical to main-generator.cpp)
Precompute at Load time:
* `mass_k = (sum over ia of D[k][ia]) * (ehigh[k]-ehigh[k-1])`; for k == ie511 (only if ie511 > 0) add `flux511_line` to `mass_k`. `etab[k]` = cumulative sum of `mass_k` normalised to 1 at k = nebin. `EnergyBinProb(k) = mass_k / sum`.
* For every k: `atab[k][ia]` = cumulative sum over ia of `D[k][ia]`, normalised to 1 at ia = nabin.
* `line_frac = flux511_line / mass_ie511` (only when ie511 > 0).
`Sample`:
1. Draw `r = rand01()`; find energy bin k = the smallest index in 1..nebin with `r <= etab[k]` (binary search; if none, k = nebin). Then `r2 = rand01()`; `e = ehigh[k-1]*r2 + ehigh[k]*(1-r2)`. If `k == ie511 && ie511 > 0`: draw `r3 = rand01()`; if `r3 < line_frac` then `e = 0.51099895` (electron mass, MeV).
2. Draw `r = rand01()`, `r2 = rand01()`; find angle bin ia = the smallest index in 1..nabin with `r <= atab[k][ia]`; `cx = ahigh[ia-1]*r2 + ahigh[ia]*(1-r2)` where ahigh as above. (`cx` = cos_down.)
3. `phi = 2*pi*(rand01()-0.5)`.
4. Position on the disc: repeat { `xd = (rand01()-0.5)*2*radius_cm; yd = (rand01()-0.5)*2*radius_cm;` } while `sqrt(xd*xd+yd*yd) > radius_cm`; `zd = radius_cm`.
5. `sx = sqrt(1-cx*cx)`;
   `x = xd*cx*cos(phi) - yd*sin(phi) + zd*sx*cos(phi)`;
   `y = xd*cx*sin(phi) + yd*cos(phi) + zd*sx*sin(phi)`;
   `z = -xd*sx + zd*cx`;
   `u = -sx*cos(phi); v = -sx*sin(phi); w = -cx`.
6. Fill and return the `Event` (`cos_down = cx`).

## Requirements
* Read the file with `std::ifstream`, check every read; verify that the number of values read is exactly as declared; reject (return false, meaningful `err` text) if `nebin`, `nabin` are not positive, if any `D` value is negative or not finite, if `total_flux` is not positive and finite, or if the sum of masses is zero.
* Use `double` everywhere. Use `M_PI` only after `#define _USE_MATH_DEFINES` before the includes, or define your own `constexpr double kPi = 3.14159265358979323846;` (prefer the latter).
* `#include` everything you use (`<fstream>`, `<sstream>`, `<cmath>`, `<algorithm>`, `<stdexcept>` if needed).
* `Sample` must not allocate.
* Comments in Russian, short. Must compile with MSVC `cl /EHsc /std:c++17 /W4` without warnings and with g++ -std=c++17 -Wall.
