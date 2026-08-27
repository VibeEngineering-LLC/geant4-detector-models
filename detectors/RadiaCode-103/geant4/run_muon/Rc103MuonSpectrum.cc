#include "Rc103MuonSpectrum.hh"

#include "Randomize.hh"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

double Rc103MuonSpectrum::Gaisser(double eGeV) {
  const double t1 = 1.0 / (1.0 + 1.1 * eGeV / 115.0);
  const double t2 = 0.054 / (1.0 + 1.1 * eGeV / 850.0);
  return 0.14 * std::pow(eGeV, -2.7) * (t1 + t2);
}

Rc103MuonSpectrum::Rc103MuonSpectrum() {
  fE.resize(kNNodes);
  fCdf.resize(kNNodes);

  const double lnLo = std::log(kELoGeV);
  const double lnHi = std::log(kEHiGeV);
  for (int i = 0; i < kNNodes; ++i) {
    const double t = double(i) / double(kNNodes - 1);
    fE[static_cast<std::size_t>(i)] = std::exp(lnLo + (lnHi - lnLo) * t);
  }

  // Кумулятив по трапециям: интеграл dN/dE по каждому интервалу сетки.
  fCdf[0] = 0.0;
  for (int i = 1; i < kNNodes; ++i) {
    const double e0 = fE[static_cast<std::size_t>(i - 1)];
    const double e1 = fE[static_cast<std::size_t>(i)];
    const double area = 0.5 * (Gaisser(e0) + Gaisser(e1)) * (e1 - e0);
    fCdf[static_cast<std::size_t>(i)] =
        fCdf[static_cast<std::size_t>(i - 1)] + area;
  }
  const double total = fCdf.back();
  if (!(total > 0.0)) {
    std::fprintf(stderr,
                 "Rc103MuonSpectrum: FATAL cumulative integral is not "
                 "positive (%g)\n",
                 total);
    std::abort();
  }
  for (auto& c : fCdf) c /= total;

  std::fprintf(stdout,
               "Rc103MuonSpectrum: Gaisser E=%.3f..%.1f GeV, %d log nodes, "
               "raw integral=%.6e (arb.units)\n",
               kELoGeV, kEHiGeV, kNNodes, total);
}

double Rc103MuonSpectrum::SampleEnergyGeV() const {
  const double u = G4UniformRand();
  // upper_bound -> первый узел со значением CDF строго больше u.
  auto it = std::upper_bound(fCdf.begin(), fCdf.end(), u);
  std::size_t hi = static_cast<std::size_t>(it - fCdf.begin());
  if (hi == 0) hi = 1;
  if (hi >= fCdf.size()) hi = fCdf.size() - 1;
  const std::size_t lo = hi - 1;

  const double c0 = fCdf[lo];
  const double c1 = fCdf[hi];
  const double frac = (c1 > c0) ? (u - c0) / (c1 - c0) : 0.0;
  return fE[lo] + frac * (fE[hi] - fE[lo]);
}
