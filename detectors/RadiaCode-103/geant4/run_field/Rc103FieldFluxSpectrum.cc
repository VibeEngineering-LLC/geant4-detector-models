#include "Rc103FieldFluxSpectrum.hh"

#include "Randomize.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

bool Rc103FieldFluxSpectrum::Load(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    std::fprintf(stderr, "Rc103FieldFluxSpectrum: FAILED to open '%s'\n",
                 path.c_str());
    return false;
  }

  fEnergyKeV.clear();
  fCumulative.clear();
  fHeaderTotal = 0.0;
  fSumColumn = 0.0;
  bool headerFound = false;

  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;

    if (line[0] == '#') {
      const std::size_t key = line.find("fluence_total_cm2_s");
      if (key != std::string::npos) {
        const std::size_t eq = line.find('=', key);
        if (eq != std::string::npos) {
          fHeaderTotal = std::atof(line.substr(eq + 1).c_str());
          headerFound = true;
        }
      }
      continue;
    }

    // Строка-заголовок колонок "E_keV,fluence_cm2_s" — пропускаем.
    if (!line.empty() && (std::isalpha(static_cast<unsigned char>(line[0])) ||
                          line[0] == '_')) {
      continue;
    }

    const std::size_t comma = line.find(',');
    if (comma == std::string::npos) continue;
    const double e = std::atof(line.substr(0, comma).c_str());
    const double f = std::atof(line.substr(comma + 1).c_str());
    if (f <= 0.0) continue;  // нулевые/отрицательные бины в розыгрыш не идут

    fSumColumn += f;
    fEnergyKeV.push_back(e);
    fCumulative.push_back(fSumColumn);
  }

  if (!headerFound) {
    std::fprintf(stderr,
                 "Rc103FieldFluxSpectrum: FATAL - header line "
                 "'fluence_total_cm2_s = ...' not found in '%s'\n",
                 path.c_str());
    return false;
  }
  if (fEnergyKeV.empty()) {
    std::fprintf(stderr,
                 "Rc103FieldFluxSpectrum: FATAL - no positive data rows in "
                 "'%s'\n",
                 path.c_str());
    return false;
  }

  for (auto& c : fCumulative) c /= fSumColumn;
  fCumulative.back() = 1.0;  // защита от накопленной ошибки округления

  std::fprintf(stdout,
               "Rc103FieldFluxSpectrum: '%s' bins=%zu header_total=%.6e "
               "column_sum=%.6e\n",
               path.c_str(), fEnergyKeV.size(), fHeaderTotal, fSumColumn);

  // Самопроверка чтения входа (спека): расхождение >1% - WARNING, не падение.
  if (fHeaderTotal > 0.0) {
    const double rel = (fSumColumn - fHeaderTotal) / fHeaderTotal;
    if (std::abs(rel) > 0.01) {
      std::fprintf(stdout,
                   "Rc103FieldFluxSpectrum: WARNING column sum differs from "
                   "header by %.3f%% (sum=%.6e header=%.6e)\n",
                   rel * 100.0, fSumColumn, fHeaderTotal);
    }
  }
  return true;
}

double Rc103FieldFluxSpectrum::SampleEnergyKeV() const {
  const double u = G4UniformRand();
  const auto it = std::lower_bound(fCumulative.begin(), fCumulative.end(), u);
  const std::size_t idx =
      (it == fCumulative.end())
          ? fCumulative.size() - 1
          : static_cast<std::size_t>(it - fCumulative.begin());
  return fEnergyKeV[idx];
}
