#include "Rc103RunRunAction.hh"

#include "G4Run.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>

Rc103RunRunAction::Rc103RunRunAction(std::string outputCsv, double solidAngleFractionCone)
    : fOutputCsv(std::move(outputCsv)), fSolidAngleFractionCone(solidAngleFractionCone) {}

void Rc103RunRunAction::BeginOfRunAction(const G4Run*) {
  fNEvents = 0;
  fNHits = 0;
  fNPhotopeak = 0;
  fHistogram.fill(0);
}

void Rc103RunRunAction::RecordEvent(G4double edepMeV) {
  ++fNEvents;
  if (edepMeV > 0.0) {
    ++fNHits;
    const double edepKeV = edepMeV / keV;
    const int bin = static_cast<int>(edepKeV / kBinWidthKeV);
    if (bin >= 0 && bin < kNBins) {
      ++fHistogram[static_cast<std::size_t>(bin)];
    }
    if (std::fabs(edepKeV - kLineKeV) <= kWindowHalfKeV) {
      ++fNPhotopeak;
    }
  }
}

void Rc103RunRunAction::EndOfRunAction(const G4Run* run) {
  if (!IsMaster()) return;  // Serial run manager: всегда true, пишем явно

  const long long nEventsFromRun = run->GetNumberOfEvent();
  if (nEventsFromRun != fNEvents) {
    std::cerr << "Rc103RunRunAction: WARNING mismatch nEvents(G4Run)="
              << nEventsFromRun << " vs fNEvents(counted)=" << fNEvents
              << "\n";
  }

  const double effFullInCone =
      (fNEvents > 0) ? double(fNHits) / double(fNEvents) : 0.0;
  const double effPeakInCone =
      (fNEvents > 0) ? double(fNPhotopeak) / double(fNEvents) : 0.0;
  const double effFullAbsolute = effFullInCone * fSolidAngleFractionCone;
  const double effPeakAbsolute = effPeakInCone * fSolidAngleFractionCone;

  std::cout << "\n=== Rc103Run RESULT ===\n";
  std::cout << "N_events= " << fNEvents << "\n";
  std::cout << "N_hits(depE>0)= " << fNHits << "\n";
  std::cout << "N_photopeak(|E-" << kLineKeV << "keV|<=" << kWindowHalfKeV
             << "keV)= " << fNPhotopeak << "\n";
  std::cout << "eff_full_inCone(hits/N)= " << effFullInCone << "\n";
  std::cout << "eff_photopeak_inCone(hits/N)= " << effPeakInCone << "\n";
  std::cout << "Omega_cone_over_4pi= " << fSolidAngleFractionCone << "\n";
  std::cout << "eff_full_ABSOLUTE(vs isotropic 4pi source)= " << effFullAbsolute
             << "\n";
  std::cout << "eff_photopeak_ABSOLUTE(vs isotropic 4pi source)= "
             << effPeakAbsolute << "\n";
  std::cout << "=== END Rc103Run RESULT ===\n" << std::endl;

  namespace fs = std::filesystem;
  fs::path outPath(fOutputCsv);
  if (outPath.has_parent_path()) {
    std::error_code ec;
    fs::create_directories(outPath.parent_path(), ec);
  }

  std::ofstream csv(fOutputCsv);
  if (!csv) {
    std::cerr << "Rc103RunRunAction: FAILED to open output CSV: "
              << fOutputCsv << "\n";
    return;
  }
  csv << "metric,value\n";
  csv << "n_events," << fNEvents << "\n";
  csv << "n_hits," << fNHits << "\n";
  csv << "n_photopeak," << fNPhotopeak << "\n";
  csv << "eff_full_in_cone," << effFullInCone << "\n";
  csv << "eff_photopeak_in_cone," << effPeakInCone << "\n";
  csv << "omega_cone_over_4pi," << fSolidAngleFractionCone << "\n";
  csv << "eff_full_absolute," << effFullAbsolute << "\n";
  csv << "eff_photopeak_absolute," << effPeakAbsolute << "\n";
  csv << "\nbin_keV,counts\n";
  for (int i = 0; i < kNBins; ++i) {
    csv << i << "," << fHistogram[static_cast<std::size_t>(i)] << "\n";
  }
  csv.close();
  std::cout << "Rc103RunRunAction: wrote " << fOutputCsv << std::endl;
}
