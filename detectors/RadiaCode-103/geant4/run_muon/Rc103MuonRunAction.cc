#include "Rc103MuonRunAction.hh"

#include "G4Run.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

Rc103MuonRunAction::Rc103MuonRunAction(std::string outputCsv, double rDiskMm,
                                       double zDiskMm, double diskAreaCm2,
                                       double pdgExpectedPerS, double eLoGeV,
                                       double eHiGeV)
    : fOutputCsv(std::move(outputCsv)),
      fRDiskMm(rDiskMm),
      fZDiskMm(zDiskMm),
      fDiskAreaCm2(diskAreaCm2),
      fPdgExpectedPerS(pdgExpectedPerS),
      fELoGeV(eLoGeV),
      fEHiGeV(eHiGeV) {}

void Rc103MuonRunAction::BeginOfRunAction(const G4Run*) {
  fNEvents = 0;
  fNHits = 0;
  fNOverflow = 0;
  fMaxEdepKeV = 0.0;
  fHistogram.fill(0);
}

void Rc103MuonRunAction::RecordEvent(G4double edepMeV) {
  ++fNEvents;
  if (edepMeV <= 0.0) return;
  ++fNHits;
  const double edepKeV = edepMeV / keV;
  fMaxEdepKeV = std::max(fMaxEdepKeV, edepKeV);
  const int bin = static_cast<int>(edepKeV / kBinWidthKeV);
  if (bin >= 0 && bin < kNBins) {
    ++fHistogram[static_cast<std::size_t>(bin)];
  } else {
    // >= 3000 кэВ (или отрицательный, чего быть не должно) — переполнение.
    // НЕ дописывается в последний бин.
    ++fNOverflow;
  }
}

void Rc103MuonRunAction::EndOfRunAction(const G4Run* run) {
  if (!IsMaster()) return;  // Serial run manager: всегда true, пишем явно

  const long long nEventsFromRun = run->GetNumberOfEvent();
  if (nEventsFromRun != fNEvents) {
    std::cerr << "Rc103MuonRunAction: WARNING mismatch nEvents(G4Run)="
              << nEventsFromRun << " vs fNEvents(counted)=" << fNEvents << "\n";
  }

  long long inHist = 0;
  for (int i = 0; i < kNBins; ++i) inHist += fHistogram[static_cast<std::size_t>(i)];
  if (inHist + fNOverflow != fNHits) {
    std::cerr << "Rc103MuonRunAction: WARNING hit bookkeeping mismatch: "
              << "sum(hist)=" << inHist << " + overflow=" << fNOverflow
              << " != n_hits=" << fNHits << "\n";
  }

  std::cout << "\n=== Rc103Muon RESULT ===\n";
  std::cout << std::setprecision(10);
  std::cout << "n_events= " << fNEvents << "\n";
  std::cout << "n_hits_in_crystal= " << fNHits << "\n";
  std::cout << "n_overflow= " << fNOverflow << "\n";
  std::cout << "sum_hist_in_range= " << inHist << "\n";
  std::cout << "max_edep_keV= " << fMaxEdepKeV << "\n";
  std::cout << "r_disk_mm= " << fRDiskMm << "\n";
  std::cout << "z_disk_mm= " << fZDiskMm << "\n";
  std::cout << "disk_area_cm2= " << fDiskAreaCm2 << "\n";
  // Справочная величина порядка, НЕ ограничение (спека, раздел "Нормировка").
  std::cout << "pdg_expected_per_s= " << fPdgExpectedPerS << "\n";
  std::cout << "e_lo_gev= " << fELoGeV << "\n";
  std::cout << "e_hi_gev= " << fEHiGeV << "\n";
  std::cout << "=== END Rc103Muon RESULT ===\n" << std::endl;

  namespace fs = std::filesystem;
  fs::path outPath(fOutputCsv);
  if (outPath.has_parent_path()) {
    std::error_code ec;
    fs::create_directories(outPath.parent_path(), ec);
  }

  std::ofstream csv(fOutputCsv);
  if (!csv) {
    std::cerr << "Rc103MuonRunAction: FAILED to open output CSV: " << fOutputCsv
              << "\n";
    return;
  }
  csv << std::setprecision(10);
  csv << "metric,value\n";
  csv << "n_events," << fNEvents << "\n";
  csv << "n_hits_in_crystal," << fNHits << "\n";
  csv << "n_overflow," << fNOverflow << "\n";
  csv << "r_disk_mm," << fRDiskMm << "\n";
  csv << "z_disk_mm," << fZDiskMm << "\n";
  csv << "disk_area_cm2," << fDiskAreaCm2 << "\n";
  csv << "pdg_expected_per_s," << fPdgExpectedPerS << "\n";
  csv << "e_lo_gev," << fELoGeV << "\n";
  csv << "e_hi_gev," << fEHiGeV << "\n";
  csv << "max_edep_keV," << fMaxEdepKeV << "\n";

  csv << "\nbin_keV,counts,per_muon\n";
  for (int i = 0; i < kNBins; ++i) {
    const long long c = fHistogram[static_cast<std::size_t>(i)];
    const double perMuon = (fNEvents > 0) ? double(c) / double(fNEvents) : 0.0;
    csv << i << "," << c << "," << perMuon << "\n";
  }
  csv.close();
  std::cout << "Rc103MuonRunAction: wrote " << fOutputCsv << std::endl;
}
