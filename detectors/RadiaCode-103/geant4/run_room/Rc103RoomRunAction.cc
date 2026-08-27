#include "Rc103RoomRunAction.hh"

#include "G4Run.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>

Rc103RoomRunAction::Rc103RoomRunAction(std::string outputCsv,
                                        double wallEmissionRatePerSec)
    : fOutputCsv(std::move(outputCsv)),
      fWallEmissionRatePerSec(wallEmissionRatePerSec) {}

void Rc103RoomRunAction::BeginOfRunAction(const G4Run*) {
  fNEvents = 0;
  fNHits = 0;
  fHistogram.fill(0);
}

void Rc103RoomRunAction::RecordEvent(G4double edepMeV) {
  ++fNEvents;
  if (edepMeV > 0.0) {
    ++fNHits;
    const double edepKeV = edepMeV / keV;
    const int bin = static_cast<int>(edepKeV / kBinWidthKeV);
    if (bin >= 0 && bin < kNBins) {
      ++fHistogram[static_cast<std::size_t>(bin)];
    }
  }
}

void Rc103RoomRunAction::EndOfRunAction(const G4Run* run) {
  if (!IsMaster()) return;  // Serial run manager: всегда true, пишем явно

  const long long nEventsFromRun = run->GetNumberOfEvent();
  if (nEventsFromRun != fNEvents) {
    std::cerr << "Rc103RoomRunAction: WARNING mismatch nEvents(G4Run)="
              << nEventsFromRun << " vs fNEvents(counted)=" << fNEvents
              << "\n";
  }

  const double effInSim = (fNEvents > 0) ? double(fNHits) / double(fNEvents) : 0.0;
  // Абсолютная прогнозная скорость счёта в кристалле от ЕРН в стенах, cps —
  // единственный по-настоящему НОВЫЙ (по сравнению с wallfield.cc, который
  // считает только флюенс в пустой полости) результат этого приложения.
  const double predictedCps = effInSim * fWallEmissionRatePerSec;

  std::cout << "\n=== Rc103Room RESULT ===\n";
  std::cout << "N_events= " << fNEvents << "\n";
  std::cout << "N_hits_in_crystal(depE>0)= " << fNHits << "\n";
  std::cout << "eff_crystal_per_event(hits/N)= " << effInSim << "\n";
  std::cout << "wall_emission_rate_photon_per_s= " << fWallEmissionRatePerSec
             << "\n";
  std::cout << "predicted_crystal_count_rate_cps= " << predictedCps << "\n";
  if (fNHits < 100) {
    std::cout << "WARNING: LOW STATISTICS — N_hits=" << fNHits
              << " < 100. Без importance biasing (см. задание/P-004/P-009 —"
                 " биасинг НЕ внедрялся в этом прогоне). Спектр по бинам и"
                 " predicted_crystal_count_rate_cps ниже статистически"
                 " ненадёжны, приводятся как честное ограничение smoke-теста,"
                 " не как итоговая оценка фона.\n";
  }
  std::cout << "=== END Rc103Room RESULT ===\n" << std::endl;

  namespace fs = std::filesystem;
  fs::path outPath(fOutputCsv);
  if (outPath.has_parent_path()) {
    std::error_code ec;
    fs::create_directories(outPath.parent_path(), ec);
  }

  std::ofstream csv(fOutputCsv);
  if (!csv) {
    std::cerr << "Rc103RoomRunAction: FAILED to open output CSV: "
              << fOutputCsv << "\n";
    return;
  }
  csv << "metric,value\n";
  csv << "n_events," << fNEvents << "\n";
  csv << "n_hits_in_crystal," << fNHits << "\n";
  csv << "eff_crystal_per_event," << effInSim << "\n";
  csv << "wall_emission_rate_photon_per_s," << fWallEmissionRatePerSec << "\n";
  csv << "predicted_crystal_count_rate_cps," << predictedCps << "\n";
  csv << "low_statistics_warning," << (fNHits < 100 ? "1" : "0") << "\n";
  csv << "\nbin_keV,counts\n";
  for (int i = 0; i < kNBins; ++i) {
    if (fHistogram[static_cast<std::size_t>(i)] > 0) {
      csv << (i * kBinWidthKeV) << "," << fHistogram[static_cast<std::size_t>(i)]
          << "\n";
    }
  }
  csv.close();
  std::cout << "Rc103RoomRunAction: wrote " << fOutputCsv << std::endl;
}
