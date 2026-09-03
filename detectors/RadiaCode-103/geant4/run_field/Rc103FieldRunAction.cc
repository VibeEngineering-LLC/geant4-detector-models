#include "Rc103FieldRunAction.hh"

#include "Rc103FieldDetectorConstruction.hh"  // состояние домика для CSV
#include "Rc103FieldPhysicsList.hh"           // порог и режим деэкситации для CSV

#include "G4LogicalVolume.hh"
#include "G4Navigator.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

Rc103FieldRunAction::Rc103FieldRunAction(std::string outputCsv,
                                         double fluxTotalCm2S, double surfaceCm2,
                                         double ratePerS, double tRunS,
                                         bool checkNormMode,
                                         double checkVolumeCm3)
    : fOutputCsv(std::move(outputCsv)),
      fFluxTotalCm2S(fluxTotalCm2S),
      fSurfaceCm2(surfaceCm2),
      fRatePerS(ratePerS),
      fTRunS(tRunS),
      fCheckNormMode(checkNormMode),
      fCheckVolumeCm3(checkVolumeCm3) {}

namespace {

// Контроль ОРИЕНТАЦИИ прибора (P-005, 30.08.2026). Кристалл сидит на 0.55 мм
// от центра габарита в сторону экрана, полувысота его 5 мм. Точка на 5 мм выше
// центра габарита попадает в кристалл ТОЛЬКО когда экран смотрит вверх: при
// экране вниз кристалл занимает [z-5.55, z+4.45] и точка оказывается снаружи.
// Одна точка различает две ориентации — проверка умеет краснеть (#SA-3).
void ProbeCrystalOrientation() {
  using DC = Rc103FieldDetectorConstruction;
  if (DC::GetCrystalLogicalVolume() == nullptr) return;   // --check-norm
  auto* nav = G4TransportationManager::GetTransportationManager()
                  ->GetNavigatorForTracking();
  const G4ThreeVector p(-49.5 * mm, 0.0, (DC::gDeviceZMm + 5.0) * mm);
  auto* pv = nav->LocateGlobalPointAndSetup(p, nullptr, false, true);
  const std::string got = pv ? pv->GetLogicalVolume()->GetName() : "(none)";
  const bool inCrystal = (got == DC::GetCrystalLogicalVolume()->GetName());
  std::cout << "Rc103FieldRunAction: контроль ориентации — точка (-49.5, 0, "
            << (DC::gDeviceZMm + 5.0) << ") в объёме '" << got << "'; экран "
            << (DC::gFlipUp ? "ВВЕРХ" : "вниз") << ", ожидание "
            << (DC::gFlipUp ? "кристалл" : "НЕ кристалл") << " — "
            << ((inCrystal == DC::gFlipUp) ? "СОШЛОСЬ" : "РАСХОЖДЕНИЕ")
            << std::endl;
}

}  // namespace

void Rc103FieldRunAction::BeginOfRunAction(const G4Run*) {
  ProbeCrystalOrientation();
  fNEvents = 0;
  fNHits = 0;
  fSumTrackLenMm = 0.0;
  fHistogram.fill(0);
  for (auto& h : fHistOrigin) h.fill(0);
  fNHitsOrigin.fill(0);
}

void Rc103FieldRunAction::RecordEvent(G4double edepMeV, int category) {
  ++fNEvents;
  if (edepMeV > 0.0) {
    ++fNHits;
    const double edepKeV = edepMeV / keV;
    const int bin = static_cast<int>(edepKeV / kBinWidthKeV);
    if (bin >= 0 && bin < kNBins) {
      ++fHistogram[static_cast<std::size_t>(bin)];
      if (category >= 0 && category < kNCat) {
        ++fHistOrigin[static_cast<std::size_t>(category)]
                     [static_cast<std::size_t>(bin)];
        ++fNHitsOrigin[static_cast<std::size_t>(category)];
      }
    }
  }
}

void Rc103FieldRunAction::EndOfRunAction(const G4Run* run) {
  if (!IsMaster()) return;  // Serial run manager: всегда true, пишем явно

  const long long nEventsFromRun = run->GetNumberOfEvent();
  if (nEventsFromRun != fNEvents) {
    std::cerr << "Rc103FieldRunAction: WARNING mismatch nEvents(G4Run)="
              << nEventsFromRun << " vs fNEvents(counted)=" << fNEvents << "\n";
  }

  const double cpsTotal = (fTRunS > 0.0) ? double(fNHits) / fTRunS : 0.0;

  double fluxMeasured = 0.0;
  double fluxRatio = 0.0;
  if (fCheckNormMode) {
    const double sumLenCm = fSumTrackLenMm / 10.0;  // мм -> см
    if (fCheckVolumeCm3 > 0.0 && fTRunS > 0.0) {
      fluxMeasured = sumLenCm / (fCheckVolumeCm3 * fTRunS);
    }
    if (fFluxTotalCm2S > 0.0) fluxRatio = fluxMeasured / fFluxTotalCm2S;
  }

  std::cout << "\n=== Rc103Field RESULT ===\n";
  std::cout << std::setprecision(10);
  std::cout << "n_events= " << fNEvents << "\n";
  std::cout << "n_hits_in_crystal= " << fNHits << "\n";
  std::cout << "flux_total_cm2_s= " << fFluxTotalCm2S << "\n";
  std::cout << "surface_cm2= " << fSurfaceCm2 << "\n";
  std::cout << "rate_per_s= " << fRatePerS << "\n";
  std::cout << "t_run_s= " << fTRunS << "\n";
  std::cout << "cps_total= " << cpsTotal << "\n";
  if (fCheckNormMode) {
    std::cout << "check_volume_cm3= " << fCheckVolumeCm3 << "\n";
    std::cout << "gamma_track_len_cm= " << (fSumTrackLenMm / 10.0) << "\n";
    std::cout << "FLUX_NOMINAL_cm2_s= " << fFluxTotalCm2S << "\n";
    std::cout << "FLUX_MEASURED_cm2_s= " << fluxMeasured << "\n";
    std::cout << "FLUX_RATIO(measured/nominal)= " << fluxRatio << "\n";
  }
  std::cout << "=== END Rc103Field RESULT ===\n" << std::endl;

  namespace fs = std::filesystem;
  fs::path outPath(fOutputCsv);
  if (outPath.has_parent_path()) {
    std::error_code ec;
    fs::create_directories(outPath.parent_path(), ec);
  }

  std::ofstream csv(fOutputCsv);
  if (!csv) {
    std::cerr << "Rc103FieldRunAction: FAILED to open output CSV: "
              << fOutputCsv << "\n";
    return;
  }
  csv << std::setprecision(10);
  csv << "metric,value\n";
  csv << "n_events," << fNEvents << "\n";
  csv << "n_hits_in_crystal," << fNHits << "\n";
  // Разбиение того же счёта по происхождению: мимо свинца / рассеян в свинце /
  // рождён в свинце (флуоресценция K-серии, тормозное).
  csv << "n_hits_direct," << fNHitsOrigin[0] << "\n";
  csv << "n_hits_pb_scat," << fNHitsOrigin[1] << "\n";
  csv << "n_hits_pb_born," << fNHitsOrigin[2] << "\n";
  csv << "flux_total_cm2_s," << fFluxTotalCm2S << "\n";
  csv << "surface_cm2," << fSurfaceCm2 << "\n";
  csv << "rate_per_s," << fRatePerS << "\n";
  csv << "t_run_s," << fTRunS << "\n";
  csv << "cps_total," << cpsTotal << "\n";
  // Состояние свинцового домика — чтобы CSV сам себя описывал и пару
  // «без домика / с домиком» нельзя было перепутать постфактум.
  {
    using DC = Rc103FieldDetectorConstruction;
    const bool shieldBuilt = (DC::GetShieldLogicalVolume() != nullptr);
    csv << "shield," << (shieldBuilt ? 1 : 0) << "\n";
    csv << "shield_pb_mm," << (shieldBuilt ? DC::kShieldPbMm : 0.0) << "\n";
    csv << "shield_cavity_mm,";
    if (shieldBuilt) {
      csv << DC::kShieldCavityXMm << "x" << DC::kShieldCavityYMm << "x"
          << DC::kShieldCavityZMm << "\n";
    } else {
      csv << "none\n";
    }
    // Посадка прибора (P-005, 30.08.2026): различение постановок по имени
    // файла — та самая дыра, из-за которой сравнение вышло без оговорки.
    csv << "stand_mm," << DC::gStandMm << "\n";
    csv << "device_z_mm," << DC::gDeviceZMm << "\n";
    csv << "screen_up," << (DC::gFlipUp ? 1 : 0) << "\n";
    // Физика (03.09.2026): порог продукции и режим деэкситации — иначе пара
    // прогонов std/deex/max различима только по имени файла.
    csv << "em_cut_mm," << Rc103FieldPhysicsList::gCutMm << "\n";
    csv << "em_deex," << Rc103FieldPhysicsList::gDeexMode << "\n";
  }
  if (fCheckNormMode) {
    csv << "check_volume_cm3," << fCheckVolumeCm3 << "\n";
    csv << "gamma_track_len_cm," << (fSumTrackLenMm / 10.0) << "\n";
    csv << "flux_measured_cm2_s," << fluxMeasured << "\n";
    csv << "flux_ratio," << fluxRatio << "\n";
  }
  // ⚠ Три первые колонки обязаны остаться на местах и с прежним смыслом: весь
  // разбор в analysis/ читает их позиционно.
  csv << "\nbin_keV,counts,cps,counts_direct,counts_pb_scat,counts_pb_born\n";
  for (int i = 0; i < kNBins; ++i) {
    const std::size_t b = static_cast<std::size_t>(i);
    const long long c = fHistogram[b];
    const double cps = (fTRunS > 0.0) ? double(c) / fTRunS : 0.0;
    csv << i << "," << c << "," << cps << "," << fHistOrigin[0][b] << ","
        << fHistOrigin[1][b] << "," << fHistOrigin[2][b] << "\n";
  }
  csv.close();
  std::cout << "Rc103FieldRunAction: wrote " << fOutputCsv << std::endl;
}
