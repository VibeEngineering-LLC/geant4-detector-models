// RadiaCode-103 — рабочий физический прогон (BeamOn) на провалидированной
// GDML-геометрии (RC103_detector.gdml, SSOT), не только рендер. Точечный
// источник Cs-137 (661.657 кэВ) в 100 мм от лицевой грани корпуса, узкий
// конус на кристалл, полная и фотопиковая эффективность per BeamOn.
//
// НЕ использует и НЕ читает detectors/RadiaCode-103/geometry/{main.cc,
// RCDetector.cc,.hh} (другая, незакоммиченная физическая линия того же
// прибора — оператором запрещена как образец/источник кода). Вся геометрия
// грузится напрямую из SSOT через G4GDMLParser, паттерн
// common/tools/gdml_vis.cc (GdmlDetector).
//
// usage: rc103_run.exe [gdml_path] [n_events] [out_csv]
#include "Rc103RunDetectorConstruction.hh"
#include "Rc103RunEventAction.hh"
#include "Rc103RunPhysicsList.hh"
#include "Rc103RunPrimaryGeneratorAction.hh"
#include "Rc103RunRunAction.hh"
#include "Rc103RunSteppingAction.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#ifndef RC103_GDML_DEFAULT_PATH
#error "RC103_GDML_DEFAULT_PATH must be defined by CMake (see run/CMakeLists.txt)"
#endif
#ifndef RC103_OUTPUT_DEFAULT_PATH
#error "RC103_OUTPUT_DEFAULT_PATH must be defined by CMake (see run/CMakeLists.txt)"
#endif

int main(int argc, char** argv) {
  const std::string gdmlPath = (argc > 1) ? argv[1] : RC103_GDML_DEFAULT_PATH;
  const long long nEvents = (argc > 2) ? std::atoll(argv[2]) : 1000000LL;
  const std::string outCsv = (argc > 3) ? argv[3] : RC103_OUTPUT_DEFAULT_PATH;

  // Omega_cone/4pi конуса источника (mintheta=175 maxtheta=180, полураствор
  // 5°) — см. полное обоснование в Rc103RunPrimaryGeneratorAction.cc.
  // Omega_cone = 2*pi*(1-cos(halfAngle)).
  const double halfAngleDeg = Rc103RunPrimaryGeneratorAction::kConeHalfAngleDeg;
  const double halfAngleRad = halfAngleDeg * M_PI / 180.0;
  const double omegaCone = 2.0 * M_PI * (1.0 - std::cos(halfAngleRad));
  const double omegaFraction = omegaCone / (4.0 * M_PI);

  std::fprintf(stdout,
               "rc103_run: gdml=%s n_events=%lld out_csv=%s "
               "cone_half_angle_deg=%.3f omega_cone_over_4pi=%.8f\n",
               gdmlPath.c_str(), nEvents, outCsv.c_str(), halfAngleDeg,
               omegaFraction);

  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  runManager->SetVerboseLevel(0);

  auto* detector = new Rc103RunDetectorConstruction(gdmlPath);
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Rc103RunPhysicsList());

  auto* runAction = new Rc103RunRunAction(outCsv, omegaFraction);
  auto* eventAction = new Rc103RunEventAction(runAction);
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(new Rc103RunPrimaryGeneratorAction());
  runManager->SetUserAction(new Rc103RunSteppingAction(eventAction));

  runManager->Initialize();

  if (!Rc103RunDetectorConstruction::GetCrystalLogicalVolume()) {
    std::fprintf(stderr,
                 "rc103_run: FATAL crystal logical volume not resolved "
                 "after Initialize() — aborting before BeamOn.\n");
    delete runManager;
    return 3;
  }

  G4UImanager* ui = G4UImanager::GetUIpointer();
  ui->ApplyCommand("/run/verbose 1");
  ui->ApplyCommand("/event/verbose 0");
  ui->ApplyCommand("/tracking/verbose 0");
  const long long progressEvery = (nEvents / 10 > 0) ? (nEvents / 10) : 1;
  ui->ApplyCommand("/run/printProgress " + std::to_string(progressEvery));

  runManager->BeamOn(static_cast<G4int>(nEvents));

  delete runManager;
  std::fprintf(stdout, "EXITCODE=0\n");
  return 0;
}
