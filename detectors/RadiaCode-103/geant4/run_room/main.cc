// RadiaCode-103 — «поле помещения»: реальный отклик (edep в Crystal_log)
// провалидированного GDML-детектора (SSOT RC103_detector.gdml) на равномерно
// распределённые в бетонных стенах ЕРН (K-40/Ra-226 ряд/Th-232 ряд, UNSCEAR
// активности), а не флюенс в пустой полости (тот считает wallfield.cc —
// другая, старая физическая линия того же прибора, оператором запрещённая
// как код-образец; здесь переиспользована только физика/данные, §33).
//
// БЕЗ importance biasing (сознательно — задание): геометрическая
// эффективность маленького кристалла (10x10x10мм) в центре 80-см бетонной
// сферы крайне мала, статистика в кристалле честно может оказаться < 100
// hits даже на несколько миллионов первичных — предупреждение печатается,
// не скрывается (см. Rc103RoomRunAction.cc).
//
// usage: rc103_room.exe [gdml_path] [n_events] [out_csv]
#include "Rc103RoomDetectorConstruction.hh"
#include "Rc103RoomEventAction.hh"
#include "Rc103RoomPhysicsList.hh"
#include "Rc103RoomPrimaryGeneratorAction.hh"
#include "Rc103RoomRunAction.hh"
#include "Rc103RoomSteppingAction.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"

#include <cstdio>
#include <cstdlib>
#include <string>

#ifndef RC103_ROOM_GDML_DEFAULT_PATH
#error "RC103_ROOM_GDML_DEFAULT_PATH must be defined by CMake (see run_room/CMakeLists.txt)"
#endif
#ifndef RC103_ROOM_OUTPUT_DEFAULT_PATH
#error "RC103_ROOM_OUTPUT_DEFAULT_PATH must be defined by CMake (see run_room/CMakeLists.txt)"
#endif

int main(int argc, char** argv) {
  const std::string gdmlPath = (argc > 1) ? argv[1] : RC103_ROOM_GDML_DEFAULT_PATH;
  const long long nEvents = (argc > 2) ? std::atoll(argv[2]) : 3000000LL;
  const std::string outCsv = (argc > 3) ? argv[3] : RC103_ROOM_OUTPUT_DEFAULT_PATH;

  std::fprintf(stdout,
               "rc103_room: gdml=%s n_events=%lld out_csv=%s "
               "(sphere R_wall=80cm concrete / R_cav=20cm air, K-40+Ra226+Th232 "
               "UNSCEAR, NO importance biasing)\n",
               gdmlPath.c_str(), nEvents, outCsv.c_str());

  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  runManager->SetVerboseLevel(0);

  auto* detector = new Rc103RoomDetectorConstruction(gdmlPath);
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Rc103RoomPhysicsList());

  auto* primaryGen = new Rc103RoomPrimaryGeneratorAction();
  auto* runAction = new Rc103RoomRunAction(
      outCsv, Rc103RoomPrimaryGeneratorAction::GetWallEmissionRatePerSec());
  auto* eventAction = new Rc103RoomEventAction(runAction);
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(primaryGen);
  runManager->SetUserAction(new Rc103RoomSteppingAction(eventAction));

  runManager->Initialize();

  if (!Rc103RoomDetectorConstruction::GetCrystalLogicalVolume()) {
    std::fprintf(stderr,
                 "rc103_room: FATAL crystal logical volume not resolved "
                 "after Initialize() — aborting before BeamOn.\n");
    delete runManager;
    return 3;
  }
  if (!Rc103RoomDetectorConstruction::GetCavityLogicalVolume()) {
    std::fprintf(stderr,
                 "rc103_room: FATAL cavity logical volume not resolved "
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
