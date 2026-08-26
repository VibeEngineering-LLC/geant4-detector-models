// mini Geant4 vis — универсальный headless-рендер ЛЮБОГО GDML-файла через
// драйвер TOOLSSG_OFFSCREEN (не RayTracer — тот крашится 0xC0000005 на этой
// сборке, см. GEANT4\SESSION-STATE.md, линия B/#AUD-3, все 4 гипотезы
// отвергнуты фактом). TSG_OFFSCREEN не тестировался ДО этого инструмента —
// синтаксис взят из эталонного examples/basic/B5/tsg_offscreen.mac
// (WebFetch, не по памяти — тот же класс риска, что уже дважды подвёл
// с RayTracer-командами в этой сессии).
//
// Общий инструмент, НЕ привязан к RC-103/RC-110 — принимает путь к GDML
// аргументом, поэтому живёт в common/tools/, не в detectors/<Х>/geometry/.
//
// Запуск:  gdml_vis.exe <файл.gdml> <вых.png> [ширина] [высота] [theta] [phi]
#include "G4GDMLParser.hh"
#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserDetectorConstruction.hh"

#include <cstdio>
#include <cstdlib>
#include <string>

namespace {

class GdmlDetector : public G4VUserDetectorConstruction {
 public:
  explicit GdmlDetector(const G4String& path) : fPath(path) {}
  G4VPhysicalVolume* Construct() override {
    fParser.SetOverlapCheck(true);  // не глушить P-025-подобные наложения
    fParser.Read(fPath, false);     // false = не валидировать против схемы
    return fParser.GetWorldVolume();
  }

 private:
  G4String fPath;
  G4GDMLParser fParser;
};

// Минимальный физ-лист — рендер геометрии не запускает событий (BeamOn не
// вызывается), физика нужна только чтобы Initialize() не упал без списка.
class MinimalPhysList : public G4VModularPhysicsList {};

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: gdml_vis.exe <in.gdml> <out.png> "
                          "[W] [H] [thetaDeg] [phiDeg]\n");
    return 2;
  }
  const std::string gdmlPath = argv[1];
  const std::string outPng = argv[2];
  const int w = (argc > 3) ? std::atoi(argv[3]) : 1200;
  const int h = (argc > 4) ? std::atoi(argv[4]) : 900;
  const std::string theta = (argc > 5) ? argv[5] : "60";
  const std::string phi = (argc > 6) ? argv[6] : "35";

  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  runManager->SetVerboseLevel(0);
  runManager->SetUserInitialization(new GdmlDetector(gdmlPath));
  runManager->SetUserInitialization(new MinimalPhysList());
  runManager->Initialize();

  auto* visManager = new G4VisExecutive("all");
  visManager->Initialize();

  G4UImanager* ui = G4UImanager::GetUIpointer();
  ui->ApplyCommand("/control/verbose 0");
  ui->ApplyCommand("/run/verbose 0");
  ui->ApplyCommand("/vis/open TSG_OFFSCREEN " + std::to_string(w) + "x" +
                    std::to_string(h));
  ui->ApplyCommand("/vis/tsg/offscreen/set/file " + outPng);
  ui->ApplyCommand("/vis/viewer/set/viewpointThetaPhi " + theta + " " + phi +
                    " deg");
  ui->ApplyCommand("/vis/drawVolume");
  ui->ApplyCommand("/vis/viewer/rebuild");

  std::fprintf(stdout, "EXITCODE=0 gdml=%s out=%s\n", gdmlPath.c_str(),
               outPng.c_str());

  delete visManager;
  delete runManager;
  return 0;
}
