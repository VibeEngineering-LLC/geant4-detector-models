// Headless-рендер геометрии ASN16Detector (АтомНано 16) через файловый драйвер
// TOOLSSG_OFFSCREEN — та же схема, что донорский
// geant4-detector-models/detectors/RadiaCode-110/geant4/geometry/vis_render.cc
// (RayTracer крашится 0xC0000005 на этой сборке, см. GEANT4\SESSION-STATE.md).
//
// Запуск: asn16_vis_render.exe [макрос.mac]
//   без аргумента — встроенный вид "три четверти" (см. DefaultView ниже).
#include "ASN16Detector.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4VModularPhysicsList.hh"

#include <cstdio>
#include <string>

namespace {

// Минимальный физ-лист: BeamOn не вызывается, физика нужна только чтобы
// Initialize() не упал без зарегистрированного физ-листа (как у донора).
class MinimalPhysList : public G4VModularPhysicsList {
 public:
  MinimalPhysList() {
    RegisterPhysics(new G4EmStandardPhysics_option4());
    SetDefaultCutValue(0.05 * mm);
  }
};

// Встроенный вид: три четверти, мировая система — центр кристалла (см.
// ASN16Detector.hh), прибор занимает Z от -43 до +43 мм — кадр по умолчанию
// на весь корпус (86х42х25 мм) с запасом.
void DefaultView(G4UImanager* ui, const std::string& outFile) {
  ui->ApplyCommand("/vis/open TOOLSSG_OFFSCREEN 1600x900");
  ui->ApplyCommand("/vis/tsg/offscreen/set/file " + outFile);
  ui->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 60 45 deg");
  ui->ApplyCommand("/vis/viewer/set/projection perspective 30 deg");
  ui->ApplyCommand("/vis/drawVolume");
  ui->ApplyCommand("/vis/viewer/rebuild");
}

}  // namespace

int main(int argc, char** argv) {
  const std::string macroFile = (argc > 1) ? argv[1] : "";

  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  runManager->SetVerboseLevel(0);

  runManager->SetUserInitialization(new ASN16Detector());
  runManager->SetUserInitialization(new MinimalPhysList());
  runManager->Initialize();

  auto* visManager = new G4VisExecutive("quiet");
  visManager->Initialize();

  G4UImanager* ui = G4UImanager::GetUIpointer();
  ui->ApplyCommand("/control/verbose 0");
  ui->ApplyCommand("/run/verbose 0");

  if (!macroFile.empty()) {
    ui->ExecuteMacroFile(macroFile.c_str());
  } else {
    DefaultView(ui, "asn16_geant4_overview.png");
  }

  std::fprintf(stdout, "EXITCODE=0 macro=%s\n",
               macroFile.empty() ? "(default)" : macroFile.c_str());

  delete visManager;
  delete runManager;
  return 0;
}

