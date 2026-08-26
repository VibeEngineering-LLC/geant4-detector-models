// Headless-рендер геометрии RC110Detector через файловый драйвер
// TOOLSSG_OFFSCREEN (не RayTracer — тот крашится 0xC0000005 на этой сборке
// Geant4 даже на тривиальной сцене, проверено фактом в проекте, см. бриф
// задачи / common/tools/gdml_vis.cc, откуда взят синтаксис vis-команд).
//
// Запуск:  rc110_vis_render.exe [макрос.mac]
//   макрос.mac — файл с командами /vis/... (по умолчанию — встроенный вид
//                "три четверти", см. DefaultView() ниже).
//
// Требует переменную PATH с C:\geant4\bin (см. g4setup.ps1) — без него
// EXITCODE=127/53 без диагностики Geant4, общая ловушка всех бинарников
// проекта, не специфика этого файла.
#include "RC110Detector.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4VModularPhysicsList.hh"

#include <cstdio>
#include <string>

namespace {

// Минимальный физ-лист: рендер геометрии не запускает ни одного события
// (BeamOn не вызывается), физика нужна только чтобы Initialize() не упал без
// зарегистрированного физ-листа. EM-физика — чтобы не зависеть от
// недокументированного поведения полностью пустого G4VModularPhysicsList
// (тот же выбор, что в RadiaCode-103/geometry/vis_render.cc).
class MinimalPhysList : public G4VModularPhysicsList {
 public:
  MinimalPhysList() {
    RegisterPhysics(new G4EmStandardPhysics_option4());
    SetDefaultCutValue(0.05 * mm);
  }
};

// Встроенный вид по умолчанию: три четверти, снаружи, весь прибор в кадре.
// Мировая система координат здесь МИРОВАЯ (не центр кристалла, как у
// RC-103) — прибор занимает X от -63.3 до +63.3, центр кадра в начале
// координат. Драйвер TSG_OFFSCREEN: открытие с разрешением "WxH" одной
// строкой, вывод файла ЗАДАЁТСЯ ДО /vis/drawVolume, обновление сцены —
// /vis/viewer/rebuild (НЕ /vis/viewer/refresh — это специфика
// TSG_OFFSCREEN, refresh был бы командой для RayTracer).
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

  runManager->SetUserInitialization(new RC110Detector());
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
    DefaultView(ui, "rc110_geant4_overview.png");
  }

  std::fprintf(stdout, "EXITCODE=0 macro=%s\n",
               macroFile.empty() ? "(default)" : macroFile.c_str());

  delete visManager;
  delete runManager;
  return 0;
}