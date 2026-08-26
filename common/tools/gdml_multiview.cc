// gdml_multiview.cc — headless-рендер ЛЮБОГО GDML-файла в 4 вида одним
// прогоном (спереди / сбоку / разрез Y=0 / три четверти). Общий инструмент,
// НЕ привязан к RC-103/RC-110 — путь к GDML приходит аргументом, поэтому
// живёт в common/tools/, не в detectors/<X>/geometry/ (тот же принцип, что
// у gdml_vis.cc).
//
// Склейка двух эталонов, дословно (не переизобретено):
//   - загрузка GDML (GdmlDetector/MinimalPhysList/RunManagerFactory,
//     /vis/open TSG_OFFSCREEN) — common/tools/gdml_vis.cc:13-64,70-73
//   - командная последовательность 4 видов (viewpointVector/upVector/
//     sectionPlane/iso с явным сбросом upVector) —
//     detectors/RadiaCode-110/geant4/geometry/vis_multiview.cc:41-143
//     (включая баг, исправленный 26.08.2026: без явного
//     "/vis/viewer/set/upVector 0 1 0" перед iso-видом камера наследует
//     upVector "0 0 1" от разреза (вид 3) и iso зеркалится).
//
// Запуск:  gdml_multiview.exe <in.gdml> <outDir> <baseName> [W] [H]
//   Пишет <outDir>/<baseName>_front.png, _side.png, _section.png, _iso.png
//   W,H опциональны, по умолчанию 1600x900 (как в vis_multiview.cc).
//   Разрез — плоскость Y=0, нормаль "0 1 0" (тот же выбор нормали, что
//   значением по умолчанию в vis_multiview.cc).
#include "G4GDMLParser.hh"
#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserDetectorConstruction.hh"

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>

namespace {

// Дословно gdml_vis.cc:27-39 — загрузка произвольного GDML-файла по пути.
class GdmlDetector : public G4VUserDetectorConstruction {
 public:
  explicit GdmlDetector(const G4String& path) : fPath(path) {}
  G4VPhysicalVolume* Construct() override {
    fParser.SetOverlapCheck(true);  // не глушить наложения объёмов
    fParser.Read(fPath, false);     // false = не валидировать против схемы
    return fParser.GetWorldVolume();
  }

 private:
  G4String fPath;
  G4GDMLParser fParser;
};

// Дословно gdml_vis.cc:43 — рендер геометрии не запускает событий (BeamOn
// не вызывается), физика нужна только чтобы Initialize() не упал без
// зарегистрированного физ-листа.
class MinimalPhysList : public G4VModularPhysicsList {};

// Один кадр на уже открытом TSG_OFFSCREEN-вьюере: сменить выходной файл,
// ракурс, проекцию и (опционально) режущую плоскость, перерисовать. Тот же
// порядок команд, что RenderView() в vis_multiview.cc:59-77.
// sectionPlaneArgs пустая строка -> "/vis/viewer/set/sectionPlane off"
// (важно явно выключать между кадрами — вьюер держит состояние между
// вызовами /vis/viewer/rebuild, иначе разрез "утёк" бы в следующий вид).
void RenderView(G4UImanager* ui, const std::string& outFile,
                 const std::string& viewpointVector,
                 const std::string& upVector, const std::string& projection,
                 const std::string& sectionPlaneArgs) {
  ui->ApplyCommand("/vis/tsg/offscreen/set/file " + outFile);
  ui->ApplyCommand("/vis/viewer/set/projection " + projection);
  ui->ApplyCommand("/vis/viewer/set/viewpointVector " + viewpointVector);
  ui->ApplyCommand("/vis/viewer/set/upVector " + upVector);
  if (sectionPlaneArgs.empty()) {
    ui->ApplyCommand("/vis/viewer/set/sectionPlane off");
  } else {
    ui->ApplyCommand("/vis/viewer/set/sectionPlane " + sectionPlaneArgs);
  }
  ui->ApplyCommand("/vis/drawVolume");
  ui->ApplyCommand("/vis/viewer/rebuild");
  std::fprintf(stdout, "  wrote %s (viewpoint %s, proj %s, section %s)\n",
               outFile.c_str(), viewpointVector.c_str(), projection.c_str(),
               sectionPlaneArgs.empty() ? "off" : sectionPlaneArgs.c_str());
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 4) {
    std::fprintf(stderr,
                  "usage: gdml_multiview.exe <in.gdml> <outDir> <baseName> "
                  "[W] [H]\n");
    return 2;
  }
  const std::string gdmlPath = argv[1];
  const std::string outDir = argv[2];
  const std::string baseName = argv[3];
  const int w = (argc > 4) ? std::atoi(argv[4]) : 1600;
  const int h = (argc > 5) ? std::atoi(argv[5]) : 900;

  std::filesystem::create_directories(outDir);
  const std::string prefix = outDir + "/" + baseName;

  auto* runManager =
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  runManager->SetVerboseLevel(0);
  runManager->SetUserInitialization(new GdmlDetector(gdmlPath));
  runManager->SetUserInitialization(new MinimalPhysList());
  runManager->Initialize();

  auto* visManager = new G4VisExecutive("quiet");
  visManager->Initialize();

  G4UImanager* ui = G4UImanager::GetUIpointer();
  ui->ApplyCommand("/control/verbose 0");
  ui->ApplyCommand("/run/verbose 0");
  ui->ApplyCommand("/vis/open TSG_OFFSCREEN " + std::to_string(w) + "x" +
                    std::to_string(h));

  // 1. Спереди: камера на -Z смотрит в +Z.
  RenderView(ui, prefix + "_front.png", "0 0 -1", "0 1 0", "orthogonal", "");

  // 2. Сбоку: камера на -Y смотрит в +Y.
  RenderView(ui, prefix + "_side.png", "0 -1 0", "0 0 1", "orthogonal", "");

  // 3. Разрез Y=0: тот же ракурс, что вид 2, + настоящая режущая плоскость
  //    Geant4 (sectionPlane), не транспарентность. Нормаль "0 1 0" — тот же
  //    выбор, что дефолт в vis_multiview.cc.
  RenderView(ui, prefix + "_section.png", "0 -1 0", "0 0 1", "orthogonal",
             "on 0 0 0 mm 0 1 0");

  // 4. Три четверти — не через RenderView() (та ждёт viewpointVector, тут
  //    viewpointThetaPhi), порядок команд дословно vis_multiview.cc:121-135.
  //    sectionPlane выключаем явно, иначе унаследовался бы разрез от вида 3.
  //    upVector ЯВНО сброшен на дефолт Geant4 (0 1 0) — без этой строки вид
  //    наследует "0 0 1", оставленный видом 3 (разрез), и камера
  //    зеркалится (баг, исправленный 26.08.2026 в оригинале vis_multiview.cc).
  ui->ApplyCommand("/vis/tsg/offscreen/set/file " + prefix + "_iso.png");
  ui->ApplyCommand("/vis/viewer/set/sectionPlane off");
  ui->ApplyCommand("/vis/viewer/set/upVector 0 1 0");
  ui->ApplyCommand("/vis/viewer/set/projection perspective 30 deg");
  ui->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 60 45 deg");
  ui->ApplyCommand("/vis/drawVolume");
  ui->ApplyCommand("/vis/viewer/rebuild");
  std::fprintf(stdout,
               "  wrote %s_iso.png (viewpointThetaPhi 60 45 deg, proj "
               "perspective 30 deg, section off)\n",
               prefix.c_str());

  std::fprintf(stdout, "EXITCODE=0 gdml=%s outDir=%s baseName=%s\n",
               gdmlPath.c_str(), outDir.c_str(), baseName.c_str());

  delete visManager;
  delete runManager;
  return 0;
}
