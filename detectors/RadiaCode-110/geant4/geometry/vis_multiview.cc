// Headless-рендер геометрии RC110Detector — 4 отдельных вида одним прогоном
// (спереди / сбоку / разрез Y=0 / три четверти), тот же драйвер
// TOOLSSG_OFFSCREEN, что в vis_render.cc (RayTracer крашится 0xC0000005 на
// этой сборке — не пробовать заново, факт уже проверен в проекте).
//
// Донор структуры: vis_render.cc того же каталога (run/vis manager setup,
// MinimalPhysList — дословно тот же класс). Новое здесь — ТОЛЬКО
// многокадровая логика (несколько /vis/tsg/offscreen/set/file на одном
// открытом вьюере) и sectionPlane, которых в проекте раньше не было (grep
// по репозиторию на "sectionPlane" пуст на момент написания — проверено,
// не переизобретение, а первое использование).
//
// Запуск:  rc110_vis_multiview.exe ["nx ny nz"]
//   Необязательный аргумент — нормаль плоскости разреза для вида 3
//   (разрез Y=0), по умолчанию "0 1 0". Синтаксис /vis/viewer/set/
//   sectionPlane не задокументирован достаточно однозначно, чтобы угадать
//   знак нормали заранее (какая половина остаётся) — предусмотрена
//   быстрая эмпирическая проверка БЕЗ пересборки: запустить второй раз
//   с "0 -1 0" и сравнить PNG (Pillow getcolors()/побайтно) с видом 2
//   (сбоку, без разреза) и с первым результатом.
//
// Пишет в текущий рабочий каталог (cwd процесса = build/RadiaCode-110 при
// запуске из _build.cmd/ninja-выхода):
//   rc110_view_front.png   — спереди, вдоль -Z, ortho
//   rc110_view_side.png    — сбоку, вдоль -Y, ortho
//   rc110_view_section.png — сбоку + sectionPlane on (разрез Y=0)
//   rc110_view_iso.png     — три четверти, perspective (та же логика,
//                             что DefaultView() в vis_render.cc)
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

// Дословно тот же минимальный физ-лист, что в vis_render.cc — рендер не
// запускает событий (BeamOn не вызывается), физика нужна только чтобы
// Initialize() не упал без зарегистрированного физ-листа.
class MinimalPhysList : public G4VModularPhysicsList {
 public:
  MinimalPhysList() {
    RegisterPhysics(new G4EmStandardPhysics_option4());
    SetDefaultCutValue(0.05 * mm);
  }
};

// Один кадр на уже открытом TOOLSSG_OFFSCREEN-вьюере: сменить выходной
// файл, ракурс, проекцию и (опционально) режущую плоскость, перерисовать.
// sectionPlaneArgs пустая строка -> "/vis/viewer/set/sectionPlane off"
// (важно явно выключать между кадрами — вьюер держит состояние между
// вызовами /vis/viewer/rebuild, иначе разрез "утёк" бы в вид 4).
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
  // Нормаль разреза — единственный параметр командной строки, см. шапку
  // файла. По умолчанию "0 1 0" (первая попытка по брифу задачи).
  const std::string sectionNormal = (argc > 1) ? argv[1] : "0 1 0";

  auto* runManager =
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  runManager->SetVerboseLevel(0);

  runManager->SetUserInitialization(new RC110Detector());
  runManager->SetUserInitialization(new MinimalPhysList());
  runManager->Initialize();

  auto* visManager = new G4VisExecutive("quiet");
  visManager->Initialize();

  G4UImanager* ui = G4UImanager::GetUIpointer();
  ui->ApplyCommand("/control/verbose 0");
  ui->ApplyCommand("/run/verbose 0");
  ui->ApplyCommand("/vis/open TOOLSSG_OFFSCREEN 1600x900");

  // 1. Спереди: камера на -Z смотрит в +Z, видна лицевая грань (значок
  //    радиации, -Z = лицо прибора, см. GEANT4-MODEL.md "Система координат").
  RenderView(ui, "rc110_view_front.png", "0 0 -1", "0 1 0", "orthogonal", "");

  // 2. Сбоку: камера на -Y смотрит в +Y — самый информативный ракурс по
  //    Z-раскладке (кристалл/плата/АКБ разнесены по Z, не по Y).
  RenderView(ui, "rc110_view_side.png", "0 -1 0", "0 0 1", "orthogonal", "");

  // 3. Разрез Y=0: тот же ракурс, что вид 2, + настоящая режущая плоскость
  //    Geant4 (sectionPlane), не транспарентность. Нормаль — из argv,
  //    эмпирика по брифу задачи.
  RenderView(ui, "rc110_view_section.png", "0 -1 0", "0 0 1", "orthogonal",
             "on 0 0 0 mm " + sectionNormal);

  // 4. Три четверти: та же логика, что DefaultView() в vis_render.cc
  //    (viewpointThetaPhi 60 45 deg, perspective 30 deg), отдельный файл.
  //    Ставится не через RenderView() (тот ждёт viewpointVector, у "три
  //    четверти" — viewpointThetaPhi, другая команда) — sectionPlane
  //    выключаем явно, иначе унаследовался бы разрез от вида 3.
  ui->ApplyCommand("/vis/tsg/offscreen/set/file rc110_view_iso.png");
  ui->ApplyCommand("/vis/viewer/set/sectionPlane off");
  // upVector ЯВНО сброшен на дефолт Geant4 (0 1 0) — без этой строки вид
  // наследует "0 0 1", оставленный видом 3 (разрез), и камера получается
  // НЕ той же, что у оригинального standalone vis_render.cc (несмотря на
  // формально идентичные viewpointThetaPhi/projection) — поймано сверкой
  // с уже опубликованным rc110_geant4_overview.png (кристалл оказался
  // справа вместо слева).
  ui->ApplyCommand("/vis/viewer/set/upVector 0 1 0");
  ui->ApplyCommand("/vis/viewer/set/projection perspective 30 deg");
  ui->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 60 45 deg");
  ui->ApplyCommand("/vis/drawVolume");
  ui->ApplyCommand("/vis/viewer/rebuild");
  std::fprintf(stdout, "  wrote rc110_view_iso.png (viewpointThetaPhi 60 45 "
                        "deg, proj perspective 30 deg, section off)\n");

  std::fprintf(stdout, "EXITCODE=0 section_normal=%s\n",
               sectionNormal.c_str());

  delete visManager;
  delete runManager;
  return 0;
}
