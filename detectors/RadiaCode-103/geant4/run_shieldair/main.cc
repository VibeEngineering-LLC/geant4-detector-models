// RadiaCode-103 - объёмный источник радона (Pb-214/Bi-214) в воздухе полости
// свинцового домика. Гипотеза 30.08.2026 (P-005): широкополосный недобор в
// мягкой части якоря в реальной посадке может объясняться неучтённым радоном
// в воздухе полости, а не только геометрией посадки. Проверка сканированием
// амплитуды мюонов показала устойчивый сигнал именно в сторону Ra226_chain
// (Pb-214/Bi-214) - см. scan_mu_fit.py и project-incidents.md P-005.
//
// Метод: та же геометрия домика (Rc103FieldDetectorConstruction, посадка
// прибора реальная), но источник - НЕ изотропная сфера снаружи, а объёмный
// ионный источник, равномерно распределённый по полости (метод 1: /gps/
// particle ion + RDM, nucleusLimits на одно звено - тот же приём, что во
// всех прочих шаблонах). Выход - тот же формат CSV (metric,value + bin_keV,
// counts,cps), что и у прочих шаблонов, ПОЭТОМУ ЧИТАЕТСЯ ftc.read_template()
// без изменений: cps трактуется как "на 1 Бк/м3 в воздухе полости", а не
// "на 1 Бк/кг родителя" - разница нормировки не в коде анализа, а в том,
// как вызывающий код интерпретирует cps на выходе.
//
// usage: rc103_shieldair.exe <nuclide> <n_events> <out_csv> [stand=<мм>|asbuilt]
//        [flip=up|down] [seed=<N>]
#include "Rc103FieldDetectorConstruction.hh"
#include "Rc103FieldEventAction.hh"
#include "Rc103FieldRunAction.hh"
#include "Rc103ShieldAirPhysicsList.hh"
#include "Rc103FieldSteppingAction.hh"
#include "Rc103ShieldAirPrimaryGeneratorAction.hh"

#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#ifndef RC103_FIELD_GDML_DEFAULT_PATH
#error "RC103_FIELD_GDML_DEFAULT_PATH must be defined by CMake"
#endif

namespace {
const char* kNucNames[] = {"Pb214", "Bi214"};
const int kNucZ[] = {82, 83};
const int kNucA[] = {214, 214};
const int kNNuc = 2;
}  // namespace

int main(int argc, char** argv) {
  if (argc < 4) {
    std::fprintf(stderr,
                 "usage: rc103_shieldair.exe <Pb214|Bi214> <n_events> "
                 "<out_csv> [stand=<mm>|asbuilt] [flip=up|down] [seed=<N>]\n");
    return 2;
  }
  const std::string nuc = argv[1];
  int nucIdx = -1;
  for (int i = 0; i < kNNuc; ++i) if (nuc == kNucNames[i]) nucIdx = i;
  if (nucIdx < 0) {
    std::fprintf(stderr, "rc103_shieldair: FATAL нуклид '%s' не поддержан "
                 "(только Pb214, Bi214 - гамма-активные дочери радона).\n",
                 nuc.c_str());
    return 2;
  }
  const long long nEvents = std::atoll(argv[2]);
  const std::string outCsv = argv[3];
  long seed = 0;
  for (int i = 4; i < argc; ++i) {
    if (std::strncmp(argv[i], "stand=", 6) == 0) {
      const char* v = argv[i] + 6;
      Rc103FieldDetectorConstruction::gStandMm =
          (std::strcmp(v, "asbuilt") == 0) ? -1.0 : std::strtod(v, nullptr);
    } else if (std::strncmp(argv[i], "flip=", 5) == 0) {
      Rc103FieldDetectorConstruction::gFlipUp =
          (std::strcmp(argv[i] + 5, "up") == 0);
    } else if (std::strncmp(argv[i], "seed=", 5) == 0) {
      seed = std::atol(argv[i] + 5);
    }
  }
  Rc103FieldDetectorConstruction::gShieldOn = true;  // источник ТОЛЬКО в домике

  // Объём полости: тот же бокс, что в PrimaryGeneratorAction. Одна плоскость
  // 1 Бк/м3: decay_rate = 1 * V_air[м3] распадов/с; T_run = N/decay_rate.
  using DC = Rc103FieldDetectorConstruction;
  const double cavXCm = DC::kShieldCavityXMm / 10.0;
  const double cavYCm = DC::kShieldCavityYMm / 10.0;
  const double cavZCm = (DC::kShieldOuterZMm - DC::kShieldPbMm) / 10.0;
  const double airVolumeM3 = (cavXCm * cavYCm * cavZCm) * 1.0e-6;
  const double tRunS = airVolumeM3 > 0.0 ? double(nEvents) / airVolumeM3 : 0.0;
  std::fprintf(stdout,
               "rc103_shieldair: nuclide=%s Z=%d A=%d n=%lld cav_air_m3=%.6e "
               "t_run_s(при 1 Бк/м3)=%.6e\n",
               kNucNames[nucIdx], kNucZ[nucIdx], kNucA[nucIdx], nEvents,
               airVolumeM3, tRunS);

  const std::string gdmlPath = RC103_FIELD_GDML_DEFAULT_PATH;
  auto* runManager =
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  runManager->SetVerboseLevel(0);
  if (seed != 0) G4Random::setTheSeed(seed);

  runManager->SetUserInitialization(
      new Rc103FieldDetectorConstruction(gdmlPath, false));
  runManager->SetUserInitialization(new Rc103ShieldAirPhysicsList());

  auto* runAction = new Rc103FieldRunAction(
      outCsv, 0.0, 0.0, airVolumeM3, tRunS, false, 0.0);
  auto* eventAction = new Rc103FieldEventAction(runAction);
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(new Rc103ShieldAirPrimaryGeneratorAction(
      kNucZ[nucIdx], kNucA[nucIdx]));
  runManager->SetUserAction(
      new Rc103FieldSteppingAction(eventAction, runAction));

  runManager->Initialize();
  if (!Rc103FieldDetectorConstruction::GetCrystalLogicalVolume()) {
    std::fprintf(stderr, "rc103_shieldair: FATAL Crystal_log не разрешён "
                 "после Initialize().\n");
    delete runManager;
    return 3;
  }

  auto* ui = G4UImanager::GetUIpointer();
  ui->ApplyCommand("/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+30 ns");
  char lim[128];
  std::snprintf(lim, sizeof(lim), "/process/had/rdm/nucleusLimits %d %d %d %d",
               kNucA[nucIdx], kNucA[nucIdx], kNucZ[nucIdx], kNucZ[nucIdx]);
  ui->ApplyCommand(lim);
  const long long progressEvery = (nEvents / 10 > 0) ? (nEvents / 10) : 1;
  ui->ApplyCommand("/run/printProgress " + std::to_string(progressEvery));

  runManager->BeamOn(static_cast<G4int>(nEvents));
  delete runManager;
  return 0;
}