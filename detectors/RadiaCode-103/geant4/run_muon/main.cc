// RadiaCode-103 — мюонный столбец для нуклидного разложения фона: отклик
// провалидированной GDML-модели прибора на космические мюоны.
//
// Программа НЕ читает и НЕ использует detectors/RadiaCode-103/geometry/
// {main.cc,RCDetector.cc,RCDetector.hh,cosmicmu.cc} и любые результаты,
// посчитанные той моделью (в т.ч. build/RadiaCode-103/cosmicmu.csv) — запрет
// оператора от 27.08.2026. Образец кода — соседний комплект run_field/ этой же
// сессии; физика входного спектра взята из спеки _spec_run_muon.md как ДАННЫЕ.
//
// usage: rc103_muon.exe <n_events> <out_csv> [rdisk=<мм>] [seed=<N>]
//
// seed= нужен затем, что проверка насыщения по радиусу сравнивает ДВА прогона;
// с общим сидом их выборки коррелированы, и расхождение форм нельзя честно
// сопоставить с пуассоновской ошибкой. Разные сиды делают прогоны независимыми.
#include "Rc103MuonDetectorConstruction.hh"
#include "Rc103MuonEventAction.hh"
#include "Rc103MuonPhysicsList.hh"
#include "Rc103MuonPrimaryGeneratorAction.hh"
#include "Rc103MuonRunAction.hh"
#include "Rc103MuonSpectrum.hh"
#include "Rc103MuonSteppingAction.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#ifndef RC103_MUON_GDML_DEFAULT_PATH
#error "RC103_MUON_GDML_DEFAULT_PATH must be defined by CMake (see run_muon/CMakeLists.txt)"
#endif

int main(int argc, char** argv) {
  double rDiskMm = Rc103MuonPrimaryGeneratorAction::kRDiskDefaultMm;
  double zDiskMm = Rc103MuonPrimaryGeneratorAction::kZDiskMm;
  bool shieldOn = false;
  long seed = 0;  // 0 = не трогать, оставить дефолтный сид Geant4
  std::string positional[2];
  int nPositional = 0;
  for (int i = 1; i < argc; ++i) {
    if (std::strncmp(argv[i], "rdisk=", 6) == 0) {
      rDiskMm = std::atof(argv[i] + 6);
    } else if (std::strncmp(argv[i], "zdisk=", 6) == 0) {
      zDiskMm = std::atof(argv[i] + 6);
    } else if (std::strncmp(argv[i], "shield=", 7) == 0) {
      const char* v = argv[i] + 7;
      if (std::strcmp(v, "on") == 0) {
        shieldOn = true;
      } else if (std::strcmp(v, "off") == 0) {
        shieldOn = false;
      } else {
        std::fprintf(stderr,
                     "rc103_muon: FATAL непонятный аргумент '%s'. "
                     "Формат: shield=on либо shield=off\n",
                     argv[i]);
        return 2;
      }
    } else if (std::strncmp(argv[i], "stand=", 6) == 0) {
      // stand=<мм> — опора над дном полости (реально картон 25 мм);
      // stand=asbuilt — прежнее допущение: центр габарита в (0,0,0). P-005.
      const char* v = argv[i] + 6;
      if (std::strcmp(v, "asbuilt") == 0) {
        Rc103MuonDetectorConstruction::gStandMm = -1.0;
      } else {
        char* end = nullptr;
        const double sv = std::strtod(v, &end);
        if (end == v || *end != 0 || sv < 0.0) {
          std::fprintf(stderr,
                       "rc103_muon: FATAL не разобран ключ '%s'. Формат: "
                       "stand=<мм >= 0> либо stand=asbuilt\n",
                       argv[i]);
          return 2;
        }
        Rc103MuonDetectorConstruction::gStandMm = sv;
      }
    } else if (std::strncmp(argv[i], "flip=", 5) == 0) {
      const char* v = argv[i] + 5;
      if (std::strcmp(v, "up") == 0) {
        Rc103MuonDetectorConstruction::gFlipUp = true;
      } else if (std::strcmp(v, "down") == 0) {
        Rc103MuonDetectorConstruction::gFlipUp = false;
      } else {
        std::fprintf(stderr,
                     "rc103_muon: FATAL не разобран ключ '%s'. "
                     "Формат: flip=up либо flip=down\n",
                     argv[i]);
        return 2;
      }
    } else if (std::strncmp(argv[i], "seed=", 5) == 0) {
      seed = std::atol(argv[i] + 5);
    } else if (nPositional < 2) {
      positional[nPositional++] = argv[i];
    }
  }
  if (nPositional < 2) {
    std::fprintf(stderr,
                 "usage: rc103_muon.exe <n_events> <out_csv> [rdisk=<mm>] "
                 "[zdisk=<mm>] [shield=on|off] [stand=<мм>|asbuilt] "
                 "[flip=up|down] [seed=<N>]\n");
    return 2;
  }

  const long long nEvents = std::atoll(positional[0].c_str());
  const std::string outCsv = positional[1];
  const std::string gdmlPath = RC103_MUON_GDML_DEFAULT_PATH;

  if (nEvents <= 0) {
    std::fprintf(stderr, "rc103_muon: FATAL n_events must be > 0 (got %lld)\n",
                 nEvents);
    return 2;
  }
  if (!(rDiskMm > 0.0)) {
    std::fprintf(stderr, "rc103_muon: FATAL rdisk must be > 0 (got %g)\n",
                 rDiskMm);
    return 2;
  }

  if (!(zDiskMm > 0.0)) {
    std::fprintf(stderr, "rc103_muon: FATAL zdisk must be > 0 (got %g)\n",
                 zDiskMm);
    return 2;
  }

  const double PI = 3.14159265358979323846;
  const double rDiskCm = rDiskMm / 10.0;
  const double diskAreaCm2 = PI * rDiskCm * rDiskCm;
  // Справочная проверка порядка величины: интегральный поток мюонов на уровне
  // моря ~0.0167 см^-2 c^-1 через горизонтальную площадку. Печатается и
  // пишется в CSV, но НЕ используется как ограничение — отклик нормируется на
  // один мюон, пересёкший диск.
  const double kMuonFluxCm2S = 0.0167;
  const double pdgExpectedPerS = kMuonFluxCm2S * diskAreaCm2;

  // Мир обязан вмещать диск: спека фиксирует 400 мм под rdisk=150..300, но
  // допускает третий прогон при rdisk=600 — тогда 400 мм мало. С 29.08.2026
  // мир обязан вмещать ещё и поднятый диск старта, и наружный габарит домика
  // (полувысота 217,5 мм) — иначе постановка развалится молча.
  const double worldHalfMm = std::max(
      {Rc103MuonDetectorConstruction::kWorldHalfMmDefault, rDiskMm + 100.0,
       zDiskMm + 100.0, 0.5 * Rc103MuonDetectorConstruction::kShieldOuterZMm + 100.0});

  std::fprintf(stdout,
               "rc103_muon: n_events=%lld out_csv=%s rdisk=%.1f mm zdisk=%.1f "
               "mm shield=%s world_half=%.1f mm\n",
               nEvents, outCsv.c_str(), rDiskMm, zDiskMm,
               shieldOn ? "on" : "off", worldHalfMm);
  std::fprintf(stdout,
               "rc103_muon: disk_area=%.6f cm2 pdg_expected=%.6e 1/s "
               "(reference only, NOT a constraint)\n",
               diskAreaCm2, pdgExpectedPerS);
  std::fprintf(stdout, "rc103_muon: gdml=%s\n", gdmlPath.c_str());

  if (seed != 0) {
    G4Random::setTheSeed(seed);
    std::fprintf(stdout, "rc103_muon: RNG seed set explicitly to %ld\n", seed);
  } else {
    std::fprintf(stdout, "rc103_muon: RNG seed = Geant4 default (no seed= arg)\n");
  }

  Rc103MuonSpectrum spectrum;

  auto* runManager =
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  runManager->SetVerboseLevel(0);

  runManager->SetUserInitialization(
      new Rc103MuonDetectorConstruction(gdmlPath, worldHalfMm, shieldOn,
                                        zDiskMm));
  runManager->SetUserInitialization(new Rc103MuonPhysicsList());

  auto* runAction = new Rc103MuonRunAction(
      outCsv, rDiskMm, zDiskMm, diskAreaCm2, pdgExpectedPerS,
      Rc103MuonSpectrum::kELoGeV, Rc103MuonSpectrum::kEHiGeV);
  auto* eventAction = new Rc103MuonEventAction(runAction);
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(
      new Rc103MuonPrimaryGeneratorAction(&spectrum, rDiskMm, zDiskMm));
  runManager->SetUserAction(new Rc103MuonSteppingAction(eventAction));

  runManager->Initialize();

  // Приёмка геометрии — не рендер и не отсутствие overlap-warning, а факт
  // разрешения нужного объёма ПОСЛЕ Initialize() плюс сам BeamOn ниже.
  if (!Rc103MuonDetectorConstruction::GetCrystalLogicalVolume()) {
    std::fprintf(stderr,
                 "rc103_muon: FATAL crystal logical volume not resolved after "
                 "Initialize() - aborting before BeamOn.\n");
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
