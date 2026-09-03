// RadiaCode-103 — шаг 2 двухшаговой схемы: отклик провалидированной GDML-модели
// прибора на ИЗОТРОПНОЕ поле фотонов с заданным спектром флюенса.
//
// Шаг 1 (поле в помещении) НЕ пересчитывается, берётся готовым из
// detectors/RadiaCode-103/results/wallion/wf_m1_<нуклид>.csv.
//
// Программа НЕ читает и НЕ использует detectors/RadiaCode-103/geometry/
// {main.cc,RCDetector.cc,RCDetector.hh} и любые результаты, посчитанные той
// моделью (results/bare/, results/m200/) — запрет оператора от 27.08.2026.
// Образец кода — собственный комплект run/ этой же сессии.
//
// usage: rc103_field.exe <flux_csv> <n_events> <out_csv> [--check-norm]
#include "Rc103FieldDetectorConstruction.hh"
#include "Rc103FieldEventAction.hh"
#include "Rc103FieldFluxSpectrum.hh"
#include "Rc103FieldPhysicsList.hh"
#include "Rc103FieldPrimaryGeneratorAction.hh"
#include "Rc103FieldRunAction.hh"
#include "Rc103FieldSteppingAction.hh"

#include "G4RunManagerFactory.hh"
#include "Randomize.hh"
#include "G4UImanager.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#ifndef RC103_FIELD_GDML_DEFAULT_PATH
#error "RC103_FIELD_GDML_DEFAULT_PATH must be defined by CMake (see run_field/CMakeLists.txt)"
#endif

int main(int argc, char** argv) {
  bool checkNorm = false;
  // seed=<N> — зерно ГСЧ. Нужен затем же, зачем в run_muon: прогоны длиннее
  // ~2 часов на этой машине не доживают (W-040), длинная статистика набирается
  // серией коротких порций, а складывать их можно ТОЛЬКО при разных зёрнах.
  long seed = 0;  // 0 = не трогать, оставить дефолтное зерно Geant4
  // emcut=/deex= — порог продукции и режим атомной деэкситации (03.09.2026).
  // ДЕФОЛТ deex=max по правилу #CFG-2 (оператор 03.09.2026: «счёт идёт всегда по
  // максимальным настройкам, если не сказано иначе»). Прогоны ДО 03.09.2026
  // воспроизводятся только явным deex=std — молча они больше не повторяются.
  // emcut=0.05 мм — единственное отступление от «максимума», оплачено
  // измерением: 1 мм→0,05 мм меняет полосы не больше чем на 1,7 %.
  double emCutMm = 0.05;
  std::string deexMode = "max";
  std::string positional[3];
  int nPositional = 0;
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "--check-norm") == 0) {
      checkNorm = true;
    } else if (std::strncmp(argv[i], "room=", 5) == 0) {
      // room=<X>x<Y>x<Z> — ПОЛНЫЕ габариты комнаты в мм.
      double rx = 0, ry = 0, rz = 0;
      if (std::sscanf(argv[i] + 5, "%lfx%lfx%lf", &rx, &ry, &rz) != 3 ||
          rx <= 0 || ry <= 0 || rz <= 0) {
        std::fprintf(stderr,
                     "rc103_field: FATAL не разобран ключ '%s'. "
                     "Формат: room=4000x4000x2800 (мм, полные габариты)\n",
                     argv[i]);
        return 2;
      }
      Rc103FieldDetectorConstruction::gRoomHalfXMm = 0.5 * rx;
      Rc103FieldDetectorConstruction::gRoomHalfYMm = 0.5 * ry;
      Rc103FieldDetectorConstruction::gRoomHalfZMm = 0.5 * rz;
    } else if (std::strncmp(argv[i], "shield=", 7) == 0) {
      // shield=on | shield=off — свинцовый домик вокруг прибора, дефолт off.
      const char* v = argv[i] + 7;
      if (std::strcmp(v, "on") == 0) {
        Rc103FieldDetectorConstruction::gShieldOn = true;
      } else if (std::strcmp(v, "off") == 0) {
        Rc103FieldDetectorConstruction::gShieldOn = false;
      } else {
        std::fprintf(stderr,
                     "rc103_field: FATAL не разобран ключ '%s'. "
                     "Формат: shield=on либо shield=off\n",
                     argv[i]);
        return 2;
      }
    } else if (std::strncmp(argv[i], "seed=", 5) == 0) {
      seed = std::atol(argv[i] + 5);
    } else if (std::strncmp(argv[i], "stand=", 6) == 0) {
      // stand=<мм> — высота опоры над дном полости (реально картон 25 мм);
      // stand=asbuilt — прежнее допущение: центр габарита прибора в (0,0,0).
      const char* v = argv[i] + 6;
      if (std::strcmp(v, "asbuilt") == 0) {
        Rc103FieldDetectorConstruction::gStandMm = -1.0;
      } else {
        char* end = nullptr;
        const double s = std::strtod(v, &end);
        if (end == v || *end != '\0' || s < 0.0) {
          std::fprintf(stderr,
                       "rc103_field: FATAL не разобран ключ '%s'. Формат: "
                       "stand=<мм >= 0> либо stand=asbuilt\n",
                       argv[i]);
          return 2;
        }
        Rc103FieldDetectorConstruction::gStandMm = s;
      }
    } else if (std::strncmp(argv[i], "flip=", 5) == 0) {
      // flip=up — экран вверх (реальная постановка); flip=down — как до 30.08.
      const char* v = argv[i] + 5;
      if (std::strcmp(v, "up") == 0) {
        Rc103FieldDetectorConstruction::gFlipUp = true;
      } else if (std::strcmp(v, "down") == 0) {
        Rc103FieldDetectorConstruction::gFlipUp = false;
      } else {
        std::fprintf(stderr,
                     "rc103_field: FATAL не разобран ключ '%s'. "
                     "Формат: flip=up либо flip=down\n",
                     argv[i]);
        return 2;
      }
    } else if (std::strncmp(argv[i], "emcut=", 6) == 0) {
      // emcut=<мм> — порог продукции вторичных, тот же для gamma/e-/e+.
      char* end = nullptr;
      const double c = std::strtod(argv[i] + 6, &end);
      if (end == argv[i] + 6 || *end != '\0' || !(c > 0.0)) {
        std::fprintf(stderr,
                     "rc103_field: FATAL не разобран ключ '%s'. "
                     "Формат: emcut=<мм больше 0>\n",
                     argv[i]);
        return 2;
      }
      emCutMm = c;
    } else if (std::strncmp(argv[i], "deex=", 5) == 0) {
      // deex=std — как было; deex=deex — плюс Оже и DeexcitationIgnoreCut;
      // deex=max — плюс PIXE. Неизвестное значение — отказ, не молчаливый std.
      const char* v = argv[i] + 5;
      if (std::strcmp(v, "std") == 0 || std::strcmp(v, "deex") == 0 ||
          std::strcmp(v, "max") == 0) {
        deexMode = v;
      } else {
        std::fprintf(stderr,
                     "rc103_field: FATAL не разобран ключ '%s'. "
                     "Формат: deex=std|deex|max\n",
                     argv[i]);
        return 2;
      }
    } else if (nPositional < 3) {
      positional[nPositional++] = argv[i];
    }
  }
  if (nPositional < 3) {
    std::fprintf(stderr,
                 "usage: rc103_field.exe <flux_csv> <n_events> <out_csv> "
                 "[room=<X>x<Y>x<Z>] [shield=on|off] [stand=<мм>|asbuilt] "
                 "[flip=up|down] [emcut=<мм>] [deex=std|deex|max] "
                 "[--check-norm]\n");
    return 2;
  }

  const std::string fluxCsv = positional[0];
  const long long nEvents = std::atoll(positional[1].c_str());
  const std::string outCsv = positional[2];
  const std::string gdmlPath = RC103_FIELD_GDML_DEFAULT_PATH;

  if (nEvents <= 0) {
    std::fprintf(stderr, "rc103_field: FATAL n_events must be > 0 (got %lld)\n",
                 nEvents);
    return 2;
  }

  Rc103FieldFluxSpectrum spectrum;
  if (!spectrum.Load(fluxCsv)) {
    std::fprintf(stderr, "rc103_field: FATAL could not load flux CSV.\n");
    return 2;
  }

  // --- Нормировка: тождество Ф = 4N/S -------------------------------------
  const double R_SRC_CM = Rc103FieldPrimaryGeneratorAction::kRSrcCm;  // 7.0 см
  const double PI = 3.14159265358979323846;
  const double surfaceCm2 = 4.0 * PI * R_SRC_CM * R_SRC_CM;
  const double fluxTotal = spectrum.HeaderTotalCm2S();
  const double ratePerS = fluxTotal * surfaceCm2 / 4.0;
  const double tRunS = (ratePerS > 0.0) ? double(nEvents) / ratePerS : 0.0;

  // Объём контрольного шара R=2.0 см: V = 4/3*pi*R^3 = 33.5103 см^3
  const double rCheckCm = Rc103FieldDetectorConstruction::kCheckRadiusMm / 10.0;
  const double checkVolumeCm3 = 4.0 / 3.0 * PI * rCheckCm * rCheckCm * rCheckCm;

  std::fprintf(stdout,
               "rc103_field: flux_csv=%s n_events=%lld out_csv=%s "
               "check_norm=%d\n",
               fluxCsv.c_str(), nEvents, outCsv.c_str(), checkNorm ? 1 : 0);
  std::fprintf(stdout,
               "rc103_field: R_SRC=%.1f cm surface=%.6f cm2 "
               "flux_total=%.6e 1/(cm2 s) rate=%.6e 1/s t_run=%.6e s\n",
               R_SRC_CM, surfaceCm2, fluxTotal, ratePerS, tRunS);
  if (checkNorm) {
    std::fprintf(stdout,
                 "rc103_field: CHECK-NORM check_ball R=%.1f cm V=%.6f cm3\n",
                 rCheckCm, checkVolumeCm3);
  } else {
    std::fprintf(stdout, "rc103_field: gdml=%s\n", gdmlPath.c_str());
  }

  auto* runManager =
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  runManager->SetVerboseLevel(0);
  if (seed != 0) {
    G4Random::setTheSeed(seed);
    std::fprintf(stdout, "rc103_field: seed=%ld\n", seed);
  }

  runManager->SetUserInitialization(
      new Rc103FieldDetectorConstruction(gdmlPath, checkNorm));
  runManager->SetUserInitialization(
      new Rc103FieldPhysicsList(emCutMm, deexMode));

  auto* runAction =
      new Rc103FieldRunAction(outCsv, fluxTotal, surfaceCm2, ratePerS, tRunS,
                              checkNorm, checkVolumeCm3);
  auto* eventAction = new Rc103FieldEventAction(runAction);
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(new Rc103FieldPrimaryGeneratorAction(&spectrum));
  runManager->SetUserAction(
      new Rc103FieldSteppingAction(eventAction, runAction));

  runManager->Initialize();

  // Приёмка геометрии — не рендер и не отсутствие overlap-warning, а факт
  // разрешения нужного объёма ПОСЛЕ Initialize() плюс сам BeamOn ниже.
  if (checkNorm) {
    if (!Rc103FieldDetectorConstruction::GetCheckLogicalVolume()) {
      std::fprintf(stderr,
                   "rc103_field: FATAL check ball logical volume not resolved "
                   "after Initialize() - aborting before BeamOn.\n");
      delete runManager;
      return 3;
    }
  } else {
    if (!Rc103FieldDetectorConstruction::GetCrystalLogicalVolume()) {
      std::fprintf(stderr,
                   "rc103_field: FATAL crystal logical volume not resolved "
                   "after Initialize() - aborting before BeamOn.\n");
      delete runManager;
      return 3;
    }
    // Домик приёмится тем же способом: не «рендер выглядит правильно», а факт
    // разрешённого объёма после Initialize() плюс сам BeamOn ниже.
    if (Rc103FieldDetectorConstruction::gShieldOn &&
        !Rc103FieldDetectorConstruction::GetShieldLogicalVolume()) {
      std::fprintf(stderr,
                   "rc103_field: FATAL shield=on, но объём свинцового домика не "
                   "разрешён после Initialize() - aborting before BeamOn.\n");
      delete runManager;
      return 3;
    }
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
