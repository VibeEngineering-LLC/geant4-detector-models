// run_roomfield — шаг 1 двухшаговой схемы в РЕАЛЬНОЙ ГЕОМЕТРИИ КОМНАТЫ:
// спектр флюенса фотонов от естественных радионуклидов, равномерно
// распределённых в кирпиче стен и бетоне перекрытий, в точке, где стоит прибор.
// Заменяет сферическую модель results/wallion/wf_m1_*.csv (бетонная оболочка
// R=80 см, полость R=20 см), отвергнутую оператором 27.08.2026.
//
// ⚠ ПРИВАТНОСТЬ: реальные размеры комнаты оператора в исходники НЕ ЗАШИТЫ —
// только ключами CLI. Все дефолты нейтральные (4000x4000x2800 мм, стены 500 мм).
//
// Метод источника — МЕТОД 1 (решение D-001 контура): первичная частица — САМО
// ЯДРО, схему распада даёт RDM из ENSDF, nucleusLimits режет цепочку до одного
// звена. Таблица гамма-линий — устаревший метод, здесь НЕ используется.
//
// Программа НЕ читает и НЕ использует detectors/RadiaCode-103/geometry/
// {main.cc,RCDetector.cc,RCDetector.hh,PbShield.*,shieldrun.cc} — та линия
// запрещена оператором. Образец структуры и сборки — соседний run_field/;
// образец ионного источника — geometry/wallfield.cc (модели прибора не имеет).
//
// usage:
//   rc103_roomfield.exe <нуклид> <n_events> <out_csv> [ключи]
// ключи (все размеры в мм):
//   inner=<X>x<Y>x<Z>   внутренние размеры комнаты   (дефолт 4000x4000x2800)
//   wall_xm= wall_xp= wall_ym= wall_yp=   толщины кирпичных стен (дефолт 500)
//   floor= ceil=        толщины бетонных перекрытий  (дефолт 200)
//   obs=<dx>,<dy>,<dz>  точка прибора от внутренних граней X-, Y- и от пола
//   ball=<R>            радиус воздушного шара-скорера (дефолт 300)
//   cut=<мм>            порог продукции               (дефолт 1.0)
//   rho_brick=<г/см3>   плотность кирпича             (дефолт 1.8)
//   act=<Бк/кг>         удельная активность родителя  (дефолт 1.0)
#include "Rc103RoomFieldDetectorConstruction.hh"
#include "Rc103RoomFieldGeometry.hh"
#include "Rc103RoomFieldPhysicsList.hh"
#include "Rc103RoomFieldPrimaryGeneratorAction.hh"
#include "Rc103RoomFieldRunAction.hh"
#include "Rc103RoomFieldSteppingAction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

namespace {
// Нуклиды и их (Z, A) — порядок и состав по спеке.
const char* kNucNames[] = {"K40",   "Ra226", "Pb214", "Bi214",
                           "Pb212", "Ac228", "Bi212", "Tl208"};
const int kNucZ[] = {19, 88, 82, 83, 82, 89, 83, 81};
const int kNucA[] = {40, 226, 214, 214, 212, 228, 212, 208};
const int kNNuc = 8;

bool StartsWith(const char* s, const char* pref) {
  return std::strncmp(s, pref, std::strlen(pref)) == 0;
}
double ValueAfter(const char* s, const char* pref) {
  return std::atof(s + std::strlen(pref));
}
}  // namespace

#ifdef ROOMFIELD_BIRTH_STATS
void RoomFieldBirthStatsPrint();
#endif

int main(int argc, char** argv) {
  std::string positional[3];
  int nPositional = 0;
  double specificActivityBqPerKg = 1.0;

  for (int i = 1; i < argc; ++i) {
    const char* a = argv[i];
    if (StartsWith(a, "inner=")) {
      double x = 0, y = 0, z = 0;
      if (std::sscanf(a + 6, "%lfx%lfx%lf", &x, &y, &z) != 3 || x <= 0 ||
          y <= 0 || z <= 0) {
        std::fprintf(stderr,
                     "run_roomfield: FATAL не разобран ключ '%s'. Формат: "
                     "inner=4000x4000x2800 (мм, внутренние размеры)\n",
                     a);
        return 2;
      }
      gRoom.innerXMm = x;
      gRoom.innerYMm = y;
      gRoom.innerZMm = z;
    } else if (StartsWith(a, "wall_xm=")) {
      gRoom.wallXmMm = ValueAfter(a, "wall_xm=");
    } else if (StartsWith(a, "wall_xp=")) {
      gRoom.wallXpMm = ValueAfter(a, "wall_xp=");
    } else if (StartsWith(a, "wall_ym=")) {
      gRoom.wallYmMm = ValueAfter(a, "wall_ym=");
    } else if (StartsWith(a, "wall_yp=")) {
      gRoom.wallYpMm = ValueAfter(a, "wall_yp=");
    } else if (StartsWith(a, "floor=")) {
      gRoom.floorMm = ValueAfter(a, "floor=");
    } else if (StartsWith(a, "ceil=")) {
      gRoom.ceilMm = ValueAfter(a, "ceil=");
    } else if (StartsWith(a, "ball=")) {
      gRoom.ballRMm = ValueAfter(a, "ball=");
    } else if (StartsWith(a, "cut=")) {
      gRoom.cutMm = ValueAfter(a, "cut=");
    } else if (StartsWith(a, "rho_brick=")) {
      gRoom.rhoBrickGCm3 = ValueAfter(a, "rho_brick=");
    } else if (StartsWith(a, "rho_slab=")) {
      gRoom.rhoSlabGCm3 = ValueAfter(a, "rho_slab=");
    } else if (StartsWith(a, "extend=")) {
      gRoom.extendMm = ValueAfter(a, "extend=");
      if (gRoom.extendMm < 0) {
        std::fprintf(stderr, "run_roomfield: FATAL extend= не может быть "
                             "отрицательным.\n");
        return 2;
      }
    } else if (StartsWith(a, "src=")) {
      const char* v = a + 4;
      if (std::strcmp(v, "all") == 0) {
        gRoom.srcMode = RoomConfig::kSrcAll;
      } else if (std::strcmp(v, "brick") == 0) {
        gRoom.srcMode = RoomConfig::kSrcBrick;
      } else if (std::strcmp(v, "concrete") == 0) {
        gRoom.srcMode = RoomConfig::kSrcConcrete;
      } else {
        std::fprintf(stderr,
                     "run_roomfield: FATAL не разобран ключ '%s'. Формат: "
                     "src=all|brick|concrete\n",
                     a);
        return 2;
      }
    } else if (StartsWith(a, "act=")) {
      specificActivityBqPerKg = ValueAfter(a, "act=");
    } else if (StartsWith(a, "obs=")) {
      double dx = 0, dy = 0, dz = 0;
      if (std::sscanf(a + 4, "%lf,%lf,%lf", &dx, &dy, &dz) != 3) {
        std::fprintf(stderr,
                     "run_roomfield: FATAL не разобран ключ '%s'. Формат: "
                     "obs=1000,1500,750 (мм от X-, от Y-, от пола)\n",
                     a);
        return 2;
      }
      gRoom.obsDxMm = dx;
      gRoom.obsDyMm = dy;
      gRoom.obsDzMm = dz;
    } else if (nPositional < 3) {
      positional[nPositional++] = a;
    } else {
      std::fprintf(stderr, "run_roomfield: FATAL неизвестный ключ '%s'\n", a);
      return 2;
    }
  }

  if (nPositional < 3) {
    std::fprintf(stderr,
                 "usage: rc103_roomfield.exe <нуклид> <n_events> <out_csv> "
                 "[inner=XxYxZ] [wall_xm=] [wall_xp=] [wall_ym=] [wall_yp=] "
                 "[floor=] [ceil=] [obs=dx,dy,dz] [ball=] [cut=] "
                 "[rho_brick=] [act=]\n");
    std::fprintf(stderr, "нуклиды:");
    for (int i = 0; i < kNNuc; ++i) std::fprintf(stderr, " %s", kNucNames[i]);
    std::fprintf(stderr, "\n");
    return 2;
  }

  int nucIdx = -1;
  for (int i = 0; i < kNNuc; ++i) {
    if (positional[0] == kNucNames[i]) {
      nucIdx = i;
      break;
    }
  }
  if (nucIdx < 0) {
    std::fprintf(stderr, "run_roomfield: FATAL неизвестный нуклид '%s'\n",
                 positional[0].c_str());
    return 2;
  }
  const int ionZ = kNucZ[nucIdx];
  const int ionA = kNucA[nucIdx];

  const long long nEvents = std::atoll(positional[1].c_str());
  if (nEvents <= 0) {
    std::fprintf(stderr, "run_roomfield: FATAL n_events должно быть > 0 "
                         "(получено %lld)\n",
                 nEvents);
    return 2;
  }
  const std::string outCsv = positional[2];

  if (!(gRoom.cutMm > 0) || !(gRoom.ballRMm > 0) ||
      !(gRoom.rhoBrickGCm3 > 0) || !(specificActivityBqPerKg > 0)) {
    std::fprintf(stderr,
                 "run_roomfield: FATAL неположительное значение среди cut/ball/"
                 "rho_brick/act.\n");
    return 2;
  }

  // --- геометрия и нормировка ДО создания run manager ---------------------
  BuildRoomSlabs();
  auto* nist = G4NistManager::Instance();
  G4Material* concrete = nist->FindOrBuildMaterial("G4_CONCRETE");
  if (!concrete) {
    std::fprintf(stderr, "run_roomfield: FATAL материал G4_CONCRETE не "
                         "построен.\n");
    return 2;
  }
  const double rhoConcrete = (gRoom.rhoSlabGCm3 > 0.0)
                                 ? gRoom.rhoSlabGCm3
                                 : concrete->GetDensity() / (g / cm3);
  AssignRoomDensities(gRoom.rhoBrickGCm3, rhoConcrete);

  // S_v_i = A_уд * rho_i  [расп/(см3 с)],  A_уд в Бк/кг = 1e-3 Бк/г
  // R = sum_i S_v_i * V_i = A_уд[Бк/кг] * 1e-3 * M[г]
  // M — масса ТОЛЬКО того материала, который выбран ключом src=. Взять здесь
  // полную массу значило бы выдать раздельные шаблоны в разных единицах, и
  // складывать их было бы нельзя.
  const double massTotalG = RoomBrickMassG() + RoomConcreteMassG();
  const double massSourceG = RoomSelectedMassG();
  const double ratePerS = specificActivityBqPerKg * 1e-3 * massSourceG;
  const double tRunS = double(nEvents) / ratePerS;

  const double PI = 3.14159265358979323846;
  const double rBallCm = gRoom.ballRMm / 10.0;
  const double ballVolumeCm3 = 4.0 / 3.0 * PI * rBallCm * rBallCm * rBallCm;

  std::fprintf(stdout,
               "run_roomfield: нуклид=%s Z=%d A=%d n_events=%lld out=%s\n",
               kNucNames[nucIdx], ionZ, ionA, nEvents, outCsv.c_str());
  std::fprintf(stdout,
               "run_roomfield: комната %.0fx%.0fx%.0f мм, стены X-=%.0f X+=%.0f "
               "Y-=%.0f Y+=%.0f мм кирпич rho=%.3f, пол=%.0f потолок=%.0f мм "
               "бетон rho=%.3f\n",
               gRoom.innerXMm, gRoom.innerYMm, gRoom.innerZMm, gRoom.wallXmMm,
               gRoom.wallXpMm, gRoom.wallYmMm, gRoom.wallYpMm,
               gRoom.rhoBrickGCm3, gRoom.floorMm, gRoom.ceilMm, rhoConcrete);
  for (const auto& s : gSlabs) {
    std::fprintf(stdout,
                 "run_roomfield: плита %-16s %-8s V=%.6e см3  m=%.3f кг\n",
                 s.name.c_str(), s.brick ? "brick" : "concrete", s.volumeCm3,
                 s.massG / 1000.0);
  }
  std::fprintf(stdout,
               "run_roomfield: V_кирпич=%.6e см3 m=%.3f кг | V_бетон=%.6e см3 "
               "m=%.3f кг | m_всего=%.3f кг\n",
               RoomBrickVolumeCm3(), RoomBrickMassG() / 1000.0,
               RoomConcreteVolumeCm3(), RoomConcreteMassG() / 1000.0,
               massTotalG / 1000.0);
  std::fprintf(stdout,
               "run_roomfield: src=%s -> масса источника %.3f кг из %.3f кг; "
               "extend=%.0f мм; rho_перекрытий=%.3f г/см3\n",
               RoomSrcModeName(), massSourceG / 1000.0, massTotalG / 1000.0,
               gRoom.extendMm, rhoConcrete);
  std::fprintf(stdout,
               "run_roomfield: A_уд=%.6g Бк/кг -> R=%.6e расп/с, T=%.6e с, "
               "V_шара=%.6e см3, cut=%.4f мм\n",
               specificActivityBqPerKg, ratePerS, tRunS, ballVolumeCm3,
               gRoom.cutMm);

  auto* runManager =
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  runManager->SetVerboseLevel(0);
  runManager->SetUserInitialization(new Rc103RoomFieldDetectorConstruction());
  runManager->SetUserInitialization(
      new Rc103RoomFieldPhysicsList(gRoom.cutMm));

  auto* runAction = new Rc103RoomFieldRunAction(
      outCsv, kNucNames[nucIdx], ionZ, ionA, ratePerS, tRunS, ballVolumeCm3,
      specificActivityBqPerKg);
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(
      new Rc103RoomFieldPrimaryGeneratorAction(ionZ, ionA));
  runManager->SetUserAction(new Rc103RoomFieldSteppingAction(runAction));

  runManager->Initialize();

  // Приёмка геометрии — не рендер и не отсутствие overlap-WARNING (тот молчит
  // даже при объёме, целиком вылезшем за мать), а факт разрешённого объёма
  // ПОСЛЕ Initialize() плюс сам BeamOn ниже.
  if (!Rc103RoomFieldDetectorConstruction::GetBallLogicalVolume()) {
    std::fprintf(stderr,
                 "run_roomfield: FATAL объём шара-скорера не разрешён после "
                 "Initialize() — прерываюсь до BeamOn.\n");
    delete runManager;
    return 3;
  }

  G4UImanager* ui = G4UImanager::GetUIpointer();
  ui->ApplyCommand("/run/verbose 0");
  ui->ApplyCommand("/event/verbose 0");
  ui->ApplyCommand("/tracking/verbose 0");
  const long long progressEvery = (nEvents / 10 > 0) ? (nEvents / 10) : 1;
  ui->ApplyCommand("/run/printProgress " + std::to_string(progressEvery));

  // БЕЗ ЭТОЙ СТРОКИ ДОЛГОЖИВУЩИЕ ЯДРА НЕ РАСПАДАЮТСЯ ВООБЩЕ, а прогон при
  // этом завершается кодом 0 и пишет файл с нулевым флюенсом (так и случилось
  // 21.08 с K-40 и Ra-226). Порог RDM по умолчанию отсекает распады с очень
  // большим временем жизни. Сторож FLUENCE_ZERO в RunAction ловит этот случай.
  ui->ApplyCommand("/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+30 ns");

  // ОТСЕЧЕНИЕ ЦЕПОЧКИ — вторая обязательная строка метода 1. Без неё распад
  // идёт вниз по всей ветви, и шаблон звена вбирает своих потомков.
  char lim[128];
  std::snprintf(lim, sizeof(lim),
                "/process/had/rdm/nucleusLimits %d %d %d %d", ionA, ionA, ionZ,
                ionZ);
  ui->ApplyCommand(lim);
  std::fprintf(stdout,
               "CHAIN_CUT: nucleusLimits A=%d Z=%d (метод 1: распад ТОЛЬКО "
               "этого звена)\n",
               ionA, ionZ);

  runManager->BeamOn(static_cast<G4int>(nEvents));

#ifdef ROOMFIELD_BIRTH_STATS
  // Только диагностическая сборка (см. Rc103RoomFieldPrimaryGeneratorAction.cc):
  // печатает фактические доли рождений по плитам против ожидания по массе.
  RoomFieldBirthStatsPrint();
#endif

  delete runManager;
  std::fprintf(stdout, "EXITCODE=0\n");
  return 0;
}
