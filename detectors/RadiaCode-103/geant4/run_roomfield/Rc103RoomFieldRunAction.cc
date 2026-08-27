#include "Rc103RoomFieldRunAction.hh"

#include "Rc103RoomFieldGeometry.hh"

#include "G4Run.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>

Rc103RoomFieldRunAction::Rc103RoomFieldRunAction(
    std::string outputCsv, std::string nuclide, int ionZ, int ionA,
    double ratePerS, double tRunS, double ballVolumeCm3,
    double specificActivityBqPerKg)
    : fOutputCsv(std::move(outputCsv)),
      fNuclide(std::move(nuclide)),
      fIonZ(ionZ),
      fIonA(ionA),
      fRatePerS(ratePerS),
      fTRunS(tRunS),
      fBallVolumeCm3(ballVolumeCm3),
      fSpecificActivityBqPerKg(specificActivityBqPerKg),
      fLenCm(static_cast<std::size_t>(kBins) + 1, 0.0) {}

void Rc103RoomFieldRunAction::BeginOfRunAction(const G4Run*) {
  std::fill(fLenCm.begin(), fLenCm.end(), 0.0);
}

void Rc103RoomFieldRunAction::AddTrackLengthCm(double eKeV, double lenCm) {
  int b = static_cast<int>(eKeV / kBinKeV);
  if (b > kBins) b = kBins;  // канал переполнения (E >= 3000 кэВ)
  if (b >= 0) fLenCm[static_cast<std::size_t>(b)] += lenCm;
}

void Rc103RoomFieldRunAction::EndOfRunAction(const G4Run* run) {
  if (!IsMaster()) return;

  const long long nEvents = run->GetNumberOfEvent();
  if (nEvents <= 0) {
    std::fprintf(stderr,
                 "Rc103RoomFieldRunAction: FATAL число событий = %lld.\n",
                 nEvents);
    std::exit(4);
  }

  // Эквивалентное время считается по ФАКТИЧЕСКОМУ числу событий, а не по
  // запрошенному: если прогон оборвали, деление на запрошенное N молча
  // занизило бы флюенс.
  const double tRunS = double(nEvents) / fRatePerS;
  if (fTRunS > 0.0 && std::fabs(tRunS - fTRunS) / fTRunS > 1e-9) {
    std::fprintf(stdout,
                 "Rc103RoomFieldRunAction: ВНИМАНИЕ t_run пересчитан по факту: "
                 "%.6e с (до прогона было %.6e с, N=%lld).\n",
                 tRunS, fTRunS, nEvents);
  }

  const double norm = 1.0 / (fBallVolumeCm3 * tRunS);  // см -> 1/(см2 с)
  double total = 0.0;
  for (double l : fLenCm) total += l * norm;

  namespace fs = std::filesystem;
  fs::path outPath(fOutputCsv);
  if (outPath.has_parent_path()) {
    std::error_code ec;
    fs::create_directories(outPath.parent_path(), ec);
  }
  FILE* f = std::fopen(fOutputCsv.c_str(), "w");
  if (!f) {
    std::fprintf(stderr, "Rc103RoomFieldRunAction: FATAL не открыть '%s'.\n",
                 fOutputCsv.c_str());
    std::exit(4);
  }

  std::fprintf(f, "# поле ЕРН в РЕАЛЬНОЙ КОМНАТЕ: флюенс фотонов в воздушном "
                  "шаре в точке прибора\n");
  std::fprintf(f, "# program = run_roomfield (метод 1: ионный источник + RDM, "
                  "nucleusLimits на одно звено)\n");
  std::fprintf(f, "# build = %s %s\n", __DATE__, __TIME__);
  std::fprintf(f, "# nuclide = %s (Z=%d A=%d), %.6g Bq/kg parent\n",
               fNuclide.c_str(), fIonZ, fIonA, fSpecificActivityBqPerKg);
  std::fprintf(f, "# room_inner_mm = %.1f x %.1f x %.1f\n", gRoom.innerXMm,
               gRoom.innerYMm, gRoom.innerZMm);
  std::fprintf(f,
               "# walls_mm: X-=%.1f X+=%.1f Y-=%.1f Y+=%.1f (кирпич), "
               "floor=%.1f ceiling=%.1f (перекрытия)\n",
               gRoom.wallXmMm, gRoom.wallXpMm, gRoom.wallYmMm, gRoom.wallYpMm,
               gRoom.floorMm, gRoom.ceilMm);
  std::fprintf(f, "# src_material = %s (розыгрыш точки рождения ТОЛЬКО в этом "
                  "материале; нормировка по его массе)\n",
               RoomSrcModeName());
  std::fprintf(f,
               "# extend_mm = %.1f (слой ТОГО ЖЕ материала снаружи каждой "
               "стены и каждого перекрытия; 0 = изолированная комната)\n",
               gRoom.extendMm);
  std::fprintf(f,
               "# eff_thickness_mm: X-=%.1f X+=%.1f Y-=%.1f Y+=%.1f "
               "floor=%.1f ceiling=%.1f\n",
               RoomEffWallXmMm(), RoomEffWallXpMm(), RoomEffWallYmMm(),
               RoomEffWallYpMm(), RoomEffFloorMm(), RoomEffCeilMm());
  if (gRoom.rhoSlabGCm3 > 0.0) {
    std::fprintf(f,
                 "# APPROXIMATION: перекрытия — железобетонные "
                 "КРУГЛОПУСТОТНЫЕ плиты, заданы ГОМОГЕНИЗИРОВАННЫМ слоем: "
                 "состав G4_CONCRETE при ЭФФЕКТИВНОЙ плотности %.3f г/см3 "
                 "(паспортная масса / габаритный объём). Цилиндрические "
                 "пустоты явно НЕ моделируются: излучение ВДОЛЬ оси пустот "
                 "гомогенная модель ослабляет сильнее реального, но телесный "
                 "угол таких направлений мал. Поверхностная плотность "
                 "перекрытия = %.1f г/см2.\n",
                 gRoom.rhoSlabGCm3, gRoom.rhoSlabGCm3 * gRoom.floorMm / 10.0);
  } else {
    std::fprintf(f,
                 "# ASSUMPTION: перекрытия — сплошной G4_CONCRETE плотности по "
                 "умолчанию, толщина не подтверждена оператором\n");
  }
  std::fprintf(f,
               "# APPROXIMATION: кирпич задан ЭЛЕМЕНТНЫМ СОСТАВОМ G4_CONCRETE "
               "при плотности %.3f г/см3 (NIST-состава керамического кирпича в "
               "Geant4 нет). Основание: в 100..3000 кэВ доминирует комптон, "
               "сечение задаётся электронной плотностью, Z/A силикатов бетона и "
               "керамики практически совпадает; расхождение существенно только "
               "ниже ~100 кэВ, где работает фотоэффект.\n",
               gRoom.rhoBrickGCm3);
  std::fprintf(f, "# obs_point_mm_from_inner_faces = %.1f, %.1f, %.1f (от X-, "
                  "от Y-, от пола)\n",
               gRoom.obsDxMm, gRoom.obsDyMm, gRoom.obsDzMm);
  std::fprintf(f, "# obs_point_world_mm = %.1f, %.1f, %.1f\n", RoomObsXMm(),
               RoomObsYMm(), RoomObsZMm());
  std::fprintf(f, "# score_ball_R_mm = %.1f  V_cm3 = %.6e\n", gRoom.ballRMm,
               fBallVolumeCm3);
  std::fprintf(f, "# cut_mm = %.4f\n", gRoom.cutMm);
  for (const auto& s : gSlabs) {
    std::fprintf(f,
                 "# slab %-16s %-8s V_cm3 = %.6e  rho = %.3f  m_kg = %.3f\n",
                 s.name.c_str(), s.brick ? "brick" : "concrete", s.volumeCm3,
                 s.densityGCm3, s.massG / 1000.0);
  }
  std::fprintf(f, "# V_brick_cm3 = %.6e  m_brick_kg = %.3f\n",
               RoomBrickVolumeCm3(), RoomBrickMassG() / 1000.0);
  std::fprintf(f, "# V_concrete_cm3 = %.6e  m_concrete_kg = %.3f\n",
               RoomConcreteVolumeCm3(), RoomConcreteMassG() / 1000.0);
  std::fprintf(f, "# m_total_kg = %.3f\n",
               (RoomBrickMassG() + RoomConcreteMassG()) / 1000.0);
  std::fprintf(f, "# m_source_kg = %.3f  (масса, по которой нормирован ЭТОТ "
                  "файл, src=%s)\n",
               RoomSelectedMassG() / 1000.0, RoomSrcModeName());
  std::fprintf(f, "# R_decays_per_s = %.6e\n", fRatePerS);
  std::fprintf(f, "# T_equiv_s = %.6e\n", tRunS);
  std::fprintf(f, "# N = %lld\n", nEvents);
  std::fprintf(f, "# fluence_total_cm2_s = %.6e\n", total);
  std::fprintf(f, "# bin_keV = %.3f  (последний канал = переполнение)\n",
               kBinKeV);
  std::fprintf(f, "E_keV,fluence_cm2_s\n");

  // Самопроверка идёт по ТЕМ ЖЕ строкам, что уходят в файл: значение
  // форматируется, читается обратно и суммируется. Так проверяется содержимое
  // файла, а не внутренние double, которые в файл не попадают.
  double sumPrinted = 0.0;
  char buf[64];
  for (int i = 0; i <= kBins; ++i) {
    const double v = fLenCm[static_cast<std::size_t>(i)] * norm;
    if (!(v > 0.0)) continue;
    std::snprintf(buf, sizeof(buf), "%.6e", v);
    sumPrinted += std::atof(buf);
    std::fprintf(f, "%.1f,%s\n", (i + 0.5) * kBinKeV, buf);
  }
  std::fclose(f);

  std::snprintf(buf, sizeof(buf), "%.6e", total);
  const double headerPrinted = std::atof(buf);

  std::fprintf(stdout, "\n=== Rc103RoomField RESULT ===\n");
  std::fprintf(stdout, "nuclide= %s Z= %d A= %d\n", fNuclide.c_str(), fIonZ,
               fIonA);
  std::fprintf(stdout, "n_events= %lld\n", nEvents);
  std::fprintf(stdout, "V_brick_cm3= %.6e  m_brick_kg= %.3f\n",
               RoomBrickVolumeCm3(), RoomBrickMassG() / 1000.0);
  std::fprintf(stdout, "V_concrete_cm3= %.6e  m_concrete_kg= %.3f\n",
               RoomConcreteVolumeCm3(), RoomConcreteMassG() / 1000.0);
  std::fprintf(stdout, "R_decays_per_s= %.6e\n", fRatePerS);
  std::fprintf(stdout, "T_equiv_s= %.6e\n", tRunS);
  std::fprintf(stdout, "ball_V_cm3= %.6e\n", fBallVolumeCm3);
  std::fprintf(stdout, "fluence_total_cm2_s= %.6e\n", total);
  std::fprintf(stdout, "file= %s\n", fOutputCsv.c_str());

  // СТОРОЖ 1: нулевой флюенс — это ПРОВАЛ, а не «такой результат».
  // Долгоживущие ядра без /process/had/rdm/thresholdForVeryLongDecayTime не
  // распадаются вообще, а прогон при этом завершается кодом 0 и пишет файл.
  if (!(total > 0.0)) {
    std::fprintf(stderr,
                 "FLUENCE_ZERO: прогон дал нулевой флюенс. Проверь "
                 "/process/had/rdm/thresholdForVeryLongDecayTime.\n");
    std::exit(3);
  }

  // СТОРОЖ 2 (обязательная самопроверка спеки): шапка против колонки.
  const double rel =
      (headerPrinted > 0.0) ? (sumPrinted - headerPrinted) / headerPrinted : 1.0;
  std::fprintf(stdout,
               "SELFCHECK header_total= %.9e column_sum= %.9e rel_diff= %.3e "
               "tol= 1.000e-06\n",
               headerPrinted, sumPrinted, rel);
  if (!(std::fabs(rel) <= 1e-6)) {
    std::fprintf(stderr,
                 "SELFCHECK_FAILED: сумма колонки расходится с шапкой на "
                 "%.3e (порог 1e-6).\n",
                 rel);
    std::exit(4);
  }
  std::fprintf(stdout, "SELFCHECK_OK\n");
  std::fprintf(stdout, "=== END Rc103RoomField RESULT ===\n");
}
