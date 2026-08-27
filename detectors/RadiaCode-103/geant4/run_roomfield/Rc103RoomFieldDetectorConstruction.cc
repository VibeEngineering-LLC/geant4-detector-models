#include "Rc103RoomFieldDetectorConstruction.hh"

#include "Rc103RoomFieldGeometry.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>

G4LogicalVolume* Rc103RoomFieldDetectorConstruction::fgBallLV = nullptr;
double Rc103RoomFieldDetectorConstruction::fgRhoConcrete = 0.0;

G4VPhysicalVolume* Rc103RoomFieldDetectorConstruction::Construct() {
  auto* nist = G4NistManager::Instance();

  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  // Кирпич: состав бетона при плотности кирпича — ПРИБЛИЖЕНИЕ (см. заголовок).
  G4Material* brick = nist->BuildMaterialWithNewDensity(
      "G4_BRICK_approx", "G4_CONCRETE", gRoom.rhoBrickGCm3 * g / cm3);

  // Перекрытия. При rhoSlabGCm3 > 0 строится ГОМОГЕНИЗИРОВАННЫЙ материал:
  // состав G4_CONCRETE при ЭФФЕКТИВНОЙ плотности круглопустотной плиты
  // (паспортная масса / габаритный объём). Цилиндрические пустоты явно НЕ
  // моделируются — для оценочного расчёта гомогенизации по массе достаточно,
  // но излучение ВДОЛЬ оси пустот такая модель ослабляет сильнее реального.
  G4Material* concrete = nullptr;
  if (gRoom.rhoSlabGCm3 > 0.0) {
    concrete = nist->BuildMaterialWithNewDensity(
        "G4_HOLLOWCORE_approx", "G4_CONCRETE", gRoom.rhoSlabGCm3 * g / cm3);
  } else {
    concrete = nist->FindOrBuildMaterial("G4_CONCRETE");
  }

  fgRhoConcrete = concrete->GetDensity() / (g / cm3);

  // --- мир ----------------------------------------------------------------
  const double wx = RoomWorldHalfXMm() * mm;
  const double wy = RoomWorldHalfYMm() * mm;
  const double wz = RoomWorldHalfZMm() * mm;
  auto* worldLV = new G4LogicalVolume(new G4Box("world_room", wx, wy, wz), air,
                                      "world_room");
  auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "world_room", nullptr,
                                    false, 0, true);

  // --- воздух комнаты -----------------------------------------------------
  const double hx = RoomInnerHalfXMm();
  const double hy = RoomInnerHalfYMm();
  const double hz = RoomInnerHalfZMm();
  auto* roomLV = new G4LogicalVolume(
      new G4Box("room_air", hx * mm, hy * mm, hz * mm), air, "room_air");
  new G4PVPlacement(nullptr, {}, roomLV, "room_air", worldLV, false, 0, true);

  // --- ограждающие конструкции -------------------------------------------
  if (gSlabs.empty()) {
    std::fprintf(stderr,
                 "Rc103RoomFieldDetectorConstruction: FATAL gSlabs пуст — "
                 "BuildRoomSlabs() не был вызван до Construct().\n");
    std::abort();
  }
  for (const auto& s : gSlabs) {
    auto* solid =
        new G4Box(s.name, s.hxMm * mm, s.hyMm * mm, s.hzMm * mm);
    auto* lv = new G4LogicalVolume(solid, s.brick ? brick : concrete, s.name);
    // pSurfChk=true: наложения плит между собой и с воздухом комнаты не
    // проглатываются молча. Но одного этого мало — приёмка геометрии
    // закрывается реальным BeamOn (см. main.cc), а не WARNING'ом отсюда.
    new G4PVPlacement(nullptr,
                      G4ThreeVector(s.cxMm * mm, s.cyMm * mm, s.czMm * mm), lv,
                      s.name, worldLV, false, 0, true);
    // Сверка назначенной массы с фактической плотностью материала: если
    // AssignRoomDensities() получила не ту плотность, нормировка была бы
    // неверной, а прогон выглядел бы успешным.
    const double rhoActual =
        lv->GetMaterial()->GetDensity() / (g / cm3);
    if (s.densityGCm3 <= 0.0 ||
        std::fabs(rhoActual - s.densityGCm3) / rhoActual > 1e-6) {
      std::fprintf(stderr,
                   "Rc103RoomFieldDetectorConstruction: FATAL плита '%s': "
                   "плотность материала %.6f г/см3 не совпадает с принятой в "
                   "нормировке %.6f г/см3.\n",
                   s.name.c_str(), rhoActual, s.densityGCm3);
      std::abort();
    }
  }

  // --- шар-скорер ---------------------------------------------------------
  const double ox = RoomObsXMm(), oy = RoomObsYMm(), oz = RoomObsZMm();
  const double R = gRoom.ballRMm;
  // Шар обязан целиком лежать в воздухе комнаты: иначе часть треков считалась
  // бы в кирпиче, а объём в знаменателе остался бы полным — молчаливое
  // искажение флюенса.
  const bool fits = (ox - R >= -hx) && (ox + R <= hx) && (oy - R >= -hy) &&
                    (oy + R <= hy) && (oz - R >= -hz) && (oz + R <= hz);
  if (!fits) {
    std::fprintf(stderr,
                 "Rc103RoomFieldDetectorConstruction: FATAL шар R=%.1f мм в "
                 "точке (%.1f, %.1f, %.1f) не помещается в воздух комнаты "
                 "полуразмерами (%.1f, %.1f, %.1f) мм.\n",
                 R, ox, oy, oz, hx, hy, hz);
    std::abort();
  }
  auto* ballLV = new G4LogicalVolume(new G4Orb("score_ball", R * mm), air,
                                     "score_ball");
  new G4PVPlacement(nullptr, G4ThreeVector(ox * mm, oy * mm, oz * mm), ballLV,
                    "score_ball", roomLV, false, 0, true);
  fgBallLV = ballLV;

  std::fprintf(stdout,
               "Rc103RoomFieldDetectorConstruction: комната %.0fx%.0fx%.0f мм; "
               "стены X- %.0f X+ %.0f Y- %.0f Y+ %.0f мм (кирпич, rho=%.3f); "
               "пол %.0f потолок %.0f мм (G4_CONCRETE, rho=%.3f); шар R=%.0f мм "
               "в (%.1f, %.1f, %.1f) мм.\n",
               gRoom.innerXMm, gRoom.innerYMm, gRoom.innerZMm, gRoom.wallXmMm,
               gRoom.wallXpMm, gRoom.wallYmMm, gRoom.wallYpMm,
               gRoom.rhoBrickGCm3, gRoom.floorMm, gRoom.ceilMm, fgRhoConcrete,
               R, ox, oy, oz);
  std::fprintf(stdout,
               "Rc103RoomFieldDetectorConstruction: ПРИБЛИЖЕНИЕ — кирпич задан "
               "элементным составом G4_CONCRETE при плотности %.3f г/см3.\n",
               gRoom.rhoBrickGCm3);

  return worldPV;
}
