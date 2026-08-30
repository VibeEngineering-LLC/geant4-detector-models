#include "Rc103FieldDetectorConstruction.hh"
#include "Rc103FieldPrimaryGeneratorAction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>

double Rc103FieldDetectorConstruction::gRoomHalfXMm = 2000.0;
double Rc103FieldDetectorConstruction::gRoomHalfYMm = 2000.0;
double Rc103FieldDetectorConstruction::gRoomHalfZMm = 1400.0;
bool Rc103FieldDetectorConstruction::gShieldOn = false;
// Дефолты = РЕАЛЬНАЯ постановка оператора (P-005), а не прежнее допущение.
double Rc103FieldDetectorConstruction::gStandMm = 25.0;
bool Rc103FieldDetectorConstruction::gFlipUp = true;
double Rc103FieldDetectorConstruction::gDeviceZMm = 0.0;

G4LogicalVolume* Rc103FieldDetectorConstruction::fgCrystalLV = nullptr;
G4LogicalVolume* Rc103FieldDetectorConstruction::fgCheckLV = nullptr;
G4LogicalVolume* Rc103FieldDetectorConstruction::fgShieldLV = nullptr;

Rc103FieldDetectorConstruction::Rc103FieldDetectorConstruction(
    const G4String& gdmlPath, bool checkNormMode)
    : fGdmlPath(gdmlPath), fCheckNormMode(checkNormMode) {}

namespace {
// Детерминированный поиск логического объёма по имени с фоллбэком на
// ЕДИНСТВЕННОЕ substring-совпадение — паттерн дословно из
// run/Rc103RunDetectorConstruction.cc (эталон), не гадание.
G4LogicalVolume* FindVolumeOrFallback(const char* exactName,
                                      const char* substrLower) {
  auto* lvStore = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* lv = lvStore->GetVolume(exactName, false);
  if (lv) return lv;

  std::fprintf(stderr,
               "Rc103FieldDetectorConstruction: '%s' NOT FOUND verbatim after "
               "GDML parse. Dumping ALL logical volumes in store (%zu total):\n",
               exactName, lvStore->size());
  for (auto* v : *lvStore) {
    std::fprintf(stderr, "  - '%s'\n", v->GetName().c_str());
  }

  G4LogicalVolume* candidate = nullptr;
  int nCandidates = 0;
  for (auto* v : *lvStore) {
    G4String lower = v->GetName();
    for (auto& c : lower) {
      c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    if (lower.find(substrLower) != G4String::npos) {
      candidate = v;
      ++nCandidates;
    }
  }
  if (nCandidates == 1) {
    std::fprintf(stderr,
                 "Rc103FieldDetectorConstruction: fallback matched EXACTLY one "
                 "candidate by substring '%s': '%s'\n",
                 substrLower, candidate->GetName().c_str());
    return candidate;
  }
  std::fprintf(stderr,
               "Rc103FieldDetectorConstruction: fallback found %d candidates "
               "(need exactly 1) - refusing to guess.\n",
               nCandidates);
  return nullptr;
}
}  // namespace

void Rc103FieldDetectorConstruction::BuildLeadShield(G4LogicalVolume* worldLV,
                                                    G4LogicalVolume* deviceLV,
                                                    double deviceZMm) {
  auto* nist = G4NistManager::Instance();

  // Полуразмеры в мм (наружный габарит центрирован в (0,0,0) — см. развёрнутое
  // обоснование допущения в заголовке класса).
  const double outHX = 0.5 * kShieldOuterXMm;   // 125
  const double outHY = 0.5 * kShieldOuterYMm;   // 125
  const double outHZ = 0.5 * kShieldOuterZMm;   // 217.5

  // Полость в системе координат домика (= мировой, домик в (0,0,0)):
  // по x,y симметрична; по z снизу ограничена дном 50 мм, сверху ОТКРЫТА.
  const double cavHX = 0.5 * kShieldCavityXMm;  // 75
  const double cavHY = 0.5 * kShieldCavityYMm;  // 75
  const double cavZMin = -outHZ + kShieldPbMm;  // -167.5, верх дна свинца
  const double cavZMax = +outHZ;                // +217.5, срез открытого верха

  // --- ПРОВЕРКА 1: вместимость полости. Габарит прибора берётся ФАКТОМ из
  // солида, разрешённого после разбора GDML, а не переписывается числом.
  auto* devBox = dynamic_cast<G4Box*>(deviceLV->GetSolid());
  if (!devBox) {
    std::fprintf(stderr,
                 "Rc103FieldDetectorConstruction: FATAL солид прибора '%s' не "
                 "G4Box — габарит не проверить фактом, отказываюсь строить "
                 "домик молча.\n",
                 deviceLV->GetSolid()->GetName().c_str());
    std::abort();
  }
  const double devHX = devBox->GetXHalfLength() / mm;
  const double devHY = devBox->GetYHalfLength() / mm;
  const double devHZ = devBox->GetZHalfLength() / mm;
  // Прибор стоит на z = deviceZMm (P-005: посадка — параметр, а не ноль),
  // поэтому его габарит по Z в мировых координатах — [zLo, zHi].
  // Разворот на 180° вокруг X габарит не меняет: тело симметрично по Y и Z.
  const double zLo = deviceZMm - devHZ;
  const double zHi = deviceZMm + devHZ;
  const bool fitsX = (devHX <= cavHX) && (-devHX >= -cavHX);
  const bool fitsY = (devHY <= cavHY) && (-devHY >= -cavHY);
  const bool fitsZ = (zHi <= cavZMax) && (zLo >= cavZMin);
  std::fprintf(stdout,
               "Rc103FieldDetectorConstruction: SHIELD прибор %.1fx%.1fx%.1f мм "
               "на z=%.2f мм в полость %.1fx%.1f мм (z от %.1f до %.1f мм); "
               "зазоры: x=%.1f y=%.1f вниз=%.1f вверх=%.1f мм\n",
               2 * devHX, 2 * devHY, 2 * devHZ, deviceZMm, kShieldCavityXMm,
               kShieldCavityYMm, cavZMin, cavZMax, cavHX - devHX,
               cavHY - devHY, zLo - cavZMin, cavZMax - zHi);
  if (!fitsX || !fitsY || !fitsZ) {
    std::fprintf(stderr,
                 "Rc103FieldDetectorConstruction: FATAL прибор %.1fx%.1fx%.1f мм "
                 "НЕ помещается в полость %.1fx%.1fx%.1f мм "
                 "(fitsX=%d fitsY=%d fitsZ=%d).\n",
                 2 * devHX, 2 * devHY, 2 * devHZ, kShieldCavityXMm,
                 kShieldCavityYMm, kShieldCavityZMm, int(fitsX), int(fitsY),
                 int(fitsZ));
    std::abort();
  }

  // --- ПРОВЕРКА 2: вложенность наружного габарита в сферу-источник.
  // Считаем максимум по ВСЕМ 8 углам от (0,0,0), а не «на глаз».
  double maxCornerMm = 0.0;
  for (int sx = -1; sx <= 1; sx += 2) {
    for (int sy = -1; sy <= 1; sy += 2) {
      for (int sz = -1; sz <= 1; sz += 2) {
        const double x = sx * outHX, y = sy * outHY, z = sz * outHZ;
        maxCornerMm = std::max(maxCornerMm, std::sqrt(x * x + y * y + z * z));
      }
    }
  }
  const double rSrcMm = Rc103FieldPrimaryGeneratorAction::kRSrcMm;
  std::fprintf(stdout,
               "Rc103FieldDetectorConstruction: SHIELD наружный габарит "
               "%.0fx%.0fx%.0f мм, полудиагональ от центра = %.2f мм, "
               "R_SRC = %.2f мм\n",
               kShieldOuterXMm, kShieldOuterYMm, kShieldOuterZMm, maxCornerMm,
               rSrcMm);
  if (maxCornerMm >= rSrcMm) {
    std::fprintf(stderr,
                 "Rc103FieldDetectorConstruction: FATAL домик НЕ вложен в сферу "
                 "источника: полудиагональ %.2f мм >= R_SRC %.2f мм. Часть "
                 "источника оказалась бы внутри свинца.\n",
                 maxCornerMm, rSrcMm);
    std::abort();
  }
  // Отдельно — комната обязана вмещать сам домик.
  if (outHX >= gRoomHalfXMm || outHY >= gRoomHalfYMm || outHZ >= gRoomHalfZMm) {
    std::fprintf(stderr,
                 "Rc103FieldDetectorConstruction: FATAL домик %.0fx%.0fx%.0f мм "
                 "не влезает в комнату %.0fx%.0fx%.0f мм.\n",
                 kShieldOuterXMm, kShieldOuterYMm, kShieldOuterZMm,
                 2 * gRoomHalfXMm, 2 * gRoomHalfYMm, 2 * gRoomHalfZMm);
    std::abort();
  }

  // --- Солид: наружный параллелепипед МИНУС полость. Вычитаемый бокс нарочно
  // выше полости и выступает за верхнюю грань — так верх остаётся ОТКРЫТЫМ
  // (крышки нет), а дно 50 мм сохраняется.
  const double cutHZ = 0.5 * (cavZMax - cavZMin) + kShieldPbMm;  // выступ вверх
  const double cutZ0 = cavZMin + cutHZ;  // центр выреза
  auto* outerBox = new G4Box("shield_pb_outer", outHX * mm, outHY * mm, outHZ * mm);
  auto* cutBox = new G4Box("shield_pb_cut", cavHX * mm, cavHY * mm, cutHZ * mm);
  auto* shieldSolid = new G4SubtractionSolid(
      "shield_pb_solid", outerBox, cutBox, nullptr,
      G4ThreeVector(0, 0, cutZ0 * mm));

  auto* shieldLV = new G4LogicalVolume(
      shieldSolid, nist->FindOrBuildMaterial("G4_Pb"), "shield_pb_log");
  // pSurfChk=true — обязательная проверка наложений с уже поставленным прибором.
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), shieldLV, "pv_shield_pb",
                    worldLV, false, 0, true);
  fgShieldLV = shieldLV;

  std::fprintf(stdout,
               "Rc103FieldDetectorConstruction: SHIELD ON — G4_Pb, стенки и дно "
               "%.0f мм, полость %.0fx%.0fx%.0f мм, верх ОТКРЫТ; центр полости "
               "на z=%+.1f мм, прибор в (0,0,0).\n",
               kShieldPbMm, kShieldCavityXMm, kShieldCavityYMm,
               kShieldCavityZMm, 0.5 * (cavZMin + cavZMax));
}

G4VPhysicalVolume* Rc103FieldDetectorConstruction::Construct() {
  auto* nist = G4NistManager::Instance();

  // Мир = комната (прямоугольный объём воздуха). Источник — сфера R_SRC
  // ВНУТРИ неё, поэтому комната обязана быть заметно больше сферы; проверка
  // ниже отказывается строить заведомо неверную геометрию, а не молчит.
  const double hx = gRoomHalfXMm * mm;
  const double hy = gRoomHalfYMm * mm;
  const double hz = gRoomHalfZMm * mm;
  const double rSrc = Rc103FieldPrimaryGeneratorAction::kRSrcMm * mm;
  if (hx <= rSrc || hy <= rSrc || hz <= rSrc) {
    std::fprintf(stderr,
                 "Rc103FieldDetectorConstruction: FATAL комната %.0fx%.0fx%.0f мм "
                 "не вмещает сферу источника R=%.0f мм.\n",
                 2 * gRoomHalfXMm, 2 * gRoomHalfYMm, 2 * gRoomHalfZMm,
                 Rc103FieldPrimaryGeneratorAction::kRSrcMm);
    std::abort();
  }
  auto* worldLV = new G4LogicalVolume(
      new G4Box("world_field", hx, hy, hz),
      nist->FindOrBuildMaterial("G4_AIR"), "world_field");
  auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "world_field",
                                    nullptr, false, 0, true);

  if (fCheckNormMode) {
    // Самопроверка нормировки: прибора нет вовсе, GDML не парсится.
    // Контрольный шар воздуха R=20 мм, V = 4/3*pi*2.0^3 = 33.5103 см^3.
    auto* checkLV =
        new G4LogicalVolume(new G4Orb("check_ball", kCheckRadiusMm * mm),
                            nist->FindOrBuildMaterial("G4_AIR"), "check_ball");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), checkLV, "pv_check_ball",
                      worldLV, false, 0, true);
    fgCheckLV = checkLV;
    fgCrystalLV = nullptr;
    fgShieldLV = nullptr;
    if (gShieldOn) {
      std::fprintf(stderr,
                   "Rc103FieldDetectorConstruction: FATAL shield=on несовместим "
                   "с --check-norm: домик экранирует контрольный шар и ломает "
                   "смысл самопроверки нормировки.\n");
      std::abort();
    }
    std::fprintf(stdout,
                 "Rc103FieldDetectorConstruction: CHECK-NORM mode - GDML NOT "
                 "parsed, only air check ball R=%.1f mm at origin.\n",
                 kCheckRadiusMm);
    return worldPV;
  }

  fParser.SetOverlapCheck(true);           // не глушить известные наложения
  fParser.Read(fGdmlPath, false);          // false = не валидировать против XSD

  // Берём "RC103_device_log" (сам прибор), НЕ "World" из парсера (тот 400 мм и
  // здесь не нужен: он остаётся orphan-объёмом в сторе, из Construct() не
  // возвращается).
  G4LogicalVolume* deviceLV =
      FindVolumeOrFallback("RC103_device_log", "device_log");
  if (!deviceLV) {
    std::fprintf(stderr,
                 "Rc103FieldDetectorConstruction: FATAL - could not resolve "
                 "RC103_device_log after GDML parse. Aborting.\n");
    std::abort();
  }

  G4LogicalVolume* crystalLV = FindVolumeOrFallback("Crystal_log", "crystal");
  if (!crystalLV) {
    std::fprintf(stderr,
                 "Rc103FieldDetectorConstruction: FATAL - could not resolve "
                 "Crystal_log after GDML parse. Aborting.\n");
    std::abort();
  }
  fgCrystalLV = crystalLV;
  fgCheckLV = nullptr;

  // --- Посадка прибора (P-005) -------------------------------------------
  // Толщина берётся из СОЛИДА, а не из константы: если GDML изменится, высота
  // посадки поедет вместе с ним, а не разойдётся молча.
  double devHZmm = 0.0;
  if (auto* devBox0 = dynamic_cast<G4Box*>(deviceLV->GetSolid())) {
    devHZmm = devBox0->GetZHalfLength() / mm;
  } else {
    std::fprintf(stderr,
                 "Rc103FieldDetectorConstruction: FATAL солид прибора '%s' не "
                 "G4Box — высоту посадки не вычислить фактом.\n",
                 deviceLV->GetSolid()->GetName().c_str());
    std::abort();
  }

  double zDevMm = 0.0;   // режим «как считалось до 30.08»
  if (gShieldOn && gStandMm >= 0.0) {
    const double cavZMinMm = -0.5 * kShieldOuterZMm + kShieldPbMm;  // верх дна
    zDevMm = cavZMinMm + gStandMm + devHZmm;
  } else if (!gShieldOn && gStandMm >= 0.0) {
    // Без домика дна нет — опору отсчитывать не от чего. Это не «поправим
    // молча»: пара «без домика / с домиком» обязана различаться ровно домиком.
    zDevMm = 0.0;
  }
  gDeviceZMm = zDevMm;

  // Разворот на 180° вокруг длинной оси X: экран (в GDML при -Z) уходит вверх.
  G4RotationMatrix* rot = nullptr;
  if (gFlipUp) {
    rot = new G4RotationMatrix();
    rot->rotateX(180.0 * deg);
  }
  std::fprintf(stdout,
               "Rc103FieldDetectorConstruction: посадка прибора z=%.2f мм "
               "(stand=%.1f, экран %s, полутолщина %.2f мм)\n",
               zDevMm, gStandMm, gFlipUp ? "ВВЕРХ" : "вниз", devHZmm);

  // pSurfChk=true последним аргументом — обязательная проверка наложений.
  new G4PVPlacement(rot, G4ThreeVector(0, 0, zDevMm * mm), deviceLV,
                    "pv_rc103_in_field", worldLV, false, 0, true);

  // Свинцовый домик — строго ПОСЛЕ прибора: так CheckOverlaps домика сверяется
  // с уже поставленным прибором, а не наоборот.
  fgShieldLV = nullptr;
  if (gShieldOn) {
    BuildLeadShield(worldLV, deviceLV, zDevMm);
  } else {
    std::fprintf(stdout,
                 "Rc103FieldDetectorConstruction: SHIELD OFF — свинцовый домик "
                 "не строится.\n");
  }

  return worldPV;
}
