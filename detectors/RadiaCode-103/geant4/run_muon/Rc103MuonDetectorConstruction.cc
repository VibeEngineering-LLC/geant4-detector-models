// Rc103MuonDetectorConstruction.cc
//
// Provenance:
// Shield geometry ported on 29.08.2026 from run_field/Rc103FieldDetectorConstruction.cc (lines 80-194).
// Numbers reported by operator on 15.08.2026. Copy made because classes live in different build projects,
// and extracting to a common module would require rebuilding 'run_field' during ongoing runs.
// If one copy is modified, the other MUST be modified too.
// Files geometry/RCDetector.* are NOT sources — forbidden by operator on 27.08.2026.

#include "Rc103MuonDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <cstdlib>

G4LogicalVolume* Rc103MuonDetectorConstruction::fgCrystalLV = nullptr;
G4LogicalVolume* Rc103MuonDetectorConstruction::fgShieldLV = nullptr;
// Дефолты = РЕАЛЬНАЯ постановка оператора (P-005), как и в run_field.
double Rc103MuonDetectorConstruction::gStandMm = 25.0;
bool Rc103MuonDetectorConstruction::gFlipUp = true;
double Rc103MuonDetectorConstruction::gDeviceZMm = 0.0;

Rc103MuonDetectorConstruction::Rc103MuonDetectorConstruction(
    const G4String& gdmlPath, double worldHalfMm, bool shieldOn, double zDiskMm)
    : fGdmlPath(gdmlPath), fWorldHalfMm(worldHalfMm), fShieldOn(shieldOn),
      fZDiskMm(zDiskMm) {}

namespace {
// Детерминированный поиск логического объёма по имени с фоллбэком на
// ЕДИНСТВЕННОЕ substring-совпадение — паттерн дословно из
// run_field/Rc103FieldDetectorConstruction.cc (эталон), не гадание.
G4LogicalVolume* FindVolumeOrFallback(const char* exactName,
                                      const char* substrLower) {
  auto* lvStore = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* lv = lvStore->GetVolume(exactName, false);
  if (lv) return lv;

  std::fprintf(stderr,
               "Rc103MuonDetectorConstruction: '%s' NOT FOUND verbatim after "
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
                 "Rc103MuonDetectorConstruction: fallback matched EXACTLY one "
                 "candidate by substring '%s': '%s'\n",
                 substrLower, candidate->GetName().c_str());
    return candidate;
  }
  std::fprintf(stderr,
               "Rc103MuonDetectorConstruction: fallback found %d candidates "
               "(need exactly 1) - refusing to guess.\n",
               nCandidates);
  return nullptr;
}
}  // namespace

void Rc103MuonDetectorConstruction::BuildLeadShield(G4LogicalVolume* worldLV,
                                                    G4LogicalVolume* deviceLV,
                                                    double deviceZMm) {
  auto* nist = G4NistManager::Instance();

  double outHX = (kShieldCavityXMm + 2.0 * kShieldPbMm) / 2.0;
  double outHY = (kShieldCavityYMm + 2.0 * kShieldPbMm) / 2.0;
  double outHZ = (kShieldCavityZMm + kShieldPbMm) / 2.0;

  double cavHX = kShieldCavityXMm / 2.0;
  double cavHY = kShieldCavityYMm / 2.0;

  double cavZMin = -outHZ + kShieldPbMm;
  double cavZMax = outHZ;

  // Check device fit
  auto* devBox = dynamic_cast<G4Box*>(deviceLV->GetSolid());
  if (!devBox) {
    std::fprintf(stderr, "Rc103MuonDetectorConstruction: SHIELD FATAL - Device solid '%s' is not a G4Box. Cannot check dimensions.\n", deviceLV->GetSolid()->GetName().c_str());
    std::abort();
  }

  double devHX = devBox->GetXHalfLength() / mm;
  double devHY = devBox->GetYHalfLength() / mm;
  double devHZ = devBox->GetZHalfLength() / mm;

  std::fprintf(stdout, "Rc103MuonDetectorConstruction: SHIELD Device dims (mm): X=%.2f Y=%.2f Z=%.2f\n", devHX*2, devHY*2, devHZ*2);
  std::fprintf(stdout, "Rc103MuonDetectorConstruction: SHIELD Cavity dims (mm): X=%.2f Y=%.2f Z=%.2f (range [%.2f, %.2f])\n", 
               kShieldCavityXMm, kShieldCavityYMm, kShieldCavityZMm, cavZMin, cavZMax);

  double gapX = cavHX - devHX;
  double gapY = cavHY - devHY;
  // Device Z range: [-devHZ, +devHZ]. Cavity Z range: [cavZMin, cavZMax].
  // Зазор вниз = (-devHZ) - cavZMin, зазор вверх = cavZMax - devHZ.
  // Генерация давала cavZMin + devHZ — то же по модулю, но с обратным знаком
  // (печаталось -158.75 вместо 158.75). На геометрию не влияло, только на печать.
  // P-005: прибор стоит на deviceZMm, зазоры считаются от фактического места.
  const double devZLo = deviceZMm - devHZ;
  const double devZHi = deviceZMm + devHZ;
  double gapZBottom = devZLo - cavZMin;
  double gapZTop = cavZMax - devZHi;

  std::fprintf(stdout, "Rc103MuonDetectorConstruction: SHIELD Gaps (mm): X=%.2f Y=%.2f Z_Bottom=%.2f Z_Top=%.2f\n", 
               gapX, gapY, gapZBottom, gapZTop);

  if (devHX > cavHX || devHY > cavHY || devZLo < cavZMin || devZHi > cavZMax) {
    std::fprintf(stderr, "Rc103MuonDetectorConstruction: SHIELD FATAL - Device does not fit in cavity.\n");
    std::abort();
  }

  // Check world fit
  double maxDist = std::sqrt(outHX*outHX + outHY*outHY + outHZ*outHZ);
  std::fprintf(stdout, "Rc103MuonDetectorConstruction: SHIELD Max corner distance from origin: %.2f mm\n", maxDist);

  if (outHX >= fWorldHalfMm || outHY >= fWorldHalfMm || outHZ >= fWorldHalfMm) {
    std::fprintf(stderr, "Rc103MuonDetectorConstruction: SHIELD FATAL - Shield (%.2f x %.2f x %.2f) exceeds world (%.2f).\n", 
                 outHX*2, outHY*2, outHZ*2, fWorldHalfMm);
    std::abort();
  }

  // Check muon disk position
  if (fZDiskMm <= cavZMax) {
    std::fprintf(stderr, "Rc103MuonDetectorConstruction: SHIELD FATAL - Muon start disk Z=%.2f is inside or below cavity top Z=%.2f. Muons would be generated behind shield.\n", 
                 fZDiskMm, cavZMax);
    std::abort();
  }

  // Build solid
  double cutHZ = 0.5 * (cavZMax - cavZMin) + kShieldPbMm;
  double cutZ0 = cavZMin + cutHZ;

  G4Box* outerBox = new G4Box("shield_pb_outer", outHX*mm, outHY*mm, outHZ*mm);
  G4Box* cutBox = new G4Box("shield_pb_cut", cavHX*mm, cavHY*mm, cutHZ*mm);
  
  // Сигнатура именно такая: (имя, A, B, поворот, сдвиг). Вариант без
  // указателя на поворот не существует — здесь генерация ошиблась.
  G4SubtractionSolid* shieldSolid = new G4SubtractionSolid(
      "shield_pb_solid", outerBox, cutBox, nullptr,
      G4ThreeVector(0, 0, cutZ0 * mm));

  G4LogicalVolume* shieldLV = new G4LogicalVolume(shieldSolid, nist->FindOrBuildMaterial("G4_Pb"), "shield_pb_log");
  
  new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), shieldLV, "pv_shield_pb", worldLV, false, 0, true);
  
  fgShieldLV = shieldLV;

  std::fprintf(stdout, 
               "Rc103MuonDetectorConstruction: SHIELD Built. Material=G4_Pb, Wall/Bottom=%.2f mm, Cavity=%.2fx%.2fx%.2f mm, Top OPEN, Center Z=%.2f, Device at (0,0,0).\n",
               kShieldPbMm, kShieldCavityXMm, kShieldCavityYMm, kShieldCavityZMm, (cavZMin+cavZMax)/2.0);
}

G4VPhysicalVolume* Rc103MuonDetectorConstruction::Construct() {
  auto* nist = G4NistManager::Instance();

  const double halfWorld = fWorldHalfMm * mm;
  auto* worldLV = new G4LogicalVolume(
      new G4Box("world_muon", halfWorld, halfWorld, halfWorld),
      nist->FindOrBuildMaterial("G4_AIR"), "world_muon");
  auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "world_muon", nullptr,
                                    false, 0, true);

  fParser.SetOverlapCheck(true);   // не глушить известные наложения
  fParser.Read(fGdmlPath, false);  // false = не валидировать против XSD

  // Берём "RC103_device_log" (сам прибор), НЕ "World" из парсера — тот
  // остаётся orphan-объёмом в сторе и из Construct() не возвращается.
  G4LogicalVolume* deviceLV =
      FindVolumeOrFallback("RC103_device_log", "device_log");
  if (!deviceLV) {
    std::fprintf(stderr,
                 "Rc103MuonDetectorConstruction: FATAL - could not resolve "
                 "RC103_device_log after GDML parse. Aborting.\n");
    std::abort();
  }

  G4LogicalVolume* crystalLV = FindVolumeOrFallback("Crystal_log", "crystal");
  if (!crystalLV) {
    std::fprintf(stderr,
                 "Rc103MuonDetectorConstruction: FATAL - could not resolve "
                 "Crystal_log after GDML parse. Aborting.\n");
    std::abort();
  }
  fgCrystalLV = crystalLV;

  // Посадка прибора (P-005) — та же арифметика, что в run_field: толщина из
  // СОЛИДА, высота от верха дна полости, разворот на 180° вокруг длинной оси X.
  double devHZmm = 0.0;
  if (auto* devBox0 = dynamic_cast<G4Box*>(deviceLV->GetSolid())) {
    devHZmm = devBox0->GetZHalfLength() / mm;
  } else {
    std::fprintf(stderr,
                 "Rc103MuonDetectorConstruction: FATAL солид прибора '%s' не "
                 "G4Box — высоту посадки не вычислить фактом.\n",
                 deviceLV->GetSolid()->GetName().c_str());
    std::abort();
  }
  double zDevMm = 0.0;
  if (fShieldOn && gStandMm >= 0.0) {
    zDevMm = -0.5 * kShieldOuterZMm + kShieldPbMm + gStandMm + devHZmm;
  }
  gDeviceZMm = zDevMm;

  G4RotationMatrix* rot = nullptr;
  if (gFlipUp) {
    rot = new G4RotationMatrix();
    rot->rotateX(180.0 * deg);
  }
  std::fprintf(stdout,
               "Rc103MuonDetectorConstruction: посадка прибора z=%.2f мм "
               "(stand=%.1f, экран %s)\n",
               zDevMm, gStandMm, gFlipUp ? "ВВЕРХ" : "вниз");

  // pSurfChk=true последним аргументом — проверка наложений с миром.
  new G4PVPlacement(rot, G4ThreeVector(0, 0, zDevMm * mm), deviceLV,
                    "pv_rc103_in_muon", worldLV, false, 0, true);

  if (fShieldOn) {
    BuildLeadShield(worldLV, deviceLV, zDevMm);
  }

  std::fprintf(stdout,
               "Rc103MuonDetectorConstruction: world half=%.1f mm, device "
               "placed at origin, Crystal_log resolved. Shield=%s, ZDisk=%.2f mm.\n",
               fWorldHalfMm, fShieldOn ? "ON" : "OFF", fZDiskMm);
  return worldPV;
}
