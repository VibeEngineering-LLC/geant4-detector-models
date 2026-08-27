#include "Rc103RoomDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"

#include <cctype>
#include <cstdio>
#include <cstdlib>

G4LogicalVolume* Rc103RoomDetectorConstruction::fgCrystalLV = nullptr;
G4LogicalVolume* Rc103RoomDetectorConstruction::fgCavityLV = nullptr;

Rc103RoomDetectorConstruction::Rc103RoomDetectorConstruction(
    const G4String& detectorGdmlPath)
    : fDetectorGdmlPath(detectorGdmlPath) {}

namespace {
// Детерминированный поиск логического объёма по имени с фоллбэком на
// единственную substring-подстроку — тот же паттерн, что в
// Rc103RunDetectorConstruction.cc (эталон), не гадание.
G4LogicalVolume* FindVolumeOrFallback(const char* exactName,
                                       const char* substrLower) {
  auto* lvStore = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* lv = lvStore->GetVolume(exactName, false);
  if (lv) return lv;

  std::fprintf(stderr,
               "Rc103RoomDetectorConstruction: '%s' NOT FOUND verbatim after "
               "GDML parse. Dumping ALL logical volumes in store (%zu "
               "total):\n",
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
                 "Rc103RoomDetectorConstruction: fallback matched EXACTLY "
                 "one candidate by substring '%s': '%s'\n",
                 substrLower, candidate->GetName().c_str());
    return candidate;
  }
  std::fprintf(stderr,
               "Rc103RoomDetectorConstruction: fallback found %d candidates "
               "(need exactly 1) — refusing to guess.\n",
               nCandidates);
  return nullptr;
}
}  // namespace

G4VPhysicalVolume* Rc103RoomDetectorConstruction::Construct() {
  auto* nist = G4NistManager::Instance();

  // --- «Поле помещения»: сфера R_WALL бетон, полость R_CAV воздух ----------
  // (переиспользовано из wallfield.cc::WallGeom, режим по умолчанию БЕЗ
  // corner=, см. комментарий в .hh). Имена отличны от объёмов детекторного
  // GDML ("World"/"Case_outer"/...), коллизий в сторах нет.
  const double rWall = kRWallCm * cm;
  const double rCav = kRCavCm * cm;

  auto* worldLV = new G4LogicalVolume(
      new G4Box("world_room", 1.1 * rWall, 1.1 * rWall, 1.1 * rWall),
      nist->FindOrBuildMaterial("G4_AIR"), "world_room");
  auto* worldPV =
      new G4PVPlacement(nullptr, {}, worldLV, "world_room", nullptr, false, 0, true);

  auto* wallLV = new G4LogicalVolume(new G4Orb("wall_concrete", rWall),
                                      nist->FindOrBuildMaterial("G4_CONCRETE"),
                                      "wall_concrete");
  new G4PVPlacement(nullptr, {}, wallLV, "wall_concrete", worldLV, false, 0, true);

  auto* cavLV = new G4LogicalVolume(new G4Orb("cav_air", rCav),
                                     nist->FindOrBuildMaterial("G4_AIR"), "cav_air");
  new G4PVPlacement(nullptr, {}, cavLV, "cav_air", wallLV, false, 0, true);
  fgCavityLV = cavLV;

  // --- GDML-детектор: отдельный парсер, ТОЛЬКО детекторный файл ------------
  fParser.SetOverlapCheck(true);  // не глушить P-004-подобные наложения
  fParser.Read(fDetectorGdmlPath, false);  // false = не валидировать против XSD

  // Берём "RC103_device_log" (сам прибор, солид Case_outer 123x34x17.5мм —
  // см. .hh), НЕ "World" из парсера (тот 400x400x400мм, больше R_CAV=200мм
  // диаметром бы конфликтовал; но мы его и не используем вовсе — просто не
  // возвращаем из Construct(), он остаётся orphan-физобъёмом в сторе).
  G4LogicalVolume* deviceLV =
      FindVolumeOrFallback("RC103_device_log", "device_log");
  if (!deviceLV) {
    std::fprintf(stderr,
                 "Rc103RoomDetectorConstruction: FATAL — could not resolve "
                 "RC103_device_log after GDML parse. Aborting.\n");
    std::abort();
  }

  G4LogicalVolume* crystalLV = FindVolumeOrFallback("Crystal_log", "crystal");
  if (!crystalLV) {
    std::fprintf(stderr,
                 "Rc103RoomDetectorConstruction: FATAL — could not resolve "
                 "Crystal_log after GDML parse. Aborting.\n");
    std::abort();
  }
  fgCrystalLV = crystalLV;

  // Прибор — в истинном центре полости (0,0,0), без поворота: та же система
  // координат, что использует rc103_run (Rc103RunPrimaryGeneratorAction.cc),
  // единица правды по геометрии прибора не дублируется. pSurfChk=true —
  // последним аргументом, ОБЯЗАТЕЛЬНАЯ проверка новых наложений между этим
  // physvol и полостью (задание, п.7): геометрически невозможны (полудиагональ
  // прибора ~64мм << R_CAV=200мм), но проверяется фактом, не по расчёту.
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), deviceLV,
                     "pv_rc103_in_room", cavLV, false, 0, true);

  return worldPV;
}
