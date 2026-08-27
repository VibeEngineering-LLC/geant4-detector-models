#include "Rc103MuonDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"

#include <cctype>
#include <cstdio>
#include <cstdlib>

G4LogicalVolume* Rc103MuonDetectorConstruction::fgCrystalLV = nullptr;

Rc103MuonDetectorConstruction::Rc103MuonDetectorConstruction(
    const G4String& gdmlPath, double worldHalfMm)
    : fGdmlPath(gdmlPath), fWorldHalfMm(worldHalfMm) {}

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

  // Прибор в (0,0,0) без поворота: корпус Case_outer 123 x 34 x 17.5 мм,
  // длинная ось по X, по Z всего +-8.75 мм. pSurfChk=true последним
  // аргументом — обязательная проверка наложений с миром.
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), deviceLV,
                    "pv_rc103_in_muon", worldLV, false, 0, true);

  std::fprintf(stdout,
               "Rc103MuonDetectorConstruction: world half=%.1f mm, device "
               "placed at origin, Crystal_log resolved.\n",
               fWorldHalfMm);
  return worldPV;
}
