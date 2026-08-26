#include "Rc103RunDetectorConstruction.hh"

#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"

#include <cctype>
#include <cstdio>
#include <cstdlib>

G4LogicalVolume* Rc103RunDetectorConstruction::fgCrystalLV = nullptr;

Rc103RunDetectorConstruction::Rc103RunDetectorConstruction(const G4String& gdmlPath)
    : fGdmlPath(gdmlPath) {}

G4VPhysicalVolume* Rc103RunDetectorConstruction::Construct() {
  fParser.SetOverlapCheck(true);  // не глушить P-025-подобные наложения
  fParser.Read(fGdmlPath, false);  // false = не валидировать против XSD-схемы
  G4VPhysicalVolume* world = fParser.GetWorldVolume();

  auto* lvStore = G4LogicalVolumeStore::GetInstance();

  // Основной кандидат — буквальное имя <volume name="Crystal_log"> из SSOT
  // (RC103_detector.gdml:175). GetVolume(..., false) не печатает встроенное
  // предупреждение Geant4 — диагностику делаем сами ниже, если не нашли.
  G4LogicalVolume* crystal = lvStore->GetVolume("Crystal_log", false);

  if (!crystal) {
    std::fprintf(stderr,
                 "Rc103RunDetectorConstruction: 'Crystal_log' NOT FOUND "
                 "verbatim after GDML parse. Dumping ALL logical volumes in "
                 "store (%zu total) for diagnosis:\n",
                 lvStore->size());
    for (auto* lv : *lvStore) {
      std::fprintf(stderr, "  - '%s'\n", lv->GetName().c_str());
    }

    // Фоллбэк — детерминированный, не гадание: единственный объём, чьё имя
    // содержит подстроку "crystal" (без учёта регистра). Если кандидатов не
    // ровно один — не выбираем наугад, падаем ниже.
    G4LogicalVolume* candidate = nullptr;
    int nCandidates = 0;
    for (auto* lv : *lvStore) {
      G4String lower = lv->GetName();
      for (auto& c : lower) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
      }
      if (lower.find("crystal") != G4String::npos) {
        candidate = lv;
        ++nCandidates;
      }
    }
    if (nCandidates == 1) {
      crystal = candidate;
      std::fprintf(stderr,
                   "Rc103RunDetectorConstruction: fallback matched EXACTLY "
                   "one candidate by substring 'crystal': '%s'\n",
                   crystal->GetName().c_str());
    } else {
      std::fprintf(stderr,
                   "Rc103RunDetectorConstruction: fallback found %d "
                   "candidates (need exactly 1) — refusing to guess.\n",
                   nCandidates);
    }
  }

  if (!crystal) {
    std::fprintf(stderr,
                 "Rc103RunDetectorConstruction: FATAL — could not resolve "
                 "crystal logical volume by any method. Aborting.\n");
    std::abort();
  }

  fgCrystalLV = crystal;
  return world;
}
