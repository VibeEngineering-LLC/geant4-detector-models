// Rc103RunDetectorConstruction — геометрия ТОЛЬКО из провалидированного GDML
// SSOT (detectors/RadiaCode-103/geant4/gdml/detector/RC103_detector.gdml),
// загружаемого через G4GDMLParser. Паттерн дословно как GdmlDetector в
// common/tools/gdml_vis.cc:27-39 (эталон), НЕ hardcoded C++ геометрия.
//
// НЕ имеет отношения к detectors/RadiaCode-103/geometry/RCDetector.cc/.hh
// (другая, старая физическая линия того же прибора — оператором запрещена
// как образец/источник кода для этого приложения).
#pragma once

#include "G4VUserDetectorConstruction.hh"
#include "G4GDMLParser.hh"
#include "G4String.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

class Rc103RunDetectorConstruction : public G4VUserDetectorConstruction {
 public:
  explicit Rc103RunDetectorConstruction(const G4String& gdmlPath);
  G4VPhysicalVolume* Construct() override;

  // Логический объём кристалла CsI(Tl), найденный ПОСЛЕ парсинга GDML по
  // имени, проверенному фактом через G4LogicalVolumeStore (см. .cc) — не
  // взято по памяти из текста .gdml. nullptr, пока Construct() не выполнен
  // либо если объём не удалось разрешить (тогда main.cc обязан упасть).
  static G4LogicalVolume* GetCrystalLogicalVolume() { return fgCrystalLV; }

 private:
  G4String fGdmlPath;
  G4GDMLParser fParser;
  static G4LogicalVolume* fgCrystalLV;
};
