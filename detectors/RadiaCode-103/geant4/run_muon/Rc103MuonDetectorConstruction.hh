// Rc103MuonDetectorConstruction — мир из G4_AIR, в центре прибор RC-103,
// загруженный ТОЛЬКО из провалидированного GDML SSOT
// (detectors/RadiaCode-103/geant4/gdml/detector/RC103_detector.gdml).
//
// Паттерн поиска логического объёма и обработки отказа взят дословно из
// run_field/Rc103FieldDetectorConstruction.cc (соседний комплект той же
// сессии, эталон).
//
// НЕ имеет отношения к detectors/RadiaCode-103/geometry/{main.cc,RCDetector.cc,
// RCDetector.hh,cosmicmu.cc} — та физическая линия прибора запрещена оператором
// (27.08.2026) и как источник чисел, и как образец кода.
//
// Полуразмер мира задаётся конструктором, а НЕ константой: спека фиксирует
// 400 мм под R_DISK=150..300, но обязательная проверка насыщения допускает
// третий прогон при R_DISK=600 — источник обязан оставаться внутри мира.
// Правило: halfWorld = max(400, R_DISK + 100) мм (см. main.cc).
#pragma once

#include "G4GDMLParser.hh"
#include "G4String.hh"
#include "G4VUserDetectorConstruction.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

class Rc103MuonDetectorConstruction : public G4VUserDetectorConstruction {
 public:
  // Базовый полуразмер мира по спеке. Фактический = max(этого, R_DISK+100).
  static constexpr double kWorldHalfMmDefault = 400.0;

  Rc103MuonDetectorConstruction(const G4String& gdmlPath, double worldHalfMm);
  G4VPhysicalVolume* Construct() override;

  // Кристалл CsI(Tl). До Initialize() — nullptr.
  static G4LogicalVolume* GetCrystalLogicalVolume() { return fgCrystalLV; }

 private:
  G4String fGdmlPath;
  double fWorldHalfMm;
  G4GDMLParser fParser;
  static G4LogicalVolume* fgCrystalLV;
};
