// Rc103RoomDetectorConstruction — «поле помещения»: сферическая полость
// R_CAV=20 см воздух внутри бетонной (G4_CONCRETE, NIST) сферы R_WALL=80 см,
// а В ЦЕНТРЕ полости — провалидированный GDML-детектор RadiaCode-103
// (SSOT: RC103_detector.gdml, P-001/P-004 закрыты, эталон интеграции —
// detectors/RadiaCode-103/geant4/run/Rc103RunDetectorConstruction.cc).
//
// Геометрия сферы/полости и активности ЕРН — ПЕРЕИСПОЛЬЗОВАНЫ КАК ФИЗИКА/ДАННЫЕ
// из detectors/RadiaCode-103/geometry/wallfield.cc (§33 доктрины: не
// изобретать заново готовый метод) — это ДРУГАЯ, старая физическая линия
// прибора (RCDetector.cc), оператором запрещённая как ОБРАЗЕЦ КОДА; отсюда
// новый файл, а не include/копия wallfield.cc. Взяты дословно: R_CAV=20см,
// R_WALL=80см (сферический дефолтный режим, "открыт по фону 0,88..1,09"),
// G4_CONCRETE NIST, таблица линий K-40/Ra-226/Th-232 (ENSDF/LNHB) и
// активности UNSCEAR (K-40=400, Ra-226=40, Th-232=30 Бк/кг) — см.
// Rc103RoomPrimaryGeneratorAction.cc.
//
// ИНТЕГРАЦИЯ GDML-ДЕТЕКТОРА (путь «а» из задания, сработал с первой попытки):
// один G4GDMLParser парсит ТОЛЬКО детекторный файл (парсер комнаты не нужен —
// комната строится программно, как WallGeom в wallfield.cc). Из
// G4LogicalVolumeStore после парсинга берётся логический объём
// "RC103_device_log" (НЕ "World" из GDML) — это сам прибор, солид Case_outer
// 123x34x17.5мм (RC103_detector.gdml:127,296-298), многократно меньше
// R_CAV=200мм полости, поэтому вкладывается как обычный physvol внутрь
// gCavLV БЕЗ конфликта с GDML-World-боксом 400x400x400мм (тот просто не
// используется — не возвращается из Construct(), остаётся orphan в
// G4PhysicalVolumeStore, что безвредно: навигация Geant4 идёт только по
// дереву, реально возвращённому из Construct()). Путь «б» (генерировать GDML
// комнаты и сливать секции regex'ом) не потребовался и не пробовался — путь
// «а» единственным парсером сработал без ошибок с первой попытки.
#pragma once

#include "G4VUserDetectorConstruction.hh"
#include "G4GDMLParser.hh"
#include "G4String.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

class Rc103RoomDetectorConstruction : public G4VUserDetectorConstruction {
 public:
  explicit Rc103RoomDetectorConstruction(const G4String& detectorGdmlPath);
  G4VPhysicalVolume* Construct() override;

  static G4LogicalVolume* GetCrystalLogicalVolume() { return fgCrystalLV; }
  static G4LogicalVolume* GetCavityLogicalVolume() { return fgCavityLV; }

  // Геометрия «поля помещения» — переиспользована из wallfield.cc (§33,
  // см. комментарий выше файла). см/дальнейший расчёт в PrimaryGeneratorAction.
  static constexpr double kRWallCm = 80.0;
  static constexpr double kRCavCm = 20.0;

 private:
  G4String fDetectorGdmlPath;
  G4GDMLParser fParser;
  static G4LogicalVolume* fgCrystalLV;
  static G4LogicalVolume* fgCavityLV;
};
