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

  // Свинцовый домик — те же числа, что в run_field (оператор, 15.08.2026):
  // стенки и дно 50 мм, полость 150 x 150 x 385 мм, ВЕРХ ОТКРЫТ (крышки нет),
  // поэтому наружная высота = полость + одно дно, а не + два.
  static constexpr double kShieldPbMm = 50.0;
  static constexpr double kShieldCavityXMm = 150.0;
  static constexpr double kShieldCavityYMm = 150.0;
  static constexpr double kShieldCavityZMm = 385.0;
  static constexpr double kShieldOuterXMm =
      kShieldCavityXMm + 2.0 * kShieldPbMm;         // 250
  static constexpr double kShieldOuterYMm =
      kShieldCavityYMm + 2.0 * kShieldPbMm;         // 250
  static constexpr double kShieldOuterZMm =
      kShieldCavityZMm + kShieldPbMm;               // 435: дно есть, крышки нет

  // --- Посадка прибора в полости (P-005, 30.08.2026) ----------------------
  // Те же параметры и та же семантика, что в run_field: реальная постановка —
  // прибор ПЛАШМЯ, ЭКРАНОМ ВВЕРХ, на картонной коробке 25 мм над дном полости.
  // gStandMm < 0 — прежнее допущение «центр габарита в (0,0,0)».
  // Мюонный расчёт обязан жить в той же посадке, что гамма: иначе складываются
  // две компоненты, посчитанные для разных геометрий, и сравнение с измерением
  // перестаёт означать то, что написано на этикетке.
  static double gStandMm;                           // = 25.0 (картон)
  static bool gFlipUp;                              // = true (экран вверх)
  static double gDeviceZMm;                         // фактическое z центра

  Rc103MuonDetectorConstruction(const G4String& gdmlPath, double worldHalfMm,
                                bool shieldOn = false, double zDiskMm = 100.0);
  G4VPhysicalVolume* Construct() override;

  // Кристалл CsI(Tl). До Initialize() — nullptr.
  static G4LogicalVolume* GetCrystalLogicalVolume() { return fgCrystalLV; }
  // Свинец домика. nullptr, если домик выключен либо до Initialize().
  static G4LogicalVolume* GetShieldLogicalVolume() { return fgShieldLV; }

 private:
  // deviceZMm — фактическое положение центра габарита прибора по Z: проверка
  // вместимости обязана считать зазоры от него, а не от нуля (P-005).
  void BuildLeadShield(G4LogicalVolume* worldLV, G4LogicalVolume* deviceLV,
                       double deviceZMm);

  G4String fGdmlPath;
  double fWorldHalfMm;
  bool fShieldOn = false;
  // Высота диска старта мюонов. Нужна здесь ради проверки: при включённом
  // домике диск обязан быть ВЫШЕ его открытого верха, иначе мюоны рождались
  // бы внутри полости, за защитой.
  double fZDiskMm = 100.0;
  G4GDMLParser fParser;
  static G4LogicalVolume* fgCrystalLV;
  static G4LogicalVolume* fgShieldLV;
};
