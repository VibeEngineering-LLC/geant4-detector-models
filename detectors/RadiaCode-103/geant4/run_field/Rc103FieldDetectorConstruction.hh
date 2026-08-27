// Rc103FieldDetectorConstruction — мир 300x300x300 мм из G4_AIR, в центре
// прибор RC-103, загруженный ТОЛЬКО из провалидированного GDML SSOT
// (detectors/RadiaCode-103/geant4/gdml/detector/RC103_detector.gdml).
// Паттерн поиска логического объёма и обработки отказа взят дословно из
// run/Rc103RunDetectorConstruction.cc и run_room/Rc103RoomDetectorConstruction.cc
// (собственные боевые комплекты этой же сессии, эталон).
//
// НЕ имеет отношения к detectors/RadiaCode-103/geometry/{main.cc,RCDetector.cc,
// RCDetector.hh} — та физическая линия прибора запрещена оператором (27.08.2026)
// и как источник чисел, и как образец кода.
//
// Режим kCheckNorm: прибор НЕ грузится (GDML вообще не парсится), в центре мира
// стоит контрольный шар G4_AIR радиусом 20 мм для track-length оценки флюенса.
#pragma once

#include "G4VUserDetectorConstruction.hh"
#include "G4GDMLParser.hh"
#include "G4String.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

class Rc103FieldDetectorConstruction : public G4VUserDetectorConstruction {
 public:
  // Мир — КОМНАТА (прямоугольная), полуразмеры по осям в мм. Задаются
  // ключами CLI room=<X>x<Y>x<Z> (полные габариты), поэтому не constexpr.
  // Дефолт помечен как ПРЕДПОЛОЖЕНИЕ: это округлённые габариты комнаты
  // разумного размера, а не реальное помещение (те данные локальные).
  // Обязан вмещать: сферу источника R_SRC = 300 мм и свинцовый домик
  // 250x250x435 мм (наружный габарит, полудиагональ 280.3 мм).
  static double gRoomHalfXMm;   // = 2000.0 по умолчанию
  static double gRoomHalfYMm;   // = 2000.0
  static double gRoomHalfZMm;   // = 1400.0
  // Радиус контрольного шара режима --check-norm.
  static constexpr double kCheckRadiusMm = 20.0;

  // --- Свинцовый домик (опция, ключ CLI shield=on|off, дефолт off) --------
  // РЕАЛЬНЫЕ размеры домика оператора (сообщены 15.08.2026): полость
  // 150x150 мм в плане, высота 385 мм; свинец 50 мм — стенки И дно; ВЕРХ
  // ОТКРЫТ (крышки нет); меди и кадмия нет; материал G4_Pb (NIST).
  // Наружный габарит 250x250x435 мм.
  //
  // ДОПУЩЕНИЕ оценочного расчёта (отступление от натуры): наружный габарит
  // домика центрируется на (0,0,0), где стоит прибор, а НЕ ставится дном на
  // пол. Причина строгая: сфера-источник имеет R_SRC = 300 мм и центрирована
  // в (0,0,0). При посадке прибора на дно полости верхний угол домика ушёл бы
  // на sqrt(125^2+125^2+376^2) = 415.7 мм от центра, то есть ЗА сферу —
  // источник оказался бы внутри свинца и расчёт стал бы бессмысленным. При
  // центрировании наружного габарита полудиагональ равна
  // sqrt(125^2+125^2+217.5^2) = 280.3 мм < 300 мм, вложенность корректна.
  //
  // Следствие допущения, которое надо знать при чтении результата: прибор в
  // (0,0,0) стоит НЕ точно в геометрическом центре полости, а на 25 мм ниже
  // него (центр полости при таком центрировании лежит на z = +25 мм, потому
  // что дно съедает 50 мм снизу, а сверху крышки нет). Проверено счётом:
  // вариант «центр ПОЛОСТИ в (0,0,0)» уводит дальний нижний угол на
  // sqrt(125^2+125^2+242.5^2) = 300.09 мм, то есть ЗА сферу — он упирается в
  // проверку вложенности ниже и потому отвергнут. При принятом варианте
  // прибор всё равно висит в объёме полости: 167.5 мм над дном свинца и
  // 217.5 мм под открытым верхом.
  static bool gShieldOn;                            // = false по умолчанию
  static constexpr double kShieldPbMm = 50.0;       // стенки и дно
  static constexpr double kShieldCavityXMm = 150.0;
  static constexpr double kShieldCavityYMm = 150.0;
  static constexpr double kShieldCavityZMm = 385.0;
  static constexpr double kShieldOuterXMm =
      kShieldCavityXMm + 2.0 * kShieldPbMm;         // 250
  static constexpr double kShieldOuterYMm =
      kShieldCavityYMm + 2.0 * kShieldPbMm;         // 250
  static constexpr double kShieldOuterZMm =
      kShieldCavityZMm + kShieldPbMm;               // 435: дно есть, крышки нет

  Rc103FieldDetectorConstruction(const G4String& gdmlPath, bool checkNormMode);
  G4VPhysicalVolume* Construct() override;

  // Кристалл CsI(Tl) — только в обычном режиме; в --check-norm остаётся nullptr.
  static G4LogicalVolume* GetCrystalLogicalVolume() { return fgCrystalLV; }
  // Контрольный шар — только в --check-norm; в обычном режиме nullptr.
  static G4LogicalVolume* GetCheckLogicalVolume() { return fgCheckLV; }
  // Свинец домика — ненулевой только если shield=on и домик реально построен.
  static G4LogicalVolume* GetShieldLogicalVolume() { return fgShieldLV; }

 private:
  // Строит свинцовый домик в worldLV вокруг уже поставленного deviceLV.
  // Обе обязательные проверки (вместимость полости, вложенность в сферу)
  // делаются здесь ФАКТОМ габаритов, отказ — FATAL + abort, не молчание.
  void BuildLeadShield(G4LogicalVolume* worldLV, G4LogicalVolume* deviceLV);

  G4String fGdmlPath;
  bool fCheckNormMode;
  G4GDMLParser fParser;
  static G4LogicalVolume* fgCrystalLV;
  static G4LogicalVolume* fgCheckLV;
  static G4LogicalVolume* fgShieldLV;
};
