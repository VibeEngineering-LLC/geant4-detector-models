// Геометрия RadiaCode-110: НАТИВНО в C++ (G4Box + G4SubtractionSolid), без
// GDML-парсера (в этой сборке Geant4 11.2.1 нет Xerces-C, G4GDMLParser
// недоступен — проверено фактом, см. бриф задачи/GEANT4-MODEL.md).
//
// ПЛОСКОЕ РАЗМЕЩЕНИЕ: все компоненты подвешены напрямую к мировому
// логическому объёму (мать = world), без вложенной иерархии. Координаты уже
// пересчитаны в мировой кадр в источнике — GDML-модели прибора
// (RC110_detector.gdml, дважды исправлена и независимо перепроверена
// 26.08.2026, числа финальные). Единицы — мм, начало координат мировое
// (не привязано к центру кристалла, в отличие от RadiaCode-103/RCDetector.hh).
//
// Донор паттерна кода (хелпер Put(), стиль класса, наследование от
// G4VUserDetectorConstruction): RadiaCode-103/geometry/RCDetector.{hh,cc}.
// Сосуда/матрицы пробы здесь нет — это vis-only геометрия для рендера,
// BeamOn не вызывается, метрологическая точность материалов не гарантируется
// (см. комментарии DefineMaterials() в .cc, там же — какие NIST-имена не
// нашлись в этой сборке и чем заменены).
#pragma once
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4Material;

class RC110Detector : public G4VUserDetectorConstruction {
public:
  RC110Detector() = default;
  G4VPhysicalVolume* Construct() override;

  // Полуразмеры мира, мм. Прибор — 63.3 x 17.05 x 10.85 (половины), запас
  // по X больше (вдоль длинной оси корпуса), см. бриф задачи.
  double fWorldHalfX = 150.0;
  double fWorldHalfY = 80.0;
  double fWorldHalfZ = 80.0;

protected:
  void DefineMaterials();
  void BuildDevice(G4LogicalVolume* world);

  // Материалы, которые могут отсутствовать в этой NIST DB и требуют
  // запасного варианта — собираются один раз в DefineMaterials() и
  // используются в BuildDevice().
  G4Material* fMatPcb     = nullptr;  // G4_G10 либо запасной FR4
  G4Material* fMatBattery = nullptr;  // G4_POLYVINYL_ACETATE либо запасной
};