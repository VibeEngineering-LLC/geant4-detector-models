// Источник поля помещения: фотон из случайной точки в толще бетона
// (равномерно по объёму сферического слоя R_CAV..R_WALL, отбраковка точек
// внутри полости), энергия и направление — розыгрыш по таблице линий
// ЕРН K-40/Ra-226-ряд/Th-232-ряд, изотропно (истинное направление вылета
// кванта из точки распада, НЕ "на детектор" — иначе поле перестало бы быть
// полем). Таблица линий и активности — ДАННЫЕ, переиспользованные из
// detectors/RadiaCode-103/geometry/wallfield.cc (см. .cc, §33 доктрины).
#pragma once

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

class G4Event;

class Rc103RoomPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  Rc103RoomPrimaryGeneratorAction();
  ~Rc103RoomPrimaryGeneratorAction() override = default;

  void GeneratePrimaries(G4Event* event) override;

  // Полная объёмная скорость испускания фотонов во всей толще бетона,
  // фотон/с — для пересчёта hits/N_events в АБСОЛЮТНУЮ прогнозную скорость
  // счёта в кристалле (см. Rc103RoomRunAction). Считается один раз в
  // конструкторе (BuildWeights), детерминированно по таблице линий и
  // объёму стены; не зависит от розыгрыша.
  static double GetWallEmissionRatePerSec() { return sWallEmissionRatePerSec; }

 private:
  G4ParticleGun fGun;
  static double sWallEmissionRatePerSec;
};
