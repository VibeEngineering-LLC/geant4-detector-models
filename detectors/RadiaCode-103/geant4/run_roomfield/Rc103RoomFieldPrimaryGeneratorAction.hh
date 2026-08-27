// Источник — МЕТОД 1 (решение D-001 контура, 21.08.2026): первичная частица
// это САМО ЯДРО в покое, а не список гамма-линий. Схему распада со всеми
// гаммами, бета-спектрами, конверсией, характеристическим рентгеном и
// тормозным даёт G4RadioactiveDecay из ENSDF; /process/had/rdm/nucleusLimits
// (выставляется в main.cc) ограничивает распад ОДНИМ звеном, поэтому шаблон
// звена остаётся шаблоном звена, а не куском цепочки.
//
// РОЗЫГРЫШ ТОЧКИ РОЖДЕНИЯ — ПО МАССЕ, НЕ ПО ОБЪЁМУ. Удельная активность
// одинакова (1 Бк/кг) для кирпича и бетона, а плотности разные, поэтому плита
// выбирается с вероятностью, пропорциональной её МАССЕ, и только внутри
// выбранной плиты точка равномерна по объёму. Розыгрыш по объёму дал бы
// систематически неверную долю распадов в бетоне (rho 2,3 против 1,8).
#pragma once

#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include <vector>

class G4Event;
class G4ParticleDefinition;

class Rc103RoomFieldPrimaryGeneratorAction
    : public G4VUserPrimaryGeneratorAction {
 public:
  Rc103RoomFieldPrimaryGeneratorAction(int ionZ, int ionA);
  ~Rc103RoomFieldPrimaryGeneratorAction() override = default;

  void GeneratePrimaries(G4Event* event) override;

 private:
  int fIonZ;
  int fIonA;
  G4ParticleDefinition* fIon = nullptr;  // определяется лениво: таблица ионов
                                         // существует только после Initialize()
  std::vector<double> fCumMass;          // нормированные кумулятивные массы
  std::vector<std::size_t> fSlabIndex;   // индексы ВЫБРАННЫХ плит в gSlabs
  G4ParticleGun fGun;
};
