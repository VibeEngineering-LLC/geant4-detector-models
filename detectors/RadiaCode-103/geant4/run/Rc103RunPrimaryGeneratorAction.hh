// Точечный источник Cs-137 (гамма 661.657 кэВ) через G4GeneralParticleSource,
// узкий конус на кристалл — обоснование числами (телесный угол, координаты
// из RC103_detector.gdml) см. в .cc.
#pragma once

#include "G4VUserPrimaryGeneratorAction.hh"

class G4GeneralParticleSource;
class G4Event;

class Rc103RunPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  Rc103RunPrimaryGeneratorAction();
  ~Rc103RunPrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event* event) override;

  // Полураствор конуса и телесный угол — экспортируются, чтобы main.cc мог
  // передать в RunAction тот же Omega_cone/4pi без дублирования числа.
  static constexpr double kConeHalfAngleDeg = 5.0;

 private:
  G4GeneralParticleSource* fGPS;
};
