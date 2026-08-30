// Объёмный источник: ион (Pb-214 или Bi-214, продукты радона) рождается
// равномерно по КАРКАСУ полости домика (тот же бокс, что kShieldCavityXMm/
// YMm и высота между полом и открытым верхом, см. Rc103FieldDetectorConstruction).
// Прибор внутри полости занимает ~1% её объёма (123x34x17.5 мм из
// 150x150x335 мм) — отбраковка точек внутри прибора СОЗНАТЕЛЬНО не сделана:
// смещение оценки первого порядка ниже погрешности, которую несёт сама
// гипотеза "радон в воздухе" (30.08.2026, см. project-incidents.md P-005).
#pragma once

#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4Event;
class G4ParticleDefinition;

class Rc103ShieldAirPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  Rc103ShieldAirPrimaryGeneratorAction(int ionZ, int ionA);
  ~Rc103ShieldAirPrimaryGeneratorAction() override = default;

  void GeneratePrimaries(G4Event* event) override;

 private:
  int fIonZ;
  int fIonA;
  G4ParticleGun fGun;
  G4ParticleDefinition* fIon = nullptr;
};