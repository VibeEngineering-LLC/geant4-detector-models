#pragma once
#include "parma_source.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleTable.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include <string>
#include <stdexcept>

class ParmaPrimary : public G4VUserPrimaryGeneratorAction {
 public:
  // src: не владеет, должен жить дольше этого объекта; particle_name: имя частицы в Geant4 (например, "neutron","proton","mu+","mu-","e-","e+","gamma");
  // radius_cm: радиус цели в см, передается в ParmaSource::Sample
  ParmaPrimary(const ParmaSource* src, const std::string& particle_name, double radius_cm)
    : fSrc(src), fGun(new G4ParticleGun(1)), fRadiusCm(radius_cm) {
    // Ищем частицу по имени
    G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle(particle_name);
    if (!particle) {
      throw std::runtime_error("неизвестная частица: " + particle_name);
    }
    fGun->SetParticleDefinition(particle);
  }

  ~ParmaPrimary() override { delete fGun; }

  void GeneratePrimaries(G4Event* event) override {
    // Получаем событие от источника
    ParmaSource::Event ev = fSrc->Sample(fRadiusCm, []{ return G4UniformRand(); });
    
    // Устанавливаем энергию
    fGun->SetParticleEnergy(ev.e_MeV * MeV);
    
    // Устанавливаем позицию
    fGun->SetParticlePosition(G4ThreeVector(ev.x_cm, ev.y_cm, ev.z_cm) * cm);
    
    // Устанавливаем направление импульса
    fGun->SetParticleMomentumDirection(G4ThreeVector(ev.u, ev.v, ev.w));
    
    // Генерируем первичную вершину
    fGun->GeneratePrimaryVertex(event);
  }

 private:
  const ParmaSource* fSrc;
  G4ParticleGun* fGun;
  double fRadiusCm;  // радиус целевой сферы, см
};
