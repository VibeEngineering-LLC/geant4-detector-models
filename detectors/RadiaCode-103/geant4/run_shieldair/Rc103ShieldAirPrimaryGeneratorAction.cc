#include "Rc103ShieldAirPrimaryGeneratorAction.hh"

#include "Rc103FieldDetectorConstruction.hh"

#include "G4Event.hh"
#include "G4IonTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <cstdio>
#include <cstdlib>

Rc103ShieldAirPrimaryGeneratorAction::Rc103ShieldAirPrimaryGeneratorAction(
    int ionZ, int ionA)
    : fIonZ(ionZ), fIonA(ionA), fGun(1) {}

void Rc103ShieldAirPrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  using DC = Rc103FieldDetectorConstruction;
  // Тот же бокс, что и полость домика в Rc103FieldDetectorConstruction
  // (см. BuildLeadShield): x,y симметричны, z от верха дна до открытого среза.
  const double hx = 0.5 * DC::kShieldCavityXMm;
  const double hy = 0.5 * DC::kShieldCavityYMm;
  const double zLo = -0.5 * DC::kShieldOuterZMm + DC::kShieldPbMm;
  const double zHi = 0.5 * DC::kShieldOuterZMm;

  const double x = (2.0 * G4UniformRand() - 1.0) * hx;
  const double y = (2.0 * G4UniformRand() - 1.0) * hy;
  const double z = zLo + G4UniformRand() * (zHi - zLo);
  fGun.SetParticlePosition(G4ThreeVector(x * mm, y * mm, z * mm));

  if (!fIon) {
    fIon = G4IonTable::GetIonTable()->GetIon(fIonZ, fIonA, 0.0);
    if (!fIon) {
      std::fprintf(stderr,
                   "Rc103ShieldAirPrimaryGeneratorAction: FATAL ион Z=%d A=%d "
                   "не найден в таблице ионов.\n",
                   fIonZ, fIonA);
      std::abort();
    }
  }
  fGun.SetParticleDefinition(fIon);
  fGun.SetParticleEnergy(0.0);
  fGun.SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));  // покоится
  fGun.GeneratePrimaryVertex(event);
}