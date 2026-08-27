#include "Rc103FieldPrimaryGeneratorAction.hh"

#include "Rc103FieldFluxSpectrum.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <cmath>

Rc103FieldPrimaryGeneratorAction::Rc103FieldPrimaryGeneratorAction(
    const Rc103FieldFluxSpectrum* spectrum)
    : fSpectrum(spectrum), fGun(1) {
  fGun.SetParticleDefinition(
      G4ParticleTable::GetParticleTable()->FindParticle("gamma"));
}

void Rc103FieldPrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  const double R_SRC = kRSrcMm * mm;

  // 1) точка равномерно по сфере радиуса R_SRC
  const double cosT = 2.0 * G4UniformRand() - 1.0;
  const double sinT = std::sqrt(1.0 - cosT * cosT);
  const double phi = twopi * G4UniformRand();
  const G4ThreeVector n(sinT * std::cos(phi), sinT * std::sin(phi),
                        cosT);  // наружная нормаль
  const G4ThreeVector pos = R_SRC * n;

  // 2) направление ВНУТРЬ, косинусное относительно -n
  const G4ThreeVector e3 = -n;               // внутренняя нормаль
  G4ThreeVector e1 = e3.orthogonal().unit();  // любой перпендикуляр
  const G4ThreeVector e2 = e3.cross(e1).unit();
  const double ct = std::sqrt(G4UniformRand());  // cos отн. e3, косинусный закон
  const double st = std::sqrt(1.0 - ct * ct);
  const double psi = twopi * G4UniformRand();
  const G4ThreeVector dir =
      st * std::cos(psi) * e1 + st * std::sin(psi) * e2 + ct * e3;

  fGun.SetParticlePosition(pos);
  fGun.SetParticleMomentumDirection(dir);
  fGun.SetParticleEnergy(fSpectrum->SampleEnergyKeV() * keV);
  fGun.GeneratePrimaryVertex(event);
}
