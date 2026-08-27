#include "Rc103MuonPrimaryGeneratorAction.hh"

#include "Rc103MuonSpectrum.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

Rc103MuonPrimaryGeneratorAction::Rc103MuonPrimaryGeneratorAction(
    const Rc103MuonSpectrum* spectrum, double rDiskMm)
    : fSpectrum(spectrum), fRDiskMm(rDiskMm), fGun(1) {
  auto* mu = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
  if (!mu) {
    std::fprintf(stderr,
                 "Rc103MuonPrimaryGeneratorAction: FATAL particle 'mu-' not "
                 "found in particle table.\n");
    std::abort();
  }
  fGun.SetParticleDefinition(mu);
}

void Rc103MuonPrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  // 1) точка равномерно ПО ПЛОЩАДИ диска
  const double r = fRDiskMm * std::sqrt(G4UniformRand()) * mm;
  const double ph = twopi * G4UniformRand();
  const G4ThreeVector pos(r * std::cos(ph), r * std::sin(ph), kZDiskMm * mm);

  // 2) направление вниз; p(cosT) ~ cosT^3  =>  cosT = U^(1/4).
  //    Именно 0.25, НЕ std::cbrt — см. развёрнутое обоснование в .hh.
  const double cosT = std::pow(G4UniformRand(), 0.25);  // p(cosT) ~ cosT^3
  const double sinT = std::sqrt(std::max(0.0, 1.0 - cosT * cosT));
  const double phd = twopi * G4UniformRand();
  const G4ThreeVector dir(sinT * std::cos(phd), sinT * std::sin(phd), -cosT);

  fGun.SetParticlePosition(pos);
  fGun.SetParticleMomentumDirection(dir);
  fGun.SetParticleEnergy(fSpectrum->SampleEnergyGeV() * GeV);
  fGun.GeneratePrimaryVertex(event);
}
