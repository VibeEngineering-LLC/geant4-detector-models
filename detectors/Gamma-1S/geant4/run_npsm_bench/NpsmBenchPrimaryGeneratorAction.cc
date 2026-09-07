#include "NpsmBenchPrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4Event.hh"

#include "NpsmBenchDetectorConstruction.hh"

#include <cstdio>
#include <cstdlib>

// Статические поля
std::string NpsmBenchPrimaryGeneratorAction::gParticle = "gamma";
bool NpsmBenchPrimaryGeneratorAction::gPencilBeam = false;

NpsmBenchPrimaryGeneratorAction::NpsmBenchPrimaryGeneratorAction(double energyKeV)
    : fGun(1), fEnergyKeV(energyKeV) {
    G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle(gParticle);
    if (!particle) {
        std::fprintf(stderr, "Ошибка: неизвестная частица '%s'\n", gParticle.c_str());
        std::abort();
    }
    fGun.SetParticleDefinition(particle);
}

void NpsmBenchPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    if (gPencilBeam) {
        // Карандашный пучок: центр кристалла, направление вверх
        const G4ThreeVector pos(0., 0., 0.);
        const G4ThreeVector dir(0., 0., 1.);

        fGun.SetParticlePosition(pos);
        fGun.SetParticleMomentumDirection(dir);
        fGun.SetParticleEnergy(fEnergyKeV * keV);

        fGun.GeneratePrimaryVertex(anEvent);
        return;
    }

    // Известный изотропный розыгрыш по сфере
    const double R_SRC = NpsmBenchDetectorConstruction::SourceRadiusMm() * mm;

    // 1) точка равномерно по сфере радиуса R_SRC
    const double cosT = 2.0 * G4UniformRand() - 1.0;
    const double sinT = std::sqrt(1.0 - cosT * cosT);
    const double phi  = twopi * G4UniformRand();
    const G4ThreeVector n(sinT * std::cos(phi), sinT * std::sin(phi), cosT); // наружная нормаль
    const G4ThreeVector pos = R_SRC * n;

    // 2) направление ВНУТРЬ, косинусное относительно -n
    const G4ThreeVector e3 = -n;
    G4ThreeVector e1 = e3.orthogonal().unit();
    const G4ThreeVector e2 = e3.cross(e1).unit();
    const double ct  = std::sqrt(G4UniformRand());
    const double st  = std::sqrt(1.0 - ct * ct);
    const double psi = twopi * G4UniformRand();
    const G4ThreeVector dir = st * std::cos(psi) * e1 + st * std::sin(psi) * e2 + ct * e3;

    fGun.SetParticlePosition(pos);
    fGun.SetParticleMomentumDirection(dir);
    fGun.SetParticleEnergy(fEnergyKeV * keV);

    fGun.GeneratePrimaryVertex(anEvent);
}
