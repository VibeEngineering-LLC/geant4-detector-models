#pragma once

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include <string>

class NpsmBenchPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    // Постановка задаётся ДО построения объекта (main разбирает ключи раньше),
    // поэтому поля статические: "gamma" либо "e-"; iso-поток либо пучок из центра.
    static std::string gParticle;
    static bool gPencilBeam;

    explicit NpsmBenchPrimaryGeneratorAction(double energyKeV);
    void GeneratePrimaries(G4Event* anEvent) override;

private:
    G4ParticleGun fGun;
    double fEnergyKeV;
};
