#pragma once

#include "G4VUserActionInitialization.hh"
#include <string>

class NpsmLightYield;
class NpsmBenchRunAction;
class NpsmBenchPrimaryGeneratorAction;
class NpsmBenchEventAction;
class NpsmBenchSteppingAction;

class NpsmBenchActionInitialization : public G4VUserActionInitialization {
public:
    NpsmBenchActionInitialization(std::string outCsv,
                                  double energyKeV,
                                  long long nEventsRequested,
                                  long seed,
                                  const NpsmLightYield* lightYield,
                                  double crystalHalfZmm);

    void BuildForMaster() const override;
    void Build() const override;

private:
    std::string fOutCsv;
    double fEnergyKeV;
    long long fNEventsRequested;
    long fSeed;
    const NpsmLightYield* fLightYield;
    double fCrystalHalfZmm;
};
