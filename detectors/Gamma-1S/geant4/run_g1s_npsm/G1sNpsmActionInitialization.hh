#pragma once

#include <G4VUserActionInitialization.hh>
#include <string>

class NpsmLightYield;

class G1sNpsmActionInitialization : public G4VUserActionInitialization {
public:
    G1sNpsmActionInitialization(std::string outCsv,
                                double energyKeV,
                                long long nEventsRequested,
                                long seed,
                                const NpsmLightYield* lightYield);

    void BuildForMaster() const override;
    void Build() const override;

private:
    std::string fOutCsv;
    double fEnergyKeV;
    long long fNEventsRequested;
    long fSeed;
    const NpsmLightYield* fLightYield;
};
