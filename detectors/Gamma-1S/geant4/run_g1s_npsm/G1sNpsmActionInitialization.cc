#include "G1sNpsmActionInitialization.hh"
#include "NpsmBenchRunAction.hh"
#include "G1sNpsmPrimaryGeneratorAction.hh"
#include "NpsmBenchEventAction.hh"
#include "NpsmBenchSteppingAction.hh"
#include <cstdio>
#include <cstdlib>

G1sNpsmActionInitialization::G1sNpsmActionInitialization(
    std::string outCsv,
    double energyKeV,
    long long nEventsRequested,
    long seed,
    const NpsmLightYield* lightYield)
    : fOutCsv(std::move(outCsv)),
      fEnergyKeV(energyKeV),
      fNEventsRequested(nEventsRequested),
      fSeed(seed),
      fLightYield(lightYield) {
    if (fLightYield == nullptr) {
        std::fprintf(stderr, "G1sNpsmActionInitialization: FATAL lightYield == nullptr\n");
        std::abort();
    }
}

void G1sNpsmActionInitialization::BuildForMaster() const {
    auto* runAction = new NpsmBenchRunAction(fOutCsv, fEnergyKeV, fNEventsRequested, fSeed, fLightYield);
    SetUserAction(runAction);
}

void G1sNpsmActionInitialization::Build() const {
    // Создаем RunAction
    auto* runAction = new NpsmBenchRunAction(fOutCsv, fEnergyKeV, fNEventsRequested, fSeed, fLightYield);
    SetUserAction(runAction);

    // Создаем PrimaryGeneratorAction
    auto* generatorAction = new G1sNpsmPrimaryGeneratorAction();
    SetUserAction(generatorAction);

    // Создаем EventAction
    auto* eventAction = new NpsmBenchEventAction(runAction);
    SetUserAction(eventAction);

    // Создаем SteppingAction
    auto* steppingAction = new NpsmBenchSteppingAction(eventAction, fLightYield);
    SetUserAction(steppingAction);
}
