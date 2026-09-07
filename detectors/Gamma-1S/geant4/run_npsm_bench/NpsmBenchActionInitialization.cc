#include "NpsmBenchActionInitialization.hh"
#include "NpsmBenchRunAction.hh"
#include "NpsmBenchPrimaryGeneratorAction.hh"
#include "NpsmBenchEventAction.hh"
#include "NpsmBenchSteppingAction.hh"
#include "NpsmBenchStackingAction.hh"

#include <cstdio>
#include <cstdlib>

// Предварительное объявление класса
class NpsmLightYield;

NpsmBenchActionInitialization::NpsmBenchActionInitialization(
    std::string outCsv,
    double energyKeV,
    long long nEventsRequested,
    long seed,
    const NpsmLightYield* lightYield,
    double crystalHalfZmm)
    : fOutCsv(std::move(outCsv)),
      fEnergyKeV(energyKeV),
      fNEventsRequested(nEventsRequested),
      fSeed(seed),
      fLightYield(lightYield),
      fCrystalHalfZmm(crystalHalfZmm) {
    if (fLightYield == nullptr) {
        std::fprintf(stderr, "NpsmBenchActionInitialization: FATAL lightYield == nullptr\n");
        std::abort();
    }
}

void NpsmBenchActionInitialization::BuildForMaster() const {
    auto* runAction = new NpsmBenchRunAction(fOutCsv, fEnergyKeV, fNEventsRequested, fSeed,
                                            fLightYield);
    SetUserAction(runAction);
}

void NpsmBenchActionInitialization::Build() const {
    // Создаем и регистрируем действия в правильном порядке
    auto* runAction = new NpsmBenchRunAction(fOutCsv, fEnergyKeV, fNEventsRequested, fSeed,
                                            fLightYield);
    SetUserAction(runAction);

    auto* primaryGeneratorAction = new NpsmBenchPrimaryGeneratorAction(fEnergyKeV);
    SetUserAction(primaryGeneratorAction);

    auto* eventAction = new NpsmBenchEventAction(runAction);
    SetUserAction(eventAction);

    auto* steppingAction = new NpsmBenchSteppingAction(eventAction, fLightYield);
    SetUserAction(steppingAction);

    // Счётчик оптических фотонов (линия Г). Регистрируется всегда: без
    // сцинтилляции оптических фотонов просто не возникает, и действие ни
    // разу не срабатывает. Условная регистрация добавила бы ветку, которую
    // нечем проверить.
    SetUserAction(new NpsmBenchStackingAction());
}
