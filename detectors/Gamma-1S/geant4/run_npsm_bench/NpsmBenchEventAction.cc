#include "NpsmBenchEventAction.hh"
#include "NpsmBenchRunAction.hh"

#include <cstdio>
#include <cstdlib>

NpsmBenchEventAction::NpsmBenchEventAction(NpsmBenchRunAction* runAction)
    : fNCompt(0), fNRayl(0), fPhotAbsorbed(false), fConv(false), fEscaped(false), fEdepMeV(0.0), fEdepLightMeV(0.0),
      fRunAction(runAction) {
    // Отказ громкий, а не тихий: без RunAction прогон не имеет смысла, и лучше
    // не стартовать вовсе, чем отработать 10 000 событий с нулевой статистикой.
    if (!fRunAction) {
        std::fprintf(stderr, "NpsmBenchEventAction: FATAL runAction == nullptr\n");
        std::abort();
    }
}

void NpsmBenchEventAction::BeginOfEventAction(const G4Event*) {
    fNCompt = 0;
    fNRayl = 0;
    fPhotAbsorbed = false;
    fConv = false;
    fEscaped = false;
    fEdepMeV = 0.0;
    fEdepLightMeV = 0.0;
}

void NpsmBenchEventAction::EndOfEventAction(const G4Event*) {
    // Проверки на nullptr здесь намеренно нет: он исключён в конструкторе, а
    // тихая проверка в этом месте и была причиной нулевой статистики.
    fRunAction->RecordEvent(fNCompt, fNRayl, fPhotAbsorbed, fConv, fEscaped, fEdepMeV, fEdepLightMeV);
}
