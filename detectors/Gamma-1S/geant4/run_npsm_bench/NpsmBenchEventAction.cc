#include "NpsmBenchEventAction.hh"
#include "NpsmBenchRunAction.hh"

#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <charconv>
#include <cstdio>
#include <cstdlib>

double NpsmBenchEventAction::gPulseWindowS = -1.0;

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
    if (gPulseWindowS >= 0.0) fDeposits.clear();
}

void NpsmBenchEventAction::EndOfEventAction(const G4Event*) {
    // Проверки на nullptr здесь намеренно нет: он исключён в конструкторе, а
    // тихая проверка в этом месте и была причиной нулевой статистики.
    if (gPulseWindowS < 0.0) {
        fRunAction->RecordEvent(fNCompt, fNRayl, fPhotAbsorbed, fConv, fEscaped, fEdepMeV, fEdepLightMeV);
        return;
    }
    // Режим импульсов. Сортировка устойчивая: при равном времени сохраняется
    // порядок шагов. Новый импульс — когда t − t_начала > τ.
    std::stable_sort(fDeposits.begin(), fDeposits.end(),
                     [](const Deposit& a, const Deposit& b) { return a.t < b.t; });
    const double tau = gPulseWindowS * CLHEP::s;
    fPulses.clear();
    double t0 = 0.0;
    for (const Deposit& d : fDeposits) {
        if (fPulses.empty() || d.t - t0 > tau) {
            fPulses.emplace_back(0.0, 0.0);
            t0 = d.t;
        }
        fPulses.back().first += d.e;
        fPulses.back().second += d.l;
    }
    // Один импульс (или ни одного депозита) — это всё событие: в гистограмму
    // идут накопленные в порядке шагов суммы fEdepMeV/fEdepLightMeV, как без
    // ключа. Сумма, пересчитанная в порядке времени, могла бы отличаться в
    // последнем разряде и сдвинуть бин — тогда мутант τ=1e60 не совпал бы
    // побайтно с прогоном без ключа.
    if (fPulses.size() <= 1) {
        fRunAction->RecordEvent(fNCompt, fNRayl, fPhotAbsorbed, fConv, fEscaped, fEdepMeV, fEdepLightMeV);
    } else {
        fRunAction->RecordEvent(fNCompt, fNRayl, fPhotAbsorbed, fConv, fEscaped, fEdepMeV, fEdepLightMeV,
                                &fPulses);
    }
}

std::string NpsmBenchEventAction::PulseWindowStr() {
    if (gPulseWindowS < 0.0) return "none";
    char buf[64];
    const auto r = std::to_chars(buf, buf + sizeof(buf), gPulseWindowS);   // кратчайшая точная запись
    return std::string(buf, r.ptr);
}
