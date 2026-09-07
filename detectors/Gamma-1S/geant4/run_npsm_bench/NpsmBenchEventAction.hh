#pragma once

#include "G4UserEventAction.hh"
#include "G4Event.hh"

class NpsmBenchRunAction;

class NpsmBenchEventAction : public G4UserEventAction {
public:
    // Указатель на RunAction ОБЯЗАТЕЛЕН. В первой версии конструктор был без
    // аргументов, fRunAction оставался nullptr, а EndOfEventAction молча
    // пропускал каждое событие через `if (fRunAction)` — прогон 10 000 событий
    // отработал успешно и записал CSV со всеми нулями (05.09).
    explicit NpsmBenchEventAction(NpsmBenchRunAction* runAction);
    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

    // Мутаторы
    inline void AddCompt() { fNCompt++; }
    inline void AddRayl() { fNRayl++; }
    inline void SetPhotAbsorbed() { fPhotAbsorbed = true; }
    inline void SetConv() { fConv = true; }
    inline void SetEscaped() { fEscaped = true; }
    inline void AddEdep(double edepMeV) { fEdepMeV += edepMeV; }
    inline void AddEdepLight(double v) { fEdepLightMeV += v; }

private:
    int fNCompt;
    int fNRayl;
    bool fPhotAbsorbed;
    bool fConv;
    bool fEscaped;
    double fEdepMeV;
    double fEdepLightMeV;

    NpsmBenchRunAction* fRunAction;
};
