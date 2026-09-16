#pragma once

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include <string>
#include <utility>
#include <vector>

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

    // Разбиение события на импульсы по глобальному времени (16.09.2026, P-033):
    // ion-шаблон проводит всю цепочку в одном событии, и кванты дочерних ядер
    // с T½ в минуты складывались как одновременные. Окно τ в секундах; <0 —
    // ключ pulse_window_s не задан, и путь кода СТРОГО прежний (буфер не
    // заполняется, в гистограмму идёт одна запись на событие).
    static double gPulseWindowS;
    static std::string PulseWindowStr();   // "none" или τ кратчайшей записью
    inline void AddDeposit(double tGlobal, double edepMeV, double lightMeV) {
        fDeposits.push_back({tGlobal, edepMeV, lightMeV});
    }

private:
    struct Deposit { double t; double e; double l; };
    std::vector<Deposit> fDeposits;                     // буфер события, свой у потока
    std::vector<std::pair<double, double>> fPulses;     // (edep, light) по импульсам
    int fNCompt;
    int fNRayl;
    bool fPhotAbsorbed;
    bool fConv;
    bool fEscaped;
    double fEdepMeV;
    double fEdepLightMeV;

    NpsmBenchRunAction* fRunAction;
};
