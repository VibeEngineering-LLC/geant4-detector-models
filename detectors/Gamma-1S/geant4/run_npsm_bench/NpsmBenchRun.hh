#pragma once

#include "G4Run.hh"
#include <array>
#include <cstdlib>
#include <cstdio>
#include <cmath>

class NpsmBenchRun : public G4Run {
public:
    static constexpr int kMaxCompt = 32;
    static constexpr int kNBins = 3200;
    static constexpr double kBinKeV = 1.0;

    // Накопители событий
    std::array<long long, kMaxCompt + 1> fHistAbsorbed{};
    std::array<long long, kMaxCompt + 1> fHistAll{};
    std::array<long long, kNBins> fSpecEdep{};
    std::array<long long, kNBins> fSpecLight{};

    long long fNEvents = 0;
    long long fNAbsorbed = 0;
    long long fNConv = 0;
    long long fNEscaped = 0;
    long long fNOther = 0;
    long long fNFullEdep = 0;
    long long fNOverflowEdep = 0;
    long long fNOverflowLight = 0;
    long long fNWithEdep = 0;

    double fSumComptAbsorbed = 0.0;
    double fSumCompt2Absorbed = 0.0;
    double fSumEdepMeV = 0.0;
    double fSumLightMeV = 0.0;
    double fSumLight2MeV2 = 0.0;
    double fSumLightFullMeV = 0.0;
    double fSumLight2FullMeV2 = 0.0;

    // Оптические фотоны, рождённые ШТАТНЫМ G4Scintillation (режим
    // light=scint, линия Г). Считает NpsmBenchStackingAction, он же их
    // убивает: транспорт не нужен, нужно только число. Поле живёт здесь,
    // а не в действии, потому что действие своё у каждого потока и его
    // счётчик не попал бы в слияние.
    long long fNScintPhotons = 0;
    void AddScintPhoton() { ++fNScintPhotons; }

    NpsmBenchRun() = default;
    ~NpsmBenchRun() override = default;

    void RecordEvent(int nCompt, int nRayl, bool phot, bool conv, bool escaped,
                     double edepMeV, double edepLightMeV, double energyKeV);

    void Merge(const G4Run* run) override;
};
