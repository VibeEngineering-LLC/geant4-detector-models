#pragma once

#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include "NpsmBenchRun.hh"
#include <array>
#include <string>

class NpsmLightYield;

class NpsmBenchRunAction : public G4UserRunAction {
public:
    explicit NpsmBenchRunAction(std::string outCsv, double energyKeV, long long nEventsRequested, long seed, const NpsmLightYield* lightYield);
    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

    // Накопители переехали в NpsmBenchRun (наследник G4Run) — 07.09.2026,
    // линия И. В однопоточном режиме поведение прежнее: те же величины,
    // просто живут в объекте прогона, а не в полях этого класса. В MT
    // Geant4 создаёт такой объект на каждый поток и сливает через Merge;
    // без этого мастер собирал бы неполный результат, а прогон выглядел бы
    // успешным — ровно тот класс отказа, который мы вычищаем (W-067).
    G4Run* GenerateRun() override;

    void RecordEvent(int nCompt, int nRayl, bool phot, bool conv, bool escaped, double edepMeV, double edepLightMeV);

    // Постановка первички ДЛЯ ШАПКИ CSV. Заполняет main того проекта, который
    // подменил генератор стенда своим (run_g1s_npsm, режим primary=ion):
    // иначе шапка печатала `particle,gamma` и для ионного прогона — файлы
    // corr_gamma=0 и corr_gamma=1 были неразличимы по содержимому (07.09.2026).
    // Умолчания воспроизводят стенд: primary=gamma, ион не задан.
    static std::string gPrimaryKind;   // "gamma" | "ion"
    static int gIonZ;                  // только при gPrimaryKind == "ion"
    static int gIonA;

    // Режим счёта: 0 — однопоточный, >0 — число рабочих потоков. Обязан
    // лежать В ФАЙЛЕ: Serial и MT дают РАЗНЫЕ выборки одних и тех же
    // событий (разброс отдельной пары до 0,6 %, систематики нет — проверено
    // на 10 зёрнах). Сравнивать между собой можно только прогоны одного
    // режима, а для этого режим должен быть виден в результате, а не
    // помниться. Тот же класс, что W-068.
    static int gThreads;

private:
    std::string fOutCsv;
    double fEnergyKeV;
    long long fNEventsRequested;
    long fSeed;
    const NpsmLightYield* fLightYield;

    // Полей-накопителей здесь БОЛЬШЕ НЕТ: все они в NpsmBenchRun. Именно
    // хранение статистики в объекте действия делало код непригодным для MT —
    // каждый поток копил бы своё, и мастер об этом не узнал бы. Смысл каждой
    // величины (включая оговорку про свет в событиях полного поглощения:
    // непропорциональность разбивает пик на группы по числу комптоновских
    // рассеяний, поэтому σ берётся разбросом, а не подгонкой гауссианы)
    // описан в NpsmBenchRun.hh — в одном месте, чтобы не разошлось.

    void WriteCSV(const NpsmBenchRun& run);
};
