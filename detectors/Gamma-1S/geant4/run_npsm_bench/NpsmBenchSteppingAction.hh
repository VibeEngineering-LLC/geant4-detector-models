#pragma once

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "G4EmCalculator.hh"

class NpsmBenchEventAction;
class NpsmLightYield;

class NpsmBenchSteppingAction : public G4UserSteppingAction {
public:
    NpsmBenchSteppingAction(NpsmBenchEventAction* eventAction,
                            const NpsmLightYield* lightYield);
    void UserSteppingAction(const G4Step*) override;

    // Граница ветвей continuous/local, кэВ. Публичная: задаётся из main по
    // ключу ДО построения объекта, как и постановка генератора.
    static double gLocalEnergyKeV;

    // Печать таблицы dE/dx при первом шаге в кристалле (ключ dedx_table=1).
    // Умолчание false: обычные прогоны ничего лишнего не печатают.
    static bool gPrintDedxTable;

    // Измерение цены отступления Д-6 (ключ restricted_dedx=1): в локальной
    // ветви берётся ОГРАНИЧЕННАЯ порогом GetDEDX вместо полной. Умолчание
    // false — полная, как согласовано.
    static bool gRestrictedDedx;

    // Кристалл, в котором считаются депозиты. Введено 06.09.2026 (этап 5):
    // этот же класс работает и в стенде run_npsm_bench (вакуумный куб), и в
    // проекте run_g1s_npsm с полной геометрией Гамма-1С, где кристалл — лишь
    // один объём из многих (свинец, банка, отражатель, ФЭУ), и спрашивать его
    // у геометрии стенда нельзя. Нулевой указатель = ПРЕЖНЕЕ поведение:
    // кристалл берётся у NpsmBenchDetectorConstruction, поэтому прогоны стенда
    // воспроизводятся бит в бит.
    static void SetCrystalLogicalVolume(G4LogicalVolume* lv);
    static G4LogicalVolume* gCrystalLV;

private:
    NpsmBenchEventAction* fEventAction;
    const NpsmLightYield* fLightYield;
    G4EmCalculator fEmCalc;
    // Ниже этой длины шага деление edep/stepLength теряет смысл.
    static constexpr double kMinStepMm = 1e-4;
    // Энергетическая граница ветвей из Algorithm S1 первоисточника: событие
    // ниже порога транспорта относится к "local" и весится по таблице dE/dx.
    // У авторов это 1 кэВ (порог FLUKA); у нас перенос идёт до 0,1 кэВ, но
    // ГРАНИЦА ВЕТВЕЙ берётся авторская — меняется способ взвешивания, а не
    // физика переноса.
};
