#include "NpsmBenchSteppingAction.hh"
#include "NpsmBenchEventAction.hh"
#include "NpsmBenchDetectorConstruction.hh"
#include "NpsmLightYield.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHandle.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include "G4Material.hh"
#include <cmath>
#include "G4ParticleDefinition.hh"

// Умолчание — 1 кэВ, порог из Algorithm S1 первоисточника.
double NpsmBenchSteppingAction::gLocalEnergyKeV = 1.0;
bool NpsmBenchSteppingAction::gPrintDedxTable = false;
bool NpsmBenchSteppingAction::gRestrictedDedx = false;

// Нулевой указатель = кристалл спрашивается у геометрии стенда (прежнее
// поведение). Задаётся из main проекта с полной геометрией.
G4LogicalVolume* NpsmBenchSteppingAction::gCrystalLV = nullptr;

void NpsmBenchSteppingAction::SetCrystalLogicalVolume(G4LogicalVolume* lv) {
    gCrystalLV = lv;
}

NpsmBenchSteppingAction::NpsmBenchSteppingAction(NpsmBenchEventAction* eventAction,
                                                   const NpsmLightYield* lightYield)
    : fEventAction(eventAction), fLightYield(lightYield), fEmCalc() {
    if (!fEventAction) {
        std::fprintf(stderr, "Ошибка: fEventAction == nullptr\n");
        std::abort();
    }
    if (!fLightYield) {
        std::fprintf(stderr, "Ошибка: fLightYield == nullptr\n");
        std::abort();
    }
}

void NpsmBenchSteppingAction::UserSteppingAction(const G4Step* step) {
    // Получаем предшествующий объем
    const G4VPhysicalVolume* preVol = step->GetPreStepPoint()->GetPhysicalVolume();
    if (!preVol) return;

    // Получаем логический объем кристалла: заданный извне (полная геометрия)
    // либо, если не задан, из геометрии стенда — прежнее поведение.
    G4LogicalVolume* crystalLV = gCrystalLV
        ? gCrystalLV
        : NpsmBenchDetectorConstruction::GetCrystalLogicalVolume();
    if (!crystalLV) return;

    // Одноразовая печать таблицы тормозной способности при ПЕРВОМ шаге в
    // кристалле: здесь таблицы заведомо построены, тогда как отдельный режим
    // в main падал даже после холостого BeamOn(0). Нужна для эталона-2 —
    // кривой отклика, построенной на ТОЙ ЖЕ S(E), что и стенд.
    if (gPrintDedxTable && preVol->GetLogicalVolume() == crystalLV) {
        gPrintDedxTable = false;
        const G4Material* mat = step->GetPreStepPoint()->GetMaterial();
        const G4ParticleDefinition* e = G4Electron::Definition();
        std::printf("DEDX_TABLE_BEGIN,material=%s\n", mat->GetName().c_str());
        for (int i = 0; i <= 400; ++i) {
            const double ek = 0.1 * std::pow(10.0, i * 4.5 / 400.0);  // 0,1 кэВ … 3,16 МэВ
            const double d = fEmCalc.ComputeElectronicDEDX(ek * keV, e, mat, DBL_MAX);
            std::printf("DEDX,%.6f,%.6f\n", ek, d / (MeV / cm));
        }
        std::printf("DEDX_TABLE_END\n");
    }

    // Энергия, депонированная в кристалле
    double edep = step->GetTotalEnergyDeposit();
    // У G4StepPoint нет GetLogicalVolume(); логический объём берётся через
    // физический объём пре-шага (preVol получен выше).
    if (edep > 0 && preVol->GetLogicalVolume() == crystalLV) {
        double edepMeV = edep / MeV;
        double w = 1.0;

        const G4Track* track = step->GetTrack();
        const G4ParticleDefinition* particleDef = track->GetDefinition();

        if (fLightYield->IsEnabled() &&
            (particleDef == G4Electron::Definition() ||
             particleDef == G4Positron::Definition())) {
            // Вычисление S — удельные потери энергии, МэВ/см.
            //
            // Деление на две ветви следует Algorithm S1 первоисточника, где
            // событие относится к "local" по ЭНЕРГИИ (порог транспорта 1 кэВ),
            // а не только по длине шага. Первая редакция делила лишь по длине
            // (10⁻⁴ мм); признаки не эквивалентны, и приёмка этапа 2 это
            // показала (ACCEPTANCE-stage2.md). Геометрический признак оставлен
            // как второй: у локального депозита длина шага стремится к нулю,
            // и edep/stepLength дало бы деление на ноль.
            double S;
            const double stepLength = step->GetStepLength();
            const G4double Ek = step->GetPreStepPoint()->GetKineticEnergy();
            const bool isLocal = (stepLength <= kMinStepMm * mm) ||
                                 (Ek < gLocalEnergyKeV * keV);
            if (!isLocal) {
                // continuous: потери размазаны вдоль криволинейного пути
                S = edepMeV / (stepLength / cm);
            } else {
                // local: энергия оставлена на месте, S берётся из таблицы
                const G4Material* material = step->GetPreStepPoint()->GetMaterial();
                // Д-6: по умолчанию ПОЛНАЯ тормозная способность. Ограниченная
                // (GetDEDX) включается ключом только для измерения цены этого
                // отступления от разбора первоисточника.
                const G4double dedx = gRestrictedDedx
                    ? fEmCalc.GetDEDX(Ek, particleDef, material)
                    : fEmCalc.ComputeElectronicDEDX(Ek, particleDef, material, DBL_MAX);
                S = dedx / (MeV / cm);
            }
            w = fLightYield->Weight(S);
        }

        fEventAction->AddEdep(edepMeV);
        fEventAction->AddEdepLight(edepMeV * w);
    }

    // Только для первичного трека (parent ID == 0)
    if (step->GetTrack()->GetParentID() != 0) return;

    // Получаем процесс, определивший шаг
    const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
    std::string name = "";
    if (process) {
        name = process->GetProcessName();
    }

    // Проверяем, находится ли постшаг внутри кристалла
    const G4VPhysicalVolume* postVol = step->GetPostStepPoint()->GetPhysicalVolume();
    bool inCrystal = (postVol && postVol->GetLogicalVolume() == crystalLV);

    if (inCrystal) {
        // Если шаг завершен внутри кристалла, обновляем статистику
        if (name == "compt") {
            fEventAction->AddCompt();
        } else if (name == "Rayl") {
            fEventAction->AddRayl();  // Почему Rayleigh считается отдельно: это упругий процесс, меняющий только направление, а не количество Compton-взаимодействий
        } else if (name == "phot") {
            fEventAction->SetPhotAbsorbed();
        } else if (name == "conv") {
            fEventAction->SetConv();
        }
    } else if (!postVol) {
        // Если трек покинул мир, помечаем как ушедший
        fEventAction->SetEscaped();
    }
}
