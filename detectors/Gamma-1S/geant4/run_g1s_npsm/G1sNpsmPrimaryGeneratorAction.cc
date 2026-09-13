#include "G1sNpsmPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4IonTable.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Exception.hh"
#include <cstdio>
#include <cstdlib>
#include <cmath>

double G1sNpsmPrimaryGeneratorAction::gSourceZmm = 91.0;
double G1sNpsmPrimaryGeneratorAction::gEnergyKeV = 661.657;
std::string G1sNpsmPrimaryGeneratorAction::gPrimary = "gamma";
int G1sNpsmPrimaryGeneratorAction::gIonZ = 27;   // Co
int G1sNpsmPrimaryGeneratorAction::gIonA = 60;   // Co-60
std::string G1sNpsmPrimaryGeneratorAction::gSourceMode = "point";
double G1sNpsmPrimaryGeneratorAction::gSrcZFrac = 1.0;
double G1sNpsmPrimaryGeneratorAction::gSrcRFrac = 1.0;

G4ThreeVector G1sNpsmPrimaryGeneratorAction::SamplePointInSample() {
    if (!fSampleSolid) {
        G4LogicalVolume* lv =
            G4LogicalVolumeStore::GetInstance()->GetVolume("Sample", false);
        if (!lv) {
            // Громкий отказ, а не тихий возврат нуля: молча вылетающая из
            // центра мира первичка дала бы правдоподобный, но неверный спектр.
            G4Exception("G1sNpsmPrimaryGeneratorAction::SamplePointInSample",
                        "NoSample", FatalException,
                        "src=sample, но логического объёма Sample нет: "
                        "сосуд не построен (vessel=none?)");
        }
        fSampleSolid = lv->GetSolid();
        fSampleSolid->BoundingLimits(fSampleMin, fSampleMax);
    }
    // Отбор с отклонением: тело пробы — G4Polycone (кольцо вокруг колодца плюс
    // слой над ним), аналитического равномерного розыгрыша для него нет, а
    // Inside() у солида точен. Счётчик срыва защищает от бесконечного цикла.
    for (int i = 0; i < 10000; ++i) {
        const G4ThreeVector p(
            fSampleMin.x() + (fSampleMax.x() - fSampleMin.x()) * G4UniformRand(),
            fSampleMin.y() + (fSampleMax.y() - fSampleMin.y()) * G4UniformRand(),
            fSampleMin.z() + (fSampleMax.z() - fSampleMin.z()) * G4UniformRand());
        if (fSampleSolid->Inside(p) != kInside) continue;
        // Ограничение области: снизу по высоте и снаружи по радиусу.
        if (gSrcZFrac < 1.0) {
            const double zTop = fSampleMin.z()
                + (fSampleMax.z() - fSampleMin.z()) * gSrcZFrac;
            if (p.z() > zTop) continue;
        }
        if (gSrcRFrac < 1.0) {
            const double rMax = std::max(std::abs(fSampleMax.x()),
                                         std::abs(fSampleMin.x()));
            const double rMin = rMax * (1.0 - gSrcRFrac);
            if (p.perp() < rMin) continue;
        }
        return p;
    }
    G4Exception("G1sNpsmPrimaryGeneratorAction::SamplePointInSample",
                "RejectionFailed", FatalException,
                "10000 попыток отбора не дали точки внутри пробы");
    return G4ThreeVector();
}

G1sNpsmPrimaryGeneratorAction::G1sNpsmPrimaryGeneratorAction()
    : fGun(1) {
    // В режиме "ion" частицу здесь НЕ задаём: таблица ионов доступна только
    // после инициализации ядра, и обращение к ней из конструктора действия
    // вернёт нуль. Ион создаётся при первом событии, см. GeneratePrimaries.
    if (gPrimary != "ion") {
        G4ParticleDefinition* gamma = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
        fGun.SetParticleDefinition(gamma);
    }
}

void G1sNpsmPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    // Позиция источника: точка на оси либо равномерно по объёму пробы.
    const G4ThreeVector position = (gSourceMode == "sample")
        ? SamplePointInSample()
        : G4ThreeVector(0., 0., gSourceZmm * mm);
    fGun.SetParticlePosition(position);

    if (gPrimary == "eplus" || gPrimary == "eplus_gamma") {
        // β⁺ (Sc-44): позитрон покоится в точке рождения и аннигилирует сам —
        // две встречные 511 кэВ; при eplus_gamma из той же точки в том же
        // событии — изотропный квант gEnergyKeV. Пробег позитрона и аннигиляция
        // на лету не учитываются — приближение, названное в отчёте.
        static G4ParticleDefinition* ep = G4ParticleTable::GetParticleTable()->FindParticle("e+");
        static G4ParticleDefinition* gm = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
        fGun.SetParticleDefinition(ep);
        fGun.SetParticleEnergy(0.0);
        fGun.SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
        fGun.GeneratePrimaryVertex(anEvent);
        if (gPrimary == "eplus") return;
        fGun.SetParticleDefinition(gm);
        fGun.SetParticleEnergy(gEnergyKeV * keV);
        const G4double ct = 1.0 - 2.0 * G4UniformRand(), ph = 2.0 * M_PI * G4UniformRand();
        const G4double st = std::sqrt(1.0 - ct * ct);
        fGun.SetParticleMomentumDirection(G4ThreeVector(st * std::cos(ph), st * std::sin(ph), ct));
        fGun.GeneratePrimaryVertex(anEvent);
        return;
    }

    if (gPrimary == "ion") {
        // Ядро создаётся один раз, при первом событии: таблица ионов готова
        // только после Initialize. Энергия нулевая — ядро покоится, распад
        // разыгрывает G4RadioactiveDecay, он же испускает весь каскад.
        if (!fIon) {
            fIon = G4IonTable::GetIonTable()->GetIon(gIonZ, gIonA, 0.0);
            if (!fIon) {
                std::fprintf(stderr, "Не найден ион Z=%d A=%d\n", gIonZ, gIonA);
                std::abort();
            }
            fGun.SetParticleDefinition(fIon);
            fGun.SetParticleCharge(0.0);
        }
        fGun.SetParticleEnergy(0.0);
    } else {
        // Задаем энергию
        fGun.SetParticleEnergy(gEnergyKeV * keV);
    }

    // Генерируем изотропное направление
    // Разыгрываем косинус полярного угла для равномерного распределения по телесному углу
    G4double cosTheta = 1.0 - 2.0 * G4UniformRand();
    G4double phi = 2.0 * M_PI * G4UniformRand();
    G4double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);

    G4ThreeVector direction(sinTheta * std::cos(phi),
                            sinTheta * std::sin(phi),
                            cosTheta);
    fGun.SetParticleMomentumDirection(direction);

    // Генерируем первичную вершину
    fGun.GeneratePrimaryVertex(anEvent);
}
