#include "G1sNpsmPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4IonTable.hh"
#include <cstdio>
#include <cstdlib>

double G1sNpsmPrimaryGeneratorAction::gSourceZmm = 91.0;
double G1sNpsmPrimaryGeneratorAction::gEnergyKeV = 661.657;
std::string G1sNpsmPrimaryGeneratorAction::gPrimary = "gamma";
int G1sNpsmPrimaryGeneratorAction::gIonZ = 27;   // Co
int G1sNpsmPrimaryGeneratorAction::gIonA = 60;   // Co-60

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
    // Задаем позицию источника
    G4ThreeVector position(0., 0., gSourceZmm * mm);
    fGun.SetParticlePosition(position);

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
