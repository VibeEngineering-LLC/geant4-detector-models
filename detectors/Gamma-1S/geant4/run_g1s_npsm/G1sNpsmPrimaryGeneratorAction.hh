#ifndef G1S_NPSM_PRIMARY_GENERATOR_ACTION_HH
#define G1S_NPSM_PRIMARY_GENERATOR_ACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ThreeVector.hh"

#include <string>

class G1sNpsmPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    static double gSourceZmm;
    static double gEnergyKeV;

    // Что вылетает из источника. "gamma" — одиночный квант gEnergyKeV, как было
    // всегда. "ion" — покоящееся ядро (gIonZ, gIonA), которое распадается сам
    // Geant4: тогда каскад испускается целиком и в спектре появляется каскадное
    // суммирование, которого при одиночном кванте нет и быть не может.
    // Режим введён 07.09.2026: расчёты этапа 4 по Co-60 шли одиночным квантом
    // 1173,23 кэВ, а измерялся реальный источник, где каскад 1173+1332 есть
    // всегда — эта разница нигде не была учтена и даже не записана.
    static std::string gPrimary;
    static int gIonZ;
    static int gIonA;

public:
    G1sNpsmPrimaryGeneratorAction();
    ~G1sNpsmPrimaryGeneratorAction() override = default;

    void GeneratePrimaries(G4Event* anEvent) override;

private:
    G4ParticleGun fGun;
    G4ParticleDefinition* fIon = nullptr;  // создаётся при первом событии
};

#endif // G1S_NPSM_PRIMARY_GENERATOR_ACTION_HH
