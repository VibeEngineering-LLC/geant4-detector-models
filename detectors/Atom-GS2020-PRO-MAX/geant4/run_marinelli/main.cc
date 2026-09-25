#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "G4RunManagerFactory.hh"
#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmParameters.hh"
#include "G4ThreeVector.hh"
#include "G4IonTable.hh"
#include "G4UImanager.hh"

const std::string kGdmlPath = "C:/g4work/gs2020/GS2020_marinelli_th232.gdml";
static G4VSolid* gSourceSolid = nullptr;
static G4LogicalVolume* gSourceLV = nullptr;
static G4LogicalVolume* gCrystalLV = nullptr;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    G4VPhysicalVolume* Construct() override {
        G4GDMLParser parser;
        parser.SetOverlapCheck(false);
        parser.Read(kGdmlPath, false);

        G4VPhysicalVolume* world = parser.GetWorldVolume();

        G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();
        gSourceLV = store->GetVolume("LV_SourceMatrix");
        gCrystalLV = store->GetVolume("LV_NaI_crystal");

        if (!gSourceLV || !gCrystalLV) {
            std::fprintf(stderr, "FATAL: volume not found\n");
            std::exit(3);
        }

        gSourceSolid = gSourceLV->GetSolid();

        return world;
    }
};

class GS2020PhysicsList : public G4VModularPhysicsList {
public:
    GS2020PhysicsList() {
        RegisterPhysics(new G4EmStandardPhysics_option4());
        RegisterPhysics(new G4DecayPhysics());
        RegisterPhysics(new G4RadioactiveDecayPhysics());

        SetDefaultCutValue(0.05*mm);
        G4EmParameters* p = G4EmParameters::Instance();
        p->SetFluo(true);
        p->SetAuger(true);
        p->SetPixe(false);
        p->SetDeexcitationIgnoreCut(true);
    }

    void SetCuts() override {
        G4VUserPhysicsList::SetCuts();
    }
};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
private:
    G4ParticleGun* fGun;
    double fEnergyKeV;
    int fIonZ;
    int fIonA;
    G4ParticleDefinition* fIon;

public:
    PrimaryGeneratorAction(double energyKeV, int ionZ, int ionA) : fEnergyKeV(energyKeV), fIonZ(ionZ), fIonA(ionA), fIon(nullptr) {
        fGun = new G4ParticleGun(1);
        if (fIonZ > 0) {
            // ион берётся лениво в GeneratePrimaries: до Initialize() таблица ионов не готова
        } else {
            G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
            G4String particleName = "gamma";
            G4ParticleDefinition* particleDef = particleTable->FindParticle(particleName);
            fGun->SetParticleDefinition(particleDef);
        }
    }

    void GeneratePrimaries(G4Event* event) override {
        const double R = 70.5 * mm;
        const double HALF_Z = 52.5 * mm;

        G4ThreeVector position;
        for (;;) {
            double r = R * std::sqrt(G4UniformRand());
            double phi = CLHEP::twopi * G4UniformRand();
            double x = r * std::cos(phi);
            double y = r * std::sin(phi);
            double z = (2.0 * G4UniformRand() - 1.0) * HALF_Z;
            G4ThreeVector pLocal(x, y, z);

            if (gSourceSolid->Inside(pLocal) == kInside) {
                position = pLocal;
                break;
            }
        }

        // Translate to world coordinates
        G4ThreeVector worldPos = position + G4ThreeVector(0, 0, 40.5 * mm);

        double ct = 2.0 * G4UniformRand() - 1.0;
        double st = std::sqrt(1.0 - ct * ct);
        double ph2 = CLHEP::twopi * G4UniformRand();
        G4ThreeVector dir(st * std::cos(ph2), st * std::sin(ph2), ct);

        fGun->SetParticlePosition(worldPos);
        fGun->SetParticleMomentumDirection(dir);
        if (fIonZ > 0) {
            if (!fIon) {
                fIon = G4IonTable::GetIonTable()->GetIon(fIonZ, fIonA, 0.0);
                if (!fIon) { std::fprintf(stderr, "FATAL: ion Z=%d A=%d not found\n", fIonZ, fIonA); std::abort(); }
                fGun->SetParticleDefinition(fIon);
            }
            fGun->SetParticleEnergy(0.0 * keV);
        } else {
            fGun->SetParticleEnergy(fEnergyKeV * keV);
        }
        fGun->GeneratePrimaryVertex(event);
    }
};

class SteppingAction : public G4UserSteppingAction {
public:
    static double gEdepThisEvent;

    void UserSteppingAction(const G4Step* step) override {
        G4LogicalVolume* lv = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
        if (lv == gCrystalLV) {
            gEdepThisEvent += step->GetTotalEnergyDeposit() / keV;
        }
    }
};

double SteppingAction::gEdepThisEvent = 0.0;

class RunAction : public G4UserRunAction {
private:
    std::vector<long long> fHist;
    std::string fOutCsv;
    double fEnergyKeV;
    long long fSeed;
    long long fNEvents;
    int fIonZ;
    int fIonA;
    long long fNWithEdep;

public:
    RunAction(const std::string& outCsv, double energyKeV, long long seed, long long nEvents, int ionZ, int ionA)
        : fOutCsv(outCsv), fEnergyKeV(energyKeV), fSeed(seed), fNEvents(nEvents), fIonZ(ionZ), fIonA(ionA) {
        fHist.resize(3700, 0);
        fNWithEdep = 0;
    }

    void Fill(double edepKeV) {
        if (edepKeV > 0.5) {
            fNWithEdep++;
            int ch = (int)std::floor(edepKeV);
            if (ch >= 0 && ch <= 3699) {
                fHist[ch]++;
            }
        }
    }

    void EndOfRunAction(const G4Run* run) override {
        {
            std::ofstream f(fOutCsv);
            f << "# gs2020_marinelli volume source in Marinelli, GDML " << kGdmlPath << "\n";
            f << "energy_keV," << fEnergyKeV << "\n";
            f << "n_events_requested," << fNEvents << "\n";
            f << "n_events_processed," << run->GetNumberOfEvent() << "\n";
            f << "seed," << fSeed << "\n";
            f << "em_cut_mm,0.05\n";
            f << "em_deex,deex\n";
            f << "em_option,opt4\n";
            f << "particle," << (fIonZ > 0 ? "ion" : "gamma") << "\n";
            f << "ion_z," << fIonZ << "\n";
            f << "ion_a," << fIonA << "\n";
            f << "decay," << (fIonZ > 0 ? 1 : 0) << "\n";
            f << "vessel,marinelli\n";
            f << "src_mode,sample\n";
            f << "sample_matrix,Epoxy_crumb\n";
            f << "sample_rho_g_cm3,0.7987\n";
            f << "npsm_enabled,0\n";
            f << "n_with_edep," << fNWithEdep << "\n";
            f << "bin_keV,count_edep,count_light\n";

            for (int i = 0; i < 3700; ++i) {
                f << i + 0.5 << "," << fHist[i] << "," << fHist[i] << "\n";
            }
            f.close();

            long long total = 0;
            for (long long x : fHist) total += x;
            std::cout << "[gs2020_marinelli] written " << fOutCsv << " total_hits=" << total << "\n";
        }
    }
};

class EventAction : public G4UserEventAction {
private:
    RunAction* fRunAction;

public:
    EventAction(RunAction* runAction) : fRunAction(runAction) {}

    void BeginOfEventAction(const G4Event*) override {
        SteppingAction::gEdepThisEvent = 0.0;
    }

    void EndOfEventAction(const G4Event*) override {
        fRunAction->Fill(SteppingAction::gEdepThisEvent);
    }
};

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "usage: gs2020_marinelli.exe <energy_keV|ion:Z:A> <n_events> <out_csv> [seed]\n");
        return 2;
    }

    double energyKeV = 0.0;
    int ionZ = 0;
    int ionA = 0;

    std::string arg1(argv[1]);
    if (arg1.rfind("ion:", 0) == 0) {
        size_t pos1 = arg1.find(':', 4);
        if (pos1 != std::string::npos) {
            ionZ = std::stoi(arg1.substr(4, pos1 - 4));
            ionA = std::stoi(arg1.substr(pos1 + 1));
        } else {
            std::fprintf(stderr, "FATAL: invalid ion format\n");
            return 3;
        }
    } else {
        energyKeV = std::stod(argv[1]);
    }

    long long nEvents = std::stoll(argv[2]);
    std::string outCsv = std::string(argv[3]);
    unsigned long seed = (argc > 4) ? std::stoul(argv[4]) : (unsigned long)std::time(nullptr);

    G4Random::setTheSeed(seed);

    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
    runManager->SetUserInitialization(new DetectorConstruction());
    runManager->SetUserInitialization(new GS2020PhysicsList());

    std::cout << "cut_mm=0.0500 deex=deex fluo=1 auger=1 pixe=0 ignore_cut=1\n";

    runManager->SetUserAction(new PrimaryGeneratorAction(energyKeV, ionZ, ionA));
    RunAction* runAction = new RunAction(outCsv, energyKeV, seed, nEvents, ionZ, ionA);
    runManager->SetUserAction(runAction);
    runManager->SetUserAction(new EventAction(runAction));
    runManager->SetUserAction(new SteppingAction());

    runManager->Initialize();

    if (ionZ > 0) {
        G4UImanager* ui = G4UImanager::GetUIpointer();
        ui->ApplyCommand("/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+30 ns");
        ui->ApplyCommand("/process/had/rdm/nucleusLimits " + std::to_string(ionA) + " " + std::to_string(ionA) + " " + std::to_string(ionZ) + " " + std::to_string(ionZ));
        std::cout << "CHAIN_CUT: nucleusLimits A=" << ionA << " Z=" << ionZ << "\n";
    }

    runManager->BeamOn(nEvents);

    delete runManager;
    return 0;
}
