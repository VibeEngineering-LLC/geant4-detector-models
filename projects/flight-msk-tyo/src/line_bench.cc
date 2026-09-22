#include "G4RunManagerFactory.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4ParticleHPManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "capture_step.hh"   // подмена гамма захвата каскадом из базы IAEA/ENSDF (необязательный 8-й аргумент)
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

struct Ctx {
    std::string material;
    double neutron_energy_MeV;
    int N_events;
    int hp_mode;
    long long seed;
    std::string out_csv;
    double radius_cm;
    std::vector<int> H_capture;
    std::vector<int> H_inelastic;
    int n_capture_gamma;
    int n_inelastic_gamma;
    int out_of_range;
};

Ctx g_ctx;

class Det : public G4VUserDetectorConstruction {
public:
    G4VPhysicalVolume* Construct() override {
        G4NistManager* nist = G4NistManager::Instance();
        G4Material* target_material = nullptr;
        if (g_ctx.material == "NYLON") {
            target_material = nist->FindOrBuildMaterial("G4_NYLON-6-6");
        } else if (g_ctx.material == "WATER") {
            target_material = nist->FindOrBuildMaterial("G4_WATER");
        } else if (g_ctx.material == "AIR") {
            target_material = nist->FindOrBuildMaterial("G4_AIR");
        } else {
            target_material = nist->FindOrBuildMaterial(g_ctx.material.c_str());
        }
        if (!target_material) {
            std::cerr << "Ошибка: неизвестный материал '" << g_ctx.material << "'\n";
            exit(2);
        }

        G4double world_size = g_ctx.radius_cm + 20.0;
        G4Box* world_solid = new G4Box("World", world_size * cm, world_size * cm, world_size * cm);
        G4LogicalVolume* world_logic = new G4LogicalVolume(world_solid, nist->FindOrBuildMaterial("G4_Galactic"), "World");
        G4VPhysicalVolume* world_phys = new G4PVPlacement(0, G4ThreeVector(), world_logic, "World", 0, false, 0);

        G4Orb* target_solid = new G4Orb("Target", g_ctx.radius_cm * cm);
        G4LogicalVolume* target_logic = new G4LogicalVolume(target_solid, target_material, "Target");
        new G4PVPlacement(0, G4ThreeVector(), target_logic, "Target", world_logic, false, 0);

        return world_phys;
    }
};

class Gun : public G4VUserPrimaryGeneratorAction {
public:
    void GeneratePrimaries(G4Event* anEvent) override {
        G4ParticleGun* gun = new G4ParticleGun(1);
        gun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("neutron"));
        gun->SetParticleEnergy(g_ctx.neutron_energy_MeV * MeV);
        gun->SetParticlePosition(G4ThreeVector());
        gun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));
        gun->GeneratePrimaryVertex(anEvent);
        delete gun;
    }
};

class Trk : public G4UserTrackingAction {
public:
    void PreUserTrackingAction(const G4Track* t) override {
        if (t->GetDefinition() != G4Gamma::Gamma()) return;
        const G4VProcess* creator = t->GetCreatorProcess();
        if (!creator) return;
        std::string proc_name = creator->GetProcessName();
        double E_keV = t->GetKineticEnergy() / keV;
        int bin = static_cast<int>(floor(E_keV));
        if (bin < 0 || bin >= 12000) {
            g_ctx.out_of_range++;
            return;
        }
        if (proc_name == "nCapture") {
            g_ctx.H_capture[bin]++;
            g_ctx.n_capture_gamma++;
        } else if (proc_name == "neutronInelastic") {
            g_ctx.H_inelastic[bin]++;
            g_ctx.n_inelastic_gamma++;
        }
    }
};

int main(int argc, char** argv) {
    if (argc < 7 || argc > 9) {
        std::cerr << "Использование: line_bench <material> <neutron_energy_MeV> <N_events> <hp_mode> <seed> <out_csv> [radius_cm [capture_db]]\n";
        exit(2);
    }
    static CaptureEmitter emitter;   // база каскадов; загружается только если задан 8-й аргумент
    if (argc == 9) {
        std::string emitErr;
        if (!emitter.Load(argv[8], emitErr)) { std::cerr << "capture_db: " << emitErr << "\n"; exit(2); }
    }

    g_ctx.material = argv[1];
    g_ctx.neutron_energy_MeV = std::stod(argv[2]);
    g_ctx.N_events = std::stoi(argv[3]);
    g_ctx.hp_mode = std::stoi(argv[4]);
    g_ctx.seed = std::stoll(argv[5]);
    g_ctx.out_csv = argv[6];
    g_ctx.radius_cm = (argc >= 8) ? std::stod(argv[7]) : 3.0;
    g_ctx.H_capture.assign(12000, 0);    // гистограммы 1 кэВ/бин обязаны иметь размер до записи
    g_ctx.H_inelastic.assign(12000, 0);

    G4RunManager* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
    G4PhysListFactory factory;
    // Режим -> физ-лист: 0,1 = FTFP_BERT_HP (1 = + фотоиспарение); 2 = FTFP_BERT без HP (захват через G4NeutronRadCapture);
    // 3 = QGSP_BIC_AllHP (как у Roy et al. 2021); 4 = Shielding.
    const char* listName = "FTFP_BERT_HP";
    if (g_ctx.hp_mode == 2) listName = "FTFP_BERT";
    else if (g_ctx.hp_mode == 3) listName = "QGSP_BIC_AllHP";
    else if (g_ctx.hp_mode == 4) listName = "Shielding";
    G4VModularPhysicsList* physList = factory.GetReferencePhysList(listName);
    if (!physList) { std::cerr << "Ошибка: нет физ-листа " << listName << "\n"; return 2; }
    physList->ReplacePhysics(new G4EmStandardPhysics_option4());
    physList->RegisterPhysics(new G4RadioactiveDecayPhysics());
    if (g_ctx.hp_mode == 1) {
        G4ParticleHPManager::GetInstance()->SetUseOnlyPhotoEvaporation(true);
        // Программный вызов один не менял гистограмму (проверено 20.09: файлы режимов 0 и 1 побайтово равны) —
        // дублируем штатной командой мессенджера HP; ответ команды печатаем, чтобы видеть, принята ли она.
        G4UImanager::GetUIpointer()->ApplyCommand("/process/had/particle_hp/use_photo_evaporation true");
    }
    runManager->SetUserInitialization(physList);

    Det* detector = new Det();
    runManager->SetUserInitialization(detector);

    G4Random::setTheSeed(g_ctx.seed);
    runManager->Initialize();

    runManager->SetUserAction(new Gun());
    runManager->SetUserAction(new Trk());
    CaptureStep* captureStep = nullptr;
    if (argc == 9) {
        captureStep = new CaptureStep(&emitter);
        runManager->SetUserAction(captureStep);
    }

    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    UImanager->ApplyCommand("/run/verbose 0");
    UImanager->ApplyCommand("/event/verbose 0");
    UImanager->ApplyCommand("/tracking/verbose 0");

    runManager->BeamOn(g_ctx.N_events);

    if (captureStep) {
        std::cout << "CAPTURE_STEP replaced=" << captureStep->nReplaced << " fallback=" << captureStep->nFallback
                  << " no_ion=" << captureStep->nNoIon << " ion_ground=" << captureStep->nIonGround << " gammas=" << captureStep->nGammasEmitted
                  << " sumE_keV=" << std::setprecision(10) << captureStep->sumEmittedKeV << std::endl;
    }
    if (captureStep) for (auto& kv : captureStep->keyCount) std::cout << "CAPTURE_KEY " << kv.first << " " << kv.second << std::endl;
    std::ofstream out(g_ctx.out_csv, std::ios::out);
    out << "# material=" << g_ctx.material
        << " E_MeV=" << std::setprecision(10) << g_ctx.neutron_energy_MeV
        << " N=" << g_ctx.N_events
        << " hp_mode=" << g_ctx.hp_mode
        << " seed=" << g_ctx.seed
        << " radius_cm=" << std::setprecision(10) << g_ctx.radius_cm
        << " n_capture_gamma=" << g_ctx.n_capture_gamma
        << " n_inelastic_gamma=" << g_ctx.n_inelastic_gamma
        << " out_of_range=" << g_ctx.out_of_range
        << "\n";

    out << "capture";
    for (int i = 0; i < 12000; ++i) {
        out << "," << g_ctx.H_capture[i];
    }
    out << "\n";

    out << "inelastic";
    for (int i = 0; i < 12000; ++i) {
        out << "," << g_ctx.H_inelastic[i];
    }
    out << "\n";

    out.close();
    return 0;
}
