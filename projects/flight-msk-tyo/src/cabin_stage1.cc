#include "parma_source.hh"
#include "parma_primary.hh"
#include "cabin_geom.hh"
#include "psp_record.hh"
#include "capture_step.hh"
#include "capture_emitter.hh"

#include <G4RunManagerFactory.hh>
#include <G4PhysListFactory.hh>
#include <G4VModularPhysicsList.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4RadioactiveDecayPhysics.hh>
#include <G4EmParameters.hh>
#include <G4IonTable.hh>
#include <G4ParticleGun.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4UserTrackingAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4Track.hh>
#include <G4EventManager.hh>
#include <G4Event.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VProcess.hh>
#include <G4UImanager.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>

static const double kSphereR_cm = std::getenv("CABIN_SPHERE_R") ? std::atof(std::getenv("CABIN_SPHERE_R")) : 480.0;   // 480: сектор ±4 м; больше — чтобы источник захватывал крылья и двигатели (вариант CABIN_FUEL)

struct Ctx {
    ParmaSource src;
    int ip;
    PspWriter writer;
    const CabinDet* det;
    double Ep_MeV, wp;
    long long nCross;
    bool k40 = false;      // режим CABIN_K40: распады K-40 в слое людей (естественная активность тела)
    double t0 = 0.0;       // момент распада K-40 (единицы G4): время записей отсчитывается от него
    bool t0set = false;
    double k40ActBq = 0.0;
};

static Ctx g_ctx;

class Trk : public G4UserTrackingAction {
public:
    void PreUserTrackingAction(const G4Track* t) override {
        if (g_ctx.k40) {   // время отсчитываем от момента распада, а не от рождения ядра (оно на ~1e25 нс раньше)
            if (t->GetParentID() == 0) g_ctx.t0set = false;
            else if (!g_ctx.t0set) { g_ctx.t0 = t->GetGlobalTime(); g_ctx.t0set = true; }
        }
        if (t->GetParentID() == 0) {
            g_ctx.Ep_MeV = t->GetKineticEnergy()/MeV;
            g_ctx.wp = t->GetMomentumDirection().z();
        }
    }
};

class Stp : public G4UserSteppingAction {
private:
    CaptureStep* cap;

public:
    explicit Stp(CaptureStep* c = nullptr) : cap(c) {}

    void SetSteppingManagerPointer(G4SteppingManager* p) override {
        G4UserSteppingAction::SetSteppingManagerPointer(p);
        if (cap) cap->SetSteppingManagerPointer(p);
    }

    void UserSteppingAction(const G4Step* s) override {
        if (cap) cap->UserSteppingAction(s);

        const auto pre = s->GetPreStepPoint();
        if (pre->GetStepStatus() != fGeomBoundary ||
            !pre->GetPhysicalVolume() ||
            pre->GetPhysicalVolume()->GetLogicalVolume() != g_ctx.det->TubeLV()) {
            return;
        }

        const auto trk = s->GetTrack();
        const auto& P = g_ctx.det->Params();

        const auto name = trk->GetDefinition()->GetParticleName();
        if (name.rfind("nu_", 0) == 0 || name.rfind("anti_nu_", 0) == 0) {
            return;
        }

        if (name == "gamma" && pre->GetKineticEnergy() < 0.010*MeV) {
            return;
        }

        PspRec r{};
        const auto evt = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
        r.evt = (int32_t)evt;
        r.pdg = trk->GetDefinition()->GetPDGEncoding();
        r.E_MeV = pre->GetKineticEnergy()/MeV;

        const auto pos = pre->GetPosition()/cm;
        // Входы через торцы трубки не пишем (дизайн: только боковая поверхность; торцевая точка лежит строго внутри радиуса)
        if (std::hypot(pos.x() - P.tubeX, pos.z() - P.tubeZ) < P.tubeR - 0.01) return;
        r.phi = std::atan2(pos.z() - P.tubeZ, pos.x() - P.tubeX);
        r.y_cm = pos.y() - P.tubeY;

        const auto dir = pre->GetMomentumDirection();
        r.ux = dir.x();
        r.uy = dir.y();
        r.uz = dir.z();

        r.t_ns = (pre->GetGlobalTime() - (g_ctx.k40 ? g_ctx.t0 : 0.0))/ns;
        r.Ep_MeV = g_ctx.Ep_MeV;
        r.wp = g_ctx.wp;

        const auto cp = trk->GetCreatorProcess();
        r.proc = ProcCode(cp ? cp->GetProcessName() : "");
        const auto lv = trk->GetLogicalVolumeAtVertex();
        r.vol = VolCode(lv ? lv->GetName() : "");
        r.ip = (int16_t) g_ctx.ip;
        r.pad = 0;

        g_ctx.writer.Write(r);
        ++g_ctx.nCross;
    }
};

class K40Gen : public G4VUserPrimaryGeneratorAction {   // ядро K-40 покоится в случайной точке слоя людей (вне трубки)
    G4ParticleGun gun{1};
public:
    K40Gen() { gun.SetParticleDefinition(G4IonTable::GetIonTable()->GetIon(19, 40, 0.0)); gun.SetParticleEnergy(0.0);
               gun.SetParticleMomentumDirection(G4ThreeVector(0, 0, 1)); }
    void GeneratePrimaries(G4Event* e) override {
        const auto& P = g_ctx.det->Params();
        for (;;) {
            const double x = (2 * G4UniformRand() - 1) * P.paxHalfX, y = (2 * G4UniformRand() - 1) * P.halfLen;
            const double z = P.floorTop + G4UniformRand() * (P.paxTopZ - P.floorTop);
            if (std::hypot(x - P.tubeX, z - P.tubeZ) < P.tubeR && std::fabs(y - P.tubeY) < P.tubeHalfL) continue;
            gun.SetParticlePosition(G4ThreeVector(x, y, z) * cm);
            break;
        }
        gun.GeneratePrimaryVertex(e);
    }
};

int main(int argc, char* argv[]) {
    if (argc < 7 || argc > 8) {
        std::cerr << "Usage: cabin_stage1 <table.tab> <particle> <ip> <N> <seed> <out.psp> [capture_db]\n";
        return 2;
    }

    const char* table = argv[1];
    const char* particle = argv[2];
    int ip = std::stoi(argv[3]);
    long long N = std::stoll(argv[4]);
    long long seed = std::stoll(argv[5]);
    const char* out = argv[6];

    std::string err;
    if (!g_ctx.src.Load(table, err)) {
        std::cerr << "Failed to load table: " << err << "\n";
        return 2;
    }

    if (!g_ctx.writer.Open(out)) {
        std::cerr << "Failed to open output file: " << out << "\n";
        return 2;
    }

    g_ctx.ip = ip;
    g_ctx.nCross = 0;

    G4PhysListFactory factory;
    auto* pl = factory.GetReferencePhysList(std::getenv("CABIN_PHYS") ? std::getenv("CABIN_PHYS") : "FTFP_BERT_HPT");   // HPT = HP + тепловое рассеяние S(alpha,beta)
    if (!pl) {
        std::cerr << "Failed to create physics list FTFP_BERT_HP\n";
        return 2;
    }

    pl->ReplacePhysics(new G4EmStandardPhysics_option4());
    G4EmParameters::Instance()->SetFluo(true);   // флуоресценция в конструкции (K-рентген Cu/Fe/Zn ниже порога записи 10 кэВ, но физика включена)
    pl->RegisterPhysics(new G4RadioactiveDecayPhysics());

    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
    CabinParams cp;
    cp.worldHalf = 1.46 * kSphereR_cm; cp.fuel = std::getenv("CABIN_FUEL") != nullptr;
    if (std::getenv("CABIN_EMPTY")) { cp.empty = true; cp.airOutRho = 1e-10; }   // замыкающий тест: пустой мир, воздух ~вакуум
    auto* det = new CabinDet(cp);
    g_ctx.det = det;
    runManager->SetUserInitialization(det);
    runManager->SetUserInitialization(pl);

    G4Random::setTheSeed(seed);
    runManager->Initialize();

    g_ctx.k40 = std::getenv("CABIN_K40") != nullptr;
    if (g_ctx.k40) {
        const auto& P = det->Params();
        const double vol_cm3 = 2 * P.paxHalfX * 2 * P.halfLen * (P.paxTopZ - P.floorTop) - 3.14159265358979323846 * P.tubeR * P.tubeR * 2 * P.tubeHalfL;
        g_ctx.k40ActBq = 54.8 * P.paxRho * vol_cm3 / 1000.0;   // 54,8 Бк/кг (UNSCEAR 2000) × масса слоя, кг
        runManager->SetUserAction(new K40Gen());
    } else
    try {
        auto* primary = new ParmaPrimary(&g_ctx.src, particle, kSphereR_cm);
        runManager->SetUserAction(primary);
    } catch (const std::runtime_error& e) {
        std::cerr << "ParmaPrimary error: " << e.what() << "\n";
        return 2;
    }

    runManager->SetUserAction(new Trk());

    CaptureStep* cap = nullptr;
    if (argc == 8) {
        static CaptureEmitter emitter;
        if (!emitter.Load(argv[7], err)) {
            std::cerr << "Failed to load capture DB: " << err << "\n";
            return 2;
        }
        cap = new CaptureStep(&emitter);
    }

    runManager->SetUserAction(new Stp(cap));

    const char* uiCmds[] = {
        "/run/setCut 1 mm",   // порог продукции 1 мм в салоне (CFG-1 п. 3); до Initialize нет региона по умолчанию
        g_ctx.k40 ? "/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 s" : "/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+5 s",   // K-40 живёт 1,8e9 лет
        "/run/verbose 0",
        "/event/verbose 0",
        "/tracking/verbose 0"
    };

    auto* UImanager = G4UImanager::GetUIpointer();
    for (const auto& cmd : uiCmds) {
        int status = UImanager->ApplyCommand(cmd);
        std::cout << "UI " << cmd << " -> " << status << "\n";
    }

    runManager->BeamOn((G4int) N);

    g_ctx.writer.Close();

    const auto& P = g_ctx.det->Params();
    double T_sim_s = g_ctx.k40 ? N / g_ctx.k40ActBq : N / (g_ctx.src.TotalFlux() * 3.14159265358979323846 * kSphereR_cm * kSphereR_cm);
    std::ofstream meta(std::string(out) + ".meta");   // метаданные — рядом, а не поверх .psp
    meta << "particle=" << particle << "\n";
    meta << "ip=" << ip << "\n";
    meta << "N=" << N << "\n";
    meta << "seed=" << seed << "\n";
    meta << "R_cm=" << kSphereR_cm << "\n";
    meta << "k40_activity_Bq=" << g_ctx.k40ActBq << std::endl;
    meta << "flux_cm2s=" << g_ctx.src.TotalFlux() << "\n";
    meta << "T_sim_s=" << T_sim_s << "\n";
    meta << "crossings=" << g_ctx.nCross << "\n";
    meta << "psp_records=" << g_ctx.writer.Count() << "\n";

    if (cap) {
        meta << "cap_replaced=" << cap->nReplaced << "\n";
        meta << "cap_fallback=" << cap->nFallback << "\n";
        meta << "cap_no_ion=" << cap->nNoIon << "\n";
        meta << "cap_ion_ground=" << cap->nIonGround << "\n";
    }

    std::cout << "STAGE1 events=" << N << " crossings=" << g_ctx.nCross << " T_sim_s=" << T_sim_s << "\n";

    if (g_ctx.writer.Count() != g_ctx.nCross) {
        std::cerr << "Mismatch between written records and crossings\n";
        return 3;
    }

    return 0;
}
