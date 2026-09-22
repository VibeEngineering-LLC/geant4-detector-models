#include "G4RunManagerFactory.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4UserSteppingAction.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "Randomize.hh"
#include "G4VPhysicalVolume.hh"
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "parma_primary.hh"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <cstdlib>

struct Ctx {
    ParmaSource src;
    std::vector<long long> E_hist, C_hist;
    long long underflow = 0, overflow = 0, hits = 0;
    double R_cm = 0, probe_r_cm = 0;
    G4LogicalVolume* probe_lv = nullptr;
};

static Ctx g_ctx;

class Det : public G4VUserDetectorConstruction {
public:
    G4VPhysicalVolume* Construct() override {
        auto world_solid = new G4Box("World", g_ctx.R_cm * cm + 50 * cm, g_ctx.R_cm * cm + 50 * cm, g_ctx.R_cm * cm + 50 * cm);
        auto world_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
        auto world_log = new G4LogicalVolume(world_solid, world_mat, "World");
        auto world_phys = new G4PVPlacement(nullptr, G4ThreeVector(), world_log, "World", nullptr, false, 0);

        auto probe_solid = new G4Orb("Probe", g_ctx.probe_r_cm * cm);
        auto probe_log = new G4LogicalVolume(probe_solid, world_mat, "Probe");
        new G4PVPlacement(nullptr, G4ThreeVector(), probe_log, "Probe", world_log, false, 0);

        g_ctx.probe_lv = probe_log;

        return world_phys;
    }
};

class Step : public G4UserSteppingAction {
public:
    void UserSteppingAction(const G4Step* step) override {
        auto pre = step->GetPreStepPoint();
        if (pre->GetPhysicalVolume() == nullptr) return;
        if (pre->GetStepStatus() != fGeomBoundary) return;
        if (step->GetTrack()->GetParentID() != 0) return;
        if (pre->GetPhysicalVolume()->GetLogicalVolume() != g_ctx.probe_lv) return;

        double E = pre->GetKineticEnergy() / MeV;
        double cos_down = -pre->GetMomentumDirection().z();

        auto& edges = g_ctx.src.EdgesMeV();
        auto nebin = g_ctx.src.NEBin();
        auto nabin = g_ctx.src.NABin();

        int k = static_cast<int>(std::lower_bound(edges.begin() + 1, edges.end(), E) - edges.begin());
        if (E <= edges[0]) {
            ++g_ctx.underflow;
        } else if (E > edges[nebin]) {
            ++g_ctx.overflow;
        } else {
            g_ctx.E_hist[k - 1]++;
        }

        int ia = std::min(nabin, std::max(1, static_cast<int>(std::floor((cos_down + 1) / 2 * nabin) + 1)));
        g_ctx.C_hist[ia - 1]++;

        ++g_ctx.hits;
    }
};

int main(int argc, char** argv) {
    if (argc != 8) {
        std::cerr << "Использование: src_check <table_file> <particle> <N_events> <R_cm> <probe_r_cm> <seed> <out_csv>\n";
        return 2;
    }

    const char* table_file = argv[1];
    const char* particle = argv[2];
    long long N_events = std::atoll(argv[3]);
    double R_cm = std::atof(argv[4]);
    double probe_r_cm = std::atof(argv[5]);
    long long seed = std::atoll(argv[6]);
    const char* out_csv = argv[7];

    std::string err;
    if (!g_ctx.src.Load(table_file, err)) {
        std::cerr << err << "\n";
        return 2;
    }

    g_ctx.R_cm = R_cm;
    g_ctx.probe_r_cm = probe_r_cm;

    auto nebin = g_ctx.src.NEBin();
    auto nabin = g_ctx.src.NABin();
    g_ctx.E_hist.resize(nebin, 0);
    g_ctx.C_hist.resize(nabin, 0);

    auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
    runManager->SetUserInitialization(new Det());
    runManager->SetUserInitialization(new FTFP_BERT(0));
    G4Random::setTheSeed(seed);
    runManager->Initialize();

    runManager->SetUserAction(new Step());
    runManager->SetUserAction(new ParmaPrimary(&g_ctx.src, particle, R_cm));

    G4UImanager::GetUIpointer()->ApplyCommand("/run/verbose 0");
    G4UImanager::GetUIpointer()->ApplyCommand("/event/verbose 0");
    G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 0");

    runManager->BeamOn(N_events);

    std::ofstream f(out_csv, std::ios::out | std::ios::binary);
    // Локаль не задаём: в русской локали Windows числа получили бы разделители тысяч.
    f << std::setprecision(10);

    double TotalFlux = g_ctx.src.TotalFlux();
    f << "# N=" << N_events << " hits=" << g_ctx.hits
      << " underflow=" << g_ctx.underflow << " overflow=" << g_ctx.overflow
      << " TotalFlux=" << TotalFlux << " R_cm=" << R_cm
      << " probe_r_cm=" << probe_r_cm << " seed=" << seed
      << " particle=" << particle << "\n";

    f << "E_hist";
    for (int i = 0; i < nebin; ++i) {
        f << "," << g_ctx.E_hist[i];
    }
    f << "\n";

    f << "C_hist";
    for (int i = 0; i < nabin; ++i) {
        f << "," << g_ctx.C_hist[i];
    }
    f << "\n";

    return 0;
}
