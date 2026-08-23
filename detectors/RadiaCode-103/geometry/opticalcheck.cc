// opticalcheck.cc - Stage-1 light collection efficiency check for CsI(Tl) crystal
// #FIT-1 follow-up, D-006, 2026-08-22
// Simplified geometry: single crystal with mirror wraps, gel layer, SiPM detector.
// Monochromatic 550nm (2.254 eV) optical photons fired from fixed Y position,
// isotropic emission direction, measure fraction reaching SiPM volume.
// Known simplifications:
// - Flat (non-dispersive) refractive index
// - No SiPM quantum efficiency or optical coupling modeled
// - No real scintillation generation — direct optical photon source only
// - Crystal material is fresh (not NIST CsI) to allow custom properties

#include "G4RunManagerFactory.hh"
#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Event.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4RandomDirection.hh"
#include "G4UIManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <iostream>
#include <vector>
#include <sstream>
#include <cstdlib>

using namespace std;

namespace {
    double gSourceY_mm = 0.0;
    long gDetected = 0;

    class OpticalDetectorConstruction : public G4VUserDetectorConstruction {
    private:
        G4VPhysicalVolume* fWorldPV = nullptr;
        G4VPhysicalVolume* fCrystalPV = nullptr;
        G4VPhysicalVolume* fSiPMPV = nullptr;

    public:
        G4VPhysicalVolume* Construct() override {
            auto* nist = G4NistManager::Instance();

            // World
            auto* worldBox = new G4Box("world", 50.*mm, 50.*mm, 50.*mm);
            auto* worldLV = new G4LogicalVolume(worldBox, nist->FindOrBuildMaterial("G4_AIR"), "worldLV");
            fWorldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "world", nullptr, false, 0);

            // CsI(Tl) crystal
            auto* csi = new G4Material("CsITl_opt", 4.51*g/cm3, 2);
            csi->AddElement(nist->FindOrBuildElement("Cs"), 1);
            csi->AddElement(nist->FindOrBuildElement("I"), 1);

            auto* crystalBox = new G4Box("crystal", 5.*mm, 5.*mm, 5.*mm);
            auto* crystalLV = new G4LogicalVolume(crystalBox, csi, "crystalLV");
            fCrystalPV = new G4PVPlacement(nullptr, G4ThreeVector(), crystalLV, "crystal", worldLV, false, 0);

            // Air wrap plates (mirror surfaces)
            auto* air = nist->FindOrBuildMaterial("G4_AIR");
            auto* wrapBox = new G4Box("wrapBox", 0.5*mm, 6.*mm, 6.*mm);
            auto* wrapLV = new G4LogicalVolume(wrapBox, air, "wrapLV");

            // +X
            auto* wrapPX = new G4PVPlacement(nullptr, G4ThreeVector(5.5*mm, 0., 0.), wrapLV, "wrapPX", worldLV, false, 0);
            // -X
            auto* wrapNX = new G4PVPlacement(nullptr, G4ThreeVector(-5.5*mm, 0., 0.), wrapLV, "wrapNX", worldLV, false, 0);
            // +Z
            auto* wrapPZ = new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 5.5*mm), wrapLV, "wrapPZ", worldLV, false, 0);
            // -Z
            auto* wrapNZ = new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -5.5*mm), wrapLV, "wrapNZ", worldLV, false, 0);
            // +Y
            auto* wrapPY = new G4PVPlacement(nullptr, G4ThreeVector(0., 5.5*mm, 0.), wrapLV, "wrapPY", worldLV, false, 0);

            // -Y with hole (SiPM window)
            auto* wrapNYBox = new G4Box("wrapNYBox", 6.*mm, 0.5*mm, 6.*mm);
            auto* wrapNYHole = new G4Box("wrapNYHole", 3.*mm, 1.5*mm, 3.*mm);
            auto* wrapNYSolid = new G4SubtractionSolid("wrapNYSolid", wrapNYBox, wrapNYHole);
            auto* wrapNYLV = new G4LogicalVolume(wrapNYSolid, air, "wrapNYLV");
            auto* wrapNY = new G4PVPlacement(nullptr, G4ThreeVector(0., -5.5*mm, 0.), wrapNYLV, "wrapNY", worldLV, false, 0);

            // Optical gel layer
            auto* gel = new G4Material("OpticalGel", 1.0*g/cm3, 1);
            gel->AddElement(nist->FindOrBuildElement("Si"), 1); // trivial composition

            auto* gelBox = new G4Box("gel", 3.*mm, 0.025*mm, 3.*mm);
            auto* gelLV = new G4LogicalVolume(gelBox, gel, "gelLV");
            auto* gelPV = new G4PVPlacement(nullptr, G4ThreeVector(0., -5.025*mm, 0.), gelLV, "gel", worldLV, false, 0);

            // SiPM detector
            auto* sipm = nist->FindOrBuildMaterial("G4_Si");
            auto* sipmBox = new G4Box("sipm", 3.*mm, 0.2*mm, 3.*mm);
            auto* sipmLV = new G4LogicalVolume(sipmBox, sipm, "sipmLV");
            fSiPMPV = new G4PVPlacement(nullptr, G4ThreeVector(0., -5.25*mm, 0.), sipmLV, "sipm", worldLV, false, 0);

            // Material properties
            G4double energies[2] = {1.5*eV, 3.5*eV};
            G4double rindex_csi[2] = {1.79, 1.79};
            G4double abs_csi[2] = {300*mm, 300*mm};

            auto* csiMPT = new G4MaterialPropertiesTable();
            csiMPT->AddProperty("RINDEX", energies, rindex_csi, 2);
            csiMPT->AddProperty("ABSLENGTH", energies, abs_csi, 2);
            csi->SetMaterialPropertiesTable(csiMPT);

            G4double rindex_gel[2] = {1.46, 1.46};
            auto* gelMPT = new G4MaterialPropertiesTable();
            gelMPT->AddProperty("RINDEX", energies, rindex_gel, 2);
            gel->SetMaterialPropertiesTable(gelMPT);

            // BUGFIX 22.08: spec required RINDEX=1.0 for air (world + wrap
            // plates share the same G4_AIR pointer) but the first generation
            // silently omitted it. Without RINDEX on the neighbor material,
            // G4OpBoundaryProcess does not recognize the boundary as optical
            // and photons pass through wrap plates untouched (confirmed:
            // REFLECTIVITY 0.98 vs 1.0 gave byte-identical detected counts,
            // meaning the mirror border surface was never actually invoked).
            G4double rindex_air[2] = {1.0, 1.0};
            auto* airMPT = new G4MaterialPropertiesTable();
            airMPT->AddProperty("RINDEX", energies, rindex_air, 2);
            air->SetMaterialPropertiesTable(airMPT);

            // Optical surface
            auto* mirrorSurface = new G4OpticalSurface("mirrorSurface");
            mirrorSurface->SetType(dielectric_metal);
            mirrorSurface->SetModel(unified);
            mirrorSurface->SetFinish(ground);
            mirrorSurface->SetSigmaAlpha(0.3); // ENGINEERING ASSUMPTION: roughness
            G4double reflectivity[2] = {0.98, 0.98};
            auto* surfMPT = new G4MaterialPropertiesTable();
            surfMPT->AddProperty("REFLECTIVITY", energies, reflectivity, 2);
            mirrorSurface->SetMaterialPropertiesTable(surfMPT);

            // Border surfaces - G4LogicalBorderSurface takes PHYSICAL volumes
            // (G4VPhysicalVolume*), not logical ones - crystalLV would not compile.
            new G4LogicalBorderSurface("crystal_wrapPX", fCrystalPV, wrapPX, mirrorSurface);
            new G4LogicalBorderSurface("crystal_wrapNX", fCrystalPV, wrapNX, mirrorSurface);
            new G4LogicalBorderSurface("crystal_wrapPZ", fCrystalPV, wrapPZ, mirrorSurface);
            new G4LogicalBorderSurface("crystal_wrapNZ", fCrystalPV, wrapNZ, mirrorSurface);
            new G4LogicalBorderSurface("crystal_wrapPY", fCrystalPV, wrapPY, mirrorSurface);
            new G4LogicalBorderSurface("crystal_wrapNY", fCrystalPV, wrapNY, mirrorSurface);

            return fWorldPV;
        }

        G4VPhysicalVolume* GetCrystalPV() const { return fCrystalPV; }
        G4VPhysicalVolume* GetSiPMPV() const { return fSiPMPV; }
    };

    class OpticalPhysList : public G4VModularPhysicsList {
    public:
        OpticalPhysList() {
            // G4OpticalPhysics registers Scintillation/Cerenkov processes that
            // attach to charged particles (e-) - without base EM physics those
            // particles have no ProcessManager and G4OpticalPhysics::
            // ConstructProcess() aborts with a fatal exception. We never fire
            // e-/gamma primaries here (only opticalphoton), so this EM physics
            // is never actually exercised - it only exists to satisfy
            // G4OpticalPhysics's particle-construction requirement.
            RegisterPhysics(new G4EmStandardPhysics());
            RegisterPhysics(new G4OpticalPhysics());
        }
    };

    class OpticalPrimaryGenerator : public G4VUserPrimaryGeneratorAction {
    private:
        G4ParticleGun* fParticleGun = nullptr;

    public:
        OpticalPrimaryGenerator() {
            auto* particleTable = G4ParticleTable::GetParticleTable();
            auto* opticalPhoton = particleTable->FindParticle("opticalphoton");
            if (!opticalPhoton) {
                G4Exception("OpticalPrimaryGenerator", "NoOpticalPhoton", FatalException, "opticalphoton not found");
            }

            fParticleGun = new G4ParticleGun(1);
            fParticleGun->SetParticleDefinition(opticalPhoton);
            fParticleGun->SetParticleEnergy(2.254*eV); // 550nm peak
        }

        void GeneratePrimaries(G4Event* evt) override {
            G4ThreeVector dir = G4RandomDirection();
            G4ThreeVector pol = dir.orthogonal().unit(); // perpendicular to momentum
            fParticleGun->SetParticleMomentumDirection(dir);
            fParticleGun->SetParticlePolarization(pol);
            fParticleGun->SetParticlePosition(G4ThreeVector(0., gSourceY_mm, 0.));
            fParticleGun->GeneratePrimaryVertex(evt);
        }

        ~OpticalPrimaryGenerator() override {
            delete fParticleGun;
        }
    };

    class OpticalSteppingAction : public G4UserSteppingAction {
    public:
        void UserSteppingAction(const G4Step* step) override {
            // NON-const: G4Step::GetTrack() returns non-const, and
            // SetTrackStatus() below is a non-const method - const here would
            // not compile.
            G4Track* track = step->GetTrack();
            if (track->GetDefinition()->GetParticleName() != "opticalphoton") return;

            const G4VPhysicalVolume* vol = step->GetPostStepPoint()->GetPhysicalVolume();
            if (!vol || vol->GetName() != "sipm") return;

            gDetected++;
            track->SetTrackStatus(fStopAndKill);
        }
    };

} // anonymous namespace

int main(int argc, char** argv) {
    long N = (argc > 1) ? atol(argv[1]) : 200000L;
    vector<double> y_positions;
    if (argc > 2) {
        string s(argv[2]);
        size_t pos = 0;
        while ((pos = s.find(',')) != string::npos) {
            y_positions.push_back(stod(s.substr(0, pos)));
            s.erase(0, pos + 1);
        }
        y_positions.push_back(stod(s));
    } else {
        y_positions = {-4.5, -3., -1.5, 0., 1.5, 3., 4.5};
    }

    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
    runManager->SetVerboseLevel(0);

    runManager->SetUserInitialization(new OpticalDetectorConstruction());
    runManager->SetUserInitialization(new OpticalPhysList());
    runManager->SetUserAction(new OpticalPrimaryGenerator());
    runManager->SetUserAction(new OpticalSteppingAction());

    runManager->Initialize();

    for (double y : y_positions) {
        gSourceY_mm = y;
        gDetected = 0;
        runManager->BeamOn(N);
        printf("Y_mm=%.2f  N=%ld  detected=%ld  LCE=%.4f\n", y, N, gDetected, (double)gDetected/N);
    }

    delete runManager;
    return 0;
}
