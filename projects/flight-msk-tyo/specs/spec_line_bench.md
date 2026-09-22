You are a C++ engineer working with Geant4 11.4.2 on Windows (MSVC). Write ONE complete source file `line_bench.cc` (a Geant4 application) and output ONLY the code (no markdown fences, no prose).

## Purpose
A "line bench": shoot mono-energetic neutrons at the centre of a sphere of a chosen material and record every GAMMA created by neutron-induced reactions (radiative capture and inelastic scattering), so that the energies and yields of the produced gamma lines can be compared with nuclear data tables.

## Command line
    line_bench <material> <neutron_energy_MeV> <N_events> <hp_mode> <seed> <out_csv> [radius_cm]
* `material`: a Geant4 NIST material name (for example `G4_Al`, `G4_Fe`, `G4_POLYETHYLENE`, `G4_CESIUM_IODIDE`) OR one of the special names `NYLON` = `G4_NYLON-6-6`, `WATER` = `G4_WATER`, `AIR` = `G4_AIR`. Just pass the string to `G4NistManager::FindOrBuildMaterial`; if it returns null print an error to stderr and exit 2.
* `neutron_energy_MeV`: kinetic energy in MeV (thermal neutrons: `2.53e-8`).
* `hp_mode`: `0` = default high-precision neutron models; `1` = call `G4ParticleHPManager::GetInstance()->SetUseOnlyPhotoEvaporation(true)` BEFORE `runManager->Initialize()` (gamma cascades from the photon-evaporation database).
* `radius_cm` optional, default 3.0: radius of the target sphere.
* Wrong argument count (fewer than 6 or more than 7 user arguments, i.e. `argc` outside 7..8) -> usage to stderr, exit 2.

## Application
* Run manager `G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial)`. Physics list `new FTFP_BERT_HP` hmm — use `#include "G4PhysListFactory.hh"` and `G4PhysListFactory().GetReferencePhysList("FTFP_BERT_HP")`; then replace EM physics: `physList->ReplacePhysics(new G4EmStandardPhysics_option4());` and register `new G4RadioactiveDecayPhysics()` with `physList->RegisterPhysics(...)`.
* Register the detector construction, then the physics list, then `SetUserInitialization(...)` in that order; call `G4Random::setTheSeed(seed)`; after the physics list is registered and the HP switch is applied, `runManager->Initialize()`.
* Detector (class `Det`): world = `G4Box` half-size `radius_cm + 20` cm of `G4_Galactic`; target = `G4Orb` "Target" of radius `radius_cm` cm at the origin made of the chosen material. Read parameters from a single global struct `g_ctx` (defined before the classes): material name, radius, histograms, counters. Do not assign private members from `main`.
* Primary generator (class `Gun : G4VUserPrimaryGeneratorAction`, created AFTER `Initialize()` and registered with `runManager->SetUserAction`): a `G4ParticleGun(1)` with particle `neutron`, the requested kinetic energy, position exactly at the origin, direction `+z`. (The neutron starts inside the target; this is intended.)
* Tracking action (class `Trk : G4UserTrackingAction`, `PreUserTrackingAction(const G4Track* t)`): if the track is a gamma (`t->GetDefinition() == G4Gamma::Gamma()`), its creator process exists (`t->GetCreatorProcess() != nullptr`), and the creator process name is either `nCapture` or `neutronInelastic` (compare `GetProcessName()`), then record the gamma kinetic energy in keV:
  * 1-keV-wide histogram bins from 0 to 12000 keV (12000 bins; energies outside are counted in `out_of_range`), separately for the two processes (arrays `H_capture[12000]`, `H_inelastic[12000]`), bin index = `floor(E_keV)`;
  * also count the total number of gammas per process.
  Register it with `runManager->SetUserAction(new Trk)` (after `Initialize()`).
* After `BeamOn(N_events)` write `out_csv` (UTF-8, `\n` newlines, use `std::ofstream`, DO NOT call `imbue`, use `std::setprecision(10)` for doubles):
  line 1: `# material=<material> E_MeV=<..> N=<N> hp_mode=<0|1> seed=<seed> radius_cm=<..> n_capture_gamma=<..> n_inelastic_gamma=<..> out_of_range=<..>`
  line 2: `capture` followed by `,` and the 12000 counts comma-separated
  line 3: `inelastic` followed by `,` and the 12000 counts comma-separated
* Exit code 0 on success. Use `G4UImanager` to set `/run/verbose 0`, `/event/verbose 0`, `/tracking/verbose 0` before `BeamOn`. Print progress nothing else.

## Requirements
* Include every Geant4 header you use (`G4RunManagerFactory.hh`, `G4PhysListFactory.hh`, `G4VModularPhysicsList.hh`, `G4EmStandardPhysics_option4.hh`, `G4RadioactiveDecayPhysics.hh`, `G4ParticleHPManager.hh`, `G4VUserDetectorConstruction.hh`, `G4VUserPrimaryGeneratorAction.hh`, `G4UserTrackingAction.hh`, `G4ParticleGun.hh`, `G4ParticleTable.hh`, `G4Gamma.hh`, `G4Track.hh`, `G4VProcess.hh`, `G4Box.hh`, `G4Orb.hh`, `G4LogicalVolume.hh`, `G4PVPlacement.hh`, `G4NistManager.hh`, `G4SystemOfUnits.hh`, `G4UImanager.hh`, `Randomize.hh`, `<fstream>`, `<iomanip>`, `<cmath>`, `<vector>`, `<string>`).
* `Trk::PreUserTrackingAction` must skip tracks whose creator process is null (primaries).
* `Det::Construct()` must return the world `G4VPhysicalVolume*`.
* Comments in Russian, short. Must compile with MSVC `cl /EHsc /std:c++17 /utf-8` against Geant4 11.4.2 headers.
