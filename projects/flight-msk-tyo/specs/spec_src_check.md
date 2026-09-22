You are a C++ engineer working with Geant4 11.4.2 on Windows (MSVC). Write ONE complete source file `src_check.cc` (a Geant4 application) and output ONLY the code (no markdown fences, no prose).

## Purpose
Check the cosmic-ray source in an EMPTY world: primaries drawn by `ParmaPrimary` (header `parma_primary.hh`, already written) fly through vacuum; every primary that ENTERS a small probe sphere at the origin is recorded (kinetic energy and the downward cosine). The recorded histograms are later compared with the source table.

## Command line
    src_check <table_file> <particle> <N_events> <R_cm> <probe_r_cm> <seed> <out_csv>
* `particle` is a Geant4 name ("neutron", "proton", "mu+", "mu-", "e-", "e+", "gamma").
* Wrong argument count -> usage to stderr, exit code 2. Table load failure -> print `err`, exit 2.

## Application
* `#include "parma_primary.hh"` (it includes `parma_source.hh`).
* Run manager: `G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial)`.
* Physics list: `new FTFP_BERT(0)` (`#include "FTFP_BERT.hh"`); set seed with `G4Random::setTheSeed(seed)` before `BeamOn`.
* Detector construction (a small class `Det : G4VUserDetectorConstruction`): world = `G4Box` of half-size `R_cm + 50` cm in x,y,z, material `G4_Galactic` (NIST), logical volume name "World"; inside it, at the origin, a `G4Orb` "Probe" of radius `probe_r_cm` cm, material `G4_Galactic`. Keep a pointer to the probe logical volume for the stepping action. Use `G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic")`.
* Primary generator: `new ParmaPrimary(&src, particle, R_cm)`.
* Stepping action (class `Step : G4UserSteppingAction`): for a step where `step->GetPreStepPoint()->GetPhysicalVolume()` has logical volume equal to the probe's logical volume AND `step->GetPreStepPoint()->GetStepStatus() == fGeomBoundary` AND `step->GetTrack()->GetParentID() == 0`, record one hit: kinetic energy `step->GetPreStepPoint()->GetKineticEnergy()/MeV` and `cos_down = -step->GetPreStepPoint()->GetMomentumDirection().z()`.
* Histograms filled at each hit (all in plain `std::vector<long long>`, owned by a global struct is acceptable since the run is serial):
  * energy: `nebin = src.NEBin()` bins with the edges `src.EdgesMeV()`; bin index = the smallest `k` in 1..nebin with `E <= edges[k]` (use `std::lower_bound` on the edges), hits with `E` below `edges[0]` or above `edges[nebin]` go to counters `underflow` / `overflow`;
  * cosine: `nabin = src.NABin()` bins, index `min(nabin, max(1, floor((cos_down+1)/2*nabin)+1))`.
* After `BeamOn(N_events)` write the CSV `out_csv` (UTF-8, `\n` newlines):
  line 1: `# N=<N> hits=<total hits> underflow=<..> overflow=<..> TotalFlux=<%.9e> R_cm=<%.9e> probe_r_cm=<%.9e> seed=<seed> particle=<particle>`
  line 2: `E_hist` followed by `,` and the `nebin` energy counts comma-separated
  line 3: `C_hist` followed by `,` and the `nabin` cosine counts comma-separated
* Exit code 0 on success.

## Mandatory structure (the first attempt failed to compile and was logically wrong; follow this literally)
1. Define ONE global context struct BEFORE all classes and one global instance `g_ctx`:
   `struct Ctx { ParmaSource src; std::vector<long long> E_hist, C_hist; long long underflow=0, overflow=0, hits=0; double R_cm=0, probe_r_cm=0; G4LogicalVolume* probe_lv=nullptr; };  static Ctx g_ctx;`
   `Det` reads `g_ctx.R_cm`, `g_ctx.probe_r_cm` in `Construct()` and stores the probe logical volume into `g_ctx.probe_lv`; `Step` reads `g_ctx.probe_lv` and `g_ctx.src`. No private data members are ever assigned from `main`.
2. `main` receives EXACTLY 7 user arguments, so the check is `if (argc != 8)`; the arguments are `argv[1]..argv[7]` in the order of the command line above (table, particle, N, R_cm, probe_r_cm, seed, out_csv). `N_events` must be parsed with `std::atoll` into a `long long`.
3. Load the table with `std::string err; if (!g_ctx.src.Load(table_file, err)) { std::cerr << err << "\n"; return 2; }` (signature: `bool Load(const std::string& path, std::string& err)`). Fill `g_ctx.R_cm`, `g_ctx.probe_r_cm` and size the histograms (`E_hist` size `nebin`, `C_hist` size `nabin`) BEFORE `runManager->Initialize()`.
4. Energy bin convention: 1-based bin `k` covers `(edges[k-1], edges[k]]`; `k = static_cast<int>(std::lower_bound(edges.begin()+1, edges.end(), E) - edges.begin())`; if `E <= edges[0]` count `underflow`; else if `E > edges[nebin]` count `overflow`; else `E_hist[k-1]++`. The CSV column j (0-based) of `E_hist` is therefore bin j+1.
5. Cosine convention: `ia = min(nabin, max(1, (int)std::floor((cos_down+1)/2*nabin) + 1))`, `C_hist[ia-1]++`.
6. `hits` is a SEPARATE counter incremented exactly once per recorded hit. Never compute it by summing the histograms. Write `hits=<g_ctx.hits>` in the CSV header line.
7. The ParmaPrimary must be created AFTER `runManager->SetUserInitialization(new FTFP_BERT(0))` and AFTER `runManager->Initialize()` (so that the particle table is populated), then registered with `runManager->SetUserAction(new ParmaPrimary(&g_ctx.src, particle, R_cm))`, then `BeamOn`. Also register `new Step` with `SetUserAction` (both actions after `Initialize()`).
8. In the stepping action test `pre->GetPhysicalVolume() != nullptr` before dereferencing it.

## Requirements
* Geant4 headers (`G4RunManagerFactory.hh`, `G4VUserDetectorConstruction.hh`, `G4UserSteppingAction.hh`, `G4Box.hh`, `G4Orb.hh`, `G4LogicalVolume.hh`, `G4PVPlacement.hh`, `G4NistManager.hh`, `G4SystemOfUnits.hh`, `G4Step.hh`, `G4Track.hh`, `G4StepPoint.hh`, `Randomize.hh`, `G4VPhysicalVolume.hh`) — include every one you use. Register the actions with `runManager->SetUserAction(new Step)` and `runManager->SetUserAction(new ParmaPrimary(...))` — the primary generator is registered as a user action of the (serial) run manager.
* `DetectorConstruction::Construct()` must return the world `G4VPhysicalVolume*`.
* Use `G4cout` sparingly; set `/run/verbose 0`, `/event/verbose 0`, `/tracking/verbose 0` through `G4UImanager::GetUIpointer()->ApplyCommand(...)` before `BeamOn`.
* Comments in Russian, short. Must compile with MSVC `cl /EHsc /std:c++17 /utf-8` against Geant4 11.4.2 headers.
