Write a complete Geant4 11.2 application in C++17 for Windows/MSVC. Output the
files as separate fenced blocks, each preceded by a line `=== FILE: <name> ===`.
No prose outside the blocks. Files to produce, exactly these names:

1. `CMakeLists.txt`
2. `main.cc`
3. `NpsmBenchDetectorConstruction.hh` / `.cc`
4. `NpsmBenchPrimaryGeneratorAction.hh` / `.cc`
5. `NpsmBenchEventAction.hh` / `.cc`
6. `NpsmBenchRunAction.hh` / `.cc`
7. `NpsmBenchSteppingAction.hh` / `.cc`

All comments in Russian. Identifiers in English. Every header has `#pragma once`.

# Purpose

A benchmark of GEOMETRY AND TRANSPORT ONLY (no light yield, no non-proportionality):
a single NaI crystal in vacuum is irradiated by an isotropic monoenergetic photon
flux, and the program counts, for each primary photon that is fully absorbed in
the crystal, how many Compton scatterings it underwent before absorption. The mean
of that count as a function of energy is compared to a published reference.

# Reused donor (do NOT rewrite it, do NOT copy it)

The physics list is taken from an existing project by relative path. In
`CMakeLists.txt` add these two donor files to the executable sources:

```
../../../RadiaCode-103/geant4/run_field/Rc103FieldPhysicsList.cc
```
and add `../../../RadiaCode-103/geant4/run_field` to the include directories so
that `#include "Rc103FieldPhysicsList.hh"` resolves. The donor class is:

```cpp
class Rc103FieldPhysicsList : public G4VModularPhysicsList {
 public:
  explicit Rc103FieldPhysicsList(double cutMm = 0.05, const std::string& deexMode = "std");
  static double gCutMm;          // for the CSV header
  static std::string gDeexMode;  // for the CSV header
};
```
`main.cc` constructs it as `new Rc103FieldPhysicsList(emCutMm, deexMode)` with
default `emCutMm = 0.05` and default `deexMode = "max"`.

# CMakeLists.txt

Model it on this recipe (direct link against import libraries, NOT
`find_package(Geant4)`):

```cmake
cmake_minimum_required(VERSION 3.16)
project(npsm_bench CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(DEFINED ENV{GEANT4_ROOT} AND NOT DEFINED G4ROOT)
  set(G4ROOT "$ENV{GEANT4_ROOT}" CACHE PATH "Geant4 root")
endif()
if(NOT DEFINED G4ROOT)
  message(FATAL_ERROR "Не задан корень Geant4. Укажите GEANT4_ROOT или -DG4ROOT=<путь>.")
endif()
file(GLOB G4_IMPORT_LIBS "${G4ROOT}/lib/*.lib")
set(DONOR_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../RadiaCode-103/geant4/run_field")
add_executable(npsm_bench
  main.cc
  NpsmBenchDetectorConstruction.cc
  NpsmBenchPrimaryGeneratorAction.cc
  NpsmBenchEventAction.cc
  NpsmBenchRunAction.cc
  NpsmBenchSteppingAction.cc
  "${DONOR_DIR}/Rc103FieldPhysicsList.cc"
)
target_include_directories(npsm_bench PRIVATE "${G4ROOT}/include/Geant4" "${DONOR_DIR}")
target_link_libraries(npsm_bench PRIVATE ${G4_IMPORT_LIBS})
```
No Xerces, no GDML.

# main.cc

Usage: `npsm_bench.exe <energy_keV> <n_events> <out_csv> [seed=<N>] [emcut=<mm>]
[deex=std|deex|max] [crystal=<X>x<Y>x<Z>]`

- `energy_keV` — a positive double, the primary photon energy.
- `n_events` — positive long long.
- `out_csv` — output path.
- `seed=` — if nonzero, `G4Random::setTheSeed(seed)` BEFORE `Initialize()`.
- `emcut=` — production cut in mm, default 0.05; must be > 0.
- `deex=` — one of `std`, `deex`, `max`; default `max`; unknown value → print
  error to stderr and return 2 (never silently fall back).
- `crystal=` — full crystal dimensions in mm as `X x Y x Z` parsed with
  `sscanf("%lfx%lfx%lf")`; default `102x102x406`. Stored into static members
  of the detector construction class BEFORE construction.

Parse with `std::strncmp` on prefixes as shown. Any parse failure → message to
stderr with the expected format and return 2.

Use `G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial)`,
`SetVerboseLevel(0)`. Register: detector construction, physics list, run action,
event action (takes run action pointer), primary generator (takes energy in keV),
stepping action (takes event action pointer). Then `Initialize()`.

After `Initialize()` and BEFORE `BeamOn`: check that
`NpsmBenchDetectorConstruction::GetCrystalLogicalVolume()` is non-null;
if null → stderr message and return 3.

Apply UI commands `/run/verbose 1`, `/event/verbose 0`, `/tracking/verbose 0`,
`/run/printProgress <n_events/10, at least 1>`.

Print to stdout before BeamOn one line:
`npsm_bench: energy_keV=%.3f n_events=%lld out_csv=%s crystal_mm=%gx%gx%g r_src_mm=%g seed=%ld`

`BeamOn(n_events)`, delete the run manager, print `EXITCODE=0`, return 0.

# NpsmBenchDetectorConstruction

Static members (doubles, mm): `gCrystalXMm`, `gCrystalYMm`, `gCrystalZMm`
(defaults 102.0, 102.0, 406.0). A static method `double SourceRadiusMm()` that
returns the half-diagonal of the crystal box plus 30 mm margin:
`0.5 * sqrt(X*X + Y*Y + Z*Z) + 30.0`. Document in a comment WHY: the virtual
source sphere must enclose the crystal, otherwise the isotropic-flux identity is
broken.

`Construct()`:
- world: `G4Box`, cube with half-size `SourceRadiusMm() + 50 mm`, material
  `G4_Galactic` (vacuum) from `G4NistManager`;
- crystal: `G4Box` with half-sizes X/2, Y/2, Z/2 in mm, material
  `G4_SODIUM_IODIDE` from `G4NistManager`, placed at origin without rotation,
  `pSurfChk = true`;
- print to stdout one line with the crystal dimensions, its volume in cm³
  (`X*Y*Z/1000`), the material density read back from the material object
  (`GetDensity()/(g/cm3)`), and the world half-size — so the log carries the
  actual constructed geometry, not the intended one;
- store the crystal logical volume in a static pointer, exposed by
  `static G4LogicalVolume* GetCrystalLogicalVolume()`.

# NpsmBenchPrimaryGeneratorAction

Constructor takes `double energyKeV`. Holds a `G4ParticleGun fGun{1}` set to
`gamma`. Implements the isotropic-flux sampling EXACTLY as below (this is the
key physics of the benchmark; copy it verbatim):

```cpp
const double R_SRC = NpsmBenchDetectorConstruction::SourceRadiusMm() * mm;
// 1) точка равномерно по сфере радиуса R_SRC
const double cosT = 2.0 * G4UniformRand() - 1.0;
const double sinT = std::sqrt(1.0 - cosT * cosT);
const double phi  = twopi * G4UniformRand();
const G4ThreeVector n(sinT * std::cos(phi), sinT * std::sin(phi), cosT); // наружная нормаль
const G4ThreeVector pos = R_SRC * n;
// 2) направление ВНУТРЬ, косинусное относительно -n
const G4ThreeVector e3 = -n;
G4ThreeVector e1 = e3.orthogonal().unit();
const G4ThreeVector e2 = e3.cross(e1).unit();
const double ct  = std::sqrt(G4UniformRand());
const double st  = std::sqrt(1.0 - ct * ct);
const double psi = twopi * G4UniformRand();
const G4ThreeVector dir = st * std::cos(psi) * e1 + st * std::sin(psi) * e2 + ct * e3;
```
Comment in Russian explaining: uniform point on the sphere + COSINE law relative
to the inward normal gives an isotropic field inside the sphere; sampling the
direction isotropically instead would NOT. Then set position, direction, energy
(`energyKeV * keV`) and `GeneratePrimaryVertex`.

Include `G4PhysicalConstants.hh` (for `twopi`), `G4SystemOfUnits.hh`,
`Randomize.hh`, `G4ParticleTable.hh`.

# NpsmBenchEventAction

Per-event state, reset in `BeginOfEventAction`:
- `int fNCompt` — number of Compton scatterings of the PRIMARY photon in the crystal;
- `int fNRayl` — number of Rayleigh scatterings of the primary photon in the crystal
  (counted separately, NOT added to Compton);
- `bool fPhotAbsorbed` — primary photon ended by photoelectric effect inside the crystal;
- `bool fConv` — primary photon ended by pair production inside the crystal;
- `bool fEscaped` — primary photon left the world;
- `double fEdepMeV` — total energy deposited in the crystal by ALL particles.

Public inline mutators for each. In `EndOfEventAction` call
`fRunAction->RecordEvent(fNCompt, fNRayl, fPhotAbsorbed, fConv, fEscaped, fEdepMeV)`.

# NpsmBenchSteppingAction

Constructor takes the event action pointer. In `UserSteppingAction(const G4Step*)`:

1. Get pre-step volume; if null return. Get the crystal LV via the static getter;
   if null return.
2. **Energy deposit (all particles):** if pre-step LV == crystal LV and
   `GetTotalEnergyDeposit() > 0`, add it to the event action (in MeV: divide by `MeV`).
3. **Primary photon bookkeeping — only for the track with `GetParentID() == 0`:**
   - get the post-step process: `step->GetPostStepPoint()->GetProcessDefinedStep()`;
     it may be null (transportation-only steps) — guard it;
   - let `name = process->GetProcessName()` when non-null;
   - if the POST-step point is inside the crystal (post-step volume non-null and its
     LV == crystal LV):
       - `name == "compt"` → `AddCompt()`;
       - `name == "Rayl"`  → `AddRayl()`;
       - `name == "phot"`  → `SetPhotAbsorbed()`;
       - `name == "conv"`  → `SetConv()`;
   - if the post-step volume is null (track left the world) → `SetEscaped()`.

   Comment in Russian: why Rayleigh is counted separately (it is elastic, changes
   only direction, and the reference quantity is the number of COMPTON events),
   and why we key on the post-step point (the process that defined the step
   happens at its end).

Include `G4Step.hh`, `G4Track.hh`, `G4VProcess.hh`, `G4LogicalVolume.hh`,
`G4VPhysicalVolume.hh`, `G4TouchableHandle.hh`, `G4SystemOfUnits.hh`.

# NpsmBenchRunAction

Constructor takes `std::string outCsv, double energyKeV, long long nEventsRequested, long seed`.

Accumulators:
- `static constexpr int kMaxCompt = 32;`
- `std::array<long long, kMaxCompt+1> fHistAbsorbed{}` — histogram of `nCompt`
  over events with `photAbsorbed == true` (values > kMaxCompt go into the last bin);
- `std::array<long long, kMaxCompt+1> fHistAll{}` — same over ALL events;
- counters: `fNEvents, fNAbsorbed, fNConv, fNEscaped, fNOther` (other = none of the
  three flags set, e.g. photon stopped by some other process or still alive);
- `fNFullEdep` — events where `fEdepMeV` is within 1 keV of the primary energy
  (an INDEPENDENT second criterion of full absorption, to be compared with the
  process-based flag in the CSV);
- `fSumComptAbsorbed`, `fSumCompt2Absorbed` — for mean and standard error.

`RecordEvent(int nCompt, int nRayl, bool phot, bool conv, bool escaped, double edepMeV)`
updates all of the above.

`EndOfRunAction`: write the CSV with `std::ofstream` opened in text mode. Format:

```
# npsm_bench geometry/transport benchmark
energy_keV,<value>
n_events_requested,<value>
n_events_processed,<value>
seed,<value>
em_cut_mm,<Rc103FieldPhysicsList::gCutMm>
em_deex,<Rc103FieldPhysicsList::gDeexMode>
lowest_electron_energy_keV,<G4EmParameters::Instance()->LowestElectronEnergy()/keV>
crystal_mm,<X>x<Y>x<Z>
crystal_volume_cm3,<value>
r_src_mm,<value>
n_absorbed_phot,<value>
n_full_edep_1keV,<value>
n_conv,<value>
n_escaped,<value>
n_other,<value>
mean_ncompt_absorbed,<mean with 4 decimals>
sem_ncompt_absorbed,<standard error of the mean, 4 decimals>
n_compt,count_absorbed,count_all
0,<..>,<..>
1,<..>,<..>
...
32,<..>,<..>
```

`mean = sum/n`, `sem = sqrt((sum2/n - mean*mean)/n)`; if `n == 0` write `nan`
for both. Include `G4EmParameters.hh` and `G4SystemOfUnits.hh` for the header.

Also print the mean and sem to stdout at end of run as
`npsm_bench: E=%.1f keV  N_absorbed=%lld  mean_ncompt=%.4f +- %.4f  escaped=%lld`.

# Requirements checklist the code must satisfy

- Compiles under MSVC with `/W3` without warnings about signed/unsigned in loops
  (use `std::size_t` or explicit casts).
- No `using namespace std;`.
- `std::fprintf(stdout, ...)` for diagnostics as in the donor project.
- Every static data member defined exactly once in its `.cc`.
- Never silently fall back on a bad CLI value.
