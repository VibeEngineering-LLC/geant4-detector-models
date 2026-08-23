OLLAMA TASK SPEC — generate ONE self-contained C++ file for Geant4 11.02, output ONLY
raw C++ code, no markdown fences, no explanations before/after.

# Target file: opticalcheck.cc (new standalone Geant4 executable, same folder as
# wallfield.cc / cosmicem.cc in this project — follow their general style: single
# .cc file, no separate .hh, everything in one anonymous namespace + main()).

## Purpose (#FIT-1 follow-up, D-006, 2026-08-22)
Cheap Stage-1 check of the "incomplete light collection" hypothesis for a CsI(Tl)
scintillator crystal: measure the light collection efficiency (LCE) as a function of
the position (along Y) where scintillation photons originate inside the crystal.
Does NOT simulate real gamma events or scintillation generation from energy deposit —
it directly fires optical photons from a point source and measures what fraction
reaches the photo-sensor (SiPM) face. If LCE varies a lot with position, that
supports the light-collection-tailing hypothesis; if it's nearly flat, it does not.

## Physical geometry (simplified — NOT the full RCDetector, this is a standalone
## minimal setup for THIS test only, do not include RCDetector.hh/.cc)

World: G4Box, air, half-size 50mm cube (plenty of margin), centered at origin.

Crystal: G4Box "crystal", CsI(Tl) material, HALF-SIZE 5.0mm cube (i.e. full 10x10x10mm),
centered at world origin (0,0,0). Use G4Material* csi = new G4Material("CsITl_opt",
density 4.51*g/cm3, ncomponents=2); csi->AddElement(nist->FindOrBuildElement("Cs"), 1);
csi->AddElement(nist->FindOrBuildElement("I"), 1); (a fresh material, NOT
G4_CESIUM_IODIDE from NIST, because we need to attach a G4MaterialPropertiesTable to
it and creating our own avoids touching any shared NIST material instance).

FIVE plain mirror plates, one per face +X, -X, +Z, -Z, +Y (everything except -Y).
Each plate: G4Box, half-thickness 0.5mm, material = same air as world (the plate's
own material doesn't matter physically — its only purpose is to give the crystal a
distinct physical-volume neighbor on that face so a G4LogicalBorderSurface can be
attached). Place each plate flush against its crystal face with a small lateral
margin over the 10x10mm face (e.g. the +X plate: half-size (0.5mm, 6mm, 6mm) in
(X,Y,Z), centered at x=5.0+0.5=5.5mm, y=0, z=0). Name them "wrapPX", "wrapNX",
"wrapPZ", "wrapNZ", "wrapPY". Each needs its own distinct G4VPhysicalVolume* (5
separate G4PVPlacement calls; they may share one G4LogicalVolume/G4Box shape
definition if convenient, doesn't matter).

ONE special mirror plate on -Y, "wrapNY", covering the full -Y face (half-size
(6mm, 0.5mm, 6mm), centered at x=0, y=-(5.0+0.5)=-5.5mm, z=0) BUT with a 6x6mm
CENTERED through-hole, built as a G4SubtractionSolid: take the full plate box
(half-size 6mm x 0.5mm x 6mm) and subtract a hole box of half-size 3mm x 1.5mm x 3mm
(Y half-thickness larger than the plate's own 0.5mm so the hole fully perforates it,
no floor left), with ZERO relative translation (hole centered on the plate — use
G4SubtractionSolid(name, solidA, solidB) with no placement argument, or an explicit
G4ThreeVector(0,0,0) if your Geant4 version's constructor requires one). Comment:
"SiPM lateral offset from crystal face center is UNKNOWN for RC103 specifically — a
teardown photo exists only for RC101, whose SiPM chip is ~3x3mm (operator estimate,
2026-08-22), roughly 1/4 the AREA of RC103's 6x6mm chip on the same ~10mm round PCB —
RC101's much smaller chip had far more room to sit off-center than RC103's chip does,
so that photo's apparent offset is not transferable. Operator decision 2026-08-22:
use CENTERED (offset=0) for this Stage-1 check — the main question (does LCE vary
along Y) doesn't depend on the lateral offset; a lateral-offset sensitivity scan is a
separate follow-up if needed."

Optical gel layer: G4Box "gel", half-size (3mm, 0.025mm, 3mm) in (X,Y,Z) — 0.05mm
full thickness (reuses this project's existing RCDetector.hh "window" thickness
convention, comment this provenance), material = a new G4Material "OpticalGel"
(silicone optical coupling gel — elemental composition is physically irrelevant here,
build it trivially, e.g. density 1.0*g/cm3 with a single element; ONLY its RINDEX
property matters). Center at x=0, y=-(5.0+0.025)=-5.025mm, z=0 (centered, see the
offset note above).

SiPM detector: G4Box "sipm", half-size (3mm, 0.2mm, 3mm), material = G4_Si (from
G4NistManager). Centered at x=0, z=0 (same note), y=-(5.0+0.05+0.2)=-5.25mm
(immediately behind the gel layer, zero gap).

Border surfaces: crystal<->wrapPX/NX/PZ/NZ/PY/NY (all 6) use the shared
mirrorSurface (see Optical surfaces section below). Border surface crystal<->gel:
NONE — leave it as a plain dielectric_dielectric boundary so ordinary Fresnel
refraction/reflection applies using the gel's RINDEX (index-matching between crystal
1.79 and gel 1.46 reduces reflection losses versus a bare air gap — this is the whole
physical point of a coupling gel). Border surface gel<->sipm: also NONE (plain
dielectric_dielectric; G4_Si has no RINDEX property attached in this simplified test,
see the SiPM-optics-out-of-scope note further below).

All 8 objects (5 plain mirror plates + 1 subtraction-solid wrapNY + 1 gel + 1 sipm)
are direct daughters of World, not nested. No rotation needed anywhere (everything
axis-aligned).

## Optical material properties (G4MaterialPropertiesTable)

Use a photon energy range 1.5 eV to 3.5 eV (covers ~354-827nm), with just 2 tabulated
points at the range endpoints (constant/flat properties — no wavelength dispersion,
this is an intentional simplification for a Stage-1 order-of-magnitude check, note
this in a comment).

- CsI(Tl) (the csi material): RINDEX = {1.79, 1.79} at {1.5eV, 3.5eV}. ABSLENGTH =
  {300*mm, 300*mm} (absorption length ~30cm at peak wavelength per project's sourced
  literature value, sci-search-corpus 2026-08-22 — comment this provenance).
- World air: RINDEX = {1.0, 1.0} at {1.5eV, 3.5eV}.
- The 6 wrap-plate air volumes: same air material/properties as world (no need for a
  separate G4MaterialPropertiesTable instance if you reuse the same G4Material
  pointer as world's air).
- OpticalGel (the gel material): RINDEX = {1.46, 1.46} at {1.5eV, 3.5eV} — typical
  silicone optical coupling gel, engineering assumption (not measured for this
  device), comment this explicitly as an assumption with the value's source being a
  generic silicone gel literature range, not a project-specific measurement.
  Index-matching between crystal (1.79) and gel (1.46) reduces — but does not
  eliminate — Fresnel reflection losses at that boundary compared to a bare
  crystal/air gap (n=1.0); this is the whole physical point of a coupling gel.
- G4_Si (SiPM material, from NIST): do NOT attach any material properties table (no
  RINDEX) — see the earlier note on SiPM optics being out of scope for this Stage-1
  check; the gel/SiPM boundary is left as an undefined optical surface too, so
  photons entering "sipm" volume are simply tracked with default physics until the
  stepping action detects and kills them.

Attach RINDEX via G4Material::SetMaterialPropertiesTable AFTER building each
G4MaterialPropertiesTable with ->AddProperty("RINDEX", energies, values, 2) (Geant4
11.x G4MaterialPropertiesTable::AddProperty signature takes const G4double* energies,
const G4double* values, G4int num — use this classic array-based form, not the
G4PhysicsOrderedFreeVector overload, for maximum compatibility).

## Optical surfaces (G4OpticalSurface + G4LogicalBorderSurface)

For EACH of the 6 wrap plates (mirror faces: +X,-X,+Y,+Z,-Z, and -Y-with-hole
"wrapNY"): create ONE shared G4OpticalSurface "mirrorSurface" (reused for all 6,
that's fine — G4LogicalBorderSurface takes the same G4OpticalSurface* for multiple
physical-volume pairs):
  auto* mirrorSurface = new G4OpticalSurface("mirrorSurface");
  mirrorSurface->SetType(dielectric_metal);
  mirrorSurface->SetModel(unified);
  mirrorSurface->SetFinish(ground);
  mirrorSurface->SetSigmaAlpha(0.3);
  // Crystal is NOT polished (operator, 2026-08-22) — "ground" finish with the
  // Unified model's Gaussian facet-tilt smearing (sigma_alpha in radians) gives a
  // rough/matte reflecting surface instead of pure specular. sigma_alpha=0.3 rad is
  // an ENGINEERING ASSUMPTION for "moderately matte", not a measured value — no
  // roughness measurement exists for this device. Comment this clearly as OPEN
  // UNCERTAINTY: results may be sensitive to this parameter, a follow-up sensitivity
  // scan (varying sigma_alpha) may be needed if Stage-1 results are inconclusive.
  A G4MaterialPropertiesTable attached to mirrorSurface with
  ->AddProperty("REFLECTIVITY", energies, {0.98, 0.98}, 2) (same 2-point energy array
  as above). mirrorSurface->SetMaterialPropertiesTable(thatTable).

Then for EACH of the 6 crystal/wrap-plate physical-volume pairs, create a
G4LogicalBorderSurface: new G4LogicalBorderSurface("crystal_wrapPX", crystalPV,
wrapPXPV, mirrorSurface); (and similarly for NX, PZ, NZ, PY, NY) — border surfaces are
DIRECTIONAL in principle but for a fully reflective mirror it's fine to only define
crystal->wrap (photons going crystal-to-wrap are what we care about; do NOT bother
defining the reverse wrap->crystal direction, it doesn't matter here since photons
never originate outside the crystal in this simulation).

Reminder (already stated above, do not add a border surface here again): the
gel/sipm boundary has no G4OpticalSurface either. G4_Si has no RINDEX property
attached at all in this Stage-1 test — a photon entering the "sipm" volume is just
tracked with default physics until the stepping action (below) detects and kills it.
Add a one-line comment on the SiPM material definition: real SiPM optical coupling
and quantum efficiency are NOT modeled here, we only count geometric arrival at the
SiPM volume as "collected" (a Stage-1 simplification).

## Physics list

Minimal: a G4VModularPhysicsList subclass (or just use G4VUserPhysicsList minimally)
that registers ONLY G4OpticalPhysics (from "G4OpticalPhysics.hh"). No EM physics, no
hadronic physics — we only ever fire opticalphoton primaries, nothing else needs
tracking. Something like:
  class OpticalPhysList : public G4VModularPhysicsList {
  public: OpticalPhysList() { RegisterPhysics(new G4OpticalPhysics()); }
  };
(follow standard Geant4 11.02 G4OpticalPhysics API — no constructor arguments needed
beyond the default.)

## Detector construction

A G4VUserDetectorConstruction subclass whose Construct() builds everything described
above (world, crystal, 6 wrap plates including wrapNY, gel, sipm, materials, optical
surfaces) and returns
the world physical volume. Store the crystal G4VPhysicalVolume* and the sipm
G4VPhysicalVolume* as public or accessible members (needed later to identify volumes
by pointer or by name in stepping).

## Primary generator

A G4VUserPrimaryGeneratorAction subclass using G4ParticleGun. In its constructor,
get "opticalphoton" via G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton").
Set particle energy to 2.254*eV (550nm peak, cite sci-search-corpus 2026-08-22 in a
comment). In GeneratePrimaryVertex(G4Event* evt): set an ISOTROPIC random direction
(use G4RandomDirection() from "Randomize.hh" / "G4RandomDirection.hh"), set gun
position to a FIXED Y coordinate that is READ FROM A GLOBAL/STATIC DOUBLE VARIABLE
(e.g. a namespace-level `double gSourceY_mm = 0.0;` that main() sets from argv BEFORE
initializing the run manager), with X=0, Z=0 (source always on the central axis, only
Y varies — matches the geometry's -Y-facing SiPM window). GeneratePrimaryGunPosition
sets G4ThreeVector(0, gSourceY_mm, 0) each event (same position every event within one
run, direction randomized isotropically per event). Also randomize optical photon
POLARIZATION per Geant4 convention (set gun's particle polarization to
G4RandomDirection() perpendicular-ish is not required for opticalphoton in G4ParticleGun
— G4ParticleGun does not require explicit polarization for basic photon tracking to
work; you MAY set fParticleGun->SetParticlePolarization(...) to a random perpendicular
vector for physical correctness, but if that adds complexity, a simple approach is
acceptable: set polarization to a fixed vector like G4ThreeVector(1,0,0) is NOT
correct for isotropic photons since polarization must be perpendicular to momentum —
simplest correct approach: after setting momentum direction dir = G4RandomDirection(),
compute any vector perpendicular to dir (e.g. dir.orthogonal().unit()) and use that as
polarization. Use dir.orthogonal() (G4ThreeVector has this method) for a guaranteed
perpendicular vector.

## Stepping action (detection + tail counting)

A G4UserSteppingAction subclass. In UserSteppingAction(const G4Step* step): check if
the current track is an opticalphoton AND the step's post-step point is in the volume
named "sipm" (compare step->GetPostStepPoint()->GetPhysicalVolume()->GetName() to
"sipm", guarding against a null physical volume when the track has left the world).
If so: increment a namespace-level global counter `long gDetected = 0;` (use a plain
global, not thread-local — this program will run single-threaded, no need for MT
safety) and call step->GetTrack()->SetTrackStatus(fStopAndKill) to stop tracking that
photon (it's been "collected"). Also handle photon LOSS bookkeeping is not required —
we only need the detected count vs total fired count (known from the run's event
count), so LCE = gDetected / N_events at the end.

## Run action / end-of-run reporting

Either a G4UserRunAction subclass with EndOfRunAction, OR simply read the global
gDetected counter directly in main() after run/beamOn returns (simpler — prefer this,
since we're single-threaded and don't need the full run-action machinery). Reset
gDetected to 0 before each new source-position run.

## main()

1. Parse argv: argv[1] = number of photons per source point (default 200000 if not
   given), argv[2] = comma-separated list of source Y positions in mm (default
   "-4.5,-3,-1.5,0,1.5,3,4.5" if not given — 7 points spanning the crystal from near
   the SiPM window (-Y) to the far side (+Y), staying safely inside the +-5mm crystal
   half-size).
2. Set G4UIManager verbosity to quiet (like the project's other executables:
   /run/verbose 0 equivalent via G4RunManager::SetPrintProgress(0) or similar — a
   simple `runManager->SetVerboseLevel(0)` call is sufficient, don't overengineer).
3. Construct G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial) (this
   project's other executables use MT — check none needed here; SERIAL is simpler and
   sufficient since this is a small single-threaded diagnostic, and avoids any
   thread-local global-counter complications from the stepping action design above).
4. SetUserInitialization(new OpticalDetectorConstruction()),
   SetUserInitialization(new OpticalPhysList()),
   SetUserAction(new OpticalPrimaryGenerator()),
   SetUserAction(new OpticalSteppingAction()).
5. runManager->Initialize().
6. For each source Y position parsed from argv[2]: set gSourceY_mm = thatValue, reset
   gDetected = 0, call runManager->BeamOn(N), then print one line to stdout:
   `printf("Y_mm=%.2f  N=%ld  detected=%ld  LCE=%.4f\n", y, N, gDetected, (double)gDetected/N);`
   (plain printf/std::printf is fine, this project's other executables use
   std::fprintf/printf freely — follow that convention, no G4cout needed for the
   final result line, though G4cout for progress messages during setup is fine and
   matches project style).
7. delete runManager; return 0;

## Constraints
- Single self-contained .cc file, C++17, Geant4 11.02 API (matches this project's
  other executables — same #include conventions you'd expect from wallfield.cc-style
  code: G4Box, G4PVPlacement, G4NistManager, G4RunManagerFactory, etc. from standard
  Geant4 headers).
- No G4VisManager / visualization — headless batch executable only, matches project
  convention (wallfield.cc/cosmicem.cc are headless).
- Include a file-header comment block explaining the purpose (#FIT-1 D-006 Stage-1
  light-collection-efficiency check, 2026-08-22) and the known simplifications listed
  above (flat/non-dispersive RINDEX, no SiPM quantum efficiency, no real scintillation
  spectrum — monochromatic 550nm photons only).
- Keep it as flat and readable as this project's style (see the reasoning in the
  purpose section above) — small focused classes, no over-engineering, no unused
  includes.
