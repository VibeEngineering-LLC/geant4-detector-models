// Example: RC103 STL case (CADMesh) + full internals from GDML.
// Paths relative to RC103/geant4/cadmesh/:
//   stl:  ../../Rc-103.stl
//   gdml: ../gdml/detector/RC103_detector.gdml
//
// For MC physics prefer RC103_detector.gdml alone (hollow box case).
// Do NOT duplicate case: either CADMesh STL *or* GDML Case_shell, not both.
#include "RC103CADMeshCase.hh"

#include "G4GDMLParser.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

G4VPhysicalVolume* BuildRC103Hybrid(const G4String& stlPath,
                                    const G4String& detectorGdmlPath)
{
  auto* nist = G4NistManager::Instance();
  auto* worldMat = nist->FindOrBuildMaterial("G4_AIR");
  auto* absMat = nist->FindOrBuildMaterial("G4_POLYCARBONATE");  // or custom ABS

  auto* worldBox = new G4Box("WorldBox", 200 * mm, 200 * mm, 200 * mm);
  auto* worldLog = new G4LogicalVolume(worldBox, worldMat, "World");
  auto* worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLog, "World", nullptr,
                                    false, 0);

  // 1) Case from STL (102252 triangles for Rc-103.stl, or 18516 for the .bak
  //    fallback — see RC103CADMeshCase.hh; 110666 was RC-110's count, wrongly
  //    copied here before 26.08.2026, fixed) — physics solid via CADMesh
  auto* caseLog = RC103CADMeshCase::BuildLogicalVolume(stlPath, absMat);
  new G4PVPlacement(nullptr, G4ThreeVector(), caseLog, "Case_STL", worldLog, false, 0);

  // 2) Detector primitives from GDML (fast to iterate)
  G4GDMLParser parser;
  parser.Read(detectorGdmlPath);
  auto* gdmlWorld = parser.GetWorldVolume();
  if (gdmlWorld && gdmlWorld->GetLogicalVolume()) {
    // Import only DetectorModule_log if exported as separate world — or Read partial.
    // Typical pattern: detector GDML world is a small box; place at origin.
    new G4PVPlacement(nullptr, G4ThreeVector(), gdmlWorld->GetLogicalVolume(),
                      "RC103_detector_gdml", worldLog, false, 1);
  }

  return worldPV;
}
