// RC110 case from Rc-110.stl via CADMesh (runtime tessellated solid).
// Requires: https://github.com/christopherpoole/CADMesh
//
// CMake (after find_package or add_subdirectory CADMesh):
//   target_link_libraries(your_target PRIVATE CADMesh::CADMesh)
//
// STL frame -> device frame (same as stl_to_gdml.py):
//   device_x = stl_z - 63.3 mm
//   device_y = stl_x - 122.5 mm
//   device_z = stl_y - 122.5 mm

#ifndef RC110_CADMESH_CASE_HH
#define RC110_CADMESH_CASE_HH

#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class RC110CADMeshCase
{
public:
  // stlPath: absolute path to Rc-110.stl
  // caseABS: e.g. NIST manager material
  static G4LogicalVolume* BuildLogicalVolume(const G4String& stlPath,
                                             G4Material* caseABS,
                                             const G4ThreeVector& extraOffset =
                                               G4ThreeVector(0, 0, 0));

  // Optional vertex remap hook applied after CADMesh read, before solid build.
  // Default applies RC110 device-axis rotation (length along +X).
  static void ApplyRC110VertexRemap(G4VSolid* tessSolid);
};

#endif
