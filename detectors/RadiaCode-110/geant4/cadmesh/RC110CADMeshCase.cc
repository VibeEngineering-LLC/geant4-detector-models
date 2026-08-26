// RC110 CADMesh loader — see RC110CADMeshCase.hh

#include "RC110CADMeshCase.hh"

#include "CADMesh/CADMeshTemplate.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"

namespace
{
// Remap STL scan coordinates to RC110 device frame (mm -> mm, then *CLHEP::mm in caller).
G4ThreeVector RemapVertex(G4double sx, G4double sy, G4double sz)
{
  constexpr G4double xc = 122.5;
  constexpr G4double yc = 122.5;
  constexpr G4double zc = 63.3;
  return G4ThreeVector(sz - zc, sx - xc, sy - yc);
}
}  // namespace

G4LogicalVolume* RC110CADMeshCase::BuildLogicalVolume(const G4String& stlPath,
                                                      G4Material* caseABS,
                                                      const G4ThreeVector& extraOffset)
{
  auto mesh = CADMesh::TessellatedMesh::FromSTL(stlPath);
  if (!mesh) {
    G4Exception("RC110CADMeshCase::BuildLogicalVolume", "CADM001",
                FatalException, "CADMesh returned null mesh");
    return nullptr;
  }

  // CADMesh builds G4TessellatedSolid in mesh; fetch and rebuild with remap.
  G4TessellatedSolid* raw = mesh->GetSolid();
  if (!raw) {
    G4Exception("RC110CADMeshCase::BuildLogicalVolume", "CADM002",
                FatalException, "CADMesh solid is null");
    return nullptr;
  }

  auto* solid = new G4TessellatedSolid("Case_STL_CADMesh");
  const G4int nFacets = raw->GetNumberOfFacets();
  for (G4int i = 0; i < nFacets; ++i) {
    auto* facet = raw->GetFacet(i);
    if (!facet) continue;
    G4ThreeVector v0 = RemapVertex(facet->GetVertex(0).x(), facet->GetVertex(0).y(),
                                   facet->GetVertex(0).z()) + extraOffset;
    G4ThreeVector v1 = RemapVertex(facet->GetVertex(1).x(), facet->GetVertex(1).y(),
                                   facet->GetVertex(1).z()) + extraOffset;
    G4ThreeVector v2 = RemapVertex(facet->GetVertex(2).x(), facet->GetVertex(2).y(),
                                   facet->GetVertex(2).z()) + extraOffset;
    solid->AddFacet(new G4TriangularFacet(v0, v1, v2, ABSOLUTE));
  }
  solid->SetSolidClosed(true);
  solid->ComputeBoundingBox();

  return new G4LogicalVolume(solid, caseABS, "Case_STL_log");
}

void RC110CADMeshCase::ApplyRC110VertexRemap(G4VSolid*)
{
  // Reserved: if CADMesh API exposes pre-build vertex hook, wire here.
}
