// RC103 CADMesh loader — see RC103CADMeshCase.hh

#include "RC103CADMeshCase.hh"

#include "CADMesh/CADMeshTemplate.hh"

#include "G4LogicalVolume.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4ThreeVector.hh"

G4LogicalVolume* RC103CADMeshCase::BuildLogicalVolume(const G4String& stlPath,
                                                      G4Material* caseABS,
                                                      const G4ThreeVector& extraOffset)
{
  auto mesh = CADMesh::TessellatedMesh::FromSTL(stlPath);
  if (!mesh) return nullptr;
  G4TessellatedSolid* raw = mesh->GetSolid();
  if (!raw) return nullptr;
  auto* solid = new G4TessellatedSolid("Case_STL_CADMesh");
  for (G4int i = 0; i < raw->GetNumberOfFacets(); ++i) {
    auto* facet = raw->GetFacet(i);
    if (!facet) continue;
    G4ThreeVector v0(facet->GetVertex(0).x(), facet->GetVertex(0).y(), facet->GetVertex(0).z());
    G4ThreeVector v1(facet->GetVertex(1).x(), facet->GetVertex(1).y(), facet->GetVertex(1).z());
    G4ThreeVector v2(facet->GetVertex(2).x(), facet->GetVertex(2).y(), facet->GetVertex(2).z());
    solid->AddFacet(new G4TriangularFacet(v0+extraOffset, v1+extraOffset, v2+extraOffset, ABSOLUTE));
  }
  solid->SetSolidClosed(true);
  return new G4LogicalVolume(solid, caseABS, "Case_STL_log");
}
