#include "Rc103RoomFieldSteppingAction.hh"

#include "Rc103RoomFieldDetectorConstruction.hh"
#include "Rc103RoomFieldRunAction.hh"

#include "G4Gamma.hh"
#include "G4LogicalVolume.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4TouchableHandle.hh"
#include "G4VPhysicalVolume.hh"

Rc103RoomFieldSteppingAction::Rc103RoomFieldSteppingAction(
    Rc103RoomFieldRunAction* runAction)
    : fRunAction(runAction) {}

void Rc103RoomFieldSteppingAction::UserSteppingAction(const G4Step* step) {
  if (step->GetTrack()->GetParticleDefinition() != G4Gamma::Definition()) return;

  G4LogicalVolume* ballLV =
      Rc103RoomFieldDetectorConstruction::GetBallLogicalVolume();
  if (!ballLV) return;

  G4VPhysicalVolume* prePV =
      step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  if (!prePV || prePV->GetLogicalVolume() != ballLV) return;

  fRunAction->AddTrackLengthCm(
      step->GetPreStepPoint()->GetKineticEnergy() / keV,
      step->GetStepLength() / cm);
}
