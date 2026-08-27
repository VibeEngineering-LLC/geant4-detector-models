#include "Rc103MuonSteppingAction.hh"

#include "Rc103MuonDetectorConstruction.hh"
#include "Rc103MuonEventAction.hh"

#include "G4LogicalVolume.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TouchableHandle.hh"
#include "G4VPhysicalVolume.hh"

Rc103MuonSteppingAction::Rc103MuonSteppingAction(
    Rc103MuonEventAction* eventAction)
    : fEventAction(eventAction) {}

void Rc103MuonSteppingAction::UserSteppingAction(const G4Step* step) {
  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep <= 0.0) return;

  G4VPhysicalVolume* prePV =
      step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  if (!prePV) return;

  G4LogicalVolume* crystalLV =
      Rc103MuonDetectorConstruction::GetCrystalLogicalVolume();
  if (!crystalLV) return;
  if (prePV->GetLogicalVolume() == crystalLV) fEventAction->AddEdep(edep);
}
