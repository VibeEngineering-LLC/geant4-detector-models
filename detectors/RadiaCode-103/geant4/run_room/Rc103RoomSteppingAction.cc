#include "Rc103RoomSteppingAction.hh"

#include "Rc103RoomDetectorConstruction.hh"
#include "Rc103RoomEventAction.hh"

#include "G4LogicalVolume.hh"
#include "G4Step.hh"
#include "G4TouchableHandle.hh"
#include "G4VPhysicalVolume.hh"

Rc103RoomSteppingAction::Rc103RoomSteppingAction(Rc103RoomEventAction* eventAction)
    : fEventAction(eventAction) {}

void Rc103RoomSteppingAction::UserSteppingAction(const G4Step* step) {
  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep <= 0.0) return;

  G4LogicalVolume* crystalLV = Rc103RoomDetectorConstruction::GetCrystalLogicalVolume();
  if (!crystalLV) return;

  G4LogicalVolume* preLV = step->GetPreStepPoint()
                                ->GetTouchableHandle()
                                ->GetVolume()
                                ->GetLogicalVolume();
  if (preLV == crystalLV) {
    fEventAction->AddEdep(edep);
  }
}
