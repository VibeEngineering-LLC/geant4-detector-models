#include "Rc103RunSteppingAction.hh"

#include "Rc103RunDetectorConstruction.hh"
#include "Rc103RunEventAction.hh"

#include "G4LogicalVolume.hh"
#include "G4Step.hh"
#include "G4TouchableHandle.hh"
#include "G4VPhysicalVolume.hh"

Rc103RunSteppingAction::Rc103RunSteppingAction(Rc103RunEventAction* eventAction)
    : fEventAction(eventAction) {}

void Rc103RunSteppingAction::UserSteppingAction(const G4Step* step) {
  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep <= 0.0) return;

  G4LogicalVolume* crystalLV = Rc103RunDetectorConstruction::GetCrystalLogicalVolume();
  if (!crystalLV) return;

  G4LogicalVolume* preLV = step->GetPreStepPoint()
                                ->GetTouchableHandle()
                                ->GetVolume()
                                ->GetLogicalVolume();
  if (preLV == crystalLV) {
    fEventAction->AddEdep(edep);
  }
}
