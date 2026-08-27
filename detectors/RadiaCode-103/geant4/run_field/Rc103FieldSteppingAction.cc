#include "Rc103FieldSteppingAction.hh"

#include "Rc103FieldDetectorConstruction.hh"
#include "Rc103FieldEventAction.hh"
#include "Rc103FieldRunAction.hh"

#include "G4Gamma.hh"
#include "G4LogicalVolume.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TouchableHandle.hh"
#include "G4VPhysicalVolume.hh"

Rc103FieldSteppingAction::Rc103FieldSteppingAction(
    Rc103FieldEventAction* eventAction, Rc103FieldRunAction* runAction)
    : fEventAction(eventAction), fRunAction(runAction) {}

void Rc103FieldSteppingAction::UserSteppingAction(const G4Step* step) {
  G4VPhysicalVolume* prePV =
      step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  if (!prePV) return;
  G4LogicalVolume* preLV = prePV->GetLogicalVolume();

  // --- режим --check-norm: track-length оценка флюенса в контрольном шаре ---
  G4LogicalVolume* checkLV =
      Rc103FieldDetectorConstruction::GetCheckLogicalVolume();
  if (checkLV) {
    if (preLV == checkLV &&
        step->GetTrack()->GetParticleDefinition() == G4Gamma::Definition()) {
      fRunAction->AddGammaTrackLengthMm(step->GetStepLength());
    }
    return;
  }

  // --- обычный режим: энерговыделение в кристалле --------------------------
  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep <= 0.0) return;

  G4LogicalVolume* crystalLV =
      Rc103FieldDetectorConstruction::GetCrystalLogicalVolume();
  if (!crystalLV) return;
  if (preLV == crystalLV) fEventAction->AddEdep(edep);
}
