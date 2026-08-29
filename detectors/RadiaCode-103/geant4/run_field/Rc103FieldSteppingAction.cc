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

  // --- происхождение трека относительно свинца домика -----------------------
  // Категория присваивается ТРЕКУ и наследуется потомками: энергию в кристалле
  // оставляет не фотон, а рождённый им здесь же электрон, и судить о свинце по
  // самому электрону нельзя. Geant4 доводит родительский трек до конца прежде
  // чем взять вторичные, поэтому к рождению потомка категория родителя уже есть.
  G4LogicalVolume* shieldLV =
      Rc103FieldDetectorConstruction::GetShieldLogicalVolume();
  const G4Track* trk = step->GetTrack();
  int cat = Rc103FieldEventAction::kCatDirect;

  if (shieldLV) {
    const G4int tid = trk->GetTrackID();
    cat = fEventAction->GetCategory(tid);
    if (cat < 0) {  // первый шаг этого трека — определяем категорию
      if (trk->GetLogicalVolumeAtVertex() == shieldLV) {
        cat = Rc103FieldEventAction::kCatPbBorn;   // флуоресценция, тормозное
      } else {
        cat = fEventAction->GetCategory(trk->GetParentID());  // наследование
        if (cat < 0) cat = Rc103FieldEventAction::kCatDirect;
      }
      fEventAction->SetCategory(tid, cat);
    }
    // Вошёл в свинец по дороге — становится рассеянным. Рождённый в свинце
    // остаётся рождённым: его происхождение важнее пути.
    if (preLV == shieldLV && cat == Rc103FieldEventAction::kCatDirect) {
      cat = Rc103FieldEventAction::kCatPbScat;
      fEventAction->SetCategory(tid, cat);
    }
  }

  // --- обычный режим: энерговыделение в кристалле --------------------------
  const G4double edep = step->GetTotalEnergyDeposit();
  if (edep <= 0.0) return;

  G4LogicalVolume* crystalLV =
      Rc103FieldDetectorConstruction::GetCrystalLogicalVolume();
  if (!crystalLV) return;
  if (preLV != crystalLV) return;

  fEventAction->AddEdep(edep);
  fEventAction->AddEdepByOrigin(edep, cat);
}
