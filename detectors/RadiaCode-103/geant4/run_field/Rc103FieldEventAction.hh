// Копит депонированную в кристалле энергию за одно событие
// (SteppingAction::AddEdep), в конце события сдаёт сумму в RunAction.
#pragma once

#include "G4UserEventAction.hh"
#include "globals.hh"

class Rc103FieldRunAction;
class G4Event;

class Rc103FieldEventAction : public G4UserEventAction {
 public:
  explicit Rc103FieldEventAction(Rc103FieldRunAction* runAction);
  ~Rc103FieldEventAction() override = default;

  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

  void AddEdep(G4double edep) { fEdep += edep; }

 private:
  Rc103FieldRunAction* fRunAction;
  G4double fEdep = 0.0;
};
