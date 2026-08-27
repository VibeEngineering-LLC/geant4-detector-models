// Копит депонированную в кристалле энергию за одно событие
// (SteppingAction::AddEdep), в конце события сдаёт сумму в RunAction.
#pragma once

#include "G4UserEventAction.hh"
#include "globals.hh"

class Rc103MuonRunAction;
class G4Event;

class Rc103MuonEventAction : public G4UserEventAction {
 public:
  explicit Rc103MuonEventAction(Rc103MuonRunAction* runAction);
  ~Rc103MuonEventAction() override = default;

  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

  void AddEdep(G4double edep) { fEdep += edep; }

 private:
  Rc103MuonRunAction* fRunAction;
  G4double fEdep = 0.0;
};
