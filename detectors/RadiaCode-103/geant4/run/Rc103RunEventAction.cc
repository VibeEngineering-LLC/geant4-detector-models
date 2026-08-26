#include "Rc103RunEventAction.hh"

#include "Rc103RunRunAction.hh"

Rc103RunEventAction::Rc103RunEventAction(Rc103RunRunAction* runAction)
    : fRunAction(runAction) {}

void Rc103RunEventAction::BeginOfEventAction(const G4Event*) { fEdep = 0.0; }

void Rc103RunEventAction::EndOfEventAction(const G4Event*) {
  fRunAction->RecordEvent(fEdep);
}
