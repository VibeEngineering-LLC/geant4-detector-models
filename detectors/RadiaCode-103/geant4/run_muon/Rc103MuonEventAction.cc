#include "Rc103MuonEventAction.hh"

#include "Rc103MuonRunAction.hh"

Rc103MuonEventAction::Rc103MuonEventAction(Rc103MuonRunAction* runAction)
    : fRunAction(runAction) {}

void Rc103MuonEventAction::BeginOfEventAction(const G4Event*) { fEdep = 0.0; }

void Rc103MuonEventAction::EndOfEventAction(const G4Event*) {
  fRunAction->RecordEvent(fEdep);
}
