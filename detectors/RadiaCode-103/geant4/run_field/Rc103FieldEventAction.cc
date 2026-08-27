#include "Rc103FieldEventAction.hh"

#include "Rc103FieldRunAction.hh"

Rc103FieldEventAction::Rc103FieldEventAction(Rc103FieldRunAction* runAction)
    : fRunAction(runAction) {}

void Rc103FieldEventAction::BeginOfEventAction(const G4Event*) { fEdep = 0.0; }

void Rc103FieldEventAction::EndOfEventAction(const G4Event*) {
  fRunAction->RecordEvent(fEdep);
}
