#include "Rc103RoomEventAction.hh"

#include "Rc103RoomRunAction.hh"

Rc103RoomEventAction::Rc103RoomEventAction(Rc103RoomRunAction* runAction)
    : fRunAction(runAction) {}

void Rc103RoomEventAction::BeginOfEventAction(const G4Event*) { fEdep = 0.0; }

void Rc103RoomEventAction::EndOfEventAction(const G4Event*) {
  fRunAction->RecordEvent(fEdep);
}
