// Копит депонированную в кристалле энергию за одно событие, в конце события
// сдаёт сумму в RunAction. Паттерн дословно как Rc103RunEventAction (эталон).
#pragma once

#include "G4UserEventAction.hh"
#include "globals.hh"

class Rc103RoomRunAction;
class G4Event;

class Rc103RoomEventAction : public G4UserEventAction {
 public:
  explicit Rc103RoomEventAction(Rc103RoomRunAction* runAction);
  ~Rc103RoomEventAction() override = default;

  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

  void AddEdep(G4double edep) { fEdep += edep; }

 private:
  Rc103RoomRunAction* fRunAction;
  G4double fEdep = 0.0;
};
