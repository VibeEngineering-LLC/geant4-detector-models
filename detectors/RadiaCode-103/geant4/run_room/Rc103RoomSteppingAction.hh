// Фильтрует шаги по объёму Crystal_log (указатель логического объёма из
// Rc103RoomDetectorConstruction после парсинга GDML) и суммирует
// энергодепонирование в текущее событие через EventAction. Паттерн дословно
// как Rc103RunSteppingAction (эталон).
#pragma once

#include "G4UserSteppingAction.hh"

class Rc103RoomEventAction;
class G4Step;

class Rc103RoomSteppingAction : public G4UserSteppingAction {
 public:
  explicit Rc103RoomSteppingAction(Rc103RoomEventAction* eventAction);
  ~Rc103RoomSteppingAction() override = default;

  void UserSteppingAction(const G4Step* step) override;

 private:
  Rc103RoomEventAction* fEventAction;
};
