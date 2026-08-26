// Фильтрует шаги по объёму Crystal_log (сравнение указателя логического
// объёма, полученного из Rc103RunDetectorConstruction после парсинга GDML)
// и суммирует энергодепонирование в текущее событие через EventAction.
#pragma once

#include "G4UserSteppingAction.hh"

class Rc103RunEventAction;
class G4Step;

class Rc103RunSteppingAction : public G4UserSteppingAction {
 public:
  explicit Rc103RunSteppingAction(Rc103RunEventAction* eventAction);
  ~Rc103RunSteppingAction() override = default;

  void UserSteppingAction(const G4Step* step) override;

 private:
  Rc103RunEventAction* fEventAction;
};
