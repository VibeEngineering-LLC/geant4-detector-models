// Фильтрует шаги по объёму Crystal_log (сравнение указателя логического объёма
// из Rc103MuonDetectorConstruction) и суммирует энерговыделение в текущее
// событие через EventAction.
//
// Тот же паттерн, что run_field/Rc103FieldSteppingAction, минус ветка
// --check-norm: в мюонном прогоне контрольного шара нет.
#pragma once

#include "G4UserSteppingAction.hh"

class Rc103MuonEventAction;
class G4Step;

class Rc103MuonSteppingAction : public G4UserSteppingAction {
 public:
  explicit Rc103MuonSteppingAction(Rc103MuonEventAction* eventAction);
  ~Rc103MuonSteppingAction() override = default;

  void UserSteppingAction(const G4Step* step) override;

 private:
  Rc103MuonEventAction* fEventAction;
};
