// Track-length оценка флюенса: суммирует длины шагов ФОТОНОВ внутри
// шара-скорера, раскладывая их по энергии фотона на входе в шаг.
#pragma once

#include "G4UserSteppingAction.hh"

class Rc103RoomFieldRunAction;
class G4Step;

class Rc103RoomFieldSteppingAction : public G4UserSteppingAction {
 public:
  explicit Rc103RoomFieldSteppingAction(Rc103RoomFieldRunAction* runAction);
  ~Rc103RoomFieldSteppingAction() override = default;

  void UserSteppingAction(const G4Step* step) override;

 private:
  Rc103RoomFieldRunAction* fRunAction;
};
