// Два разных дела в зависимости от режима:
//  - обычный режим: фильтрует шаги по объёму Crystal_log (сравнение указателя
//    логического объёма из Rc103FieldDetectorConstruction) и суммирует
//    энерговыделение в текущее событие через EventAction;
//  - --check-norm: копит суммарную длину треков ФОТОНОВ внутри контрольного
//    шара для track-length оценки флюенса.
#pragma once

#include "G4UserSteppingAction.hh"

class Rc103FieldEventAction;
class Rc103FieldRunAction;
class G4Step;

class Rc103FieldSteppingAction : public G4UserSteppingAction {
 public:
  Rc103FieldSteppingAction(Rc103FieldEventAction* eventAction,
                           Rc103FieldRunAction* runAction);
  ~Rc103FieldSteppingAction() override = default;

  void UserSteppingAction(const G4Step* step) override;

 private:
  Rc103FieldEventAction* fEventAction;
  Rc103FieldRunAction* fRunAction;
};
