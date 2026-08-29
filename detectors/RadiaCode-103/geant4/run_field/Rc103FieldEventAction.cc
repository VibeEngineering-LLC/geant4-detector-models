#include "Rc103FieldEventAction.hh"

#include "Rc103FieldRunAction.hh"

Rc103FieldEventAction::Rc103FieldEventAction(Rc103FieldRunAction* runAction)
    : fRunAction(runAction) {}

void Rc103FieldEventAction::BeginOfEventAction(const G4Event*) {
  fEdep = 0.0;
  fEdepCat.fill(0.0);
  fTrackCat.clear();
}

void Rc103FieldEventAction::EndOfEventAction(const G4Event*) {
  // Событие может смешивать категории. Полное энерговыделение кладём в ту, что
  // внесла БОЛЬШЕ энергии: так спектр каждой категории остаётся спектром
  // полного энерговыделения, а не обрезком.
  int cat = kCatDirect;
  for (int i = 1; i < kNCat; ++i) {
    if (fEdepCat[i] > fEdepCat[cat]) cat = i;
  }
  fRunAction->RecordEvent(fEdep, cat);
}
