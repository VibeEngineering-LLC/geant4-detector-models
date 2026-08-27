#include "Rc103FieldPhysicsList.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4SystemOfUnits.hh"

Rc103FieldPhysicsList::Rc103FieldPhysicsList() {
  SetDefaultCutValue(0.05 * mm);
  RegisterPhysics(new G4EmStandardPhysics_option4());
}
