#include "Rc103RoomPhysicsList.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4SystemOfUnits.hh"

Rc103RoomPhysicsList::Rc103RoomPhysicsList() {
  RegisterPhysics(new G4EmStandardPhysics_option4());
  SetDefaultCutValue(1.0 * mm);
}
