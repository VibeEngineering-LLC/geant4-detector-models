#include "Rc103RoomFieldPhysicsList.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4SystemOfUnits.hh"

Rc103RoomFieldPhysicsList::Rc103RoomFieldPhysicsList(double cutMm) {
  SetDefaultCutValue(cutMm * mm);
  RegisterPhysics(new G4EmStandardPhysics_option4());
  RegisterPhysics(new G4DecayPhysics());
  RegisterPhysics(new G4RadioactiveDecayPhysics());
}
