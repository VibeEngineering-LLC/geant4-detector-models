#include "Rc103RunPhysicsList.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4SystemOfUnits.hh"

Rc103RunPhysicsList::Rc103RunPhysicsList() {
  SetDefaultCutValue(0.3 * mm);  // см. geant4-win-native SKILL.md, «cut 0.3 мм» для дозиметрии
  RegisterPhysics(new G4EmStandardPhysics_option4());
  RegisterPhysics(new G4DecayPhysics());
  RegisterPhysics(new G4RadioactiveDecayPhysics());
}
