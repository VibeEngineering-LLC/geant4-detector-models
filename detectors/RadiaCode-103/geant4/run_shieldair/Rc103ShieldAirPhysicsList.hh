// Физ-лист run_shieldair: G4EmStandardPhysics_option4 (порог 0.05 мм, как в
// Rc103FieldPhysicsList - источник у самого кристалла, грубый порог обрежет
// низкоэнергичные электроны/тормозное) + G4DecayPhysics + G4RadioactiveDecayPhysics
// (без них ион-первичка не распадается вовсе - тот же принцип, что в
// Rc103RoomFieldPhysicsList; отдельный класс, а не переиспользование поля,
// потому что Rc103FieldPhysicsList рассчитан на готовый фотонный источник и
// декей-физику сознательно не тащит).
#pragma once

#include "G4VModularPhysicsList.hh"

class Rc103ShieldAirPhysicsList : public G4VModularPhysicsList {
 public:
  Rc103ShieldAirPhysicsList();
  ~Rc103ShieldAirPhysicsList() override = default;
};