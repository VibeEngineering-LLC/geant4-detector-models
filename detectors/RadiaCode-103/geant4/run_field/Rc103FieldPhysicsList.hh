// Физ-лист прогонов по полю: только ЭМ-физика.
//
// G4EmStandardPhysics_option4 — самый точный EM-конструктор Geant4 для задач
// низко-/среднеэнергетической гамма-дозиметрии.
//
// SetDefaultCutValue(0.05*mm) — здесь важна точность отклика в кристалле 10 мм
// и в тонких слоях обвязки (ESR-плёнка, SiPM, окно дисплея). В run_room/ стоял
// 1 мм, но там задача была транспорт через 80 см бетона — другая.
//
// Декей-физика НЕ регистрируется: первичные частицы этой программы — фотоны с
// разыгранной по CSV энергией, никаких ионов и распадов здесь нет.
#pragma once

#include "G4VModularPhysicsList.hh"

class Rc103FieldPhysicsList : public G4VModularPhysicsList {
 public:
  Rc103FieldPhysicsList();
  ~Rc103FieldPhysicsList() override = default;
};
