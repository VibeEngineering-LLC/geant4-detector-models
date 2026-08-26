// Физ-лист: G4EmStandardPhysics_option4 (точность на низких энергиях, важно
// для CsI) + G4DecayPhysics + G4RadioactiveDecayPhysics — тот же набор, что
// использует методика skills/geant4-spectrum-pipeline/SKILL.md (раздел
// «Архитектура: два метода»), чтобы это приложение было физически
// совместимо со следующим этапом конвейера (МК-шаблоны распадов) без смены
// физ-листа. Для ТЕКУЩЕГО smoke-теста (одиночный гамма-квант 661.657 кэВ из
// GPS, не распад ядра из /gps/ion) декей-физика не задействуется активно —
// директива /process/had/rdm/thresholdForVeryLongDecayTime из методики
// пропущена намеренно (актуальна только при генерации радиоактивного
// источника через /gps/ion, см. задание).
#pragma once

#include "G4VModularPhysicsList.hh"

class Rc103RunPhysicsList : public G4VModularPhysicsList {
 public:
  Rc103RunPhysicsList();
  ~Rc103RunPhysicsList() override = default;
};
