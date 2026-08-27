// Физ-лист прогонов поля ЕРН: G4EmStandardPhysics_option4 + G4DecayPhysics +
// G4RadioactiveDecayPhysics. Без двух последних ион-первичка просто стоит на
// месте и не распадается (метод 1 без них не работает вовсе).
//
// Порог продукции — ключ CLI cut=<мм>, дефолт 1.0 мм (спека). ⚠ ЗНАЙ ЦЕНУ:
// порог 1 мм в бетоне/кирпиче отсекает электроны примерно ниже 350 кэВ — они
// не трекаются, энергия кладётся локально, и ТОРМОЗНОГО ИЗЛУЧЕНИЯ ОТ НИХ НЕ
// БУДЕТ. Существующая серия results/wallion/wf_m1_*.csv считалась при 0,05 мм
// (geometry/wallfield.cc, ионная ветка Phys), поэтому прямое сравнение с ней
// корректно только если разница по порогу измерена, а не предположена.
#pragma once

#include "G4VModularPhysicsList.hh"

class Rc103RoomFieldPhysicsList : public G4VModularPhysicsList {
 public:
  explicit Rc103RoomFieldPhysicsList(double cutMm);
  ~Rc103RoomFieldPhysicsList() override = default;
};
