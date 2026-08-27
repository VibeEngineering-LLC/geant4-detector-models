// Физ-лист поля помещения: G4EmStandardPhysics_option4, БЕЗ декей-физики —
// источник уже задан таблицей линий (готовые энергии/выходы, wallfield.cc
// линейный режим по умолчанию), не /gps/ion, декей-физика не нужна (сравни
// с Rc103RunPhysicsList, который её включает про запас совместимости с
// будущим МК-шаблонным этапом — здесь такой связи нет).
// Cut 1мм — как в wallfield.cc линейном режиме (важен транспорт гамма через
// толщу бетона, не точность по электронам).
#pragma once

#include "G4VModularPhysicsList.hh"

class Rc103RoomPhysicsList : public G4VModularPhysicsList {
 public:
  Rc103RoomPhysicsList();
  ~Rc103RoomPhysicsList() override = default;
};
