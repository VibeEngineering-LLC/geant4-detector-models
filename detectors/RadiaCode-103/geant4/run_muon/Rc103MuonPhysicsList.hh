// Физ-лист мюонного прогона: FTFP_BERT, а НЕ G4EmStandardPhysics_option4.
//
// Обоснование (спека _spec_run_muon.md, раздел "Физ-лист"): мюону нужны
// ионизация, тормозное излучение, образование пар и распад, плюс адронное
// взаимодействие. Чисто ЭМ-конструктор run_field/ этого не даёт. Это
// единственное место, где run_muon расходится с run_field по физике.
//
// SetDefaultCutValue(0.05*mm) — как в run_field: важна точность отклика в
// кристалле 10 мм и в тонких слоях обвязки (ESR-плёнка, SiPM, окно дисплея).
#pragma once

#include "FTFP_BERT.hh"

class Rc103MuonPhysicsList : public FTFP_BERT {
 public:
  Rc103MuonPhysicsList();
  ~Rc103MuonPhysicsList() override = default;
};
