#include "Rc103MuonPhysicsList.hh"

#include "G4SystemOfUnits.hh"

// FTFP_BERT(ver=0) — молчаливая инициализация, иначе на каждый прогон
// печатается полная таблица конструкторов физики.
Rc103MuonPhysicsList::Rc103MuonPhysicsList() : FTFP_BERT(0) {
  SetDefaultCutValue(0.05 * mm);
}
