#include "Rc103FieldPhysicsList.hh"

#include "G4EmParameters.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4SystemOfUnits.hh"

#include <cstdio>

double Rc103FieldPhysicsList::gCutMm = 0.05;
std::string Rc103FieldPhysicsList::gDeexMode = "std";

Rc103FieldPhysicsList::Rc103FieldPhysicsList(double cutMm,
                                             const std::string& deexMode) {
  SetDefaultCutValue(cutMm * mm);
  RegisterPhysics(new G4EmStandardPhysics_option4());

  // Строго ПОСЛЕ конструктора option4 — см. предупреждение в заголовке.
  auto* p = G4EmParameters::Instance();
  if (deexMode == "deex" || deexMode == "max") {
    p->SetFluo(true);
    p->SetAuger(true);
    p->SetDeexcitationIgnoreCut(true);
  }
  if (deexMode == "max") {
    p->SetPixe(true);
  }

  gCutMm = cutMm;
  gDeexMode = deexMode;

  // Диагностика ASCII-строкой: читается из лога любой кодировкой консоли и
  // служит проверяемым артефактом того, что флаги ДЕЙСТВИТЕЛЬНО выставлены.
  std::fprintf(stdout,
               "Rc103FieldPhysicsList: cut_mm=%.4f deex=%s fluo=%d auger=%d "
               "pixe=%d ignore_cut=%d\n",
               cutMm, deexMode.c_str(), static_cast<int>(p->Fluo()),
               static_cast<int>(p->Auger()), static_cast<int>(p->Pixe()),
               static_cast<int>(p->DeexcitationIgnoreCut()));
}