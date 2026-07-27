// Отклик спектрометра ГАММА-1С: спектр энерговыделения в кристалле NaI(Tl)
// 63x63 внутри свинцового экрана-защиты.
//
// Запуск:  g1s <макрос> [режим] [плотность]
//   shield — экран собран, крышка закрыта (по умолчанию)
//   open   — экран собран, крышка снята (точечная геометрия 25 см)
//   bare   — устройство детектирования в воздухе, без экрана
//   vessel[:сосуд] — экран закрыт, на головке сосуд комплекта:
//            marinelli (по умолчанию) | denta | petri
//
// Позиционные аргументы: 3 — плотность матрицы, 4 — матрица
// (OISN16 | water | risn379), 5 — объём засыпки мл, 6 — плотность MgO,
// 7 — глубина колодца маринелли мм.
#include "G1SDetector.hh"

#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4GeneralParticleSource.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4Run.hh"
#include "G4RunManagerFactory.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4UImanager.hh"
#include "G4UImessenger.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

class PhysList : public G4VModularPhysicsList {
public:
  PhysList() {
    RegisterPhysics(new G4EmStandardPhysics_option4());   // с флуоресценцией:
    RegisterPhysics(new G4DecayPhysics());                // нужна для ХРИ Pb/Cd/Cu
    RegisterPhysics(new G4RadioactiveDecayPhysics());
    SetDefaultCutValue(0.05 * mm);
  }
};

class Primary : public G4VUserPrimaryGeneratorAction {
  G4GeneralParticleSource fGPS;
public:
  void GeneratePrimaries(G4Event* e) override { fGPS.GeneratePrimaryVertex(e); }
};

// --- Накопление спектра -----------------------------------------------------
// 1024 канала — как у самого спектрометра (паспорт, п. 2.12), но шкала здесь
// линейная по энергии и без уширения: разрешение навешивается в постобработке.
class RunAct : public G4UserRunAction {
public:
  static constexpr int kBins = 3200;      // 1 кэВ на канал, диапазон 50–3000
  static constexpr double kBinKeV = 1.0;

  std::vector<long> fHist;
  // Спектр ИСПУЩЕННЫХ при распаде гамма-квантов: сколько квантов каждой
  // энергии рождается на один распад. Нужен для поправки на каскадное
  // суммирование — чтобы выход линии p_gamma брать из ТОЙ ЖЕ базы
  // PhotonEvaporation, что и транспорт, а не вписывать из справочника руками.
  std::vector<long> fEmit;
  long fWithSignal = 0;
  double fSumEprim = 0;
  G4String fPart = "?";
  G4String fMode = "shield";
  G4String fOut = "spectrum.csv";

  RunAct() : fHist(kBins + 1, 0), fEmit(kBins + 1, 0) {}

  void BeginOfRunAction(const G4Run*) override {
    std::fill(fHist.begin(), fHist.end(), 0L);
    std::fill(fEmit.begin(), fEmit.end(), 0L);
    fWithSignal = 0;
    fSumEprim = 0;
    fPart = "?";
  }

  void FillEmit(double eKeV) {
    if (eKeV <= 0) return;
    int b = static_cast<int>(eKeV / kBinKeV);
    if (b > kBins) b = kBins;
    ++fEmit[b];
  }

  void Fill(double edepKeV, double eprim) {
    fSumEprim += eprim;
    if (edepKeV <= 0) return;
    ++fWithSignal;
    int b = static_cast<int>(edepKeV / kBinKeV);
    if (b > kBins) b = kBins;
    ++fHist[b];
  }

  void EndOfRunAction(const G4Run* run) override {
    const long N = run->GetNumberOfEvent();
    if (N == 0) return;
    FILE* f = std::fopen(fOut.c_str(), "w");
    if (!f) {
      G4cerr << "!! не открыть " << fOut << G4endl;
      return;
    }
    std::fprintf(f, "# GAMMA-1S, UDS-GC-63x63-USB, NaI(Tl) 63x63 mm\n");
    std::fprintf(f, "# mode = %s\n", fMode.c_str());
    std::fprintf(f, "# particle = %s\n", fPart.c_str());
    std::fprintf(f, "# E_prim_keV = %.4f\n", fSumEprim / N);
    std::fprintf(f, "# N_primaries = %ld\n", N);
    std::fprintf(f, "# N_with_signal = %ld\n", fWithSignal);
    std::fprintf(f, "# bin_keV = %.3f  (последний канал = переполнение)\n", kBinKeV);
    std::fprintf(f, "E_keV,counts\n");
    for (int i = 0; i <= kBins; ++i)
      if (fHist[i]) std::fprintf(f, "%.1f,%ld\n", (i + 0.5) * kBinKeV, fHist[i]);
    std::fclose(f);

    // Спектр испускания — отдельным файлом, если он не пуст (прогон распада)
    long emitted = 0;
    for (long c : fEmit) emitted += c;
    if (emitted > 0) {
      G4String en = fOut;
      const size_t dot = en.rfind('.');
      en = (dot == G4String::npos ? en : en.substr(0, dot)) + "_emit.csv";
      FILE* g = std::fopen(en.c_str(), "w");
      if (g) {
        std::fprintf(g, "# гамма, испущенные при распаде, на %ld распадов\n", N);
        std::fprintf(g, "# N_primaries = %ld\n", N);
        std::fprintf(g, "E_keV,counts\n");
        for (int i = 0; i <= kBins; ++i)
          if (fEmit[i]) std::fprintf(g, "%.1f,%ld\n", (i + 0.5) * kBinKeV, fEmit[i]);
        std::fclose(g);
      }
      G4cout << "EMIT всего " << emitted << " квантов на " << N
             << " распадов -> " << en << G4endl;
    }

    G4cout << "RESULT E_keV= " << fSumEprim / N << " N= " << N
           << " hits= " << fWithSignal
           << " eff_total= " << double(fWithSignal) / N
           << " file= " << fOut << G4endl;
  }
};

class EventAct : public G4UserEventAction {
  RunAct* fRun;
public:
  double fEdep = 0;
  explicit EventAct(RunAct* r) : fRun(r) {}
  void BeginOfEventAction(const G4Event*) override { fEdep = 0; }
  void EndOfEventAction(const G4Event* e) override {
    double ep = 0;
    if (e->GetNumberOfPrimaryVertex() > 0) {
      auto* p = e->GetPrimaryVertex(0)->GetPrimary(0);
      ep = p->GetKineticEnergy() / keV;
      if (fRun->fPart == "?" && p->GetParticleDefinition())
        fRun->fPart = p->GetParticleDefinition()->GetParticleName();
    }
    fRun->Fill(fEdep / keV, ep);
  }
};

class Stepping : public G4UserSteppingAction {
  EventAct* fEvt;
  const G4LogicalVolume* fCry;
public:
  Stepping(EventAct* ev, const G4LogicalVolume* c) : fEvt(ev), fCry(c) {}
  void UserSteppingAction(const G4Step* s) override {
    auto* h = s->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    if (h && h->GetLogicalVolume() == fCry)
      fEvt->fEdep += s->GetTotalEnergyDeposit();
  }
};

// Учёт квантов, РОЖДЁННЫХ радиоактивным распадом (а не рассеянием). Так выход
// линии p_gamma приходит из той же базы PhotonEvaporation, что и транспорт.
class Tracking : public G4UserTrackingAction {
  RunAct* fRun;
public:
  explicit Tracking(RunAct* r) : fRun(r) {}
  void PreUserTrackingAction(const G4Track* t) override {
    if (t->GetDefinition() != G4Gamma::Definition()) return;
    const G4VProcess* p = t->GetCreatorProcess();
    if (!p) return;                       // первичная частица, не распад
    const G4String& n = p->GetProcessName();
    if (n == "RadioactiveDecay" || n == "Radioactivation")
      fRun->FillEmit(t->GetKineticEnergy() / keV);
  }
};

class OutMessenger : public G4UImessenger {
  RunAct* fRun;
  G4UIdirectory* fDir;
  G4UIcmdWithAString* fCmd;
public:
  explicit OutMessenger(RunAct* r) : fRun(r) {
    fDir = new G4UIdirectory("/g1s/");
    fDir->SetGuidance("ГАММА-1С: управление выводом");
    fCmd = new G4UIcmdWithAString("/g1s/outFile", this);
    fCmd->SetGuidance("Файл CSV для спектра следующего прогона");
    fCmd->AvailableForStates(G4State_Idle, G4State_PreInit);
  }
  ~OutMessenger() override { delete fCmd; delete fDir; }
  void SetNewValue(G4UIcommand*, G4String v) override { fRun->fOut = v; }
};

int main(int argc, char** argv) {
  const std::string mode = (argc > 2) ? argv[2] : "shield";

  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  auto* det = new G1SDetector();
  det->fWithShield = (mode != "bare");
  det->fShield.lidClosed = (mode != "open");
  // vessel[:сосуд] — marinelli (по умолчанию) | denta | petri
  det->fWithVessel = (mode.rfind("vessel", 0) == 0);
  if (det->fWithVessel) {
    const size_t c = mode.find(':');
    det->fVessel = VesselGeom::Preset(c == std::string::npos ? "marinelli"
                                                             : mode.substr(c + 1));
  }
  if (argc > 3) det->fVessel.sampleDensity = std::atof(argv[3]);
  if (argc > 4) det->fVessel.sampleMatrix = argv[4];   // OISN16 | water | risn379
  // Объём засыпки можно переопределить: у источников комплекта разные массы
  // при одном номинальном объёме кюветы, см. kit_inventory.py.
  if (argc > 5) det->fVessel.sampleCm3 = std::atof(argv[5]);
  // Два параметра, вынесенные по итогам сверки с ЛСРМ (см. STATUS.md):
  //   6-й — насыпная плотность отражателя MgO (мягкий край кривой);
  //   7-й — глубина колодца маринелли (главный подозреваемый в превышении
  //         расчёта над измерением на объёмном источнике).
  if (argc > 6) det->fHead.mgoDensity = std::atof(argv[6]);
  if (argc > 7) det->fVessel.wellDepth = std::atof(argv[7]);

  rm->SetUserInitialization(det);
  rm->SetUserInitialization(new PhysList());
  rm->SetUserAction(new Primary());

  auto* runAct = new RunAct();
  runAct->fMode = mode;
  rm->SetUserAction(runAct);
  auto* evtAct = new EventAct(runAct);
  rm->SetUserAction(evtAct);
  auto* mess = new OutMessenger(runAct);

  rm->Initialize();
  det->ReportMasses();
  rm->SetUserAction(new Stepping(evtAct, det->fCrystalLV));
  rm->SetUserAction(new Tracking(runAct));

  auto* ui = G4UImanager::GetUIpointer();
  if (argc > 1) ui->ApplyCommand(G4String("/control/execute ") + argv[1]);

  delete mess;
  delete rm;
  return 0;
}
