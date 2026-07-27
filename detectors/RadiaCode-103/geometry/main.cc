// Отклик RadiaCode 101/102/103: спектр энерговыделения в кристалле CsI(Tl).
// Из спектра в постобработке извлекаются фотопиковая эффективность, полная
// эффективность счёта и отклик-матрица.
//
// Запуск:  rc_curves <макрос>            — сосуд Маринелли построен
//          rc_curves <макрос> bare       — «голый» прибор в воздухе
#include "RCDetector.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4LogicalVolume.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4Run.hh"
#include "G4RunManagerFactory.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UImanager.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

// --- Физика -----------------------------------------------------------------
class PhysList : public G4VModularPhysicsList {
public:
  PhysList() {
    RegisterPhysics(new G4EmStandardPhysics_option4());  // с флуоресценцией
    RegisterPhysics(new G4DecayPhysics());
    RegisterPhysics(new G4RadioactiveDecayPhysics());
    SetDefaultCutValue(0.05 * mm);
  }
};

// --- Источник ---------------------------------------------------------------
class Primary : public G4VUserPrimaryGeneratorAction {
  G4GeneralParticleSource fGPS;
public:
  void GeneratePrimaries(G4Event* e) override { fGPS.GeneratePrimaryVertex(e); }
};

// --- Накопление спектра -----------------------------------------------------
class RunAct : public G4UserRunAction {
public:
  static constexpr int    kBins = 3200;     // 1 кэВ на канал
  static constexpr double kBinKeV = 1.0;

  std::vector<long> fHist;
  long   fWithSignal = 0;
  double fSumEprim = 0;
  double fSampleCm3 = 0;
  double fDensity = 0;
  G4String fMatrix = "-";
  G4String fPart = "?";
  G4String fOut = "spectrum.csv";

  RunAct() : fHist(kBins + 1, 0) {}

  void BeginOfRunAction(const G4Run*) override {
    std::fill(fHist.begin(), fHist.end(), 0L);
    fWithSignal = 0;
    fSumEprim = 0;
    // Один процесс может гонять несколько beamOn с разными первичными
    // частицами (макрос нуклидов), поэтому подпись сбрасывается: иначе все
    // файлы прогона подписываются частицей ПЕРВОГО из них.
    fPart = "?";
  }

  void Fill(double edepKeV, double eprim) {
    fSumEprim += eprim;
    if (edepKeV <= 0) return;
    ++fWithSignal;
    int b = static_cast<int>(edepKeV / kBinKeV);
    if (b > kBins) b = kBins;             // последний канал — переполнение
    ++fHist[b];
  }

  void EndOfRunAction(const G4Run* run) override {
    const long N = run->GetNumberOfEvent();
    if (N == 0) return;
    const double eMean = fSumEprim / N;   // кэВ

    FILE* f = std::fopen(fOut.c_str(), "w");
    if (!f) {
      G4cerr << "!! не открыть " << fOut << G4endl;
      return;
    }
    std::fprintf(f, "# RadiaCode 101/102/103, CsI(Tl) 10x10x10 mm\n");
    std::fprintf(f, "# particle = %s\n", fPart.c_str());
    std::fprintf(f, "# E_prim_keV = %.4f\n", eMean);
    std::fprintf(f, "# N_primaries = %ld\n", N);
    std::fprintf(f, "# N_with_signal = %ld\n", fWithSignal);
    std::fprintf(f, "# sample_cm3 = %.3f\n", fSampleCm3);
    std::fprintf(f, "# matrix = %s\n", fMatrix.c_str());
    std::fprintf(f, "# density_gcm3 = %.4f\n", fDensity);
    std::fprintf(f, "# bin_keV = %.3f  (последний канал = переполнение)\n", kBinKeV);
    std::fprintf(f, "E_keV,counts\n");
    for (int i = 0; i <= kBins; ++i)
      if (fHist[i]) std::fprintf(f, "%.1f,%ld\n", (i + 0.5) * kBinKeV, fHist[i]);
    std::fclose(f);

    G4cout << "RESULT E_keV= " << eMean << " N= " << N
           << " hits= " << fWithSignal
           << " eff_total= " << double(fWithSignal) / N
           << " file= " << fOut << G4endl;
  }
};

// --- Энерговыделение за событие ---------------------------------------------
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

// --- Команда задания выходного файла ----------------------------------------
class OutMessenger : public G4UImessenger {
  RunAct* fRun;
  G4UIdirectory* fDir;
  G4UIcmdWithAString* fCmd;
public:
  explicit OutMessenger(RunAct* r) : fRun(r) {
    fDir = new G4UIdirectory("/rc/");
    fDir->SetGuidance("RadiaCode: управление выводом");
    fCmd = new G4UIcmdWithAString("/rc/outFile", this);
    fCmd->SetGuidance("Файл CSV для спектра следующего прогона");
    fCmd->AvailableForStates(G4State_Idle, G4State_PreInit);
  }
  ~OutMessenger() override { delete fCmd; delete fDir; }
  void SetNewValue(G4UIcommand*, G4String v) override { fRun->fOut = v; }
};

// ---------------------------------------------------------------------------
// Запуск: rc_curves <макрос> [full|bare|empty] [матрица] [плотность] [сосуд]
//                            [пластик]
//   full    — сосуд с пробой (по умолчанию вода 1.0 г/см³)
//   bare    — прибор в воздухе, без сосуда
//   empty   — сосуд построен, но пустой: холостой отсчёт (фон)
//   сосуд   — m200 (по умолчанию) | m500
//   пластик — PLA (по умолчанию: из него напечатаны измеренные экземпляры)
//             | PETG (Pet-G рекомендует автор моделей)
// Плотность, матрица и сосуд задаются до построения геометрии, поэтому на один
// процесс приходится одна комбинация — сетку гоняет драйвер run_grid.py.
int main(int argc, char** argv) {
  const std::string mode = (argc > 2) ? argv[2] : "full";
  const bool bare = (mode == "bare");

  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  auto* det = new RCDetector(!bare);
  const G4String matrix = (mode == "empty") ? "air"
                        : (argc > 3 ? argv[3] : "water");
  const double rho = (argc > 4) ? std::atof(argv[4]) : 1.0;
  det->fVes = VesselGeom::Preset(argc > 5 ? argv[5] : "m200");
  det->fVes.sampleMatrix = matrix;
  det->fVes.sampleDensity = rho;
  if (argc > 6) {
    const std::string p = argv[6];
    if (p != "PLA" && p != "PETG") {
      G4cout << "пластик: ожидается PLA или PETG, получено " << p << G4endl;
      return 2;
    }
    det->fVes.plasticMat = p;
  }
  rm->SetUserInitialization(det);
  rm->SetUserInitialization(new PhysList());
  rm->SetUserAction(new Primary());

  auto* runAct = new RunAct();
  rm->SetUserAction(runAct);
  auto* evtAct = new EventAct(runAct);
  rm->SetUserAction(evtAct);
  auto* mess = new OutMessenger(runAct);

  rm->Initialize();
  runAct->fSampleCm3 = det->fSampleVolumeCm3;
  runAct->fMatrix = (mode == "bare") ? "-" : det->fVes.sampleMatrix;
  runAct->fDensity = (mode == "bare") ? 0.0 : det->fVes.sampleDensity;
  rm->SetUserAction(new Stepping(evtAct, det->fCrystalLV));

  auto* ui = G4UImanager::GetUIpointer();
  if (argc > 1) ui->ApplyCommand(G4String("/control/execute ") + argv[1]);

  delete mess;
  delete rm;
  return 0;
}
