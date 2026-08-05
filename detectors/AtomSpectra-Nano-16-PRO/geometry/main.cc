// Отклик AtomSpectra Nano 16 PRO: спектр энерговыделения в бруске CsI(Tl)
// 18 x 15 x 57 мм внутри алюминиевой экструзии.
//
// Запуск:  asn16 <макрос>
//
// Положение источника задаётся макросом через GPS. Опорные плоскости, от
// которых отмеряется расстояние, печатает сама программа при старте
// (ReportPlanes): для точечной геометрии «10 см от торца» это НАРУЖНАЯ
// поверхность передней крышки, а не грань кристалла.
#include "ASN16Detector.hh"

// Отпечаток исходников, запечённый в бинарник (provenance.cmake генерирует его
// перед каждой сборкой). Идёт в шапку каждого выходного спектра: без него
// вопрос «этот спектр посчитан ТЕКУЩЕЙ геометрией?» отвечался бы по mtime.
#if defined(__has_include)
#  if __has_include("asn16_provenance.hh")
#    include "asn16_provenance.hh"
#  endif
#endif
#ifndef ASN16_SRC_SHA1
#  define ASN16_SRC_SHA1 "БЕЗ-ШТАМПА"
#  define ASN16_GIT_DESCRIBE "БЕЗ-ШТАМПА"
#endif

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
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <system_error>
#include <vector>

class PhysList : public G4VModularPhysicsList {
public:
  PhysList() {
    RegisterPhysics(new G4EmStandardPhysics_option4());
    RegisterPhysics(new G4DecayPhysics());
    RegisterPhysics(new G4RadioactiveDecayPhysics());
    // Порог рождения вторичных. У Гамма-1С стоит 0,05 мм и проверен на мягком
    // крае (задача 100). Здесь лицевой стек тоньше (0,571 г/см² против
    // торцевой колонны Гамма-1С), поэтому тот же порог тем более достаточен.
    SetDefaultCutValue(0.05 * mm);
  }
};

class Primary : public G4VUserPrimaryGeneratorAction {
  G4GeneralParticleSource fGPS;
public:
  void GeneratePrimaries(G4Event* e) override { fGPS.GeneratePrimaryVertex(e); }

  // Доля телесного угла ФАКТИЧЕСКОГО розыгрыша, (1−cos θmax)/2, спрошенная у
  // самого генератора после исполнения макроса: прямой множитель на точечную
  // eps при конусном розыгрыше. Сообщает её тот, кто разыгрывал, а не
  // сторонняя таблица драйвера.
  double SolidAngleFrac() {
    auto* src = fGPS.GetCurrentSource();
    if (!src || !src->GetAngDist()) return 1.0;
    const double th = src->GetAngDist()->GetMaxTheta();
    return 0.5 * (1.0 - std::cos(th));
  }
};

// Время разрешения тракта: энерговыделения, разнесённые больше чем на столько,
// считаются РАЗНЫМИ срабатываниями спектрометра (обоснование — как у Гамма-1С:
// каскад приходит за наносекунды, звенья ряда разделены секундами и годами).
constexpr double kResolvingTimeNs = 1000.0;

class RunAct : public G4UserRunAction {
public:
  // 1 кэВ на канал, потолок 3700 — как у Гамма-1С, чтобы верхние узлы сетки и
  // сумм-пики не уезжали в канал переполнения. Шкала линейная и без уширения:
  // приборное разрешение навешивается в постобработке.
  static constexpr int kBins = 3700;
  static constexpr double kBinKeV = 1.0;

  std::vector<long> fHist;
  std::vector<long> fEmit;   // гамма, ИСПУЩЕННЫЕ при распаде: выход линии
  long fWithSignal = 0;
  double fSumEprim = 0;
  G4String fPart = "?";
  G4String fArgs = "?";
  double fSolidAngleFrac = 1.0;
  Primary* fPrimary = nullptr;
  ASN16Detector* fDet = nullptr;
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

  void FillPrimary(double eprim) { fSumEprim += eprim; }

  void Fill(double edepKeV) {
    if (edepKeV <= 0) return;
    ++fWithSignal;
    int b = static_cast<int>(edepKeV / kBinKeV);
    if (b > kBins) b = kBins;
    ++fHist[b];
  }

  void EndOfRunAction(const G4Run* run) override {
    const long N = run->GetNumberOfEvent();
    if (N == 0) return;
    if (fPrimary) fSolidAngleFrac = fPrimary->SolidAngleFrac();
    // Каталог вывода создаётся здесь, а не считается существующим: макросы
    // пишут в подкаталоги (spectra/, spectra_face/), и без этого прогон на
    // 25 узлов отрабатывал ЦЕЛИКОМ, возвращал ноль и не оставлял ни одного
    // файла — тихий отказ, найденный независимым аудитом 05.08.2026.
    {
      const std::size_t s = fOut.find_last_of("/\\");
      if (s != G4String::npos) {
        std::error_code ec;
        std::filesystem::create_directories(fOut.substr(0, s), ec);
      }
    }
    FILE* f = std::fopen(fOut.c_str(), "w");
    if (!f) {
      // Не предупреждение, а аварийный останов: прежде здесь стоял return, и
      // прогон завершался с кодом 0 без единого файла.
      G4Exception("RunAct::EndOfRunAction", "ASN16_OUT", FatalException,
                  ("не открыть файл вывода " + fOut).c_str());
      return;
    }
    std::fprintf(f, "# ATOMSPECTRA NANO 16 PRO, CsI(Tl) 18x15x57 mm\n");
    std::fprintf(f, "# src_sha1 = %s\n", ASN16_SRC_SHA1);
    std::fprintf(f, "# git_describe = %s\n", ASN16_GIT_DESCRIBE);
    std::fprintf(f, "# build = %s %s\n", __DATE__, __TIME__);
    std::fprintf(f, "# run_args = %s\n", fArgs.c_str());
    // Опорные плоскости — в шапку спектра. Прибор анизотропен, и «10 см» без
    // указания, от какой поверхности и к какой грани, величины не задаёт.
    if (fDet) {
      std::fprintf(f, "# front_face_z_mm = %.3f  (наружная поверхность крышки)\n",
                   fDet->FrontFaceZ());
      std::fprintf(f, "# crystal_front_z_mm = %.3f\n", fDet->CrystalFrontZ());
      std::fprintf(f, "# work_face_y_mm = %.3f  (наружная поверхность стенки)\n",
                   fDet->WorkFaceY());
      std::fprintf(f, "# crystal_top_y_mm = %.3f\n", fDet->CrystalTopY());
    }
    std::fprintf(f, "# solid_angle_frac = %.8f\n", fSolidAngleFrac);
    std::fprintf(f, "# particle = %s\n", fPart.c_str());
    std::fprintf(f, "# E_prim_keV = %.4f\n", fSumEprim / N);
    std::fprintf(f, "# N_primaries = %ld\n", N);
    std::fprintf(f, "# N_with_signal = %ld\n", fWithSignal);
    std::fprintf(f, "# resolving_time_ns = %.0f\n", kResolvingTimeNs);
    std::fprintf(f, "# bin_keV = %.3f  (последний канал = переполнение)\n",
                 kBinKeV);
    std::fprintf(f, "E_keV,counts\n");
    for (int i = 0; i <= kBins; ++i)
      if (fHist[i]) std::fprintf(f, "%.1f,%ld\n", (i + 0.5) * kBinKeV, fHist[i]);
    std::fclose(f);

    long emitted = 0;
    for (long c : fEmit) emitted += c;
    if (emitted > 0) {
      G4String en = fOut;
      const size_t dot = en.rfind('.');
      en = (dot == G4String::npos ? en : en.substr(0, dot)) + "_emit.csv";
      FILE* g = std::fopen(en.c_str(), "w");
      if (g) {
        std::fprintf(g, "# гамма, испущенные при распаде, на %ld распадов\n", N);
        std::fprintf(g, "# src_sha1 = %s\n", ASN16_SRC_SHA1);
        std::fprintf(g, "# git_describe = %s\n", ASN16_GIT_DESCRIBE);
        std::fprintf(g, "# build = %s %s\n", __DATE__, __TIME__);
        std::fprintf(g, "# N_primaries = %ld\n", N);
        std::fprintf(g, "E_keV,counts\n");
        for (int i = 0; i <= kBins; ++i)
          if (fEmit[i]) std::fprintf(g, "%.1f,%ld\n", (i + 0.5) * kBinKeV,
                                     fEmit[i]);
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

// Событие Geant4 — это НЕ событие спектрометра: энерговыделения собираются с
// отметкой глобального времени и разбиваются на группы по времени разрешения.
class EventAct : public G4UserEventAction {
  RunAct* fRun;
public:
  std::vector<std::pair<double, double>> fDep;   // (нс, кэВ)
  explicit EventAct(RunAct* r) : fRun(r) {}
  void BeginOfEventAction(const G4Event*) override { fDep.clear(); }
  void EndOfEventAction(const G4Event* e) override {
    double ep = 0;
    if (e->GetNumberOfPrimaryVertex() > 0) {
      auto* p = e->GetPrimaryVertex(0)->GetPrimary(0);
      ep = p->GetKineticEnergy() / keV;
      if (fRun->fPart == "?" && p->GetParticleDefinition())
        fRun->fPart = p->GetParticleDefinition()->GetParticleName();
    }
    fRun->FillPrimary(ep);
    if (fDep.empty()) return;

    std::sort(fDep.begin(), fDep.end());
    double sum = fDep[0].second, t0 = fDep[0].first;
    for (size_t i = 1; i < fDep.size(); ++i) {
      if (fDep[i].first - t0 > kResolvingTimeNs) {
        fRun->Fill(sum);
        sum = 0;
      }
      t0 = fDep[i].first;
      sum += fDep[i].second;
    }
    fRun->Fill(sum);
  }
};

class Stepping : public G4UserSteppingAction {
  EventAct* fEvt;
  const G4LogicalVolume* fCry;
public:
  Stepping(EventAct* ev, const G4LogicalVolume* c) : fEvt(ev), fCry(c) {}
  void UserSteppingAction(const G4Step* s) override {
    const double e = s->GetTotalEnergyDeposit();
    if (e <= 0) return;
    auto* h = s->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    if (h && h->GetLogicalVolume() == fCry)
      fEvt->fDep.emplace_back(s->GetPreStepPoint()->GetGlobalTime() / ns,
                              e / keV);
  }
};

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
    fDir = new G4UIdirectory("/asn16/");
    fDir->SetGuidance("AtomSpectra Nano 16 PRO: управление выводом");
    fCmd = new G4UIcmdWithAString("/asn16/outFile", this);
    fCmd->SetGuidance("Файл CSV для спектра следующего прогона");
    fCmd->AvailableForStates(G4State_Idle, G4State_PreInit);
  }
  ~OutMessenger() override { delete fCmd; delete fDir; }
  void SetNewValue(G4UIcommand*, G4String v) override { fRun->fOut = v; }
};

int main(int argc, char** argv) {
  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  auto* det = new ASN16Detector();
  rm->SetUserInitialization(det);
  rm->SetUserInitialization(new PhysList());

  auto* primary = new Primary();
  rm->SetUserAction(primary);

  auto* runAct = new RunAct();
  runAct->fPrimary = primary;
  runAct->fDet = det;
  {
    std::string a;
    for (int i = 1; i < argc; ++i) {
      if (i > 1) a += " ";
      // ТОЛЬКО ИМЯ ФАЙЛА, без каталогов: полный argv занёс бы путь конкретной
      // машины в шапку каждого спектра, а спектр может попасть в репозиторий.
      std::string v(argv[i]);
      const size_t s = v.find_last_of("/\\");
      a += (s == std::string::npos) ? v : v.substr(s + 1);
    }
    runAct->fArgs = a;
  }
  rm->SetUserAction(runAct);
  auto* evtAct = new EventAct(runAct);
  rm->SetUserAction(evtAct);
  auto* mess = new OutMessenger(runAct);

  rm->Initialize();
  det->ReportPlanes();
  det->ReportMasses();
  rm->SetUserAction(new Stepping(evtAct, det->fCrystalLV));
  rm->SetUserAction(new Tracking(runAct));

  auto* ui = G4UImanager::GetUIpointer();
  if (argc > 1) ui->ApplyCommand(G4String("/control/execute ") + argv[1]);

  delete mess;
  delete rm;
  return 0;
}
