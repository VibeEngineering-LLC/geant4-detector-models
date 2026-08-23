// Отклик RadiaCode 101/102/103: спектр энерговыделения в кристалле CsI(Tl).
// Из спектра в постобработке извлекаются фотопиковая эффективность, полная
// эффективность счёта и отклик-матрица.
//
// Запуск:  rc_curves <макрос>            — сосуд Маринелли построен
//          rc_curves <макрос> bare       — «голый» прибор в воздухе
#include "RCDetector.hh"

// Отпечаток исходников, запекаемый в бинарник общим генератором
// common/cmake/provenance.cmake (то же правило, что у Гамма-1С).
#if defined(__has_include)
#  if __has_include("rc_provenance.hh")
#    include "rc_provenance.hh"
#  endif
#endif
#ifndef RC_SRC_SHA1
// Сборка мимо CMake — не запрещаем, но помечаем, чтобы такой спектр нельзя
// было принять за прослеженный.
#  define RC_SRC_SHA1 "БЕЗ-ШТАМПА"
#  define RC_GIT_DESCRIBE "БЕЗ-ШТАМПА"
#endif

#include "G4Event.hh"
#include "G4Gamma.hh"
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
#include "G4Track.hh"
#include "G4VProcess.hh"
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

// --- Разложение отклика по каналам ------------------------------------------
// ПЕРЕНЕСЕНО 21.08.2026 из
//   detectors/AtomSpectra-Nano-16-PRO/geometry/main.cc  (enum Chan, Channel(),
//   пометки в Stepping, запись _chan.csv)
// по указанию оператора «посмотри эти составляющие для RC103». Импортировать
// нельзя — это другой бинарник; поэтому копия, и здесь путь оригинала и дата
// (§33). Правило приоритета менять ТОЛЬКО синхронно с донором: разложения
// сравнимы между приборами лишь при одинаковом правиле.
//
// Приоритет: 1) было рождение пар; 2) первичный квант вылетел после комптона —
// однократного или многократного; 3) иначе вылетел характеристический рентген;
// 4) иначе вылетел тормозной; 5) иначе ничего не вылетело (фотоэффект либо
// комптон с последующим поглощением).
enum Chan {
  kChPhoto = 0,     // фотоэффект, ничего не вылетело
  kChComptFull,     // комптон(ы) и поглощение, ничего не вылетело
  kChComptEsc1,     // однократный комптон, квант ушёл
  kChComptEscN,     // многократный комптон, квант ушёл
  kChXrayEsc,       // вылет характеристического рентгена
  kChBremsEsc,      // вылет тормозного кванта
  kChPairFull,      // пары, оба аннигиляционных кванта поглощены
  kChPairEsc1,      // пары, вылетел один квант 511 кэВ
  kChPairEsc2,      // пары, вылетели оба
  kChExternal,      // первичный квант в кристалле не взаимодействовал: энергию
                    // принесла вторичная частица из корпуса/отражателя/сосуда
  kChOther,         // остаточный: не должен населяться, служит сторожем
  kNChan
};

static const char* const kChanName[kNChan] = {
  "photo", "compt_full", "compt_esc1", "compt_escN", "xray_esc",
  "brems_esc", "pair_full", "pair_esc1", "pair_esc2", "external", "other"
};

// --- Накопление спектра -----------------------------------------------------
class RunAct : public G4UserRunAction {
public:
  static constexpr int    kBins = 3200;     // 1 кэВ на канал
  static constexpr double kBinKeV = 1.0;

  static constexpr int    kYBins   = 7;
  static constexpr double kYBinMM  = 1.5;     // ширина Y-бина, мм
  static constexpr double kYMinMM  = -5.25;   // нижняя граница первого бина
  // Центры бинов при этих константах: -4.5,-3.0,-1.5,0.0,1.5,3.0,4.5 мм —
  // СОВПАДАЮТ с точками карты LCE(Y) из opticalcheck.cc (D-007), чтобы
  // постобработка сопоставляла бины напрямую, без интерполяции.

  std::vector<long> fHist;
  // Каналы взаимоисключающие и в сумме обязаны давать fHist — проверяется при
  // записи. Канал ставится В МОМЕНТ СОБЫТИЯ по истории процессов: из готового
  // спектра его восстановить нельзя, форма к тому времени уже сложена.
  std::vector<std::vector<long>> fChan;
  std::vector<std::vector<long>> fYHist;   // [kYBins][kBins+1]
  long   fWithSignal = 0;
  double fSumEprim = 0;
  double fSampleCm3 = 0;
  double fDensity = 0;
  G4String fMatrix = "-";
  G4String fPart = "?";
  G4String fOut = "spectrum.csv";

  RunAct() : fHist(kBins + 1, 0),
             fChan(kNChan, std::vector<long>(kBins + 1, 0)),
             fYHist(kYBins, std::vector<long>(kBins + 1, 0)) {}

  void BeginOfRunAction(const G4Run*) override {
    std::fill(fHist.begin(), fHist.end(), 0L);
    for (auto& c : fChan) std::fill(c.begin(), c.end(), 0L);
    for (auto& y : fYHist) std::fill(y.begin(), y.end(), 0L);
    fWithSignal = 0;
    fSumEprim = 0;
    // Один процесс может гонять несколько beamOn с разными первичными
    // частицами (макрос нуклидов), поэтому подпись сбрасывается: иначе все
    // файлы прогона подписываются частицей ПЕРВОГО из них.
    fPart = "?";
  }

  void Fill(double edepKeV, double eprim, int chan = -1, double yMeanMM = 0.0) {
    fSumEprim += eprim;
    if (edepKeV <= 0) return;
    ++fWithSignal;
    int b = static_cast<int>(edepKeV / kBinKeV);
    if (b > kBins) b = kBins;
    ++fHist[b];
    if (chan >= 0 && chan < kNChan) ++fChan[chan][b];

    int iy = static_cast<int>((yMeanMM - kYMinMM) / kYBinMM);
    if (iy < 0) iy = 0;
    if (iy >= kYBins) iy = kYBins - 1;
    ++fYHist[iy][b];
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
    // Отпечаток исходников, из которых собран ЭТОТ exe. Правило одно на все
    // детекторы (common/docs/method-rules.md): вопрос «из чего получено»
    // решается бинарником, а не временем файла.
    std::fprintf(f, "# src_sha1 = %s\n", RC_SRC_SHA1);
    std::fprintf(f, "# git_describe = %s\n", RC_GIT_DESCRIBE);
    std::fprintf(f, "# build = %s %s\n", __DATE__, __TIME__);
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

    // --- разложение отклика по каналам, ОТДЕЛЬНЫМ файлом ---
    // Отдельным, а не колонками основного спектра: формат «E_keV,counts»
    // читают все разборные скрипты дерева, и менять его ради нового
    // содержимого значит чинить их все.
    {
      G4String cn = fOut;
      const size_t dot = cn.rfind('.');
      cn = (dot == G4String::npos ? cn : cn.substr(0, dot)) + "_chan.csv";
      FILE* g = std::fopen(cn.c_str(), "w");
      if (g) {
        std::fprintf(g, "# разложение отклика по каналам, канал ставится по "
                        "истории процессов события\n");
        std::fprintf(g, "# src_sha1 = %s\n", RC_SRC_SHA1);
        std::fprintf(g, "# git_describe = %s\n", RC_GIT_DESCRIBE);
        std::fprintf(g, "# particle = %s\n", fPart.c_str());
        std::fprintf(g, "# E_prim_keV = %.4f\n", eMean);
        std::fprintf(g, "# N_primaries = %ld\n", N);
        std::fprintf(g, "# N_with_signal = %ld\n", fWithSignal);
        std::fprintf(g, "# matrix = %s\n", fMatrix.c_str());
        std::fprintf(g, "# bin_keV = %.3f\n", kBinKeV);
        std::fprintf(g, "# правило приоритета: пары -> вылет первичного после "
                        "комптона -> рентген -> тормозное -> без вылета\n");
        std::fprintf(g, "E_keV");
        for (int c = 0; c < kNChan; ++c) std::fprintf(g, ",%s", kChanName[c]);
        std::fprintf(g, "\n");
        for (int i = 0; i <= kBins; ++i) {
          if (!fHist[i]) continue;
          std::fprintf(g, "%.1f", (i + 0.5) * kBinKeV);
          for (int c = 0; c < kNChan; ++c)
            std::fprintf(g, ",%ld", fChan[c][i]);
          std::fprintf(g, "\n");
        }
        std::fclose(g);
      }
      // Каналы обязаны в сумме давать полный спектр. Не дают — правило
      // приоритета пропускает случай; молча потерянные события выглядели бы
      // как «канала нет», а не как дефект разложения.
      long sumChan = 0, sumHist = 0;
      for (int i = 0; i <= kBins; ++i) {
        sumHist += fHist[i];
        for (int c = 0; c < kNChan; ++c) sumChan += fChan[c][i];
      }
      if (sumChan != sumHist)
        G4cerr << "ВНИМАНИЕ: сумма каналов " << sumChan
               << " не равна спектру " << sumHist
               << " — правило приоритета неполно" << G4endl;
      else
        G4cout << "CHAN_OK sum= " << sumChan << " file= " << cn << G4endl;
    }

    {
      G4String yn = fOut;
      const size_t dot = yn.rfind('.');
      yn = (dot == G4String::npos ? yn : yn.substr(0, dot)) + "_ypos.csv";
      FILE* g = std::fopen(yn.c_str(), "w");
      if (g) {
        std::fprintf(g, "# разложение отклика по Y-координате взаимодействия "
                        "(энерговзвешенное среднее за событие), ось как в opticalcheck.cc\n");
        std::fprintf(g, "# src_sha1 = %s\n", RC_SRC_SHA1);
        std::fprintf(g, "# git_describe = %s\n", RC_GIT_DESCRIBE);
        std::fprintf(g, "# particle = %s\n", fPart.c_str());
        std::fprintf(g, "# E_prim_keV = %.4f\n", eMean);
        std::fprintf(g, "# N_primaries = %ld\n", N);
        std::fprintf(g, "# N_with_signal = %ld\n", fWithSignal);
        std::fprintf(g, "# y_bin_mm = %.2f, y_min_mm = %.2f, y_bins = %d\n",
                     kYBinMM, kYMinMM, kYBins);
        std::fprintf(g, "E_keV");
        for (int y = 0; y < kYBins; ++y)
          std::fprintf(g, ",y%.2f", kYMinMM + (y + 0.5) * kYBinMM);
        std::fprintf(g, "\n");
        for (int i = 0; i <= kBins; ++i) {
          long rowSum = 0;
          for (int y = 0; y < kYBins; ++y) rowSum += fYHist[y][i];
          if (!rowSum) continue;
          std::fprintf(g, "%.1f", (i + 0.5) * kBinKeV);
          for (int y = 0; y < kYBins; ++y) std::fprintf(g, ",%ld", fYHist[y][i]);
          std::fprintf(g, "\n");
        }
        std::fclose(g);

        long sumY = 0, sumHist2 = 0;
        for (int i = 0; i <= kBins; ++i) {
          sumHist2 += fHist[i];
          for (int y = 0; y < kYBins; ++y) sumY += fYHist[y][i];
        }
        if (sumY != sumHist2)
          G4cerr << "ВНИМАНИЕ: сумма Y-бинов " << sumY
                 << " не равна спектру " << sumHist2 << " — дефект биннинга Y" << G4endl;
        else
          G4cout << "YPOS_OK sum= " << sumY << " file= " << yn << G4endl;
      }
    }

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
  double fEdepY = 0;   // накопитель Edep_i * Y_mid_i по шагам, для взвешенного среднего

  // Признаки истории события для разложения по каналам. Ставятся в Stepping.
  int fFirst = 0;          // 1 phot, 2 compt, 3 conv, 0 ничего неупругого
  bool fHadRayl = false;   // было упругое рассеяние (энергии не оставляет)
  int fNCompt = 0;         // сколько раз первичный квант рассеялся в кристалле
  bool fHadConv = false;   // было рождение пар в кристалле
  int fNAnnihEsc = 0;      // сколько аннигиляционных квантов покинуло кристалл
  double fEBremEsc = 0;    // энергия вылетевшего тормозного, кэВ
  double fEXrayEsc = 0;    // энергия вылетевших прочих вторичных гамма, кэВ
  bool fPrimEsc = false;   // сам первичный квант вышел из кристалла

  explicit EventAct(RunAct* r) : fRun(r) {}
  void BeginOfEventAction(const G4Event*) override {
    fEdep = 0;
    fEdepY = 0;
    fFirst = 0; fHadRayl = false; fNCompt = 0; fHadConv = false;
    fNAnnihEsc = 0; fEBremEsc = 0; fEXrayEsc = 0; fPrimEsc = false;
  }

  // Канал по правилу приоритета (см. enum Chan). Возвращает индекс канала.
  int Channel() const {
    if (fHadConv)
      return fNAnnihEsc == 0 ? kChPairFull
           : (fNAnnihEsc == 1 ? kChPairEsc1 : kChPairEsc2);
    if (fPrimEsc && fNCompt > 0)
      return fNCompt == 1 ? kChComptEsc1 : kChComptEscN;
    if (fEXrayEsc > 0) return kChXrayEsc;
    if (fEBremEsc > 0) return kChBremsEsc;
    if (fFirst == 1) return kChPhoto;
    if (fNCompt > 0) return kChComptFull;
    if (fFirst == 0) return kChExternal;
    return kChOther;
  }

  void EndOfEventAction(const G4Event* e) override {
    double ep = 0;
    if (e->GetNumberOfPrimaryVertex() > 0) {
      auto* p = e->GetPrimaryVertex(0)->GetPrimary(0);
      ep = p->GetKineticEnergy() / keV;
      if (fRun->fPart == "?" && p->GetParticleDefinition())
        fRun->fPart = p->GetParticleDefinition()->GetParticleName();
    }
    const double yMean = (fEdep > 0) ? (fEdepY / fEdep) : 0.0;
    fRun->Fill(fEdep / keV, ep, Channel(), yMean);
  }
};

class Stepping : public G4UserSteppingAction {
  EventAct* fEvt;
  const G4LogicalVolume* fCry;
public:
  Stepping(EventAct* ev, const G4LogicalVolume* c) : fEvt(ev), fCry(c) {}
  void UserSteppingAction(const G4Step* s) override {
    auto* pre = s->GetPreStepPoint();
    auto* post = s->GetPostStepPoint();
    auto* h = pre->GetTouchableHandle()->GetVolume();
    if (!h || h->GetLogicalVolume() != fCry) return;

    // --- пометка канала: что произошло с квантом ВНУТРИ кристалла ---
    const G4Track* trk = s->GetTrack();
    if (trk->GetDefinition() == G4Gamma::Gamma()) {
      const G4VProcess* pr = post->GetProcessDefinedStep();
      const G4String pn = pr ? pr->GetProcessName() : G4String();
      if (trk->GetParentID() == 0) {
        if (pn == "compt") {
          ++fEvt->fNCompt;
          if (!fEvt->fFirst) fEvt->fFirst = 2;
        } else if (pn == "phot") {
          if (!fEvt->fFirst) fEvt->fFirst = 1;
        } else if (pn == "conv") {
          fEvt->fHadConv = true;
          if (!fEvt->fFirst) fEvt->fFirst = 3;
        } else if (pn == "Rayl") {
          // Рэлеевское рассеяние УПРУГОЕ: энергии не оставляет и первым
          // взаимодействием не считается (иначе «рэлей, затем фотоэффект»
          // уходит в остаточный канал — у донора это давало 5,5 % на 180 кэВ).
          fEvt->fHadRayl = true;
        }
      }
      // Выход кванта из кристалла. Квант, вернувшийся обратно из корпуса, здесь
      // уже посчитан вылетевшим — огрубление в пользу каналов вылета; доля
      // таких возвратов НЕ измерена.
      if (post->GetStepStatus() == fGeomBoundary) {
        auto* hp = post->GetTouchableHandle()->GetVolume();
        const bool outCry = !hp || hp->GetLogicalVolume() != fCry;
        const double ek = post->GetKineticEnergy() / keV;
        if (outCry && ek > 0) {
          const G4VProcess* cp = trk->GetCreatorProcess();
          const G4String cn = cp ? cp->GetProcessName() : G4String();
          if (cn == "annihil") ++fEvt->fNAnnihEsc;
          else if (cn == "eBrem") fEvt->fEBremEsc += ek;
          else if (trk->GetParentID() == 0) fEvt->fPrimEsc = true;
          else fEvt->fEXrayEsc += ek;
        }
      }
    }

    const double de = s->GetTotalEnergyDeposit();
    fEvt->fEdep += de;
    if (de > 0) {
      const double yMid = 0.5 * (pre->GetPosition().y() + post->GetPosition().y()) / mm;
      fEvt->fEdepY += de * yMid;
    }
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
  // #SHIELD-16 диагностика: rc_curves <мак> [full|bare|empty] [матрица] [плотность]
  //   [сосуд] [пластик] [diag=<маска>] — маска бит 0=case, 1=reflector, 2=internals -> воздух.
  if (argc > 7) {
    const std::string dp = argv[7];
    const std::string pfx = "diag=";
    if (dp.rfind(pfx, 0) == 0) det->fDiagMask = static_cast<unsigned>(std::atoi(dp.c_str() + pfx.size()));
  }
  if (det->fDiagMask) G4cout << "!! diag: fDiagMask=" << det->fDiagMask << " (материалы -> воздух)" << G4endl;
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
