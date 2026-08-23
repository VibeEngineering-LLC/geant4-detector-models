// Свинцовая защита RadiaCode: прогоны и проверки.
//
// Имя файла НЕ pbshield.cc намеренно: Windows не различает регистр в именах
// файлов, и pbshield.cc затирал бы PbShield.cc. Ловушка стоила одного файла.
//
// Режимы (argv[1]):
//   geom   — построить, напечатать массы и габариты, прогнать штатную проверку
//            пересечений объёмов. Ничего не считает; нужен для приёмки
//            геометрии и для картинки.
//   resp   — ОТКЛИК ПОЛОСТИ: моноэнергетический изотропный источник на самой
//            границе полости (там же, где в дальнейшем будет граница свинца),
//            прибор+сосуд построены, свинца НЕТ. Даёт спектр энерговыделения
//            в кристалле на один пересекший границу квант данной энергии.
//            Экономия, ради которой это отдельный режим: полость ФИКСИРОВАНА
//            для всех толщин защиты, поэтому отклик считается ОДИН РАЗ на
//            сетке энергий, а не заново для каждой точки сетки толщин — эту
//            сетку сворачивают со спектром пересечений от режима trans.
//
//   trans  — стадия 1 (толща свинца). Небиасированная по умолчанию; флаг
//            bias включает G4ImportanceBiasing (importance растёт слой за
//            слоем К ПОЛОСТИ, шаг impstep=2.0 по умолчанию). Пишет CSV
//            пересечений границы полости С ВЕСОМ (10-я колонка) — без bias
//            вес всегда 1.0, с bias — вес расщеплённой/срулеченной истории.
//   replay — стадия 2: точное воспроизведение записей trans (позиция,
//            направление, энергия, ВЕС) как первичных в полной геометрии.
//   beam   — узкопучковая проверка (iii), см. ниже.
//   sample — собственная активность пробы (K-40/Cs-137) сквозь собранную
//            защиту; device+vessel+shield все построены, короткий путь,
//            биасинг не нужен. Параметр nuc=K40|Cs137.
//   muon   — космический мюон сквозь ВСЮ сборку (диск-источник над защитой,
//            cos²θ, спектр Гайссера — см. cosmicmu.cc); прямой прогон, без
//            биасинга (мюон не гасится экспоненциально, в отличие от гамма).
//
// Параметры — ключами вида имя=значение в любом порядке:
//   pb=50 cu=1.5 cd=1.2 nshell=8 rcav=50 hzcav=90 zcav=35
//   vessel=m200 matrix=organic rho=0.50 novessel noshield
//   (resp/trans/beam)  e=661.7 nprim=1000000 out=resp_661.7.csv
//   (trans, опц.)      bias impstep=2.0
//   (replay)           in=trans_output.csv
//   (sample)           nuc=K40|Cs137 nprim=... out=...
//   (muon)             nprim=... out=...
#if defined(__has_include)
#  if __has_include("rc_pbshield_provenance.hh")
#    include "rc_pbshield_provenance.hh"
#  endif
#endif
#ifndef RCPB_SRC_SHA1
#  define RCPB_SRC_SHA1 "БЕЗ-ШТАМПА"
#  define RCPB_GIT_DESCRIBE "БЕЗ-ШТАМПА"
#endif

#include "PbShield.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4StepLimiterPhysics.hh"   // P-009: исполнитель G4UserLimits
#include "G4UserLimits.hh"           // P-009: ограничение шага в тонких ячейках
#include "G4RadioactiveDecayPhysics.hh"
#include "FTFP_BERT.hh"
#include "Randomize.hh"
#include "G4LogicalVolume.hh"
#include "G4PrimaryParticle.hh"
#include "G4Run.hh"
#include "G4RunManagerFactory.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "G4Track.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4VProcess.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"

// Геометрический importance-биасинг (задача №2, реализовано 12.08.2026).
// Порядок вызовов сверен ДВАЖДЫ против реального исходника B01
// (gitlab.cern.ch/geant4/geant4, тег v11.2.1) — см. план работ, раздел
// «Статус реализации», п.1. Ключевая находка второй сверки: конструктору
// G4GeometrySampler можно и нужно передавать явный nullptr вместо мира
// (в B01 передаётся ещё-не-инициализированный указатель — латентный баг
// примера, не идиома для копирования); настоящий мир резолвится позже сам,
// при PrepareImportanceSampling()/Configure(), уже ПОСЛЕ rm->Initialize().
#include "G4GeometrySampler.hh"
#include "G4ImportanceBiasing.hh"
#include "G4IStore.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

class Phys : public G4VModularPhysicsList {
public:
  Phys() {
    RegisterPhysics(new G4EmStandardPhysics_option4());
    // Decay+RadioactiveDecay — нужны ТОЛЬКО для mode=sample (ионный источник
    // K-40/Cs-137, GPS /gps/particle ion), но регистрируются всегда: то же
    // сочетание, что в main.cc, лишние таблицы распада безвредны для гамма-
    // транспортных режимов (trans/replay/resp/beam), не портят их физику.
    RegisterPhysics(new G4DecayPhysics());
    RegisterPhysics(new G4RadioactiveDecayPhysics());
    // P-009 (17.08): ограничитель длины шага. Сам по себе НИЧЕГО не меняет —
    // работает только там, где логическому объёму явно задан G4UserLimits
    // (делается в doBiasReplay для тонких ячеек важности). Без него
    // SetMaxAllowedStep молча игнорируется: процесса, который бы его
    // исполнял, в списке нет. Для всех режимов без UserLimits поведение
    // остаётся побитово прежним.
    RegisterPhysics(new G4StepLimiterPhysics());
    // 0.05 мм: электроны локальны, важен транспорт гамма. То же значение,
    // что в основном расчёте прибора (main.cc).
    SetDefaultCutValue(0.05 * mm);
  }
};

bool KeyVal(const char* arg, const char* key, double* out) {
  const size_t n = std::strlen(key);
  if (std::strncmp(arg, key, n) != 0 || arg[n] != '=') return false;
  *out = std::atof(arg + n + 1);
  return true;
}
bool KeyStr(const char* arg, const char* key, std::string* out) {
  const size_t n = std::strlen(key);
  if (std::strncmp(arg, key, n) != 0 || arg[n] != '=') return false;
  *out = arg + n + 1;
  return true;
}

// --- источник GPS, управляется UI-командами (см. main() режима resp) --------
class Primary : public G4VUserPrimaryGeneratorAction {
  G4GeneralParticleSource fGPS;
public:
  void GeneratePrimaries(G4Event* e) override { fGPS.GeneratePrimaryVertex(e); }
};

// --- replay: точное воспроизведение записей стадии 1 ------------------------
// Родилось из провала проверки (iv): угловое распределение пересечений НЕ
// косинусное (см. план, раздел «3а») — свёртка с resp(E) незаконна. Честная
// замена: каждая запись стадии 1 стартует как первичная частица В ТОЧНОСТИ с
// тем положением, направлением и энергией, с которыми она реально пересекла
// границу полости. Циклический повтор записей допустИм и предусмотрен планом
// (раздел «Стоимость», п.3) при условии контроля корреляций — тот контроль
// делается отдельно на стороне анализа (бутстрэп по записям), не здесь.
// w — вес записи (Σ по всем записям / N_primaries_stage1 = несмещённая
// оценка вероятности пересечения). Без биасинга w=1.0 у всех записей —
// старые файлы trans без 10-й колонки читаются тем же кодом, вес по
// умолчанию 1.0 (побитово прежнее поведение, старые CSV не инвалидированы).
struct CrossRec { double eKeV, x, y, z, dx, dy, dz, w; };

std::vector<CrossRec> ReadCrossings(const std::string& path) {
  std::vector<CrossRec> v;
  std::ifstream f(path);
  if (!f) { std::cerr << "не открыть " << path << "\n"; return v; }
  std::string line;
  bool sawHeader = false;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#') continue;
    if (!sawHeader) { sawHeader = true; continue; }  // строка "E_keV,cosTheta,..."
    std::istringstream ss(line);
    std::vector<std::string> parts;
    std::string tok;
    while (std::getline(ss, tok, ',')) parts.push_back(tok);
    if (parts.size() < 9) continue;
    const double w = (parts.size() >= 10) ? std::atof(parts[9].c_str()) : 1.0;
    CrossRec r{std::atof(parts[0].c_str()), std::atof(parts[3].c_str()),
              std::atof(parts[4].c_str()), std::atof(parts[5].c_str()),
              std::atof(parts[6].c_str()), std::atof(parts[7].c_str()),
              std::atof(parts[8].c_str()), w};
    v.push_back(r);
  }
  return v;
}

class ReplayPrimary : public G4VUserPrimaryGeneratorAction {
  G4ParticleGun fGun{1};
  const std::vector<CrossRec>* fRecs;
  size_t fIdx = 0;
public:
  explicit ReplayPrimary(const std::vector<CrossRec>* recs) : fRecs(recs) {
    fGun.SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("gamma"));
  }
  void GeneratePrimaries(G4Event* e) override {
    const CrossRec& r = (*fRecs)[fIdx % fRecs->size()];
    ++fIdx;
    fGun.SetParticlePosition(G4ThreeVector(r.x, r.y, r.z) * mm);
    fGun.SetParticleMomentumDirection(G4ThreeVector(r.dx, r.dy, r.dz));
    fGun.SetParticleEnergy(r.eKeV * keV);
    fGun.SetParticleWeight(r.w);
    fGun.GeneratePrimaryVertex(e);
  }
};

// --- спектр энерговыделения в кристалле, тот же формат, что main.cc ---------
// Дословно тот же формат вывода (RunAct/EventAct/Stepping), что в основном
// расчёте прибора — analysis/rcspec.py читает оба одинаково. Не общий код с
// main.cc намеренно: shieldrun — отдельный исполняемый файл со своим
// провенансом, дублирование ~90 строк дешевле, чем общий заголовок на два
// несвязанных бинарника (см. конвенцию репо: wallfield.cc/mucalc.cc тоже
// самодостаточны).
class RunAct : public G4UserRunAction {
public:
  static constexpr int    kBins = 3200;     // 1 кэВ на канал
  static constexpr double kBinKeV = 1.0;

  // ВЗВЕШЕННЫЕ отсчёты (double, не long) — нужно для biased trans/replay:
  // одна симулированная история несёт вес w!=1 (расщепление importance-
  // биасингом), и бин обязан копить Σw, не число историй. Без биасинга
  // w всегда 1.0 (умолчание Fill), числа побитово те же, что были раньше
  // (int -> float одного и того же значения) — уже провалидированные resp/
  // trans/beam прогоны этим не затронуты.
  std::vector<double> fHist;
  // Σw² по каналам — для оценки ДОСТОВЕРНОСТИ взвешенной гистограммы.
  // Под биасингом «сумма весов» ничего не говорит о статистике: 0,73 может
  // быть и одной историей с весом 0,73, и тысячей по 0,0007. Эффективное
  // число выборок N_eff = (Σw)²/Σw², относительная погрешность канала
  // σ/значение = sqrt(Σw²)/Σw. Введено 13.08.2026 (задача №19): проверка
  // достоверности велась только в окне 662 кэВ, а низкоэнергетическая часть
  // спектра оказалась шумом (площадь группы ХРИ Pb выходила отрицательной).
  std::vector<double> fHist2;
  double fWithSignal = 0;
  double fSumEprim = 0;
  G4String fOut = "resp.csv";
  G4String fTag = "-";   // произвольная метка конфигурации в шапке файла

  RunAct() : fHist(kBins + 1, 0.0), fHist2(kBins + 1, 0.0) {}

  void BeginOfRunAction(const G4Run*) override {
    std::fill(fHist.begin(), fHist.end(), 0.0);
    std::fill(fHist2.begin(), fHist2.end(), 0.0);
    fWithSignal = 0;
    fSumEprim = 0;
  }

  // Разделено на две части (задача №13, 12.08.2026) — ОДИН раз на событие
  // (энергия первичной, диагностика) и, возможно, НЕСКОЛЬКО раз на событие
  // (депозит в кристалле): под doBiasReplay одно replay-событие может
  // содержать НЕСКОЛЬКО независимых расщеплённых клонов, каждый со своим
  // депозитом и весом — их нельзя сливать в один Fill(), иначе энергии
  // разных клонов складываются в одно фиктивное число (проверено эмпирически:
  // без разделения переполнение >3200 кэВ на 12 из 14 попаданий при
  // максимальной реальной линии 1460,8 кэВ). Fill() — обёртка для ОСТАЛЬНЫХ
  // режимов (resp/sample/muon/beam/TransStep), где в событии всегда РОВНО
  // одна история: даёт ПОБИТОВО то же поведение, что было до правки.
  void FillPrimary(double eprim) { fSumEprim += eprim; }

  void FillDeposit(double edepKeV, double w = 1.0) {
    if (edepKeV <= 0) return;
    fWithSignal += w;
    int b = static_cast<int>(edepKeV / kBinKeV);
    if (b > kBins) b = kBins;
    fHist[b] += w;
    fHist2[b] += w * w;
  }

  void Fill(double edepKeV, double eprim, double w = 1.0) {
    FillPrimary(eprim);
    FillDeposit(edepKeV, w);
  }

  void EndOfRunAction(const G4Run* run) override {
    const long N = run->GetNumberOfEvent();
    if (N == 0) return;
    const double eMean = fSumEprim / N;

    FILE* f = std::fopen(fOut.c_str(), "w");
    if (!f) {
      G4cerr << "!! не открыть " << fOut << G4endl;
      return;
    }
    std::fprintf(f, "# RadiaCode + защита: отклик полости (mode=resp)\n");
    std::fprintf(f, "# src_sha1 = %s\n", RCPB_SRC_SHA1);
    std::fprintf(f, "# git_describe = %s\n", RCPB_GIT_DESCRIBE);
    std::fprintf(f, "# build = %s %s\n", __DATE__, __TIME__);
    std::fprintf(f, "# tag = %s\n", fTag.c_str());
    std::fprintf(f, "# E_prim_keV = %.4f\n", eMean);
    std::fprintf(f, "# N_primaries = %ld\n", N);
    std::fprintf(f, "# N_with_signal = %.6g\n", fWithSignal);
    std::fprintf(f, "# bin_keV = %.3f  (последний канал = переполнение)\n", kBinKeV);
    std::fprintf(f, "E_keV,counts\n");
    for (int i = 0; i <= kBins; ++i)
      if (fHist[i] != 0.0) std::fprintf(f, "%.1f,%.6g\n", (i + 0.5) * kBinKeV, fHist[i]);
    std::fclose(f);

    // Σw² — ОТДЕЛЬНЫМ файлом, а не третьей колонкой: читатели спектров
    // (rcspec.read_spec и др.) разбирают строку как ровно две колонки
    // (`a, b = line.split(",")`), третья их сломала бы.
    const std::string out2 = std::string(fOut.c_str()) + ".sumw2.csv";
    if (FILE* f2 = std::fopen(out2.c_str(), "w")) {
      std::fprintf(f2, "# Σw² по каналам к файлу %s\n", fOut.c_str());
      std::fprintf(f2, "# N_eff канала = (Σw)²/Σw²; отн. погрешность = sqrt(Σw²)/Σw\n");
      std::fprintf(f2, "# N_primaries = %ld\n", N);
      std::fprintf(f2, "# bin_keV = %.3f\n", kBinKeV);
      std::fprintf(f2, "E_keV,sumw2\n");
      for (int i = 0; i <= kBins; ++i)
        if (fHist2[i] != 0.0)
          std::fprintf(f2, "%.1f,%.6g\n", (i + 0.5) * kBinKeV, fHist2[i]);
      std::fclose(f2);
    }

    G4cout << "RESULT resp E_keV= " << eMean << " N= " << N
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
    // Вес первичной частицы — единственное место, где считается "сколько
    // это событие стоит" для взвешенного гистограммирования. В replay-
    // режиме он приходит из CrossRec.w (стадия 1, trans, возможно biased);
    // в resp-режиме гана не биасится, вес всегда умолчательный 1.0
    // (G4ParticleGun::particle_weight по умолчанию) — числа не меняются.
    double ep = 0, w = 1.0;
    if (e->GetNumberOfPrimaryVertex() > 0) {
      auto* prim = e->GetPrimaryVertex(0)->GetPrimary(0);
      ep = prim->GetKineticEnergy() / keV;
      w  = prim->GetWeight();
    }
    fRun->Fill(fEdep / keV, ep, w);
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

// --- replay ПОД doBiasReplay: учёт ПО КЛОНУ, не по событию (задача №13) ----
// G4ImportanceProcess расщепляет трек на НЕСКОЛЬКО независимых клонов внутри
// ОДНОГО G4Event; каждый клон способен ОТДЕЛЬНО долететь до кристалла.
// Клон опознаётся по GetCreatorProcess()=="ImportanceProcess" — рождается
// НОВАЯ линия (свой бакет депозита+веса); обычные вторичные (Комптон-
// электрон и т.п.) наследуют линию родителя. Без этого разделения энергии
// разных клонов складывались бы в одно фиктивное число (см. план, задача
// №13 — эмпирически пойман переполняющий бин >3200 кэВ при максимуме 1461).
class LineageTracking : public G4UserTrackingAction {
  std::map<G4int, G4int>* fLineage;
public:
  explicit LineageTracking(std::map<G4int, G4int>* m) : fLineage(m) {}
  void PreUserTrackingAction(const G4Track* trk) override {
    const G4int id = trk->GetTrackID();
    const G4int parentId = trk->GetParentID();
    bool isClone = false;
    if (const G4VProcess* cp = trk->GetCreatorProcess())
      isClone = (cp->GetProcessName() == "ImportanceProcess");
    if (parentId == 0 || isClone) {
      (*fLineage)[id] = id;                              // новый корень линии
    } else {
      auto it = fLineage->find(parentId);
      (*fLineage)[id] = (it != fLineage->end()) ? it->second : id;
    }
  }
};

class ReplayEventAct : public G4UserEventAction {
  RunAct* fRun;
  std::map<G4int, G4int>* fLineage;
  std::map<G4int, double> fEdepByLineage;    // кэВ, накопленный депозит на линию
  std::map<G4int, double> fWeightByLineage;  // постоянен внутри линии
public:
  ReplayEventAct(RunAct* r, std::map<G4int, G4int>* lin) : fRun(r), fLineage(lin) {}
  void BeginOfEventAction(const G4Event*) override {
    fEdepByLineage.clear();
    fWeightByLineage.clear();
    fLineage->clear();
  }
  void AddDeposit(G4int trackId, double edepKeV, double w) {
    G4int lineage = trackId;
    auto it = fLineage->find(trackId);
    if (it != fLineage->end()) lineage = it->second;
    fEdepByLineage[lineage] += edepKeV;
    fWeightByLineage[lineage] = w;
  }
  void EndOfEventAction(const G4Event* e) override {
    double ep = 0;
    if (e->GetNumberOfPrimaryVertex() > 0)
      ep = e->GetPrimaryVertex(0)->GetPrimary(0)->GetKineticEnergy() / keV;
    fRun->FillPrimary(ep);   // ОДИН раз на событие, независимо от числа линий
    for (const auto& kv : fEdepByLineage)
      fRun->FillDeposit(kv.second, fWeightByLineage[kv.first]);
  }
};

class ReplayStepping : public G4UserSteppingAction {
  ReplayEventAct* fEvt;
  const G4LogicalVolume* fCry;
public:
  ReplayStepping(ReplayEventAct* ev, const G4LogicalVolume* c) : fEvt(ev), fCry(c) {}
  void UserSteppingAction(const G4Step* s) override {
    auto* h = s->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    if (!h || h->GetLogicalVolume() != fCry) return;
    const double edepKeV = s->GetTotalEnergyDeposit() / keV;
    if (edepKeV <= 0) return;
    const double w = s->GetPreStepPoint()->GetWeight();
    fEvt->AddDeposit(s->GetTrack()->GetTrackID(), edepKeV, w);
  }
};

// --- узкопучковая проверка (iii): mucalc против реального транспорта --------
// Карандашный пучок вдоль оси Z сквозь ВЕРХНЮЮ крышку защиты (толщина ровно
// pb+cd+cu, как в mucalc). Считаются только ПЕРВИЧНЫЕ кванты (trackID==1),
// дошедшие до границы полости БЕЗ единого взаимодействия — направление
// осталось точно исходным. Это узкое ("good geometry") определение и есть то,
// что считает mu в mucalc.cc (compt+phot+conv, БЕЗ Rayleigh): даже
// когерентное рассеяние меняет направление и увело бы квант из проверки,
// поэтому проверять надо направление, не только энергию (энергию Rayleigh не
// меняет вовсе — проверка по одной энергии была бы неполной).
class BeamStep : public G4UserSteppingAction {
  RunAct* fRun;
  long fSurvived = 0;
public:
  explicit BeamStep(RunAct* r) : fRun(r) {}
  long Survived() const { return fSurvived; }
  void UserSteppingAction(const G4Step* s) override {
    auto* trk = s->GetTrack();
    if (trk->GetTrackID() != 1) return;                       // только первичный
    if (trk->GetDefinition()->GetPDGEncoding() != 22) return;
    if (s->GetPostStepPoint()->GetStepStatus() != fGeomBoundary) return;
    auto* postV = s->GetPostStepPoint()->GetPhysicalVolume();
    // "cavity", не "world" — с 12.08.2026 полость отдельный физобъём (см.
    // PbShield.hh); beam-режим тоже строится с fWithDevice=false, поэтому
    // "cavity" там тоже есть.
    if (!postV || postV->GetName() != "cavity") return;
    const G4ThreeVector dir = s->GetPostStepPoint()->GetMomentumDirection();
    // Пучок идёт СВЕРХУ ВНИЗ (-z) сквозь верхнюю крышку в полость.
    if (dir.z() > -0.999999) return;   // отклонился (Rayleigh или иное) -> не считаем
    ++fSurvived;
    fRun->Fill(s->GetPostStepPoint()->GetKineticEnergy() / keV,
              s->GetPostStepPoint()->GetKineticEnergy() / keV);
    trk->SetTrackStatus(fStopAndKill);
  }
};

// --- стадия 1: запись пересечений границы полости ---------------------------
// НЕБИАСИРОВАННАЯ версия. Записывает каждый гамма-квант, переходящий из
// внутреннего слоя защиты (depth 0) в полость (мир), и сразу его убивает —
// см. разбор стыка стадий в PbShield.hh/плане: обратного хода нет по
// построению, полость в стадии 1 всегда пуста.
//
// ПОЧЕМУ ПИШЕМ УГОЛ, А НЕ ТОЛЬКО ЭНЕРГИЮ. Свёртка «отклик полости от resp(E)»
// законна ТОЛЬКО если входящие пересечения изотропно-косинусные (условие
// проговорено в плане явно, не постулируется). Значит первый прогон обязан
// это ПРОВЕРИТЬ, а не предположить: cosTheta — угол между направлением влёта
// и внутренней нормалью к поверхности полости в точке пересечения (1 = по
// нормали, 0 = по касательной). Нормаль берётся геометрически: у боковой
// стенки — радиальная, у дна/крышки — вдоль оси Z; какая из них,
// определяется по фактическому положению точки пересечения, не по имени тела
// (устойчивее к порядку слоёв).
class TransStep : public G4UserSteppingAction {
  RunAct* fRun;
  const std::vector<G4VPhysicalVolume*>* fLayerPV;
  const std::vector<int>* fLayerDepth;
  double fHxCav, fHyCav, fHzCav, fZCav;
  bool fWithLid;
  FILE* fOut;
  long fCrossings = 0;      // сырое число пересечений (после расщепления) — диагностика
  double fSumW = 0.0;       // Σвес — несмещённая оценка; ЭТО число, не fCrossings, делить на nprim
  // Отдельный учёт канала через ОТКРЫТЫЙ ВЕРХ: он не проходит сквозь свинец
  // вовсе, поэтому его доля в Σw — прямая мера того, сколько фона даёт
  // отверстие, и главный ориентир при обсуждении крышки с оператором.
  long fCrossTop = 0;
  double fSumWTop = 0.0;
public:
  TransStep(RunAct* r, const RCShieldDetector* det, FILE* out)
      : fRun(r), fLayerPV(&det->fLayerPV), fLayerDepth(&det->fLayerDepth),
        fHxCav(det->fSh.hxCav), fHyCav(det->fSh.hyCav),
        fHzCav(det->fSh.hzCav), fZCav(det->fSh.zCav),
        fWithLid(det->fWithLid), fOut(out) {}
  long Crossings() const { return fCrossings; }
  double SumWeight() const { return fSumW; }
  long CrossingsTop() const { return fCrossTop; }
  double SumWeightTop() const { return fSumWTop; }

  void UserSteppingAction(const G4Step* s) override {
    if (s->GetTrack()->GetDefinition()->GetPDGEncoding() != 22) return;  // gamma

    // ОТКРЫТЫЙ ВЕРХ: там частица не ПЕРЕСЕКАЕТ границу полости, а РОЖДАЕТСЯ
    // на ней. Верхняя грань наружного габарита при nolid совпадает с верхней
    // гранью полости (свинца над полостью нет, PlannedOuterZHi == zCav+hzCav),
    // поэтому первичные, разыгранные BoxSurfaceGun на этой грани, стартуют
    // уже ВНУТРИ cavity — шага «мир -> полость» у них не бывает, и правка,
    // разрешавшая приход из "world" (первая попытка 15.08.2026), ничего не
    // дала: проба вернула ровно те же 21049 пересечений и ни одного сверху.
    //
    // Объём рождения сам разделяет каналы и делать это вручную не нужно:
    // верхняя грань габарита 250×250 накрывает и полость 150×150 (там объём
    // cavity), и рамку над торцами свинцовых стенок (там объём слоя) — в
    // свинце трек стартует, и эта ветка не срабатывает.
    if (!fWithLid) {
      auto* trk = s->GetTrack();
      auto* v0 = s->GetPreStepPoint()->GetPhysicalVolume();
      if (trk->GetParentID() == 0 && trk->GetCurrentStepNumber() == 1 && v0 &&
          v0->GetName() == "cavity") {
        const G4ThreeVector p0 = s->GetPreStepPoint()->GetPosition() / mm;
        const G4ThreeVector d0 = s->GetPreStepPoint()->GetMomentumDirection();
        const double e0 = s->GetPreStepPoint()->GetKineticEnergy() / keV;
        const double w0 = s->GetPreStepPoint()->GetWeight();
        const double cos0 = -d0.z();   // внешняя нормаль верхней грани = +z
        if (fOut)
          std::fprintf(fOut, "%.3f,%.4f,+z,%.2f,%.2f,%.2f,%.6f,%.6f,%.6f,%.6e\n",
                       e0, cos0, p0.x(), p0.y(), p0.z(), d0.x(), d0.y(), d0.z(),
                       w0);
        fRun->Fill(e0, e0, w0);
        ++fCrossings;
        fSumW += w0;
        ++fCrossTop;
        fSumWTop += w0;
        trk->SetTrackStatus(fStopAndKill);
        return;
      }
    }

    if (s->GetPostStepPoint()->GetStepStatus() != fGeomBoundary) return;
    auto* preV = s->GetPreStepPoint()->GetPhysicalVolume();
    auto* postV = s->GetPostStepPoint()->GetPhysicalVolume();
    // "cavity", не "world" — см. PbShield.hh, правка 12.08.2026.
    if (!preV || !postV || postV->GetName() != "cavity") return;

    bool depth0 = false;
    for (size_t i = 0; i < fLayerPV->size(); ++i)
      if ((*fLayerPV)[i] == preV && (*fLayerDepth)[i] == 0) { depth0 = true; break; }
    // ОТКРЫТЫЙ ВЕРХ (nolid): над полостью свинца нет, и квант приходит туда
    // прямо из мира — слоя depth=0 на этом пути не существует. До правки
    // 15.08.2026 такие пересечения не записывались ВООБЩЕ: условие требовало
    // depth=0, и весь канал через открытое отверстие 150×150 молча выпадал
    // из стадии 1, а значит и из replay. Именно этот канал — причина, по
    // которой открытый домик вообще моделируется, так что его потеря делала
    // расчёт фона заведомо заниженным, и заметно это было только по
    // отсутствию грани +z в CSV (проба 15.08.2026: 21049 пересечений, ни
    // одного сверху).
    // Мир граничит с полостью ТОЛЬКО сверху и только при снятой крышке:
    // с боков и снизу между ними всегда есть слой защиты (pb>0 в любой
    // точке сетки). Поэтому проверки имени объёма достаточно, отдельная
    // проверка «точка лежит на верхней грани» была бы избыточной.
    const bool fromWorldTop = (!fWithLid && preV->GetName() == "world");
    if (!depth0 && !fromWorldTop) return;

    const G4ThreeVector p = s->GetPostStepPoint()->GetPosition() / mm;
    const G4ThreeVector dir = s->GetPostStepPoint()->GetMomentumDirection();
    // Грань КОРОБА выбирается по наименьшему зазору до неё — тот же принцип,
    // что был у цилиндра (сравнение расстояний до боковой и до торцов), но на
    // пяти-шести плоскостях. Имя тела по-прежнему не используется: устойчиво
    // к порядку слоёв. Верхняя грань участвует, только если крышка построена;
    // при открытом верхе кванта, влетающего сверху, снаружи не бывает — там
    // нет слоя depth=0, и до этой строки шаг просто не доходит.
    const double dx = fHxCav - std::abs(p.x());
    const double dy = fHyCav - std::abs(p.y());
    const double dzB = (p.z() - (fZCav - fHzCav));   // до дна
    const double dzT = ((fZCav + fHzCav) - p.z());   // до крышки
    double best = std::abs(dx);
    G4ThreeVector outNorm((p.x() > 0) ? 1 : -1, 0, 0);
    const char* face = (p.x() > 0) ? "+x" : "-x";
    if (std::abs(dy) < best) {
      best = std::abs(dy);
      outNorm = G4ThreeVector(0, (p.y() > 0) ? 1 : -1, 0);
      face = (p.y() > 0) ? "+y" : "-y";
    }
    if (std::abs(dzB) < best) {
      best = std::abs(dzB);
      outNorm = G4ThreeVector(0, 0, -1);
      face = "-z";
    }
    // Верхняя грань ПОЛОСТИ существует всегда — она граница объёма cavity, а
    // не свинцовой крышки. Прежнее условие «только при fWithLid» приводило к
    // тому, что кванту, влетевшему сверху, приписывалась ближайшая БОКОВАЯ
    // нормаль, и cosTheta выходил бессмысленным (вплоть до отрицательного).
    if (std::abs(dzT) < best) {
      best = std::abs(dzT);
      outNorm = G4ThreeVector(0, 0, 1);
      face = "+z";
    }
    const double cosTheta = -(dir * outNorm);   // >0 при корректном влёте внутрь
    const double eKeV = s->GetPostStepPoint()->GetKineticEnergy() / keV;
    // Вес — по конвенции B01 (см. план): PreStepPoint()->GetWeight(), не
    // GetTrack()->GetWeight(). Без биасинга всегда 1.0 (побитово прежнее
    // поведение). С биасингом — вес расщеплённой/срулеченной истории;
    // именно fSumW (не fCrossings) даёт несмещённую оценку пропускания.
    const double w = s->GetPreStepPoint()->GetWeight();

    // dx,dy,dz — ПОЛНОЕ направление, не только угол к нормали: cosTheta один
    // не восстанавливает 3-вектор (теряется азимут в касательной плоскости).
    // Нужно для честного replay — стадия 2 обязана стартовать частицу ровно
    // с тем направлением, с которым она реально пересекла границу, и С ТЕМ
    // ЖЕ ВЕСОМ (10-я колонка) — иначе стадия 2 тихо потеряет несмещённость.
    if (fOut)
      std::fprintf(fOut, "%.3f,%.4f,%s,%.2f,%.2f,%.2f,%.6f,%.6f,%.6f,%.6e\n", eKeV,
                   cosTheta, face, p.x(), p.y(), p.z(), dir.x(), dir.y(), dir.z(), w);
    fRun->Fill(eKeV, eKeV, w);   // переиспользуем гистограмму RunAct как взвешенный счётчик
    ++fCrossings;
    fSumW += w;
    if (face[0] == '+' && face[1] == 'z') { ++fCrossTop; fSumWTop += w; }

    s->GetTrack()->SetTrackStatus(fStopAndKill);
  }
};

// --- космический мюон сквозь ПОЛНУЮ сборку (прибор+сосуд+защита) ------------
// В ОТЛИЧИЕ от гамма-компонент (trans/replay), мюону НЕ нужен двухстадийный
// биасинг: 100 мм Pb не "экспоненциально гасят" ГэВ-мюон почти без остатка
// (как гамма 662 кэВ), а лишь отбирают у него несколько сотен МэВ по dE/dx —
// прямой небиасированный прогон через ВСЮ геометрию работает за разумное
// время при любой толщине (та же логика, что уже в cosmicmu.cc для голого
// прибора, тут — та же схема, но с защитой). Формулы (Гайссер, cos²θ)
// СВЕРЕНЫ построчно с cosmicmu.cc — источник истины один, здесь копия для
// самодостаточности бинарника (конвенция репо: shieldrun.cc/wallfield.cc/
// mucalc.cc не делят общий заголовок ради экономии дублирования).
double GaisserVertical(double E_GeV) {
  const double t1 = 1.0 / (1.0 + 1.1 * E_GeV / 115.0);
  const double t2 = 0.054 / (1.0 + 1.1 * E_GeV / 850.0);
  return 0.14 * std::pow(E_GeV, -2.7) * (t1 + t2);
}

const double MU_E_LO_GEV = 0.3, MU_E_HI_GEV = 100.0;
const int MU_N_EGRID = 400;
std::vector<double> gMuEGrid, gMuCumE;

void BuildMuEnergyTable() {
  gMuEGrid.resize(MU_N_EGRID);
  gMuCumE.resize(MU_N_EGRID);
  double s = 0;
  const double logLo = std::log(MU_E_LO_GEV), logHi = std::log(MU_E_HI_GEV);
  for (int i = 0; i < MU_N_EGRID; ++i)
    gMuEGrid[i] = std::exp(logLo + (logHi - logLo) * i / (MU_N_EGRID - 1));
  for (int i = 0; i < MU_N_EGRID; ++i) {
    s += GaisserVertical(gMuEGrid[i]) * gMuEGrid[i];   // якобиан d(lnE)=dE/E
    gMuCumE[i] = s;
  }
  for (int i = 0; i < MU_N_EGRID; ++i) gMuCumE[i] /= s;
}

double SampleMuEnergyGeV() {
  const double u = G4UniformRand();
  int lo = 0, hi = MU_N_EGRID - 1;
  while (lo < hi) {
    const int mid = (lo + hi) / 2;
    if (gMuCumE[mid] < u) lo = mid + 1; else hi = mid;
  }
  return gMuEGrid[lo];
}

// Диск-источник: НАД собранной защитой, радиус растёт с fOuterR (иначе на
// pb=100 мм — rOut=150 мм — диск R=70 мм из cosmicmu.cc даже вертикальные
// мюоны накрывал бы не всю крышку). Запас DISK_MARGIN_MM — тот же порядок,
// что у cosmicmu.cc (R_DISK=70 при габарите прибора ~40 мм, запас ~30 мм) и
// запас по высоте DISK_CLEAR_MM — то же соотношение (Z_DISK=130 при хвосте
// 111 мм, запас ~19 мм), округлено чуть щедрее ради косых треков.
const double DISK_MARGIN_MM = 40.0;
const double DISK_CLEAR_MM = 25.0;

// Фактические размеры диска, выбранные в main() — пишутся в шапку CSV, чтобы
// downstream-анализ брал радиус ОТТУДА, а не пересчитывал по своей копии
// формулы (run_shield_grid.muon_r_disk() именно так и делал — расхождение
// копии с исходником осталось бы незамеченным).
double gMuRDisk = 0.0, gMuZDisk = 0.0;

class MuGun : public G4VUserPrimaryGeneratorAction {
  G4ParticleGun fGun{1};
  double fRDisk, fZDisk;
public:
  MuGun(double rDisk, double zDisk) : fRDisk(rDisk), fZDisk(zDisk) {
    fGun.SetParticleDefinition(
        G4ParticleTable::GetParticleTable()->FindParticle("mu-"));
  }
  void GeneratePrimaries(G4Event* e) override {
    const double r = fRDisk * std::sqrt(G4UniformRand());
    const double ph = twopi * G4UniformRand();
    fGun.SetParticlePosition(
        G4ThreeVector(r * std::cos(ph), r * std::sin(ph), fZDisk * mm));

    // Розыгрыш угла — для ПОТОКА ЧЕРЕЗ ГОРИЗОНТАЛЬНУЮ ПЛОЩАДКУ, а не для
    // интенсивности. Интенсивность мюонов I(θ) ~ cos²θ задана на площадку,
    // ПЕРПЕНДИКУЛЯРНУЮ треку; частиц же, пересекающих горизонтальный диск,
    // приходится dN ~ I(θ)·cosθ·dΩ ~ cos³θ·sinθ·dθ, то есть p(cosθ) ~ cos³θ
    // и розыгрыш cosθ = u^(1/4). Здесь до 15.08.2026 стояло cbrt(u), что даёт
    // p(cosθ) ~ cos²θ — распределение самой интенсивности, без множителя cosθ,
    // и мюоны стартовали заметно наклоннее, чем идут на самом деле
    // (⟨cosθ⟩ 3/4 вместо 4/5, а на хвосте больших углов расхождение сильнее —
    // именно там длиннее путь сквозь свинец и корпус).
    // ⚠ ЭТО ОБЕСЦЕНИВАЕТ A_MU=315,2 из run_shield_grid.py: та величина
    // подбиралась под измеренный фон ПРИ СТАРОМ розыгрыше и часть ошибки
    // впитала в себя. Перекалибровать (или, что честнее, взять априорный
    // поток PDG ~0,0167 мюон/(см²·с) на горизонтальную поверхность) — #SHIELD-9.
    const double cosT = std::pow(G4UniformRand(), 0.25);   // p(cosT) ~ cosT^3
    const double sinT = std::sqrt(std::max(0.0, 1.0 - cosT * cosT));
    const double phd = twopi * G4UniformRand();
    fGun.SetParticleMomentumDirection(
        G4ThreeVector(sinT * std::cos(phd), sinT * std::sin(phd), -cosT));

    fGun.SetParticleEnergy(SampleMuEnergyGeV() * GeV);
    fGun.GeneratePrimaryVertex(e);
  }
};

// --- изотропное поле на поверхности КОРОБА ----------------------------------
// Заменяет /gps/pos/type Surface + shape Cylinder, которым описывалась
// цилиндрическая защита. У GPS нет формы «поверхность прямоугольного короба»,
// поэтому розыгрыш точки и направления делается здесь, а ЭНЕРГИЯ по-прежнему
// берётся штатным GPS: вызываем его GeneratePrimaryVertex() и переписываем у
// готовой вершины позицию и направление. Так весь механизм энергетических
// спектров (Mono и Arb-гистограмма из specmac=field_spectrum_*.mac) продолжает
// работать без единой правки, а геометрия розыгрыша становится верной.
//
// ФИЗИКА, НА КОТОРОЙ ЭТО ДЕРЖИТСЯ. Изотропное поле с флюенсом Ф, падающее на
// выпуклое тело с поверхностью S, даёт N = Ф·S/4 входящих частиц, причём
// распределение по углу входа — косинусное относительно внутренней нормали
// (соотношение Коши). Обратно: разыграв N частиц равномерно по площади S с
// косинусным законом внутрь, получаем внутри в точности изотропное поле
// Ф = 4N/S. Тождество то же, что в resp/run_bg.py/wallfield.cc, и оно НЕ
// зависит от формы тела — только от выпуклости. Поэтому для короба меняется
// лишь способ розыгрыша и значение S, а не сама схема нормировки.
//
// ⚠ ВЕРХНЯЯ ГРАНЬ РАЗЫГРЫВАЕТСЯ ВСЕГДА, в том числе при открытом верхе
// (nolid). Это не описка. Выпуклая оболочка защиты остаётся замкнутым
// параллелепипедом независимо от того, построена ли свинцовая крышка; поле
// помещения входит через верхнюю плоскость габарита в обоих случаях. Разница
// в том, ЧТО квант встречает за этой плоскостью: при закрытом верхе — свинец,
// при открытом — воздух полости и прямой путь к кристаллу. Убрать верхнюю
// грань из розыгрыша значило бы выбросить именно тот канал, ради которого
// открытый верх и моделируется. Соответственно S в знаменателе Ф = 4N/S —
// площадь ПОЛНОГО замкнутого габарита, тоже независимо от nolid.
class BoxSurfaceGun : public G4VUserPrimaryGeneratorAction {
  G4GeneralParticleSource fGPS;
  double fHx, fHy, fHz, fZc;
public:
  BoxSurfaceGun(double hx, double hy, double hz, double zc)
      : fHx(hx), fHy(hy), fHz(hz), fZc(zc) {}

  // Полная площадь замкнутого габарита, мм² — то самое S для Ф = 4N/S.
  // Отдаётся наружу, чтобы downstream брал число ОТСЮДА, а не пересчитывал
  // по своей копии формулы (ровно та ошибка, что уже случалась с радиусом
  // мюонного диска, см. gMuRDisk).
  static double SurfaceArea(double hx, double hy, double hz) {
    return 8.0 * (hy * hz + hx * hz + hx * hy);
  }

  // Косинусный закон вокруг направления n: cos(theta) = sqrt(u) — это и есть
  // распределение входящих частиц изотропного поля через плоскую площадку.
  static G4ThreeVector CosineAbout(const G4ThreeVector& n) {
    const double cosT = std::sqrt(G4UniformRand());
    const double sinT = std::sqrt(std::max(0.0, 1.0 - cosT * cosT));
    const double ph = twopi * G4UniformRand();
    // Касательный базис: берём орт, заведомо не параллельный n.
    const G4ThreeVector a =
        (std::abs(n.x()) < 0.9) ? G4ThreeVector(1, 0, 0) : G4ThreeVector(0, 1, 0);
    const G4ThreeVector u = (a - n * a.dot(n)).unit();
    const G4ThreeVector v = n.cross(u);
    return (u * (sinT * std::cos(ph)) + v * (sinT * std::sin(ph)) + n * cosT).unit();
  }

  void GeneratePrimaries(G4Event* e) override {
    fGPS.GeneratePrimaryVertex(e);          // частица и энергия — штатным GPS
    auto* vtx = e->GetPrimaryVertex(0);
    if (!vtx) return;

    // Грань выбирается пропорционально своей площади — только тогда точки
    // равномерны по ВСЕЙ поверхности, а не по числу граней.
    const double sx = 4.0 * fHy * fHz;      // каждая из граней +-X
    const double sy = 4.0 * fHx * fHz;      // каждая из граней +-Y
    const double sz = 4.0 * fHx * fHy;      // каждая из граней +-Z
    double t = G4UniformRand() * (2.0 * (sx + sy + sz));
    G4ThreeVector pos, inward;
    if (t < 2.0 * sx) {
      const double s = (t < sx) ? +1.0 : -1.0;
      pos = G4ThreeVector(s * fHx,
                          fHy * (2.0 * G4UniformRand() - 1.0),
                          fZc + fHz * (2.0 * G4UniformRand() - 1.0));
      inward = G4ThreeVector(-s, 0, 0);
    } else if ((t -= 2.0 * sx) < 2.0 * sy) {
      const double s = (t < sy) ? +1.0 : -1.0;
      pos = G4ThreeVector(fHx * (2.0 * G4UniformRand() - 1.0),
                          s * fHy,
                          fZc + fHz * (2.0 * G4UniformRand() - 1.0));
      inward = G4ThreeVector(0, -s, 0);
    } else {
      t -= 2.0 * sy;
      const double s = (t < sz) ? +1.0 : -1.0;
      pos = G4ThreeVector(fHx * (2.0 * G4UniformRand() - 1.0),
                          fHy * (2.0 * G4UniformRand() - 1.0),
                          fZc + s * fHz);
      inward = G4ThreeVector(0, 0, -s);
    }

    vtx->SetPosition(pos.x() * mm, pos.y() * mm, pos.z() * mm);
    if (vtx->GetPrimary(0))
      vtx->GetPrimary(0)->SetMomentumDirection(CosineAbout(inward));
  }
};

}  // namespace

int main(int argc, char** argv) {
  const std::string mode = (argc > 1) ? argv[1] : "geom";

  ShieldGeom sh;
  std::string vessel = "m200", matrix = "organic", out = "resp.csv", in;
  double rho = 0.50, nshell = 8;
  double eKeV = 661.7, nprim = 1000000;
  double impStep = 2.0;   // геометрический шаг важности между слоями (bias), как в B01
  std::string specmac;    // непрерывный спектр поля (GPS Arb) вместо моно — только mode=trans
  std::string nuc = "K40";   // K40 | Cs137 — только mode=sample
  double cryStep = 4.0;      // шаг важности case->caseAir->reflector->crystal, mode=replay bias
  // seed — зерно генератора. 0 = не трогать (прежнее поведение побитово:
  // Geant4 стартует с фиксированного зерна, и два одинаковых вызова дают
  // ОДИН И ТОТ ЖЕ ответ). Введено 15.08.2026: коробчатая полость 150×150×385
  // на два порядка просторнее прежней цилиндрической, кристалл виден из неё
  // под телесным углом ~1e-5, и одиночный прогон на 10^5 первичных даёт
  // 0,25 взвешенных попаданий в кристалл — статистики нет вовсе. Единственный
  // выход — гнать десятки процессов разом, а без своего зерна у каждого они
  // все повторяют одну и ту же историю и складывать их бессмысленно.
  double seed = 0;
  // repeat — сколько раз прогнать входной файл replay. Записи стадии 1 стоят
  // дорого (перенос сквозь 50 мм свинца), а разыгрывать по ним перенос ВНУТРИ
  // полости можно многократно с новым зерном: это независимые выборки того
  // самого редкого события, ради которого всё считается. Нормировать при этом
  // надо на N_stage1 × repeat — множитель пишется в шапку выходного файла.
  double repeat = 1;
  // rdisk — радиус диск-источника мюонов, мм. 0 = авто (rOut+DISK_MARGIN_MM,
  // прежнее поведение, побитово). Введён 13.08.2026 для проверки НАСЫЩЕНИЯ:
  // при cos²θ, θ<70° мюон приходит в кристалл со старта на расстоянии до
  // ~2,75·Z_disk от оси, тогда как авто-радиус даёт всего rOut+40 мм — то
  // есть диск ОБРЕЗАЕТ наклонные треки, и тем сильнее, чем он меньше.
  // Физический ответ (j·πR²·eff(R)) обязан выходить на полку с ростом R;
  // пока не вышел — поправка (R/70)² в run_shield_grid.mu_cps() неполна.
  double rDiskOverride = 0.0;
  double seatFloor = 1.0;   // 1 = сажать полость на дно (см. seatOnFloor)
  // P-008 (17.08): ориентация прибора в полости. Умолчание 0 — вертикально,
  // как было. Оператор сообщил, что реально прибор ЛЕЖАЛ: «горизонтально на
  // картонной коробочке, 2 см от дна домика, кристаллом вниз».
  double horizFlag = 0.0;   // 1 = прибор лежит (длинная ось поперёк)
  double liftMm = 0.0;      // ВНУТРЕННЯЯ координата: сдвиг объёма "case" по z
  // P-013 (17.08): оператор задаёт не координату кристалла, а ФИЗИЧЕСКИЙ ЗАЗОР
  // между дном полости и низом корпуса («он не прижат. зазор 25 мм, пустая
  // картонная коробочка»). Это РАЗНЫЕ числа: при горизонтальной посадке центр
  // кристалла отстоит от нижней грани корпуса на crystalToFace = 8,20 мм, а
  // пол полости при seatfloor лежит на z = -crystalZ0 = -12,00 мм, поэтому
  // lift=25 даёт зазор 28,8 мм, а не 25. Ключ gap= снимает пересчёт с
  // человека: зазор задаётся как измерен, lift вычисляется здесь.
  // gap < 0 означает «не задан» — ноль это законное значение (прибор прижат).
  double gapMm = -1.0;
  // mode=pbself: имя ячейки свинца, в которой разыгрывается Pb-210. Свинец —
  // пять тел при открытом верхе, а /gps/pos/confine принимает одно имя, поэтому
  // ячейки считаются по одной и складываются драйвером с весами по массе.
  std::string pbCell;
  bool withVessel = true, withShield = true, bias = false, withLid = true;
  bool noShieldAsked = false;   // см. разбор у mode == "replay" ниже

  for (int i = 2; i < argc; ++i) {
    const char* a = argv[i];
    if (KeyVal(a, "pb", &sh.pb)) continue;
    if (KeyVal(a, "cu", &sh.cu)) continue;
    if (KeyVal(a, "cd", &sh.cd)) continue;
    if (KeyVal(a, "hxcav", &sh.hxCav)) continue;
    if (KeyVal(a, "hycav", &sh.hyCav)) continue;
    if (KeyVal(a, "hzcav", &sh.hzCav)) continue;
    // zcav= задаётся вручную ТОЛЬКО вместе с seatfloor=0, иначе центр полости
    // всё равно пересчитается посадкой на дно и заданное значение пропадёт
    // молча. Явный ключ вместо тихой перезаписи — см. seatOnFloor в PbShield.hh.
    if (KeyVal(a, "zcav", &sh.zCav)) continue;
    if (KeyVal(a, "seatfloor", &seatFloor)) continue;
    if (KeyVal(a, "horiz", &horizFlag)) continue;   // P-008
    if (KeyVal(a, "lift", &liftMm)) continue;       // P-008, внутренняя координата
    if (KeyVal(a, "gap", &gapMm)) continue;         // P-013, физический зазор
    if (KeyVal(a, "nshell", &nshell)) continue;
    if (KeyVal(a, "rho", &rho)) continue;
    if (KeyVal(a, "e", &eKeV)) continue;
    if (KeyVal(a, "nprim", &nprim)) continue;
    if (KeyVal(a, "impstep", &impStep)) continue;
    if (KeyVal(a, "crystep", &cryStep)) continue;
    if (KeyVal(a, "seed", &seed)) continue;
    if (KeyVal(a, "repeat", &repeat)) continue;
    if (KeyVal(a, "rdisk", &rDiskOverride)) continue;
    if (KeyStr(a, "vessel", &vessel)) continue;
    if (KeyStr(a, "matrix", &matrix)) continue;
    if (KeyStr(a, "out", &out)) continue;
    if (KeyStr(a, "in", &in)) continue;
    if (KeyStr(a, "specmac", &specmac)) continue;
    if (KeyStr(a, "cell", &pbCell)) continue;   // mode=pbself: ячейка свинца
    if (KeyStr(a, "nuc", &nuc)) continue;
    if (std::strcmp(a, "novessel") == 0) { withVessel = false; continue; }
    // noShieldAsked отделяет «оператор ЯВНО попросил без свинца» от «свинца
    // нет по умолчанию режима» — без этого различия нельзя дать replay свинец
    // по умолчанию, сохранив ключ для прежнего поведения.
    if (std::strcmp(a, "noshield") == 0) {
      withShield = false;
      noShieldAsked = true;
      continue;
    }
    // bias — ТОЛЬКО mode=trans (см. ниже doBias); в остальных режимах флаг
    // молча игнорируется, а не ошибка — так проще гонять один и тот же
    // командный шаблон по режимам.
    if (std::strcmp(a, "bias") == 0) { bias = true; continue; }
    // nolid — не строить верхнюю крышку защиты (см. fWithLid в PbShield.hh).
    // Нужен для сверки с реальным домиком оператора, который собран открытым.
    if (std::strcmp(a, "nolid") == 0) { withLid = false; continue; }
    std::cerr << "неизвестный аргумент: " << a << "\n";
    return 2;
  }
  sh.nShellPb = static_cast<int>(nshell);
  // resp: отклик считается на границе полости БЕЗ свинца — источник должен
  // сидеть ровно там, где начинается защита в стадии trans. Если noshield не
  // просили явно, для resp его выключаем сами: иначе источник на границе
  // полости стартовал бы ВНУТРИ первого слоя свинца.
  if (mode == "resp") withShield = false;
  // trans: стадия 1 конвейера снижения дисперсии — полость ПУСТА по
  // построению (прибор и сосуд туда не помещаются намеренно, см. PbShield.hh),
  // а защита обязана быть.
  bool withDevice = true;
  if (mode == "trans") { withDevice = false; withVessel = false; withShield = true; }
  // replay: полная геометрия — прибор, сосуд И СВИНЕЦ.
  //
  // ИСПРАВЛЕНО 15.08.2026. Здесь стояло withShield = false с обоснованием
  // «роль свинца уже учтена тем, что до этой точки дошли только реально
  // выжившие в trans пересечения». Это верно ровно наполовину: учтён ПРЯМОЙ
  // проход сквозь толщу, но не то, что квант делает, уже оказавшись в полости.
  // А делает он вот что — летит через 385 мм пустоты, чаще всего мимо
  // кристалла, и попадает в ПРОТИВОПОЛОЖНУЮ стенку, откуда идёт обратное
  // рассеяние (альбедо свинца в области сотен кэВ немалое) и K-флуоресценция
  // 72,8-87,3 кэВ. Без свинца в стадии 2 он вместо этого улетает в мир, и обе
  // добавки теряются целиком.
  //
  // Симптом в данных (сверка с измерением 15.08.2026): выше 400 кэВ модель
  // ложится на опыт (1000-2000 кэВ — 0,985), а ниже недобирает вдвое, причём
  // ровно в области 72-88 кэВ у модели ПРОВАЛ (K-край поглощения) там, где у
  // измерения максимум. Сама флуоресценция при этом работает — в записях
  // стадии 1 линии Pb стоят на своих местах (73,5 / 75,5 / 84,5 / 85,5 / 87,5
  // кэВ), просто возбуждать их внутри полости было нечему.
  //
  // Двойного счёта прохождения не возникает: частица стартует НА внутренней
  // грани полости и направлена внутрь, сквозь записанную толщу второй раз она
  // не идёт. Прежнее поведение остаётся доступным явным ключом noshield.
  if (mode == "replay" && !noShieldAsked) withShield = true;
  // beam: та же геометрия, что trans (пустая защита, без прибора) — узкий
  // пучок сквозь неё вдоль оси.
  if (mode == "beam") { withDevice = false; withVessel = false; withShield = true; }

  std::vector<CrossRec> recs;
  if (mode == "replay") {
    if (in.empty()) {
      std::cerr << "replay: нужен аргумент in=<файл пересечений trans>\n";
      return 5;
    }
    recs = ReadCrossings(in);
    if (recs.empty()) {
      std::cerr << "replay: 0 записей прочитано из " << in << "\n";
      return 6;
    }
  }

  // Зерно ставится ДО создания RunManager — иначе часть инициализации уже
  // успевает потянуть числа из генератора, и «независимые» процессы стартуют
  // с общего куска последовательности.
  if (seed != 0) G4Random::setTheSeed(static_cast<long>(seed));

  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  auto* det = new RCShieldDetector(withVessel, withShield, withDevice);
  sh.seatOnFloor = (seatFloor != 0.0);
  det->fSh = sh;
  det->fVes = VesselGeom::Preset(vessel);
  det->fVes.sampleMatrix = matrix;
  det->fVes.sampleDensity = rho;
  det->fWithLid = withLid;
  // P-008: ориентация прибора. По умолчанию оба нуля — поведение побитово
  // прежнее (вертикально, без подъёма).
  det->fHorizontal = (horizFlag != 0.0);
  det->fLiftMm = liftMm;
  // P-013: gap= (физический зазор) -> lift= (внутренняя координата).
  // Пол полости берётся у самой геометрии (PlannedZCav — единственный носитель
  // формулы посадки), а не пересчитывается здесь копией формулы: копия
  // разъехалась бы молча при любой правке сосуда или seatOnFloor.
  if (gapMm >= 0.0) {
    if (liftMm != 0.0) {
      std::cerr << "gap= и lift= заданы одновременно — это два имени одной "
                   "величины, оставьте одно\n";
      return 9;
    }
    const double zFloor = det->PlannedZCav() - sh.hzCav;
    // Насколько НИЖНЯЯ грань корпуса ниже точки, которую двигает fLiftMm.
    //
    // Лёжа fLiftMm сдвигает НАЧАЛО ОБЪЁМА "case", а оно лежит в середине
    // корпуса по толщине — то есть низ отстоит на полтолщины caseY/2 = 8,75 мм.
    // Это НЕ crystalToFace (8,20): кристалл смещён от середины на 0,55 мм в
    // сторону фотоприёмника. Первая версия взяла 8,20 и дала зазор 24,45 вместо
    // 25,00 — поймано ЗАМЕРОМ построенной геометрии ниже, не рассуждением;
    // комментарий в RCDetector.cc, обещавший «центр кристалла на fLiftMm»,
    // оказался неверен и исправлен там же.
    //
    // Стоя fLiftMm сдвигает корпус вдоль его оси, и низ — это нос прибора,
    // отстоящий от кристалла на crystalZ0 = 12,00 мм.
    const double belowAnchor = det->fHorizontal ? 0.5 * det->fDev.caseY
                                                : det->fDev.crystalZ0;
    det->fLiftMm = zFloor + gapMm + belowAnchor;
  }
  if (det->fHorizontal || det->fLiftMm != 0.0)
    G4cout << "!! посадка прибора: " << (det->fHorizontal ? "ГОРИЗОНТАЛЬНО" : "вертикально")
           << ", lift=" << det->fLiftMm << " мм (сдвиг корпуса по z)"
           << (gapMm >= 0.0 ? " из gap=" : " задан напрямую")
           << (gapMm >= 0.0 ? std::to_string(gapMm) : std::string())
           << " (P-008/P-013)" << G4endl;
  rm->SetUserInitialization(det);

  // doBias: явный флаг bias= И (mode=trans ИЛИ mode=replay) — ДВЕ РАЗНЫЕ
  // схемы важности на ДВУХ разных наборах ячеек (слои защиты vs путь к
  // кристаллу), см. ниже. Найдено 12.08.2026 (задача №13): без биасинга
  // ВНУТРИ replay толстая защита давала <1 взвешенного попадания в кристалл
  // на миллион+ пересечений границы полости — trans-биасинг решает только
  // половину задачи (доставляет квант до полости), а "полость -> кристалл"
  // само по себе редкое событие (~0,1-1%), для которого до сих пор не было
  // ни одного уровня расщепления.
  // ⚠️ ИСТОРИЯ (задача №13, 12.08.2026): первая версия doBiasReplay давала
  // НЕВЕРНЫЙ спектр — внутри ОДНОГО replay-события ImportanceProcess может
  // расщепить трек на НЕСКОЛЬКО независимых клонов, каждый способен ОТДЕЛЬНО
  // долететь до кристалла, а старый учёт (один fEdep на G4Event) СКЛАДЫВАЛ
  // их вклады в одно фиктивное число (поймано эмпирически: переполнение
  // >3200 кэВ на 12 из 14 попаданий при максимальной реальной линии
  // 1460,8 кэВ). Исправлено: `ReplayEventAct`/`ReplayStepping`/
  // `LineageTracking` считают ПО КЛОНУ (см. класс выше), не по событию —
  // используются ТОЛЬКО для mode=replay, обычный `EventAct`/`Stepping`
  // (resp/sample/muon) не тронуты.
  const bool doBiasTrans = bias && (mode == "trans");
  const bool doBiasReplay = bias && (mode == "replay");
  const bool doBias = doBiasTrans || doBiasReplay;
  G4GeometrySampler* mgs = nullptr;
  if (doBias) mgs = new G4GeometrySampler(nullptr, "gamma");  // nullptr — см. комментарий у #include выше

  // muon: FTFP_BERT, не Phys() — мюону нужны ионизация+тормозное+пары
  // (доминируют в депозите) и мюон-ядерные (второстепенно); тот же выбор,
  // что уже провалидирован в cosmicmu.cc. Остальные режимы — Phys()
  // (G4EmStandardPhysics_option4) как раньше, ПОБИТОВО не меняются —
  // FTFP_BERT для них не подставляется (медленнее инициализация, лишние
  // адронные таблицы, и меняло бы уже провалидированные числа).
  G4VModularPhysicsList* physicsList =
      (mode == "muon") ? static_cast<G4VModularPhysicsList*>(new FTFP_BERT())
                       : static_cast<G4VModularPhysicsList*>(new Phys());
  // РегистрАция G4ImportanceBiasing строго ДО SetUserInitialization(physicsList)
  // — порядок из B01, сверено дважды (план, «Статус реализации», п.1).
  if (doBias) physicsList->RegisterPhysics(new G4ImportanceBiasing(mgs));
  rm->SetUserInitialization(physicsList);

  if (mode == "replay") rm->SetUserAction(new ReplayPrimary(&recs));
  else if (mode == "muon") {
    // Диск растёт с наружным габаритом защиты (см. разбор у MuGun) —
    // fOuterR/fOuterHz валидны ТОЛЬКО после Construct() внутри
    // rm->Initialize(), то есть позже этой точки. Но геометрия ЗНАЕТ свои
    // будущие размеры заранее и отвечает на них методами Planned*(), где
    // формула записана ОДИН раз (PbShield.cc). Прежде эти формулы дублировались
    // здесь руками; для короба с посадкой на дно такое дублирование стало
    // прямой ошибкой: zCav до Construct() ещё не пересчитан.
    //
    // Радиус диска берётся от ОПИСАННОЙ вокруг короба окружности. Полуширины
    // мало: диск, накрывающий грань, оставил бы углы защиты без облучения, а
    // мюоны идут по cos²θ и вклад углов не нулевой.
    const double rOutExp = det->PlannedOuterR();
    const double hzOutExp = det->PlannedOuterHz();
    BuildMuEnergyTable();
    gMuRDisk = (rDiskOverride > 0.0) ? rDiskOverride
                                     : rOutExp + DISK_MARGIN_MM;
    gMuZDisk = det->PlannedZCav() + hzOutExp + DISK_CLEAR_MM;
    // Мир обязан вместить ВЕСЬ диск, иначе часть первичных стартует снаружи
    // мира и молча теряется (найдено 13.08.2026, см. fWorldMinHalfXY).
    det->fWorldMinHalfXY = gMuRDisk + 20.0;
    rm->SetUserAction(new MuGun(gMuRDisk, gMuZDisk));
  } else if (mode == "resp") {
    // Поверхность ПОЛОСТИ. Розыгрыш — своим генератором по граням короба
    // (см. BoxSurfaceGun); энергию по-прежнему задаёт GPS ниже.
    rm->SetUserAction(new BoxSurfaceGun(sh.hxCav, sh.hyCav, sh.hzCav,
                                        det->PlannedZCav()));
  } else if (mode == "trans" || mode == "beam") {
    // НАРУЖНАЯ поверхность собранной защиты. beam переопределит источник на
    // точечный ниже по коду, но собственного действия ему всё равно нужно.
    rm->SetUserAction(new BoxSurfaceGun(det->PlannedOuterHx(),
                                        det->PlannedOuterHy(),
                                        det->PlannedOuterHz(),
                                        det->PlannedOuterZc()));
  } else rm->SetUserAction(new Primary());

  auto* runAct = new RunAct();
  rm->SetUserAction(runAct);
  // replay — ОТДЕЛЬНЫЙ путь (ReplayEventAct, учёт по клону, задача №13);
  // resp/sample/muon — обычный EventAct (одна история на событие, не тронут).
  EventAct* evtAct = nullptr;
  ReplayEventAct* replayEvt = nullptr;
  std::map<G4int, G4int> lineageMap;
  if (mode == "replay") {
    replayEvt = new ReplayEventAct(runAct, &lineageMap);
    rm->SetUserAction(replayEvt);
    rm->SetUserAction(new LineageTracking(&lineageMap));
  } else if (mode == "resp" || mode == "sample" || mode == "muon" ||
             mode == "pbself") {
    // pbself: см. пояснение ниже у SetUserAction(Stepping) — без EventAct
    // депозит некуда складывать, режим молча даёт нули.
    evtAct = new EventAct(runAct);
    rm->SetUserAction(evtAct);
  }

  rm->Initialize();
  // P-013: ЗАМЕР посадки по ПОСТРОЕННОЙ геометрии, а не по формуле, которой её
  // строили. Ровно этот класс ошибки уже стоил суток счёта (P-012: прогоны шли
  // вертикально, потому что ключ не передавался, — и ни одна печать этого не
  // показывала). Здесь берётся фактический физобъём "case", его габарит и
  // фактическое размещение; если знак или слагаемое в пересчёте gap->lift
  // перепутаны, зазор в этой строке разойдётся с заданным.
  if (withDevice) {
    auto* casePV = G4PhysicalVolumeStore::GetInstance()->GetVolume("case", false);
    if (casePV && casePV->GetLogicalVolume() &&
        casePV->GetLogicalVolume()->GetSolid()) {
      G4ThreeVector lo, hi;
      casePV->GetLogicalVolume()->GetSolid()->BoundingLimits(lo, hi);
      const G4RotationMatrix rot = casePV->GetObjectRotationValue();
      const G4ThreeVector tr = casePV->GetObjectTranslation();
      double zMin = 1e300, zMax = -1e300;
      for (int c = 0; c < 8; ++c) {
        const G4ThreeVector corner((c & 1) ? hi.x() : lo.x(),
                                   (c & 2) ? hi.y() : lo.y(),
                                   (c & 4) ? hi.z() : lo.z());
        const double z = (rot * corner + tr).z();
        zMin = std::min(zMin, z);
        zMax = std::max(zMax, z);
      }
      const double zFloor = (det->fSh.zCav - det->fSh.hzCav) * mm;
      G4cout << "!! ЗАМЕР посадки (построенная геометрия): пол полости z="
             << zFloor / mm << " мм, корпус z=" << zMin / mm << ".."
             << zMax / mm << " мм, ФАКТИЧЕСКИЙ ЗАЗОР=" << (zMin - zFloor) / mm
             << " мм, высота корпуса=" << (zMax - zMin) / mm << " мм" << G4endl;
      if (zMin < zFloor - 1e-6)
        G4cerr << "!! ВНИМАНИЕ: корпус УТОПЛЕН в дно полости на "
               << (zFloor - zMin) / mm << " мм — геометрия недостоверна\n";
    }
  }
  if (mode == "replay") {
    runAct->fOut = out;
    rm->SetUserAction(new ReplayStepping(replayEvt, det->fCrystalLV));
  } else if (mode == "resp" || mode == "sample" || mode == "muon" ||
             mode == "pbself") {
    // pbself добавлен 17.08. Без него режим отрабатывал с кодом 0, но БЕЗ
    // Stepping — депозит в кристалле не регистрировался вовсе (hits=0), и
    // выход уходил в умолчательный resp.csv. Тот же класс, что W-007:
    // «выход есть, результата нет». Проверять такие режимы только по наличию
    // попаданий, а не по коду возврата.
    runAct->fOut = out;
    rm->SetUserAction(new Stepping(evtAct, det->fCrystalLV));
  }

  if (doBias) {
    // G4IStore — синглтон; заполняется ПОСЛЕ Initialize() (реальная
    // геометрия построена InitializeGeometry() внутри Initialize(), СТРОГО
    // раньше InitializePhysics() — план, «Статус реализации», п.1, второй
    // разбор source).
    //
    // БЕЗ mgs->PrepareImportanceSampling()/Configure() ЗДЕСЬ НАМЕРЕННО.
    // G4ImportanceBiasing::ConstructProcess() уже вызвал их САМ, внутри
    // предыдущего rm->Initialize() (физлист строит процесс через тот же
    // G4GeometrySampler, который мы передали в его конструктор) — и
    // ImportanceProcess читает G4IStore::GetImportance() ДИНАМИЧЕСКИ на
    // каждом шаге, а не однократным снимком при Configure(). Повторный
    // ручной вызов PrepareImportanceSampling()/Configure() здесь (первая
    // попытка реализации, 12.08.2026) давал WARNING
    // "G4GeometrySampler::IsConfigured() ... use ClearSampling()" и падал
    // фатальным G4Exception GeomBias0002 при реальном прогоне — типичный
    // симптом двойной конфигурации одного сэмплера. Урок записан здесь,
    // не только в плане: следующая правка этого блока не должна возвращать
    // эти два вызова не проверив предположение заново.
    G4IStore* istore = G4IStore::GetInstance();

    if (doBiasTrans) {
      // Важность растёт К ПОЛОСТИ: слой depth=0 (fLayerDepth, ближайший к
      // полости) — самый редкий для кванта, идущего снаружи внутрь,
      // получает наибольшую важность impStep^(fNDepth-1).
      for (size_t i = 0; i < det->fLayerPV.size(); ++i) {
        const int depthFromSource = det->fNDepth - 1 - det->fLayerDepth[i];
        const double imp = std::pow(impStep, depthFromSource);
        istore->AddImportanceGeometryCell(imp, *det->fLayerPV[i], 0);
      }
      // "world" и "cavity" ОБА обязаны быть зарегистрированы: G4IStore::
      // GetImportance() бросает фатальный GeomBias0002 для ЛЮБОЙ
      // незарегистрированной ячейки (проверено на практике 12.08.2026, не
      // предположение). До правки 12.08.2026 "world" был ОДНИМ физическим
      // объёмом на обе роли (полость и внешнее пространство источника) —
      // одному указателю нельзя присвоить две разные важности сразу
      // (GeomBias1001, ipre_over_ipost=1/128 — источник, рождённый в
      // "world" с важностью полости, получал обвальный roulette на первом
      // же шаге). С отдельной "cavity" (PbShield.hh/.cc) конфликт снят:
      //   world (снаружи, где рождается источник) — важность 1, естественная;
      //   cavity (полость, куда должны дойти редкие расщеплённые истории) —
      //     важность слоя depth=0, БЕЗ разрыва на входе (trStep не рулетит
      //     трек ровно в момент, когда его нужно записать).
      G4VPhysicalVolume* worldPV = G4PhysicalVolumeStore::GetInstance()->GetVolume("world");
      if (worldPV) istore->AddImportanceGeometryCell(1.0, *worldPV, 0);
      else G4cerr << "!! doBias: физический объём \"world\" не найден для G4IStore\n";

      if (det->fCavityPV) {
        const double impCav = std::pow(impStep, det->fNDepth - 1);
        istore->AddImportanceGeometryCell(impCav, *det->fCavityPV, 0);
      } else {
        G4cerr << "!! doBias: физический объём \"cavity\" не построен (fWithDevice=true?)\n";
      }
    }

    if (doBiasReplay) {
      // Задача №13 (12.08.2026): trans-биасинг доставляет квант ДО полости,
      // но "полость -> кристалл" — редкое событие само по себе (~0,1-1%,
      // без единого уровня расщепления). Здесь — второй, независимый набор
      // ячеек НА ПУТИ к кристаллу: world -> case -> caseAir -> reflector ->
      // crystal — РЕАЛЬНАЯ вложенность G4 (case мать caseAir, caseAir мать
      // reflector, reflector мать crystal — сверено построчно с
      // RCDetector.cc::BuildDevice, не предположение). Сосуд (vessel) в эту
      // вложенность НЕ входит (его объёмы — прямые дети world, не device) —
      // не участвует в расщеплении, важность 1 как у world, безвредно.
      //
      // Регистрируем ВСЕ физобъёмы geometрии разом важностью 1 (blanket) —
      // иначе GeomBias0002 на первом же непредусмотренном имени (sipm/pcb/
      // display/batt/детали сосуда — их много, перечислять руками ловушка).
      // Затем ChangeImportance() поднимает ТОЛЬКО путь к кристаллу.
      G4PhysicalVolumeStore* pvStore = G4PhysicalVolumeStore::GetInstance();
      for (auto* pv : *pvStore) istore->AddImportanceGeometryCell(1.0, *pv, 0);

      static const char* kCrystalPath[] = {"case", "caseAir", "reflector", "crystal"};
      double imp = 1.0;
      double impCase = 1.0, impCaseAir = 1.0;
      for (const char* nm : kCrystalPath) {
        imp *= cryStep;
        G4VPhysicalVolume* pv = pvStore->GetVolume(nm);
        if (pv) istore->ChangeImportance(imp, *pv, 0);
        else G4cerr << "!! doBiasReplay: физобъём \"" << nm << "\" не найден\n";
        if (std::strcmp(nm, "case") == 0) impCase = imp;
        if (std::strcmp(nm, "caseAir") == 0) impCaseAir = imp;
      }

      // СОСЕДИ ПО caseAir — важность матери, а не blanket-1 (правка 16.08.2026).
      // Blanket выше ставит 1 ВСЕМ объёмам, а путь к кристаллу поднимает только
      // case/caseAir/reflector/crystal. Но sipm, sipmPcb, display и pcb — прямые
      // дети caseAir (RCDetector.cc:349-363), они лежат вплотную к объёмам с
      // высокой важностью, оставаясь с важностью 1. На границе получается обрыв
      // в cryStep^2 (caseAir->сосед) или cryStep^3 (reflector->sipm через окно),
      // тогда как G4ImportanceAlgorithm требует отношение соседних ячеек в
      // [0.25, 4] и на выходе за диапазон печатает GeomBias1001 и ведёт себя
      // численно неустойчиво (вплоть до порчи памяти — см. #SHIELD-5).
      //
      // Диагноз подтверждён арифметикой на ТРЁХ независимых прогонах, каждое
      // число сошлось точно, а не «похоже»:
      //   crystep=2 -> ipre_over_ipost=8   = reflector(2^3) / sipm(1)
      //   crystep=4 -> ipre_over_ipost=16  = caseAir(4^2)  / sipm|display|pcb(1)
      //   crystep=8 -> ipre_over_ipost=0.125 = world(1) / case(8)
      // Третий случай — отдельная причина (сам cryStep=8 слишком велик для
      // одного уровня); он лечится не здесь, а ограничением cryStep <= 4.
      //
      // Физика при этом не меняется: важность влияет только на расщепление и
      // рулетку, оценка остаётся несмещённой (Sum(w) сохраняется — проверено:
      // 39,35 без биасинга против 41,03 с ним на одном входе).
      // Список сверен построчно с RCDetector.cc::BuildDevice: sipm (стр. 349),
      // sipmPcb (351), display (358), pcb (361), batt (366) — все Put(..., airLV,
      // ...), то есть прямые дети caseAir. Проверено чтением файла, не по памяти:
      // batt по названию похож на деталь корпуса, но лежит именно в caseAir.
      static const char* kCaseAirSiblings[] = {"sipm", "sipmPcb", "display",
                                               "pcb", "batt"};
      for (const char* nm : kCaseAirSiblings) {
        G4VPhysicalVolume* pv = pvStore->GetVolume(nm);
        if (pv) istore->ChangeImportance(impCaseAir, *pv, 0);
        // отсутствие объёма здесь не ошибка: состав прибора зависит от сборки
      }
      (void)impCase;   // оставлен для симметрии, если появятся дети case

      // P-009 (17.08): ОГРАНИЧЕНИЕ ШАГА В ТОНКИХ ЯЧЕЙКАХ ВАЖНОСТИ.
      //
      // Симптом: при ГОРИЗОНТАЛЬНОЙ посадке прибора (P-008) прогон падал с
      // `GeomBias1001: ipre_over_ipost = 0.0625`. Число разобрано точно:
      // 0,0625 = 1/4² = 1/cryStep², то есть трек переходил из ячейки важности 1
      // сразу в `caseAir` (важность cryStep²), МИНУЯ `case`.
      //
      // Почему только при повороте: вертикально стоящий прибор квант проходил
      // вдоль оси через нос/хвост (wallNose/wallTail = 2,0 мм); лёжа он входит
      // через большую грань, где стенка `wallFace` всего 1,50 мм. Свободный
      // пробег гамма в ABS на сотнях кэВ — сантиметры, поэтому такую стенку
      // трек проскакивает ЦЕЛИКОМ за один шаг, и промежуточная ячейка `case`
      // не регистрируется вовсе. G4ImportanceAlgorithm видит скачок важности
      // через уровень и уходит за допустимый диапазон [0.25, 4].
      //
      // Тот же класс, что P-004 (там предлагалось для плёнки reflector 0,05 мм,
      // но не было применено — процесса-исполнителя в физлисте не существовало;
      // теперь `G4StepLimiterPhysics` зарегистрирован, см. класс Phys выше).
      //
      // ПЕРВАЯ ПОПЫТКА (0,25 мм на все ячейки) ПРОВАЛИЛАСЬ — приёмка 17.08
      // показала ТОТ ЖЕ `ipre_over_ipost = 0.0625`. Разбор числа вместо
      // подбора следующего значения:
      //   важности пути = case 4, caseAir 16, reflector 64, crystal 256
      //   (cryStep, cryStep², cryStep³, cryStep⁴ при cryStep=4).
      //   0,0625 = 1/16 даёт пара **caseAir(16) → crystal(256)** — то есть
      //   трек проскакивает reflector ЦЕЛИКОМ.
      // Со стороны окна фотоприёмника reflector — не 1,25 мм, а ПЛЁНКА
      // 0,05 мм (`win` в RCDetector.cc::BuildDevice): именно там кристалл
      // почти обнажён. Лимит 0,25 мм ТОЛЩЕ этой плёнки и сработать не мог.
      // Это ровно тот дефект, что описан в P-004 как гипотеза.
      //
      // Поэтому лимит РАЗНЫЙ по объёмам: reflector получает 0,01 мм (пять
      // шагов на плёнку 0,05 мм), остальные — 0,25 мм (этого хватает на
      // стенки 1,25-1,50 мм и не плодит шаги в объёмах покрупнее).
      // Цена: до ~125 шагов на пересечение толстой части чашки. На физику не
      // влияет — G4UserLimits меняет дискретизацию транспорта, не сечения.
      const double kStepThin = 0.01 * mm;   // reflector: плёнка окна 0,05 мм
      const double kStepWall = 0.25 * mm;   // прочие ячейки пути
      for (const char* nm : kCrystalPath) {
        G4VPhysicalVolume* pv = pvStore->GetVolume(nm);
        if (!pv || !pv->GetLogicalVolume()) continue;
        const double lim = (std::strcmp(nm, "reflector") == 0) ? kStepThin : kStepWall;
        pv->GetLogicalVolume()->SetUserLimits(new G4UserLimits(lim));
      }
      // Соседи caseAir тоже граничат с ячейками высокой важности — плёнки
      // между ними нет, но шаг ограничиваем той же мерой, чтобы переход
      // сосед→caseAir не проскакивался на тонких деталях (sipm 0,40 мм,
      // sipmPcb 0,60 мм, display/pcb 1,00 мм — все тоньше 0,25 мм лимита).
      for (const char* nm : kCaseAirSiblings) {
        G4VPhysicalVolume* pv = pvStore->GetVolume(nm);
        if (pv && pv->GetLogicalVolume())
          pv->GetLogicalVolume()->SetUserLimits(new G4UserLimits(kStepThin));
      }
      G4cout << "!! doBiasReplay: шаг ограничен " << kStepThin / mm
             << " мм (reflector и тонкие соседи) / " << kStepWall / mm
             << " мм (прочие ячейки пути) — P-009" << G4endl;

      if (cryStep > 4.0)
        G4cerr << "!! doBiasReplay: crystep=" << cryStep << " > 4 — отношение "
                  "важностей соседних ячеек выйдет за допустимый G4 диапазон "
                  "[0.25, 4] уже на переходе world->case\n";
    }
  }

  G4cout << "# src_sha1 = " << RCPB_SRC_SHA1 << G4endl;
  G4cout << "# git_describe = " << RCPB_GIT_DESCRIBE << G4endl;

  auto* ui = G4UImanager::GetUIpointer();

  if (mode == "geom") {
    // Проверка пересечений — штатным средством Geant4, по всей иерархии.
    // Ловит и щели между кусками слоя, и наползание защиты на прибор.
    // Самопроверка объёма в BuildShield ловит другое (арифметику стыков),
    // поэтому проверки дополняют друг друга, а не дублируют.
    ui->ApplyCommand("/geometry/test/tolerance 0.001 mm");
    ui->ApplyCommand("/geometry/test/recursion_depth 2");
    ui->ApplyCommand("/geometry/test/run");
    G4cout << "RESULT geom ok  pb=" << det->fSh.pb << " cu=" << det->fSh.cu
           << " cd=" << det->fSh.cd << " massPb=" << det->fMassPb
           << " massCd=" << det->fMassCd << " massCu=" << det->fMassCu
           << " rOut=" << det->fOuterR << " hzOut=" << det->fOuterHz
           << " ndepth=" << det->fNDepth << G4endl;

  } else if (mode == "resp") {
    // Источник: изотропный на границе полости (r=rCav, z в [zCav-hzCav,
    // zCav+hzCav]), косинусный закон внутрь — то же тождество Ф=4N/S, что
    // используется в run_bg.py и wallfield.cc, здесь не для нормировки, а
    // чтобы угловое распределение падающих квантов было физически верным
    // (то, что в реальности даст изотропное поле снаружи).
    // Позицию и направление задаёт BoxSurfaceGun (розыгрыш по граням короба),
    // GPS отвечает ТОЛЬКО за сорт частицы и энергию. Поэтому pos/type здесь
    // Point и ang/type iso: их значения всё равно перезаписываются у готовой
    // вершины, и оставлять тут Surface+Cylinder было бы обманом читателя.
    char buf[256];
    ui->ApplyCommand("/gps/particle gamma");
    ui->ApplyCommand("/gps/pos/type Point");
    ui->ApplyCommand("/gps/ene/type Mono");
    std::snprintf(buf, sizeof(buf), "/gps/ene/mono %.4f keV", eKeV);
    ui->ApplyCommand(buf);

    const double zc = det->PlannedZCav();
    std::snprintf(buf, sizeof(buf),
                  "hxcav=%.2f hycav=%.2f hzcav=%.2f zcav=%.2f S_mm2=%.4e e=%.2f",
                  sh.hxCav, sh.hyCav, sh.hzCav, zc,
                  BoxSurfaceGun::SurfaceArea(sh.hxCav, sh.hyCav, sh.hzCav), eKeV);
    runAct->fTag = buf;

    std::snprintf(buf, sizeof(buf), "/run/beamOn %ld", std::lround(nprim));
    ui->ApplyCommand(buf);

  } else if (mode == "trans") {
    // Источник: изотропный-косинусный на НАРУЖНОЙ поверхности собранной
    // защиты (то же тождество Ф=4N/S, что в resp/run_bg.py) — это то, что
    // реально видит защита снаружи: поле помещения, изотропное на её
    // габарите. Пока МОНОЭНЕРГЕТИЧЕСКИЙ источник (arg e=...): полноценный
    // спектр поля ЕРН подключается отдельным шагом, когда угловая проверка
    // ниже подтвердит саму схему.
    // Позиция и направление — BoxSurfaceGun по граням НАРУЖНОГО габарита
    // (задан выше при SetUserAction). GPS отвечает только за сорт частицы и
    // энергию, поэтому pos/type Point: значение всё равно перезаписывается.
    char buf[256];
    ui->ApplyCommand("/gps/particle gamma");
    ui->ApplyCommand("/gps/pos/type Point");
    // specmac= — непрерывный спектр поля (field_spectrum_K/Ra/Th.mac,
    // GPS Arb, тот же формат, что gen_macros() в fit_room_field.py) вместо
    // моноэнергетического источника. Макрос пишет только particle/ene/hist
    // (сверено построчно с analysis/fit_room_field.py:gen_macros()), позицию
    // и угол не трогает — а если бы и тронул, BoxSurfaceGun перезапишет их у
    // готовой вершины, так что схема розыгрыша защищена от такой правки.
    // Нормировка Ф=4N/S та же, что у mono: абсолютный масштаб (Бк/кг ->
    // имп/с) считается СНАРУЖИ, в python-драйвере, по wf_S.csv fluence_total
    // и ФАКТИЧЕСКОЙ площади поверхности защиты. ⚠ Для короба эта площадь —
    // S_mm2 из тега ниже, а НЕ пересчёт по rOut/hzOut как для цилиндра.
    if (!specmac.empty()) {
      std::snprintf(buf, sizeof(buf), "/control/execute %s", specmac.c_str());
      ui->ApplyCommand(buf);
    } else {
      ui->ApplyCommand("/gps/ene/type Mono");
      std::snprintf(buf, sizeof(buf), "/gps/ene/mono %.4f keV", eKeV);
      ui->ApplyCommand(buf);
    }

    const double sOut = BoxSurfaceGun::SurfaceArea(
        det->PlannedOuterHx(), det->PlannedOuterHy(), det->PlannedOuterHz());
    std::snprintf(buf, sizeof(buf),
                  "pb=%.2f cu=%.2f cd=%.2f e=%.2f S_mm2=%.6e specmac=%s",
                  sh.pb, sh.cu, sh.cd, eKeV, sOut,
                  specmac.empty() ? "-" : specmac.c_str());
    runAct->fTag = buf;

    FILE* raw = std::fopen(out.c_str(), "w");
    if (!raw) { std::cerr << "не открыть " << out << "\n"; delete rm; return 4; }
    std::fprintf(raw, "# пересечения границы полости, стадия 1 (%s)\n",
                 doBias ? "importance biasing" : "небиасированная");
    std::fprintf(raw, "# src_sha1 = %s\n", RCPB_SRC_SHA1);
    // Габариты пишем полуразмерами короба И готовой площадью: downstream
    // обязан брать S отсюда, а не собирать её из размеров по своей копии
    // формулы. Ключ S_mm2 — единственный источник знаменателя для Ф = 4N/S.
    std::fprintf(raw,
                 "# pb=%.2f cu=%.2f cd=%.2f e_in_keV=%.2f "
                 "hxOut=%.2f hyOut=%.2f hzOut=%.2f zCav=%.2f S_mm2=%.6e lid=%s\n",
                 sh.pb, sh.cu, sh.cd, eKeV, det->fOuterHx, det->fOuterHy,
                 det->fOuterHz, det->fSh.zCav,
                 BoxSurfaceGun::SurfaceArea(det->fOuterHx, det->fOuterHy,
                                            det->fOuterHz),
                 det->fWithLid ? "on" : "off");
    std::fprintf(raw, "# bias=%s impstep=%.4f\n", doBias ? "on" : "off", impStep);
    // N_primaries_stage1 — знаменатель для несмещённой оценки пропускания
    // (Σweight / N_primaries_stage1), нужен downstream-анализу replay:
    // recs.size() != N_primaries_stage1 при биасинге (расщепление меняет
    // число записей относительно числа реально запущенных первичных).
    std::fprintf(raw, "# N_primaries_stage1=%ld\n", std::lround(nprim));
    std::fprintf(raw, "E_keV,cosTheta,face,x_mm,y_mm,z_mm,dx,dy,dz,weight\n");

    auto* trStep = new TransStep(runAct, det, raw);
    rm->SetUserAction(trStep);
    runAct->fOut = out + ".hist.csv";

    std::snprintf(buf, sizeof(buf), "/run/beamOn %ld", std::lround(nprim));
    ui->ApplyCommand(buf);

    std::fclose(raw);
    G4cout << "RESULT trans crossings= " << trStep->Crossings()
           << " sumWeight= " << trStep->SumWeight()
           << " / " << std::lround(nprim)
           << " transmission= " << trStep->SumWeight() / nprim
           << " topCrossings= " << trStep->CrossingsTop()
           << " topSumWeight= " << trStep->SumWeightTop()
           << " topFraction= "
           << (trStep->SumWeight() > 0
                   ? trStep->SumWeightTop() / trStep->SumWeight() : 0.0)
           << " bias= " << (doBias ? "on" : "off")
           << " impstep= " << (doBias ? impStep : 0.0)
           << " file= " << out << G4endl;

  } else if (mode == "replay") {
    // Без GPS-команд: источник задаётся ReplayPrimary напрямую из recs.
    // По одной записи на событие, ровно recs.size() раз — без циклического
    // повтора здесь; повтор (если статистики resp не хватает) делается на
    // стороне драйвера явным дублированием входного файла, чтобы решение
    // «сколько раз повторять» было видно в логе прогона, а не спрятано внутри.
    char buf[256];
    const long nRep = std::max(1L, std::lround(repeat));
    const size_t nEvt = recs.size() * static_cast<size_t>(nRep);
    // repeat и seed попадают в тег, а тег — в шапку выходного спектра: без них
    // downstream не отличит «прогнали вход один раз» от «прогнали сорок тысяч»
    // и промахнётся в нормировке ровно в repeat раз.
    std::snprintf(buf, sizeof(buf), "in=%s N=%zu repeat=%ld seed=%ld",
                  in.c_str(), recs.size(), nRep, static_cast<long>(seed));
    runAct->fTag = buf;

    std::snprintf(buf, sizeof(buf), "/run/beamOn %zu", nEvt);
    ui->ApplyCommand(buf);

  } else if (mode == "beam") {
    // Точечный источник строго над верхней крышкой, направление точно (0,0,-1)
    // — не /gps/ang/type cos (там разброс направлений), а Iso с фиксированной
    // осью через /gps/direction эквивалент: GPS не даёт "точное направление"
    // напрямую, поэтому используем /gps/ang/type mono с тэта=0.
    char buf[256];
    ui->ApplyCommand("/gps/particle gamma");
    ui->ApplyCommand("/gps/pos/type Point");
    std::snprintf(buf, sizeof(buf), "/gps/pos/centre 0 0 %.4f mm",
                 sh.zCav + det->fOuterHz + 1.0);
    ui->ApplyCommand(buf);
    ui->ApplyCommand("/gps/ang/type iso");
    // Конвенция GPS: theta=0 -> направление (0,0,-1) (Pz=-cos(theta)).
    // Вырожденный диапазон mintheta=maxtheta=0 даёт ровно эту ось без разброса.
    ui->ApplyCommand("/gps/ang/mintheta 0 deg");
    ui->ApplyCommand("/gps/ang/maxtheta 0 deg");
    ui->ApplyCommand("/gps/ene/type Mono");
    std::snprintf(buf, sizeof(buf), "/gps/ene/mono %.4f keV", eKeV);
    ui->ApplyCommand(buf);

    auto* beamStep = new BeamStep(runAct);
    rm->SetUserAction(beamStep);
    runAct->fOut = out;
    std::snprintf(buf, sizeof(buf), "pb=%.2f cu=%.2f cd=%.2f e=%.2f (пучок)",
                  sh.pb, sh.cu, sh.cd, eKeV);
    runAct->fTag = buf;

    std::snprintf(buf, sizeof(buf), "/run/beamOn %ld", std::lround(nprim));
    ui->ApplyCommand(buf);

    const double t_cm = (sh.pb + sh.cu + sh.cd) / 10.0;
    G4cout << "RESULT beam survived= " << beamStep->Survived() << " / "
           << std::lround(nprim) << " T_measured= "
           << double(beamStep->Survived()) / nprim
           << "  (толщина по оси " << t_cm << " см; сравнить с exp(-mu*t) из mu_shield.csv)"
           << G4endl;

  } else if (mode == "sample") {
    // Собственная активность пробы (K-40 или Cs-137) СКВОЗЬ ПОСТРОЕННУЮ
    // ЗАЩИТУ (fWithDevice=fWithVessel=fWithShield=true — умолчания режима не
    // трогали, см. диспетчер выше). Источник — ион, полный распад со всеми
    // продуктами (гамма+бета+конверсия), а не моноэнергетическая линия —
    // тот же смысл, что у macros/nuclides.mac, но БЕЗ него: `/rc/outFile` —
    // custom UI команда main.cc, которой в shieldrun.cc нет и не нужна,
    // вывод уже идёт через RunAct/EventAct/Stepping (тот же путь, что resp).
    // Геометрия источника — m200-специфичная (радиус/полувысота/центр
    // полости сосуда), СВЕРЕНО построчно с macros/nuclides.mac.
    // P-015. Прежде здесь стояли ЧИСЛА m200 (R=33.24, halfz=33.25, z=-0.56) и
    // предупреждение в G4cerr при другом сосуде. Предупреждение в лог НЕ ПОПАЛО
    // (G4cerr мастер-потока в MT теряется), а с vessel=m500 охватывающий цилиндр
    // оказывался МЕНЬШЕ тела пробы: confine отбраковывает точки вне "sample",
    // но разыгрываются они только внутри цилиндра. Активность садилась в
    // ближнюю к прибору сердцевину сосуда -> эффективность завышена.
    // Габариты берём ЗАМЕРОМ построенного тела (тот же приём, что у посадки,
    // P-013), печатаем через G4cout — этот поток в лог доходит, проверено.
    char buf[256];
    double srcR = 33.24, srcHz = 33.25, srcZ = -0.56;
    auto* smpPV = G4PhysicalVolumeStore::GetInstance()->GetVolume("sample", false);
    if (smpPV) {
      G4ThreeVector lo, hi;
      smpPV->GetLogicalVolume()->GetSolid()->BoundingLimits(lo, hi);
      const G4ThreeVector tr = smpPV->GetObjectTranslation();
      srcR = std::max(std::max(std::fabs(lo.x() + tr.x()), std::fabs(hi.x() + tr.x())),
                      std::max(std::fabs(lo.y() + tr.y()), std::fabs(hi.y() + tr.y()))) / mm;
      srcHz = 0.5 * (hi.z() - lo.z()) / mm;
      srcZ = (tr.z() + 0.5 * (lo.z() + hi.z())) / mm;
      G4cout << "!! ЗАМЕР источника пробы (построенная геометрия): vessel="
             << vessel << " R=" << srcR << " мм, halfz=" << srcHz
             << " мм, центр z=" << srcZ << " мм, объём пробы="
             << det->fSampleVolumeCm3 << " см3, масса="
             << det->fSampleVolumeCm3 * rho << " г" << G4endl;
    } else {
      G4cout << "!! mode=sample: тело 'sample' не найдено — источник по "
                "умолчаниям m200, результат геометрически неверен" << G4endl;
    }
    ui->ApplyCommand("/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+30 ns");
    ui->ApplyCommand("/gps/particle ion");
    ui->ApplyCommand("/gps/energy 0 keV");
    ui->ApplyCommand("/gps/pos/type Volume");
    ui->ApplyCommand("/gps/pos/shape Cylinder");
    std::snprintf(buf, sizeof(buf), "/gps/pos/radius %.3f mm", srcR);
    ui->ApplyCommand(buf);
    std::snprintf(buf, sizeof(buf), "/gps/pos/halfz %.3f mm", srcHz);
    ui->ApplyCommand(buf);
    std::snprintf(buf, sizeof(buf), "/gps/pos/centre 0 0 %.3f mm", srcZ);
    ui->ApplyCommand(buf);
    ui->ApplyCommand("/gps/pos/confine sample");
    ui->ApplyCommand("/gps/ang/type iso");
    if (nuc == "K40") {
      ui->ApplyCommand("/process/had/rdm/nucleusLimits 40 40 19 19");
      ui->ApplyCommand("/gps/ion 19 40 0 0");
    } else if (nuc == "Cs137") {
      ui->ApplyCommand("/process/had/rdm/nucleusLimits 137 137 55 56");
      ui->ApplyCommand("/gps/ion 55 137 0 0");
    } else {
      G4cerr << "неизвестный nuc=" << nuc << " (K40 | Cs137)\n";
      delete rm;
      return 7;
    }

    std::snprintf(buf, sizeof(buf), "nuc=%s pb=%.2f cu=%.2f cd=%.2f",
                  nuc.c_str(), sh.pb, sh.cu, sh.cd);
    runAct->fTag = buf;

    std::snprintf(buf, sizeof(buf), "/run/beamOn %ld", std::lround(nprim));
    ui->ApplyCommand(buf);

  } else if (mode == "pbself") {
    // СОБСТВЕННАЯ АКТИВНОСТЬ СВИНЦА ЗАЩИТЫ (Pb-210). Спека:
    // docs/spec-pb210-self.md. Заведено 17.08.2026 по #SHIELD-16.
    //
    // ЗАЧЕМ. Модель фона не содержит источника ВНУТРИ самого свинца, а свинец
    // защиты практически всегда содержит Pb-210 (T½ 22,2 года, переплавкой не
    // выводится; обычный коммерческий — 10-300 Бк/кг, измеренные реперы 67-91).
    // Форма его вклада совпадает с наблюдаемым дефицитом: по замеру
    // домик/открытый = 0,40 в полосе 20-60, 0,44 в 60-100, 0,47 в 100-300 и
    // уже 0,82 в 700-1500 — то есть недостаёт именно там, куда Pb-210 и бьёт.
    //
    // ФИЗИКА — три канала, главный НЕ линия (Heusser, Annu.Rev.Nucl.Part.Sci.
    // 45 (1995) 543, стр. 551, цитата сверена с растром страницы):
    //   1. ХРИ свинца 72,8/75,0/84,9/87,4 кэВ — возбуждается ТОРМОЗНЫМ от бета
    //      Bi-210 (E_max 1161 кэВ), потому что сама гамма 46,5 кэВ лежит ниже
    //      K-края свинца (88 кэВ) и K-оболочку возбудить не может;
    //   2. тормозной континуум, максимум около 170 кэВ;
    //   3. гамма 46,5 кэВ Pb-210 — выход всего 4,25 % (переход M1, α_T=17,9).
    // Поэтому считается ПОЛНАЯ цепочка Pb-210 → Bi-210 → Po-210, а не линия:
    // расчёт «одной линии 46,5» потерял бы оба главных канала.
    //
    // ГЕОМЕТРИЯ ИСТОЧНИКА. Свинец — не одно тело, а пять при открытом верхе
    // (sh0_Pb_xhi/xlo/yhi/ylo/bot, см. PbShield.cc::BuildShield), а
    // /gps/pos/confine принимает ОДНО имя. Поэтому режим считает ОДНУ ячейку
    // за прогон (ключ cell=), а складывает их драйвер с весами по массе.
    // Розыгрыш ведётся в боксе, покрывающем весь габарит защиты; точки вне
    // указанной ячейки отбраковывает сам GPS штатным confine.
    //
    // НОРМИРОВКА. Прогон ведётся на 1 Бк/кг (как wallfield для ЕРН);
    // абсолютная активность подставляется в драйвере, чтобы сценарии
    // 25/60/100/300 Бк/кг пересчитывались без нового прогона и чтобы число
    // нельзя было незаметно подогнать под измерение.
    char buf[256];
    if (pbCell.empty()) {
      G4cerr << "!! mode=pbself: нужен ключ cell=<имя объёма свинца>, например "
                "cell=sh0_Pb_bot (см. PbShield.cc: sh<k>_Pb_{xhi,xlo,yhi,ylo,bot,top})\n";
      return 8;
    }
    const double hxO = det->PlannedOuterHx(), hyO = det->PlannedOuterHy();
    const double hzO = det->PlannedOuterHz(), zcO = det->PlannedOuterZc();
    ui->ApplyCommand("/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+30 ns");
    // Цепочка обрезана на Bi-210 (Z=82..83), Po-210 НЕ разыгрывается.
    //
    // ПОЧЕМУ ЭТО НЕ ПОТЕРЯ ФИЗИКИ. Нужны два канала: бета Bi-210 (E_max
    // 1161 кэВ) — источник тормозного и, через него, ХРИ свинца; и гамма
    // 46,5 кэВ самого Pb-210. Po-210 — альфа-излучатель (5304 кэВ), альфа из
    // металла не выходит вовсе, а его единственная гамма 803 кэВ имеет выход
    // 1,2·10⁻⁵ на распад (LNHB) — на четыре порядка ниже прочих каналов.
    //
    // ПОЧЕМУ ПРИШЛОСЬ ОБРЕЗАТЬ. С полной цепочкой (82..84) прогон ячейки
    // `sh0_Pb_bot` падал ДВАЖДЫ с разными зёрнами (9001 и 31337), не создав
    // файла, с `G4Exception TRACK001: Secondary with illegal time and/or
    // energy and/or momentum`. Дно — единственная ячейка, к которой прибор
    // прижат вплотную (стоит на полу полости), там больше всего близких
    // распадов. Стенки с полной цепочкой считались нормально, но ради
    // согласованности ВСЕ ячейки пересчитываются с обрезанной.
    ui->ApplyCommand("/process/had/rdm/nucleusLimits 210 210 82 83");
    ui->ApplyCommand("/gps/particle ion");
    ui->ApplyCommand("/gps/ion 82 210 0 0");
    ui->ApplyCommand("/gps/energy 0 keV");
    ui->ApplyCommand("/gps/pos/type Volume");
    ui->ApplyCommand("/gps/pos/shape Para");
    std::snprintf(buf, sizeof(buf), "/gps/pos/halfx %.3f mm", hxO);
    ui->ApplyCommand(buf);
    std::snprintf(buf, sizeof(buf), "/gps/pos/halfy %.3f mm", hyO);
    ui->ApplyCommand(buf);
    std::snprintf(buf, sizeof(buf), "/gps/pos/halfz %.3f mm", hzO);
    ui->ApplyCommand(buf);
    std::snprintf(buf, sizeof(buf), "/gps/pos/centre 0 0 %.3f mm", zcO);
    ui->ApplyCommand(buf);
    std::snprintf(buf, sizeof(buf), "/gps/pos/confine %s", pbCell.c_str());
    ui->ApplyCommand(buf);
    ui->ApplyCommand("/gps/ang/type iso");

    std::snprintf(buf, sizeof(buf),
                  "pbself cell=%s pb=%.2f massPb_kg=%.3f (1 Bq/kg, абсолют — в драйвере)",
                  pbCell.c_str(), sh.pb, det->fMassPb);
    runAct->fTag = buf;
    G4cout << "!! pbself: ячейка " << pbCell << ", габарит розыгрыша "
           << hxO << "x" << hyO << "x" << hzO << " мм, центр z=" << zcO
           << ", полная масса Pb " << det->fMassPb << " кг" << G4endl;

    std::snprintf(buf, sizeof(buf), "/run/beamOn %ld", std::lround(nprim));
    ui->ApplyCommand(buf);

  } else if (mode == "muon") {
    // Прямой (небиасированный) прогон сквозь ВСЮ сборку — источник (MuGun)
    // уже создан и передан в rm выше, GPS-команды не нужны (не G4GPS-based).
    // withDevice/withVessel/withShield — умолчания режима (все true).
    char buf[256];
    std::snprintf(buf, sizeof(buf),
                  "muon pb=%.2f cu=%.2f cd=%.2f rdisk_mm=%.2f zdisk_mm=%.2f",
                  sh.pb, sh.cu, sh.cd, gMuRDisk, gMuZDisk);
    runAct->fTag = buf;

    std::snprintf(buf, sizeof(buf), "/run/beamOn %ld", std::lround(nprim));
    ui->ApplyCommand(buf);

  } else {
    std::cerr << "режим '" << mode << "' ещё не реализован\n";
    delete rm;
    return 3;
  }

  // Порядок из B01 (план, п.1, шаг 7): OpenGeometry() ДО delete runManager —
  // иначе биасинг-хранилища (G4IStore и связанные с ним geometry cells)
  // подчищаются некорректно. Без биасинга безвредно, но выполняем всегда
  // для единообразия одного пути кода.
  G4GeometryManager::GetInstance()->OpenGeometry();
  delete rm;
  return 0;
}
