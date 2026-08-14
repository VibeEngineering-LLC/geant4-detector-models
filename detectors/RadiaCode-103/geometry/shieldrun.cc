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
  double fRCav, fHzCav, fZCav;
  FILE* fOut;
  long fCrossings = 0;      // сырое число пересечений (после расщепления) — диагностика
  double fSumW = 0.0;       // Σвес — несмещённая оценка; ЭТО число, не fCrossings, делить на nprim
public:
  TransStep(RunAct* r, const RCShieldDetector* det, FILE* out)
      : fRun(r), fLayerPV(&det->fLayerPV), fLayerDepth(&det->fLayerDepth),
        fRCav(det->fSh.rCav), fHzCav(det->fSh.hzCav), fZCav(det->fSh.zCav),
        fOut(out) {}
  long Crossings() const { return fCrossings; }
  double SumWeight() const { return fSumW; }

  void UserSteppingAction(const G4Step* s) override {
    if (s->GetTrack()->GetDefinition()->GetPDGEncoding() != 22) return;  // gamma
    if (s->GetPostStepPoint()->GetStepStatus() != fGeomBoundary) return;
    auto* preV = s->GetPreStepPoint()->GetPhysicalVolume();
    auto* postV = s->GetPostStepPoint()->GetPhysicalVolume();
    // "cavity", не "world" — см. PbShield.hh, правка 12.08.2026.
    if (!preV || !postV || postV->GetName() != "cavity") return;

    bool depth0 = false;
    for (size_t i = 0; i < fLayerPV->size(); ++i)
      if ((*fLayerPV)[i] == preV && (*fLayerDepth)[i] == 0) { depth0 = true; break; }
    if (!depth0) return;

    const G4ThreeVector p = s->GetPostStepPoint()->GetPosition() / mm;
    const G4ThreeVector dir = s->GetPostStepPoint()->GetMomentumDirection();
    const double r = std::sqrt(p.x() * p.x() + p.y() * p.y());
    G4ThreeVector outNorm;
    const char* face;
    if (std::abs(r - fRCav) < std::abs(std::abs(p.z() - fZCav) - fHzCav)) {
      outNorm = (r > 0) ? G4ThreeVector(p.x() / r, p.y() / r, 0) : G4ThreeVector(1, 0, 0);
      face = "side";
    } else {
      outNorm = G4ThreeVector(0, 0, (p.z() > fZCav) ? 1 : -1);
      face = "cap";
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

    const double cosT = std::cbrt(G4UniformRand());   // P(cosT) ~ cosT^2
    const double sinT = std::sqrt(std::max(0.0, 1.0 - cosT * cosT));
    const double phd = twopi * G4UniformRand();
    fGun.SetParticleMomentumDirection(
        G4ThreeVector(sinT * std::cos(phd), sinT * std::sin(phd), -cosT));

    fGun.SetParticleEnergy(SampleMuEnergyGeV() * GeV);
    fGun.GeneratePrimaryVertex(e);
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
  // rdisk — радиус диск-источника мюонов, мм. 0 = авто (rOut+DISK_MARGIN_MM,
  // прежнее поведение, побитово). Введён 13.08.2026 для проверки НАСЫЩЕНИЯ:
  // при cos²θ, θ<70° мюон приходит в кристалл со старта на расстоянии до
  // ~2,75·Z_disk от оси, тогда как авто-радиус даёт всего rOut+40 мм — то
  // есть диск ОБРЕЗАЕТ наклонные треки, и тем сильнее, чем он меньше.
  // Физический ответ (j·πR²·eff(R)) обязан выходить на полку с ростом R;
  // пока не вышел — поправка (R/70)² в run_shield_grid.mu_cps() неполна.
  double rDiskOverride = 0.0;
  bool withVessel = true, withShield = true, bias = false, withLid = true;

  for (int i = 2; i < argc; ++i) {
    const char* a = argv[i];
    if (KeyVal(a, "pb", &sh.pb)) continue;
    if (KeyVal(a, "cu", &sh.cu)) continue;
    if (KeyVal(a, "cd", &sh.cd)) continue;
    if (KeyVal(a, "rcav", &sh.rCav)) continue;
    if (KeyVal(a, "hzcav", &sh.hzCav)) continue;
    if (KeyVal(a, "zcav", &sh.zCav)) continue;
    if (KeyVal(a, "nshell", &nshell)) continue;
    if (KeyVal(a, "rho", &rho)) continue;
    if (KeyVal(a, "e", &eKeV)) continue;
    if (KeyVal(a, "nprim", &nprim)) continue;
    if (KeyVal(a, "impstep", &impStep)) continue;
    if (KeyVal(a, "crystep", &cryStep)) continue;
    if (KeyVal(a, "rdisk", &rDiskOverride)) continue;
    if (KeyStr(a, "vessel", &vessel)) continue;
    if (KeyStr(a, "matrix", &matrix)) continue;
    if (KeyStr(a, "out", &out)) continue;
    if (KeyStr(a, "in", &in)) continue;
    if (KeyStr(a, "specmac", &specmac)) continue;
    if (KeyStr(a, "nuc", &nuc)) continue;
    if (std::strcmp(a, "novessel") == 0) { withVessel = false; continue; }
    if (std::strcmp(a, "noshield") == 0) { withShield = false; continue; }
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
  // replay: полная геометрия (прибор+сосуд), БЕЗ свинца — его роль уже учтена
  // тем, что до этой точки дошли только реально выжившие в trans пересечения.
  if (mode == "replay") withShield = false;
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

  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  auto* det = new RCShieldDetector(withVessel, withShield, withDevice);
  det->fSh = sh;
  det->fVes = VesselGeom::Preset(vessel);
  det->fVes.sampleMatrix = matrix;
  det->fVes.sampleDensity = rho;
  det->fWithLid = withLid;
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
    // fOuterR/fOuterHz валидны ТОЛЬКО после rm->SetUserInitialization(det)
    // (уже случилось выше) и реально заполняются позже, при Construct()
    // внутри rm->Initialize() — но геометрия ЗНАЕТ свои будущие размеры
    // заранее: rCav/hzCav фиксированы, шаг слоёв известен из sh, поэтому
    // безопасно посчитать ожидаемый rOut/hzOut ЗДЕСЬ же по формуле
    // PbShield.cc (r+=d, h+=d по каждому слою), не дожидаясь Construct().
    const double wallT = sh.cu + sh.cd + sh.pb;
    const double rOutExp = sh.rCav + wallT, hzOutExp = sh.hzCav + wallT;
    BuildMuEnergyTable();
    gMuRDisk = (rDiskOverride > 0.0) ? rDiskOverride
                                     : rOutExp + DISK_MARGIN_MM;
    gMuZDisk = sh.zCav + hzOutExp + DISK_CLEAR_MM;
    // Мир обязан вместить ВЕСЬ диск, иначе часть первичных стартует снаружи
    // мира и молча теряется (найдено 13.08.2026, см. fWorldMinHalfXY).
    det->fWorldMinHalfXY = gMuRDisk + 20.0;
    rm->SetUserAction(new MuGun(gMuRDisk, gMuZDisk));
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
  } else if (mode == "resp" || mode == "sample" || mode == "muon") {
    evtAct = new EventAct(runAct);
    rm->SetUserAction(evtAct);
  }

  rm->Initialize();
  if (mode == "replay") {
    runAct->fOut = out;
    rm->SetUserAction(new ReplayStepping(replayEvt, det->fCrystalLV));
  } else if (mode == "resp" || mode == "sample" || mode == "muon") {
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
      for (const char* nm : kCrystalPath) {
        imp *= cryStep;
        G4VPhysicalVolume* pv = pvStore->GetVolume(nm);
        if (pv) istore->ChangeImportance(imp, *pv, 0);
        else G4cerr << "!! doBiasReplay: физобъём \"" << nm << "\" не найден\n";
      }
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
    char buf[256];
    ui->ApplyCommand("/gps/particle gamma");
    ui->ApplyCommand("/gps/pos/type Surface");
    ui->ApplyCommand("/gps/pos/shape Cylinder");
    std::snprintf(buf, sizeof(buf), "/gps/pos/radius %.4f mm", sh.rCav);
    ui->ApplyCommand(buf);
    std::snprintf(buf, sizeof(buf), "/gps/pos/halfz %.4f mm", sh.hzCav);
    ui->ApplyCommand(buf);
    std::snprintf(buf, sizeof(buf), "/gps/pos/centre 0 0 %.4f mm", sh.zCav);
    ui->ApplyCommand(buf);
    ui->ApplyCommand("/gps/ang/type cos");
    ui->ApplyCommand("/gps/ene/type Mono");
    std::snprintf(buf, sizeof(buf), "/gps/ene/mono %.4f keV", eKeV);
    ui->ApplyCommand(buf);

    std::snprintf(buf, sizeof(buf), "rcav=%.2f hzcav=%.2f zcav=%.2f e=%.2f",
                  sh.rCav, sh.hzCav, sh.zCav, eKeV);
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
    char buf[256];
    ui->ApplyCommand("/gps/particle gamma");
    ui->ApplyCommand("/gps/pos/type Surface");
    ui->ApplyCommand("/gps/pos/shape Cylinder");
    std::snprintf(buf, sizeof(buf), "/gps/pos/radius %.4f mm", det->fOuterR);
    ui->ApplyCommand(buf);
    std::snprintf(buf, sizeof(buf), "/gps/pos/halfz %.4f mm", det->fOuterHz);
    ui->ApplyCommand(buf);
    std::snprintf(buf, sizeof(buf), "/gps/pos/centre 0 0 %.4f mm", sh.zCav);
    ui->ApplyCommand(buf);
    ui->ApplyCommand("/gps/ang/type cos");
    // specmac= — непрерывный спектр поля (field_spectrum_K/Ra/Th.mac,
    // GPS Arb, тот же формат, что gen_macros() в fit_room_field.py) вместо
    // моноэнергетического источника. Позиция/угол уже заданы ВЫШЕ и НЕ
    // трогаются макросом (field_spectrum_S.mac пишет только particle/ene/
    // hist — сверено построчно с analysis/fit_room_field.py:gen_macros()).
    // Нормировка Ф=4N/S та же, что у mono: абсолютный масштаб (Бк/кг ->
    // имп/с) считается СНАРУЖИ, в python-драйвере, по wf_S.csv fluence_total
    // и ФАКТИЧЕСКОЙ площади наружной поверхности защиты (det->fOuterR/Hz на
    // этой толщине) — той же схемой, что fit_lines.model_cps_curve().
    if (!specmac.empty()) {
      std::snprintf(buf, sizeof(buf), "/control/execute %s", specmac.c_str());
      ui->ApplyCommand(buf);
    } else {
      ui->ApplyCommand("/gps/ene/type Mono");
      std::snprintf(buf, sizeof(buf), "/gps/ene/mono %.4f keV", eKeV);
      ui->ApplyCommand(buf);
    }

    std::snprintf(buf, sizeof(buf), "pb=%.2f cu=%.2f cd=%.2f e=%.2f rOut=%.2f specmac=%s",
                  sh.pb, sh.cu, sh.cd, eKeV, det->fOuterR,
                  specmac.empty() ? "-" : specmac.c_str());
    runAct->fTag = buf;

    FILE* raw = std::fopen(out.c_str(), "w");
    if (!raw) { std::cerr << "не открыть " << out << "\n"; delete rm; return 4; }
    std::fprintf(raw, "# пересечения границы полости, стадия 1 (%s)\n",
                 doBias ? "importance biasing" : "небиасированная");
    std::fprintf(raw, "# src_sha1 = %s\n", RCPB_SRC_SHA1);
    std::fprintf(raw, "# pb=%.2f cu=%.2f cd=%.2f e_in_keV=%.2f rOut=%.2f hzOut=%.2f\n",
                 sh.pb, sh.cu, sh.cd, eKeV, det->fOuterR, det->fOuterHz);
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
    std::snprintf(buf, sizeof(buf), "in=%s N=%zu", in.c_str(), recs.size());
    runAct->fTag = buf;

    std::snprintf(buf, sizeof(buf), "/run/beamOn %zu", recs.size());
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
    if (vessel != "m200") {
      G4cerr << "!! mode=sample: параметры источника (radius/halfz/centre) "
                "калиброваны ТОЛЬКО под m200, vessel=" << vessel
             << " даст геометрически неверный источник\n";
    }
    char buf[256];
    ui->ApplyCommand("/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+30 ns");
    ui->ApplyCommand("/gps/particle ion");
    ui->ApplyCommand("/gps/energy 0 keV");
    ui->ApplyCommand("/gps/pos/type Volume");
    ui->ApplyCommand("/gps/pos/shape Cylinder");
    ui->ApplyCommand("/gps/pos/radius 33.24 mm");
    ui->ApplyCommand("/gps/pos/halfz 33.25 mm");
    ui->ApplyCommand("/gps/pos/centre 0 0 -0.56 mm");
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
