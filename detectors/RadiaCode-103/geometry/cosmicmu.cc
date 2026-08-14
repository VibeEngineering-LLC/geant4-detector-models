// Космический мюонный континуум в кристалле RC-103 (голый прибор+сосуд, без
// защиты) — единичный отклик для подгонки в analysis/fit_room_field.py.
//
// ПОЧЕМУ ЭТОТ ФАЙЛ ПОЯВИЛСЯ. Подгонка поля ЕРН тремя сериями (K-40/Ra-226/
// Th-232) по семи энергетическим полосам показала системный, растущий с
// энергией недобор (0,98 на 20-100 кэВ -> 0,57-0,68 на 750-2700 кэВ),
// устойчивый к шестикратному росту статистики — не шум. Профиль «плоский
// континуум, всё заметнее там, где линейчатые спектры K/Ra/Th спадают» —
// подпись мюонного фона, не ошибка одной серии. См. план, разделы
// «Мюонный континуум — верифицированные параметры» и «Задачи №3 и №6
// оказались одной задачей».
//
// ФИЗИКА ВХОДНОГО СПЕКТРА (PDG RPP2020 §30.3.1, sci-search 11.08.2026,
// verified-facts.jsonl):
//   dNmu/dE dOmega = 0.14 E^-2.7 [1/(1+1.1E cosT/115) + 0.054/(1+1.1E cosT/850)]
//   E в ГэВ. ОГОВОРКА ИЗ ПЕРВОИСТОЧНИКА: строго применима при E>100/cosT ГэВ
//   — для интересующего нас диапазона (сотни МэВ - единицы ГэВ, то, что
//   реально пересекает кристалл 1 см³) это ЭКСТРАПОЛЯЦИЯ. Годится для формы
//   континуума, не для абсолютного числа с точностью до процентов.
//
// УПРОЩЕНИЕ, ЗАЯВЛЕННОЕ ЯВНО: энергия и угол разыгрываются НЕЗАВИСИМО —
// E из вертикальной формулы (cosT=1), угол отдельно по cos^2(T) (PDG:
// «угловое распределение мюонов у земли ~cos^2(T), характерно для Eмю~3 ГэВ»
// — приближение, а не точная связь; в самой формуле выше E и T на самом
// деле связаны через произведение E*cosT). Разъединение оправдано тем, что
// формула и так экстраполяция за пределы паспортной точности — добавлять
// точную 2D-связь означало бы ложную точность поверх уже приближённого входа.
//
// НОРМИРОВКА: НЕ добивается до абсолютного потока PDG (~1 см^-2 мин^-1).
// Отклик считается НА ОДИН МЮОН, пересёкший источник — амплитуда (сколько
// таких мюонов в секунду) ПОДГОНЯЕТСЯ NNLS вместе с K/Ra/Th в
// fit_room_field.py, как и активности. Абсолютный поток PDG используется
// только ПОСЛЕ подгонки как проверка порядка величины подобранной
// амплитуды, не как жёсткое ограничение при расчёте отклика.
//
// Запуск:  cosmicmu <N событий> <выходной .csv>
#if defined(__has_include)
#  if __has_include("rc_cosmicmu_provenance.hh")
#    include "rc_cosmicmu_provenance.hh"
#  endif
#endif
#ifndef RCCM_SRC_SHA1
#  define RCCM_SRC_SHA1 "БЕЗ-ШТАМПА"
#  define RCCM_GIT_DESCRIBE "БЕЗ-ШТАМПА"
#endif

#include "RCDetector.hh"

#include "FTFP_BERT.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4Run.hh"
#include "G4RunManagerFactory.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

// --- табличная функция Гайссера, вертикаль (cosT=1) --------------------
double GaisserVertical(double E_GeV) {
  const double t1 = 1.0 / (1.0 + 1.1 * E_GeV / 115.0);
  const double t2 = 0.054 / (1.0 + 1.1 * E_GeV / 850.0);
  return 0.14 * std::pow(E_GeV, -2.7) * (t1 + t2);
}

// Сетка 400 точек, лог-равномерно от 0.3 до 100 ГэВ — покрывает диапазон,
// где мюон ещё не сильно релятивистски однообразен (ниже 0.3 ГэВ формула
// Гайссера тем более не предназначена — там начинают работать другие
// эффекты, распад на лету и т.п., не рассматриваем) и где вклад уже
// пренебрежим (выше 100 ГэВ спектр падает как E^-2.7, вклад мал).
const double E_LO_GEV = 0.3, E_HI_GEV = 100.0;
const int N_EGRID = 400;
std::vector<double> gEGrid, gCumE;

void BuildEnergyTable() {
  gEGrid.resize(N_EGRID);
  gCumE.resize(N_EGRID);
  double s = 0;
  const double logLo = std::log(E_LO_GEV), logHi = std::log(E_HI_GEV);
  for (int i = 0; i < N_EGRID; ++i) {
    const double logE = logLo + (logHi - logLo) * i / (N_EGRID - 1);
    gEGrid[i] = std::exp(logE);
  }
  for (int i = 0; i < N_EGRID; ++i) {
    const double w = GaisserVertical(gEGrid[i]) * gEGrid[i];  // якобиан d(lnE)=dE/E
    s += w;
    gCumE[i] = s;
  }
  for (int i = 0; i < N_EGRID; ++i) gCumE[i] /= s;
}

double SampleEnergyGeV() {
  const double u = G4UniformRand();
  int lo = 0, hi = N_EGRID - 1;
  while (lo < hi) {
    const int mid = (lo + hi) / 2;
    if (gCumE[mid] < u) lo = mid + 1; else hi = mid;
  }
  return gEGrid[lo];
}

// --- источник: горизонтальный диск над сборкой, мюон mu- -----------------
// Диск радиусом R_DISK на высоте Z_DISK (в мировых координатах, центр
// кристалла = 0) — с запасом накрывает прибор+сосуд (внешний радиус сборки
// ~40 мм по CYL_M200 из fit_room_field.py, здесь R_DISK шире вдвое ради
// косых треков при не самом малом угле).
const double Z_DISK = 130.0;   // мм, выше хвоста прибора (111 мм)
const double R_DISK = 70.0;    // мм

class MuGun : public G4VUserPrimaryGeneratorAction {
  G4ParticleGun fGun{1};
public:
  MuGun() {
    fGun.SetParticleDefinition(
        G4ParticleTable::GetParticleTable()->FindParticle("mu-"));
  }
  void GeneratePrimaries(G4Event* e) override {
    const double r = R_DISK * std::sqrt(G4UniformRand());
    const double ph = twopi * G4UniformRand();
    fGun.SetParticlePosition(
        G4ThreeVector(r * std::cos(ph), r * std::sin(ph), Z_DISK * mm));

    // угол: cosT = u^(1/3) -> P(cosT) ~ cosT^2 (см. разбор в шапке файла)
    const double cosT = std::cbrt(G4UniformRand());
    const double sinT = std::sqrt(std::max(0.0, 1.0 - cosT * cosT));
    const double phd = twopi * G4UniformRand();
    // вниз: -z
    fGun.SetParticleMomentumDirection(
        G4ThreeVector(sinT * std::cos(phd), sinT * std::sin(phd), -cosT));

    fGun.SetParticleEnergy(SampleEnergyGeV() * GeV);
    fGun.GeneratePrimaryVertex(e);
  }
};

// --- скоринг: спектр энерговыделения в кристалле, тот же формат ---------
class RunAct : public G4UserRunAction {
public:
  static constexpr int kBins = 20000;   // до 20 МэВ, 1 кэВ/канал — мюон
                                        // может отдать несколько МэВ разом
  static constexpr double kBinKeV = 1.0;
  std::vector<long> fHist{std::vector<long>(kBins + 1, 0)};
  long fWithSignal = 0;
  G4String fOut = "cosmicmu.csv";

  void BeginOfRunAction(const G4Run*) override {
    std::fill(fHist.begin(), fHist.end(), 0L);
    fWithSignal = 0;
  }
  void Fill(double edepKeV) {
    if (edepKeV <= 0) return;
    ++fWithSignal;
    int b = static_cast<int>(edepKeV / kBinKeV);
    if (b > kBins) b = kBins;
    ++fHist[b];
  }
  void EndOfRunAction(const G4Run* run) override {
    const long N = run->GetNumberOfEvent();
    if (!N) return;
    FILE* f = std::fopen(fOut.c_str(), "w");
    std::fprintf(f, "# космический мюонный континуум в кристалле, ОТКЛИК НА 1 "
                   "МЮОН НА ДИСКЕ (не на абсолютный поток, см. шапку .cc)\n");
    std::fprintf(f, "# src_sha1 = %s\n", RCCM_SRC_SHA1);
    std::fprintf(f, "# git_describe = %s\n", RCCM_GIT_DESCRIBE);
    std::fprintf(f, "# build = %s %s\n", __DATE__, __TIME__);
    std::fprintf(f, "# R_disk_mm = %.1f  Z_disk_mm = %.1f  E_range_GeV = %.1f..%.1f\n",
                 R_DISK, Z_DISK, E_LO_GEV, E_HI_GEV);
    std::fprintf(f, "# N_primaries = %ld\n", N);
    std::fprintf(f, "# N_with_signal = %ld\n", fWithSignal);
    std::fprintf(f, "# bin_keV = %.3f  (последний канал = переполнение)\n", kBinKeV);
    std::fprintf(f, "E_keV,counts\n");
    for (int i = 0; i <= kBins; ++i)
      if (fHist[i]) std::fprintf(f, "%.1f,%ld\n", (i + 0.5) * kBinKeV, fHist[i]);
    std::fclose(f);
    G4cout << "RESULT cosmicmu N= " << N << " hits= " << fWithSignal
           << " eff_total= " << double(fWithSignal) / N << " file= " << fOut
           << G4endl;
  }
};

class EventAct : public G4UserEventAction {
  RunAct* fRun;
public:
  double fEdep = 0;
  explicit EventAct(RunAct* r) : fRun(r) {}
  void BeginOfEventAction(const G4Event*) override { fEdep = 0; }
  void EndOfEventAction(const G4Event*) override { fRun->Fill(fEdep / keV); }
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

}  // namespace

int main(int argc, char** argv) {
  BuildEnergyTable();
  const long n = (argc > 1) ? std::atol(argv[1]) : 2000000;

  // FTFP_BERT — не EmStandardPhysics: мюону нужны ионизация+тормозное+пары
  // (доминируют в депозите) и мюон-ядерные (второстепенно), FTFP_BERT
  // регистрирует всё разом; тот же физлист сверен по исходнику B01
  // (см. план, раздел про биасинг) как стандартный выбор для мюонов.
  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  rm->SetUserInitialization(new RCDetector(true));   // прибор+сосуд, без защиты
  rm->SetUserInitialization(new FTFP_BERT());
  rm->SetUserAction(new MuGun());

  auto* runAct = new RunAct();
  if (argc > 2) runAct->fOut = argv[2];
  rm->SetUserAction(runAct);
  auto* evtAct = new EventAct(runAct);
  rm->SetUserAction(evtAct);

  rm->Initialize();
  auto* det = dynamic_cast<RCDetector*>(
      const_cast<G4VUserDetectorConstruction*>(rm->GetUserDetectorConstruction()));
  rm->SetUserAction(new Stepping(evtAct, det->fCrystalLV));

  G4cout << "# src_sha1 = " << RCCM_SRC_SHA1 << G4endl;
  G4cout << "# git_describe = " << RCCM_GIT_DESCRIBE << G4endl;
  rm->BeamOn(n);

  delete rm;
  return 0;
}
