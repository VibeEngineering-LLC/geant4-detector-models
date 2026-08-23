// Мягкая (электронно-фотонная) компонента космического излучения в кристалле
// RC-103 на ПЕРВОМ ЭТАЖЕ двухэтажного дома: источник над ДВУМЯ бетонными
// перекрытиями, транспорт ливня через них честной физикой Geant4.
//
// ЗАЧЕМ. Открытый фон не сходится в 1500-2400 кэВ (м/и 0,58-0,68) после того,
// как мюонная ветка закрыта абсолютно (a_mu/PDG=1,25, диск R=500). Дыра того
// же порядка, что ожидаемый вклад мягкой компоненты после 2 перекрытий
// (~0,002-0,02 cps) — единственный оставшийся кандидат (#SHIELD-28).
//
// ВХОДНОЙ СПЕКТР (PDG RPP2019 §29.3.2, прочитан прямо 20.08.2026, цитата:
// "integral vertical intensity of electrons plus positrons is very
// approximately 30, 6, and 0.2 m-2 s-1 sr-1 above 10, 100, and 1000 MeV";
// "ratio of photons to electrons ... approximately 1.3 above 1 GeV and 1.7
// below the critical energy"):
//   дифференциально dN/dE ~ E^-1.70 (10..100 МэВ) с изломом на E^-2.48
//   (100..1000 МэВ) — показатели из трёх интегральных точек PDG.
//   Фотоны = 1,7 x электроны, форма спектра принята той же (ливневое
//   равновесие; ДОПУЩЕНИЕ). Розыгрыш 2..3000 МэВ лог-таблицей.
// ДОПУЩЕНИЯ, ЗАЯВЛЕННЫЕ ЯВНО:
//   - угловое распределение cos^2(T) как у мюонов (PDG: "angular dependence
//     is complex" — упрощение);
//   - продление спектра вниз до 2 МэВ тем же степенным законом (мягче 10 МэВ
//     PDG численно не даёт; частицы <2 МэВ перекрытия не проходят, а
//     равновесный мягкий хвост у прибора рождается в самих плитах);
//   - перекрытия: 2 плиты сплошного G4_CONCRETE по 0,20 м, 8x8 м, на 2,5 и
//     5,3 м над прибором (типовые ЖБ 200-220 мм; пустотность НЕ учтена);
//   - крыша/кровля НЕ учтена (лёгкая — мало вещества против 2 плит).
// НОРМИРОВКА: отклик НА ОДИН ПЕРВИЧНЫЙ с диска; абсолютная скорость — в
// анализе: I_v(e,>10 МэВ)=30 м-2с-1ср-1, поток через горизонт при cos^2 =
// pi/2*I_v; доля розыгрыша >10 МэВ (F10) печатается в шапку CSV.
//
// Запуск: cosmicem <N> <out.csv> [rdisk=<мм>] [seed=<n>]
#define RCCM_SRC_SHA1 "cosmicem-20260820"
#define RCCM_GIT_DESCRIBE "cosmicem-20260820"
#include "RCDetector.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"

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
#include <string>
#include <vector>

namespace {

// --- спектр мягкой компоненты: степенной с изломом (PDG, см. шапку) ------
double SoftDiff(double E_MeV) {
  // непрерывен в точке излома 100 МэВ
  if (E_MeV < 100.0) return std::pow(E_MeV / 100.0, -1.70);
  return std::pow(E_MeV / 100.0, -2.48);
}

const double E_LO_MEV = 2.0, E_HI_MEV = 3000.0;
const int N_EGRID = 400;
std::vector<double> gEGrid, gCumE;
double gFracAbove10 = 0.0;   // доля розыгрыша с E>10 МэВ (для нормировки)

void BuildEnergyTable() {
  gEGrid.resize(N_EGRID);
  gCumE.resize(N_EGRID);
  const double logLo = std::log(E_LO_MEV), logHi = std::log(E_HI_MEV);
  for (int i = 0; i < N_EGRID; ++i)
    gEGrid[i] = std::exp(logLo + (logHi - logLo) * i / (N_EGRID - 1));
  double s = 0, s10 = 0;
  for (int i = 0; i < N_EGRID; ++i) {
    const double w = SoftDiff(gEGrid[i]) * gEGrid[i];  // якобиан d(lnE)
    s += w;
    if (gEGrid[i] > 10.0) s10 += w;
    gCumE[i] = s;
  }
  gFracAbove10 = s10 / s;
  for (int i = 0; i < N_EGRID; ++i) gCumE[i] /= s;
}

double SampleEnergyMeV() {
  const double u = G4UniformRand();
  int lo = 0, hi = N_EGRID - 1;
  while (lo < hi) {
    const int mid = (lo + hi) / 2;
    if (gCumE[mid] < u) lo = mid + 1; else hi = mid;
  }
  return gEGrid[lo];
}

// --- геометрия: прибор + ДВА бетонных перекрытия --------------------------
const double SLAB_T = 200.0;       // мм, толщина плиты
const double SLAB_HALF_XY = 4000.0; // мм, плита 8x8 м
const double SLAB1_Z = 400.0;      // мм, низ 1-го перекрытия (воздух сжат: X0=300 м, радиографически ничтожен; реальная высота 2,5 м)
const double SLAB2_Z = 700.0;      // мм, низ 2-го (реально 5,3 м; сжато ради углового покрытия диска)
double R_DISK = 1200.0;            // мм, диск источника (rdisk= переопределяет)
const double Z_DISK_EM = 950.0;    // мм, над верхней плитой; R=1200 -> углы до 51 гр (94 % потока cos^3)

class EmGeom : public G4VUserDetectorConstruction {
public:
  RCDetector* fDet = nullptr;
  G4VPhysicalVolume* Construct() override {
    fDet = new RCDetector(true);   // прибор+сосуд, без защиты
    fDet->fWorldHalfXY = 4500.0;   // мм
    fDet->fWorldHalfZ = 1500.0;
    auto* pv = fDet->Construct();
    auto* worldLV = pv->GetLogicalVolume();
    auto* nist = G4NistManager::Instance();
    auto* concrete = nist->FindOrBuildMaterial("G4_CONCRETE");
    for (double z0 : {SLAB1_Z, SLAB2_Z}) {
      auto* s = new G4Box("slab", SLAB_HALF_XY * mm, SLAB_HALF_XY * mm,
                          0.5 * SLAB_T * mm);
      auto* lv = new G4LogicalVolume(s, concrete, "slab");
      new G4PVPlacement(nullptr, G4ThreeVector(0, 0, (z0 + 0.5 * SLAB_T) * mm),
                        lv, "slab", worldLV, false, 0, true);
    }
    return pv;
  }
};

// --- источник: диск над верхней плитой, e-/gamma 1:1.7 ---------------------
class EmGun : public G4VUserPrimaryGeneratorAction {
  G4ParticleGun fGun{1};
  G4ParticleDefinition* fE;
  G4ParticleDefinition* fG;
public:
  EmGun() {
    fE = G4ParticleTable::GetParticleTable()->FindParticle("e-");
    fG = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
  }
  void GeneratePrimaries(G4Event* e) override {
    // вид частицы: фотонов в 1,7 раза больше электронов (PDG, ниже E_c)
    fGun.SetParticleDefinition(G4UniformRand() < 1.7 / 2.7 ? fG : fE);

    const double r = R_DISK * std::sqrt(G4UniformRand());
    const double ph = twopi * G4UniformRand();
    fGun.SetParticlePosition(
        G4ThreeVector(r * std::cos(ph), r * std::sin(ph), Z_DISK_EM * mm));

    // поток через горизонтальную площадку при I(T)~cos^2: p(cosT)~cos^3
    const double cosT = std::pow(G4UniformRand(), 0.25);
    const double sinT = std::sqrt(std::max(0.0, 1.0 - cosT * cosT));
    const double phd = twopi * G4UniformRand();
    fGun.SetParticleMomentumDirection(
        G4ThreeVector(sinT * std::cos(phd), sinT * std::sin(phd), -cosT));

    fGun.SetParticleEnergy(SampleEnergyMeV() * MeV);
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
  G4String fOut = "cosmicem.csv";

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
    std::fprintf(f, "# мягкая (e/gamma) косм. компонента через 2 бетонных "
                   "перекрытия, ОТКЛИК НА 1 ПЕРВИЧНЫЙ НА ДИСКЕ\n");
    std::fprintf(f, "# src_sha1 = %s\n", RCCM_SRC_SHA1);
    std::fprintf(f, "# git_describe = %s\n", RCCM_GIT_DESCRIBE);
    std::fprintf(f, "# build = %s %s\n", __DATE__, __TIME__);
    std::fprintf(f, "# R_disk_mm = %.1f  Z_disk_mm = %.1f  E_range_MeV = %.1f..%.1f\n",
                 R_DISK, Z_DISK_EM, E_LO_MEV, E_HI_MEV);
    std::fprintf(f, "# frac_above_10MeV = %.5f (I_v e >10MeV = 30 m-2s-1sr-1, flux pi/2*I_v, photons x1.7)\n", gFracAbove10);
    std::fprintf(f, "# slabs: 2 x %.0f mm G4_CONCRETE at z=%.0f,%.0f mm\n",
                 SLAB_T, SLAB1_Z, SLAB2_Z);
    std::fprintf(f, "# N_primaries = %ld\n", N);
    std::fprintf(f, "# N_with_signal = %ld\n", fWithSignal);
    std::fprintf(f, "# bin_keV = %.3f  (последний канал = переполнение)\n", kBinKeV);
    std::fprintf(f, "E_keV,counts\n");
    for (int i = 0; i <= kBins; ++i)
      if (fHist[i]) std::fprintf(f, "%.1f,%ld\n", (i + 0.5) * kBinKeV, fHist[i]);
    std::fclose(f);
    G4cout << "RESULT cosmicem N= " << N << " hits= " << fWithSignal
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
  // rdisk=<мм> — радиус диска источника (см. комментарий у R_DISK);
  // seed=<n> — зерно ГСЧ: без него ПАРАЛЛЕЛЬНЫЕ прогоны идентичны побайтно
  // (поймано 20.08: 8 процессов дали по 14141 попаданий — один поток чисел).
  for (int i = 3; i < argc; ++i) {
    const std::string a = argv[i];
    if (a.rfind("rdisk=", 0) == 0) R_DISK = std::atof(a.c_str() + 6);
    else if (a.rfind("seed=", 0) == 0) G4Random::setTheSeed(std::atol(a.c_str() + 5));
  }

  // FTFP_BERT — не EmStandardPhysics: мюону нужны ионизация+тормозное+пары
  // (доминируют в депозите) и мюон-ядерные (второстепенно), FTFP_BERT
  // регистрирует всё разом; тот же физлист сверен по исходнику B01
  // (см. план, раздел про биасинг) как стандартный выбор для мюонов.
  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  rm->SetUserInitialization(new EmGeom());
  rm->SetUserInitialization(new FTFP_BERT());
  rm->SetUserAction(new EmGun());

  auto* runAct = new RunAct();
  if (argc > 2) runAct->fOut = argv[2];
  rm->SetUserAction(runAct);
  auto* evtAct = new EventAct(runAct);
  rm->SetUserAction(evtAct);

  rm->Initialize();
  auto* eg = dynamic_cast<EmGeom*>(
      const_cast<G4VUserDetectorConstruction*>(rm->GetUserDetectorConstruction()));
  rm->SetUserAction(new Stepping(evtAct, eg->fDet->fCrystalLV));

  G4cout << "# src_sha1 = " << RCCM_SRC_SHA1 << G4endl;
  G4cout << "# git_describe = " << RCCM_GIT_DESCRIBE << G4endl;
  rm->BeamOn(n);

  delete rm;
  return 0;
}
