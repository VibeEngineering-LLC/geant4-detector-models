// Спектр поля ЕРН в помещении: флюенс фотонов в воздушной полости внутри
// бетона с равномерно распределёнными естественными радионуклидами.
//
// ПОЧЕМУ ТАК. Точка в середине помещения окружена ограждающими конструкциями на
// полный телесный угол. У поверхности полубесконечной среды флюенс равен
// половине флюенса в бесконечной среде (источник занимает 2pi); шесть
// поверхностей дают 4pi, то есть в сумме — значение для бесконечной среды.
// Поэтому поле считается как флюенс в полости внутри толстого бетона, а не
// собирается вручную из линий: так автоматически получается и нерассеянная
// компонента, и рассеянный континуум в правильном соотношении. Это существенно:
// в области 100..400 кэВ рассеянных фотонов больше, чем первичных, а именно там
// у CsI максимум эффективности.
//
// Скоринг: сумма длин пробегов фотонов в полости, поделённая на её объём, есть
// флюенс. Абсолютная нормировка — через объёмную скорость испускания.
//
// Запуск:  wallfield [число событий] [выходной файл]
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4Run.hh"
#include "G4RunManagerFactory.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "Randomize.hh"
#include "globals.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

const double R_WALL = 80. * cm;   // бетон: 7 длин пробега даже на 2.6 МэВ
const double R_CAV  = 20. * cm;   // воздушная полость («помещение»)

// Удельные активности обычных ограждающих конструкций, Бк/кг (UNSCEAR, средние
// по миру; реальный разброс K 100..1600, Ra 10..100, Th 10..60).
const double A_K40 = 400., A_RA226 = 40., A_TH232 = 30.;
const double RHO_CONCRETE = 2.30;  // г/см³, G4_CONCRETE

struct Line { double E_keV, yield; int series; };  // 0 = K, 1 = Ra, 2 = Th

// Выходы на распад родителя при равновесии цепочки (LNHB/DDEP, основные линии).
const Line LINES[] = {
    {1460.8, 0.1055, 0},
    // U-238 -> Ra-226 -> ... (Pb-214, Bi-214)
    {186.2, 0.0359, 1},  {241.9, 0.0727, 1},  {295.2, 0.1841, 1},
    {351.9, 0.3560, 1},  {609.3, 0.4549, 1},  {665.4, 0.0153, 1},
    {768.4, 0.0489, 1},  {806.2, 0.0126, 1},  {934.1, 0.0310, 1},
    {1120.3, 0.1491, 1}, {1155.2, 0.0164, 1}, {1238.1, 0.0583, 1},
    {1280.9, 0.0143, 1}, {1377.7, 0.0397, 1}, {1401.5, 0.0133, 1},
    {1408.0, 0.0239, 1}, {1509.2, 0.0213, 1}, {1661.3, 0.0105, 1},
    {1729.6, 0.0284, 1}, {1764.5, 0.1531, 1}, {1847.4, 0.0203, 1},
    {2118.6, 0.0116, 1}, {2204.2, 0.0491, 1}, {2447.9, 0.0155, 1},
    // Th-232 -> ... (Pb-212, Ac-228, Bi-212, Tl-208 с ветвлением 35.94 %)
    {238.6, 0.4360, 2},  {240.0, 0.0410, 2},  {270.2, 0.0346, 2},
    {300.1, 0.0328, 2},  {338.3, 0.1127, 2},  {463.0, 0.0440, 2},
    {510.7, 0.0810, 2},  {583.2, 0.3055, 2},  {727.3, 0.0667, 2},
    {772.3, 0.0155, 2},  {794.9, 0.0426, 2},  {835.7, 0.0161, 2},
    {860.6, 0.0450, 2},  {911.2, 0.2580, 2},  {964.8, 0.0499, 2},
    {968.9, 0.1580, 2},  {1588.2, 0.0327, 2}, {1620.5, 0.0149, 2},
    {1630.6, 0.0170, 2}, {2614.5, 0.3585, 2},
};
const int NLINES = sizeof(LINES) / sizeof(LINES[0]);

double gCum[NLINES];      // кумулятивные веса для розыгрыша линии
double gSvTotal = 0;      // объёмная скорость испускания, фотон/(см³·с)

void BuildWeights() {
  const double a[3] = {A_K40, A_RA226, A_TH232};
  double s = 0;
  for (int i = 0; i < NLINES; ++i) {
    s += a[LINES[i].series] * LINES[i].yield * RHO_CONCRETE / 1000.0;  // 1/(см³·с)
    gCum[i] = s;
  }
  gSvTotal = s;
  for (int i = 0; i < NLINES; ++i) gCum[i] /= s;
}

G4LogicalVolume* gCavLV = nullptr;

}  // namespace

// --- геометрия ---------------------------------------------------------------
class WallGeom : public G4VUserDetectorConstruction {
public:
  G4VPhysicalVolume* Construct() override {
    auto* nist = G4NistManager::Instance();
    auto* world = new G4LogicalVolume(
        new G4Box("world", 1.1 * R_WALL, 1.1 * R_WALL, 1.1 * R_WALL),
        nist->FindOrBuildMaterial("G4_AIR"), "world");
    auto* pv = new G4PVPlacement(nullptr, {}, world, "world", nullptr, false, 0, true);

    auto* wall = new G4LogicalVolume(new G4Orb("wall", R_WALL),
                                     nist->FindOrBuildMaterial("G4_CONCRETE"), "wall");
    new G4PVPlacement(nullptr, {}, wall, "wall", world, false, 0, true);

    gCavLV = new G4LogicalVolume(new G4Orb("cav", R_CAV),
                                 nist->FindOrBuildMaterial("G4_AIR"), "cav");
    new G4PVPlacement(nullptr, {}, gCavLV, "cav", wall, false, 0, true);
    return pv;
  }
};

// --- источник ----------------------------------------------------------------
class WallSource : public G4VUserPrimaryGeneratorAction {
  G4ParticleGun fGun{1};
public:
  WallSource() {
    fGun.SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("gamma"));
  }
  void GeneratePrimaries(G4Event* e) override {
    // равномерно по бетону: розыгрыш в шаре с отбрасыванием полости
    G4ThreeVector p;
    do {
      p.set((2 * G4UniformRand() - 1) * R_WALL, (2 * G4UniformRand() - 1) * R_WALL,
            (2 * G4UniformRand() - 1) * R_WALL);
    } while (p.mag() > R_WALL || p.mag() < R_CAV);
    fGun.SetParticlePosition(p);

    const double u = G4UniformRand();
    int i = 0;
    while (i < NLINES - 1 && u > gCum[i]) ++i;
    fGun.SetParticleEnergy(LINES[i].E_keV * keV);

    const double c = 2 * G4UniformRand() - 1, s = std::sqrt(1 - c * c);
    const double ph = twopi * G4UniformRand();
    fGun.SetParticleMomentumDirection({s * std::cos(ph), s * std::sin(ph), c});
    fGun.GeneratePrimaryVertex(e);
  }
};

// --- скоринг: длина пробега фотонов в полости по энергиям --------------------
class Fluence : public G4UserRunAction {
public:
  static constexpr int kBins = 300;      // 10 кэВ на канал до 3 МэВ
  static constexpr double kBinKeV = 10.0;
  std::vector<double> fLen{std::vector<double>(kBins + 1, 0.0)};
  G4String fOut = "wallfield.csv";

  void BeginOfRunAction(const G4Run*) override {
    std::fill(fLen.begin(), fLen.end(), 0.0);
  }
  void Add(double eKeV, double lenCm) {
    int b = static_cast<int>(eKeV / kBinKeV);
    if (b > kBins) b = kBins;
    if (b >= 0) fLen[b] += lenCm;
  }
  void EndOfRunAction(const G4Run* run) override {
    const long N = run->GetNumberOfEvent();
    if (!N) return;
    const double vWall = 4. / 3. * pi * (std::pow(R_WALL / cm, 3) - std::pow(R_CAV / cm, 3));
    const double vCav = 4. / 3. * pi * std::pow(R_CAV / cm, 3);
    // Реальная скорость испускания в объёме бетона, фотон/с
    const double rate = gSvTotal * vWall;
    // Флюенс = (средняя длина пробега на фотон) * (скорость) / объём полости
    const double norm = rate / (double(N) * vCav);

    double tot = 0;
    for (int i = 0; i <= kBins; ++i) tot += fLen[i] * norm;

    FILE* f = std::fopen(fOut.c_str(), "w");
    std::fprintf(f, "# поле ЕРН в помещении: флюенс в воздушной полости\n");
    std::fprintf(f, "# concrete: K-40 %.0f, Ra-226 %.0f, Th-232 %.0f Bq/kg, rho %.2f\n",
                 A_K40, A_RA226, A_TH232, RHO_CONCRETE);
    std::fprintf(f, "# R_wall_cm = %.1f  R_cav_cm = %.1f\n", R_WALL / cm, R_CAV / cm);
    std::fprintf(f, "# N = %ld\n", N);
    std::fprintf(f, "# Sv_total_per_cm3_s = %.6e\n", gSvTotal);
    std::fprintf(f, "# fluence_total_cm2_s = %.6e\n", tot);
    std::fprintf(f, "E_keV,fluence_cm2_s\n");
    for (int i = 0; i <= kBins; ++i)
      if (fLen[i] > 0)
        std::fprintf(f, "%.1f,%.6e\n", (i + 0.5) * kBinKeV, fLen[i] * norm);
    std::fclose(f);
    G4cout << "RESULT N= " << N << " fluence_total= " << tot << " cm-2 s-1  file= "
           << fOut << G4endl;
  }
};

class Track : public G4UserSteppingAction {
  Fluence* fRun;
public:
  explicit Track(Fluence* r) : fRun(r) {}
  void UserSteppingAction(const G4Step* s) override {
    if (s->GetTrack()->GetDefinition()->GetPDGEncoding() != 22) return;
    auto* v = s->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    if (!v || v->GetLogicalVolume() != gCavLV) return;
    fRun->Add(s->GetPreStepPoint()->GetKineticEnergy() / keV, s->GetStepLength() / cm);
  }
};

class Phys : public G4VModularPhysicsList {
public:
  Phys() {
    RegisterPhysics(new G4EmStandardPhysics_option4());
    SetDefaultCutValue(1.0 * mm);   // электроны здесь не нужны, важен транспорт гамма
  }
};

int main(int argc, char** argv) {
  BuildWeights();
  const long n = (argc > 1) ? std::atol(argv[1]) : 2000000;

  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  rm->SetUserInitialization(new WallGeom());
  rm->SetUserInitialization(new Phys());
  rm->SetUserAction(new WallSource());
  auto* fl = new Fluence();
  if (argc > 2) fl->fOut = argv[2];
  rm->SetUserAction(fl);
  rm->Initialize();
  rm->SetUserAction(new Track(fl));

  auto* ui = G4UImanager::GetUIpointer();
  ui->ApplyCommand("/run/verbose 0");
  ui->ApplyCommand("/event/verbose 0");
  ui->ApplyCommand("/tracking/verbose 0");
  ui->ApplyCommand("/run/printProgress 0");
  ui->ApplyCommand("/run/beamOn " + std::to_string(n));
  delete rm;
  return 0;
}
