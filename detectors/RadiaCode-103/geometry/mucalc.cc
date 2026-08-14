// Коэффициенты ослабления матриц пробы — из той же физики, что и основной
// расчёт (EmStandardPhysics_option4), а не из посторонних таблиц.
//
// Нужны, чтобы свести поправку на самопоглощение к физической форме
//     f = (1 - exp(-mu*d)) / (mu*d)
// с ОДНИМ геометрическим параметром d (эффективная толщина пробы). Если один d
// описывает все матрицы и плотности сразу, поправка перестаёт быть таблицей и
// становится формулой, а заодно снимается шум точек на высоких энергиях, где
// эффект в 10..15 % тонет в статистике прогонов.
#include "RCDetector.hh"

// Собственный штамп провенанса (mucalc.cc + RCDetector.cc/.hh) — не общий с
// main.cc, которого в этом бинарнике нет.
#if defined(__has_include)
#  if __has_include("rc_mucalc_provenance.hh")
#    include "rc_mucalc_provenance.hh"
#  endif
#endif
#ifndef RCMU_SRC_SHA1
#  define RCMU_SRC_SHA1 "БЕЗ-ШТАМПА"
#  define RCMU_GIT_DESCRIBE "БЕЗ-ШТАМПА"
#endif

#include "G4Box.hh"
#include "G4EmCalculator.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Event.hh"

#include <cstdio>
#include <string>
#include <vector>

namespace {

struct Cfg { const char* matrix; double rho; };
const Cfg CFGS[] = {
    {"air", 0.0012}, {"organic", 0.50}, {"soil", 0.80},
    {"water", 1.00}, {"soil", 1.20}, {"soil", 1.60},
};
const int NCFG = sizeof(CFGS) / sizeof(CFGS[0]);

const double E_GRID[] = {30, 46.5, 59.5, 80, 92.6, 122, 186, 238.6, 295.2,
                         351.9, 477, 583.2, 609.3, 661.7, 795, 911.2, 968.9,
                         1120.3, 1173.2, 1332.5, 1460.8, 1764.5, 2204, 2614.5,
                         3000};
const int NE = sizeof(E_GRID) / sizeof(E_GRID[0]);

// --- материалы защиты -------------------------------------------------------
// Свинец, кадмий и медь нужны для расчёта свинцовой защиты: длина свободного
// пробега lambda(E) = 1/mu задаёт карту важностей при расщеплении в толще, а
// узкопучковое пропускание exp(-mu*t) служит проверкой геометрии.
// Считаются той же физикой, что и матрицы пробы, — это принцип файла.
const char* SHIELD_MATS[] = {"G4_Pb", "G4_Cd", "G4_Cu"};
const int NSH = sizeof(SHIELD_MATS) / sizeof(SHIELD_MATS[0]);

// СВОЯ сетка, гуще матричной внизу. Причина: у этих материалов в рабочем
// диапазоне лежат K-края, на которых mu разрывен, и одной точкой «в крае»
// коэффициент не описывается. Края взяты в вилку:
//   Pb  K 88,0 кэВ  -> 87,0 и 88,5
//   Cd  K 26,7 кэВ  -> 26,5 и 27,0
// K-край меди 8,98 кэВ лежит НИЖЕ начала сетки: на таких энергиях фотон не
// выходит ни из пробы, ни из корпуса прибора, и для задачи защиты не нужен.
// Верх сетки — 2614,5 кэВ (Tl-208), самая жёсткая линия поля ЕРН, и 3000.
const double E_SHIELD[] = {
    10, 12, 15, 20, 23, 25, 26.5, 27, 30, 40, 50, 60, 70, 80, 87, 88.5,
    92.6, 100, 122, 150, 186, 238.6, 295.2, 351.9, 477, 583.2, 609.3,
    661.7, 795, 911.2, 968.9, 1120.3, 1173.2, 1332.5, 1460.8, 1764.5,
    2204, 2614.5, 3000};
const int NESH = sizeof(E_SHIELD) / sizeof(E_SHIELD[0]);

std::vector<G4Material*> gMats;
std::vector<G4Material*> gShield;

class Geom : public G4VUserDetectorConstruction {
public:
  G4VPhysicalVolume* Construct() override {
    for (int i = 0; i < NCFG; ++i) {
      char nm[64];
      std::snprintf(nm, sizeof(nm), "M_%s_%.2f", CFGS[i].matrix, CFGS[i].rho);
      gMats.push_back(RCDetector::MakeMatrix(CFGS[i].matrix, CFGS[i].rho, nm));
    }
    auto* nist = G4NistManager::Instance();
    for (int i = 0; i < NSH; ++i)
      gShield.push_back(nist->FindOrBuildMaterial(SHIELD_MATS[i]));
    auto* lv = new G4LogicalVolume(new G4Box("w", 1 * m, 1 * m, 1 * m),
                                   gMats[0], "w");
    return new G4PVPlacement(nullptr, {}, lv, "w", nullptr, false, 0);
  }
};

class Gun : public G4VUserPrimaryGeneratorAction {
  G4ParticleGun fGun{1};
public:
  Gun() {
    fGun.SetParticleDefinition(
        G4ParticleTable::GetParticleTable()->FindParticle("gamma"));
    fGun.SetParticleEnergy(1 * MeV);
  }
  void GeneratePrimaries(G4Event* e) override { fGun.GeneratePrimaryVertex(e); }
};

class Phys : public G4VModularPhysicsList {
public:
  Phys() { RegisterPhysics(new G4EmStandardPhysics_option4()); }
};

}  // namespace

int main(int argc, char** argv) {
  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  rm->SetUserInitialization(new Geom());
  rm->SetUserInitialization(new Phys());
  rm->SetUserAction(new Gun());
  rm->Initialize();
  rm->BeamOn(1);          // инициализировать таблицы сечений

  G4EmCalculator calc;
  const char* PROC[] = {"compt", "phot", "conv"};   // без Rayleigh: он почти не
                                                    // уводит фотон из пучка
  const char* fname = (argc > 1) ? argv[1] : "mu.csv";
  FILE* f = std::fopen(fname, "w");
  std::fprintf(f, "# массовый коэффициент ослабления, см²/г (compt+phot+conv)\n");
  std::fprintf(f, "# physics: EmStandardPhysics_option4, Geant4 11.2.1\n");
  std::fprintf(f, "# src_sha1 = %s\n", RCMU_SRC_SHA1);
  std::fprintf(f, "# git_describe = %s\n", RCMU_GIT_DESCRIBE);
  std::fprintf(f, "# build = %s %s\n", __DATE__, __TIME__);
  std::fprintf(f, "matrix,rho_gcm3,E_keV,mu_over_rho_cm2_g,mu_cm-1\n");
  for (int i = 0; i < NCFG; ++i) {
    auto* m = gMats[i];
    const double rho = m->GetDensity() / (g / cm3);
    for (int j = 0; j < NE; ++j) {
      const double E = E_GRID[j] * keV;
      double mu = 0;   // 1/мм
      for (const char* p : PROC)
        mu += calc.ComputeCrossSectionPerVolume(E, "gamma", p, m->GetName());
      const double mu_cm = mu * 10.0;            // 1/см
      std::fprintf(f, "%s,%.4f,%.1f,%.6e,%.6e\n", CFGS[i].matrix, CFGS[i].rho,
                   E_GRID[j], mu_cm / rho, mu_cm);
    }
  }
  // материалы защиты — своя, более густая сетка (см. E_SHIELD)
  for (int i = 0; i < NSH; ++i) {
    auto* m = gShield[i];
    const double rho = m->GetDensity() / (g / cm3);
    for (int j = 0; j < NESH; ++j) {
      const double E = E_SHIELD[j] * keV;
      double mu = 0;   // 1/мм
      for (const char* p : PROC)
        mu += calc.ComputeCrossSectionPerVolume(E, "gamma", p, m->GetName());
      const double mu_cm = mu * 10.0;            // 1/см
      std::fprintf(f, "%s,%.4f,%.1f,%.6e,%.6e\n", SHIELD_MATS[i], rho,
                   E_SHIELD[j], mu_cm / rho, mu_cm);
    }
    // Вторая строка по тем же материалам — С КОГЕРЕНТНЫМ (Rayleigh). Это
    // РАЗНОЕ определение mu для РАЗНОЙ задачи: без Rayleigh годится для
    // самопоглощения в пробе (когерентное рассеяние на малый угол почти не
    // уводит квант от детектора, см. комментарий выше), но для узкого
    // ("good geometry") пучка сквозь толщу защиты это неверно — там любое
    // отклонение выводит квант из строгого определения пучка. Проверено
    // прямым МК-прогоном 20 мм Pb/661.7 кэВ (shieldrun beam): без Rayleigh
    // exp(-mu*t) предсказывает пропускание 0.0956, реально получено 0.0828
    // (−13.4 %) — расхождение снимается добавлением когерентного сечения.
    for (int j = 0; j < NESH; ++j) {
      const double E = E_SHIELD[j] * keV;
      double mu = calc.ComputeCrossSectionPerVolume(E, "gamma", "Rayl", m->GetName());
      for (const char* p : PROC)
        mu += calc.ComputeCrossSectionPerVolume(E, "gamma", p, m->GetName());
      const double mu_cm = mu * 10.0;
      std::fprintf(f, "%s+Rayl,%.4f,%.1f,%.6e,%.6e\n", SHIELD_MATS[i], rho,
                   E_SHIELD[j], mu_cm / rho, mu_cm);
    }
  }
  std::fclose(f);
  G4cout << "RESULT записано " << NCFG * NE + NSH * NESH << " строк в " << fname
         << G4endl;
  delete rm;
  return 0;
}
