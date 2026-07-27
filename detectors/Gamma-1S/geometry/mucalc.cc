// Коэффициент ослабления матрицы ОИСН-16 — из той же физики, что и основной
// расчёт (EmStandardPhysics_option4), а не из посторонних таблиц.
//
// Нужен, чтобы свести ход эффективности по плотности к физической форме
//     f = (1 - exp(-mu*ro*d)) / (mu*ro*d)
// с одним геометрическим параметром d — эффективной толщиной пробы. У ЛСРМ
// для Маринелли 1 л есть два независимых значения для сверки: табличное
// 26(2) мм («Прецизионные измерения», с. 11) и Thick = 31(2) мм в .efa
// этого экземпляра детектора.
//
// Без Rayleigh: когерентное рассеяние почти не уводит фотон из пучка, и
// именно так считался mu в контуре radiacode-curves — сохраняем сопоставимость.
#include "G1SDetector.hh"

#include "G4Box.hh"
#include "G4EmCalculator.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include <cstdio>

namespace {

const double E_GRID[] = {59.5, 88.0, 122.1, 165.9, 238.632, 241.995, 295.223,
                         338.32, 351.932, 463.004, 583.187, 609.32, 661.657,
                         768.36, 911.204, 1120.294, 1460.822, 1764.491,
                         2614.511, 3000.0};
const int NE = sizeof(E_GRID) / sizeof(E_GRID[0]);

G4Material* gOisn = nullptr;
G4Material* gWater = nullptr;   // для лёгких матриц источников комплекта

class Geom : public G4VUserDetectorConstruction {
public:
  G4VPhysicalVolume* Construct() override {
    gOisn = G1SDetector::MakeMatrix("OISN16", 1.0, "OISN16_unit");
    gWater = G1SDetector::MakeMatrix("water", 1.0, "Water_unit");
    auto* lv = new G4LogicalVolume(new G4Box("w", 1 * m, 1 * m, 1 * m),
                                   gOisn, "w");
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
  rm->BeamOn(1);                      // инициализировать таблицы сечений

  G4EmCalculator calc;
  const char* PROC[] = {"compt", "phot", "conv"};
  struct Out { G4Material* m; const char* fn; const char* title; };
  const Out OUTS[] = {
      {gOisn, "mu_oisn16.csv", "ОИСН-16"},
      {gWater, "mu_water.csv", "вода"},
  };
  for (const Out& o : OUTS) {
    FILE* f = std::fopen(o.fn, "w");
    std::fprintf(f, "# %s, массовый коэффициент ослабления, см²/г\n", o.title);
    std::fprintf(f, "# compt+phot+conv, EmStandardPhysics_option4, Geant4 11.2.1\n");
    std::fprintf(f, "E_keV,mu_over_rho_cm2_g\n");
    for (int j = 0; j < NE; ++j) {
      double mu = 0;                  // 1/мм при ро = 1 г/см³
      for (const char* p : PROC)
        mu += calc.ComputeCrossSectionPerVolume(E_GRID[j] * keV, "gamma", p,
                                                o.m->GetName());
      std::fprintf(f, "%.3f,%.6e\n", E_GRID[j], mu * 10.0);
    }
    std::fclose(f);
    G4cout << "RESULT записано " << NE << " строк в " << o.fn << G4endl;
  }
  delete rm;
  return 0;
}
