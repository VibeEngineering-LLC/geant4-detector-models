#include "Rc103RoomPrimaryGeneratorAction.hh"
#include "Rc103RoomDetectorConstruction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <cmath>

namespace {

// --- Таблица линий ЕРН + активности UNSCEAR ---------------------------------
// ПЕРЕИСПОЛЬЗОВАНО КАК ДАННЫЕ (§33 доктрины — не изобретать заново готовое)
// из detectors/RadiaCode-103/geometry/wallfield.cc:101,170-209 (ENSDF/LNHB,
// выгрузка 09.08.2026; см. комментарии происхождения в wallfield.cc —
// НЕ повторены здесь дословно, чтобы не дублировать источник дважды).
// Переписано в этот файл руками (не include), т.к. wallfield.cc — другая,
// запрещённая как код-образец физическая линия. Только режим "all" — три
// серии сразу со стандартными активностями UNSCEAR, без CLI-переключателя
// серии/нуклида/иона (эта программа не воспроизводит ВСЮ гибкость
// wallfield.cc, только тот срез, что нужен задаче: единый фон помещения).
const double A_K40 = 400., A_RA226 = 40., A_TH232 = 30.;   // Бк/кг
const double RHO_CONCRETE = 2.30;                          // г/см3, G4_CONCRETE

struct Line { double E_keV, yield; int series; };  // 0=K-40, 1=Ra-226 ряд, 2=Th-232 ряд

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
    // Слабые линии Bi-214 выше 2447.9 кэВ (>=0.002% выход) — источник
    // wallfield.cc:181-200 (IAEA NDS/ENSDF, 09.08.2026).
    {2472.9, 2.30e-5, 1},  {2505.5, 5.60e-5, 1},  {2694.7, 3.00e-4, 1},
    {2699.2, 2.77e-5, 1},  {2769.9, 2.45e-4, 1},  {2785.9, 5.50e-5, 1},
    {2827.0, 2.40e-5, 1},  {2880.4, 1.00e-4, 1},  {2893.6, 5.90e-5, 1},
    {2922.0, 1.36e-4, 1},  {2940.0, 3.60e-5, 1},  {2978.9, 1.36e-4, 1},
    {3000.0, 8.60e-5, 1},  {3053.9, 2.09e-4, 1},  {3081.8, 5.90e-5, 1},
    // Th-232 -> ... (Pb-212, Ac-228, Bi-212, Tl-208, ветвление 35.94%)
    {238.6, 0.4360, 2},  {240.0, 0.0410, 2},  {270.2, 0.0346, 2},
    {300.1, 0.0328, 2},  {338.3, 0.1127, 2},  {463.0, 0.0440, 2},
    {510.7, 0.0810, 2},  {583.2, 0.3055, 2},  {727.3, 0.0667, 2},
    {772.3, 0.0155, 2},  {794.9, 0.0426, 2},  {835.7, 0.0161, 2},
    {860.6, 0.0450, 2},  {911.2, 0.2580, 2},  {964.8, 0.0499, 2},
    {968.9, 0.1580, 2},  {1588.2, 0.0327, 2}, {1620.5, 0.0149, 2},
    {1630.6, 0.0170, 2}, {2614.5, 0.3585, 2},
};
const int NLINES = sizeof(LINES) / sizeof(LINES[0]);

double gCum[NLINES];
bool gWeightsBuilt = false;

void BuildWeights(double* outRatePerCm3S) {
  const double a[3] = {A_K40, A_RA226, A_TH232};
  double s = 0;
  for (int i = 0; i < NLINES; ++i) {
    s += a[LINES[i].series] * LINES[i].yield * RHO_CONCRETE / 1000.0;  // 1/(cm3 s)
    gCum[i] = s;
  }
  *outRatePerCm3S = s;
  for (int i = 0; i < NLINES; ++i) gCum[i] /= s;
  gWeightsBuilt = true;
}

}  // namespace

double Rc103RoomPrimaryGeneratorAction::sWallEmissionRatePerSec = 0.0;

Rc103RoomPrimaryGeneratorAction::Rc103RoomPrimaryGeneratorAction()
    : fGun(1) {
  fGun.SetParticleDefinition(
      G4ParticleTable::GetParticleTable()->FindParticle("gamma"));

  if (!gWeightsBuilt) {
    double ratePerCm3S = 0.0;
    BuildWeights(&ratePerCm3S);
    const double rWall = Rc103RoomDetectorConstruction::kRWallCm;  // cm
    const double rCav = Rc103RoomDetectorConstruction::kRCavCm;    // cm
    const double vWallCm3 =
        4.0 / 3.0 * pi * (std::pow(rWall, 3) - std::pow(rCav, 3));
    sWallEmissionRatePerSec = ratePerCm3S * vWallCm3;  // фотон/с во всей стене
  }
}

void Rc103RoomPrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  const double rWall = Rc103RoomDetectorConstruction::kRWallCm * cm;
  const double rCav = Rc103RoomDetectorConstruction::kRCavCm * cm;

  // Равномерно по бетону: розыгрыш в шаре с отбраковкой полости — переиспользовано
  // из wallfield.cc::WallSource (некорнерная ветка).
  G4ThreeVector p;
  do {
    p.set((2 * G4UniformRand() - 1) * rWall, (2 * G4UniformRand() - 1) * rWall,
          (2 * G4UniformRand() - 1) * rWall);
  } while (p.mag() > rWall || p.mag() < rCav);
  fGun.SetParticlePosition(p);

  const double u = G4UniformRand();
  int i = 0;
  while (i < NLINES - 1 && u > gCum[i]) ++i;
  fGun.SetParticleEnergy(LINES[i].E_keV * keV);

  // Изотропное направление вылета — истинная физика точки распада, не "на
  // детектор" (иначе поле перестало бы автоматически давать правильное
  // соотношение прямой/рассеянной компоненты, см. обоснование в wallfield.cc).
  const double c = 2 * G4UniformRand() - 1, s = std::sqrt(1 - c * c);
  const double ph = twopi * G4UniformRand();
  fGun.SetParticleMomentumDirection({s * std::cos(ph), s * std::sin(ph), c});

  fGun.GeneratePrimaryVertex(event);
}
