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
//
// Собственный штамп провенанса: нет ни main.cc, ни RCDetector.cc, поэтому
// отпечаток свой (rc_wallfield_provenance.hh), не общий с rc_curves/mucalc.
#if defined(__has_include)
#  if __has_include("rc_wallfield_provenance.hh")
#    include "rc_wallfield_provenance.hh"
#  endif
#endif
#ifndef RCWF_SRC_SHA1
#  define RCWF_SRC_SHA1 "БЕЗ-ШТАМПА"
#  define RCWF_GIT_DESCRIBE "БЕЗ-ШТАМПА"
#endif

#include "G4Box.hh"
#include "G4Event.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"                 // ионный режим D-001
#include "G4DecayPhysics.hh"             // ионный режим D-001
#include "G4RadioactiveDecayPhysics.hh"  // ионный режим D-001
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

// #SHIELD-16 (16.08.2026): не const — перекрывается ключом rwall= для проверки
// насыщения buildup, той же методикой, что уже ловила занижение 19% на
// мюонном диске (#SHIELD-8, было R=700мм). "7 длин пробега" — расчёт на
// бумаге, не проверено эмпирически до сих пор.
double R_WALL = 80. * cm;   // бетон: 7 длин пробега даже на 2.6 МэВ
const double R_CAV  = 20. * cm;   // воздушная полость («помещение»)

// #SHIELD-26 (20.08.2026): асимметричная («угловая») геометрия — домик при
// съёмке фона 23 дня стоял В УГЛУ комнаты, у пересечения двух кирпичных стен
// (факт от оператора, план 1-го этажа). Старый сферический режим (R_WALL,
// выше) НЕ трогается — он уже откалиброван по открытому фону (0,88…1,09) и
// остаётся дефолтом. Режим corner= включает отдельную несимметричную
// коробчатую геометрию ТОЛЬКО для сравнения домик/открытое место.
//
// Расстояния по осям X (лево/право от домика) и Y (низ/верх) — из ответов
// оператора: слева 1 м до стены в полкирпича (~12 см), за ней ЕЩЁ ДОМ (не
// улица) — approx достроено на ЕЩЁ 80 см того же материала (ДОПУЩЕНИЕ, не
// факт — за неимением точных данных о втором помещении); снизу 1,5 м до
// стены в 2 кирпича (~50 см), за ней лоджия и улица — материал ОБРЫВАЕТСЯ,
// добавочного источника нет. Направления без данных (+X открытая комната,
// +Y вглубь комнаты, вся ось Z пол/потолок) оставлены как в старой сфере:
// начинаются сразу от R_CAV, толщина до 80 см суммарно — совпадает со старым
// R_WALL=80 см для контрольного сравнения (regression test, см. --corner-test).
bool gCorner = false;
struct AxisSide { double innerCm, wallCm, extraCm; };  // до стены | стена | за стеной (0=обрыв)
// сторона "-" каждой оси — ближняя стена/факт; "+" — открытое (допущение = старая модель)
AxisSide gX_neg{100., 12., 80.};   AxisSide gX_pos{ 20., 60., 0.};
AxisSide gY_neg{150., 50.,  0.};   AxisSide gY_pos{ 20., 60., 0.};
AxisSide gZ_neg{ 20., 60.,  0.};   AxisSide gZ_pos{ 20., 60., 0.};
const double RHO_BRICK = 1.8;      // г/см3, силикатный/керамический кирпич (типичный диапазон 1,6-1,9)
G4Material* gWallMat = nullptr;    // назначается в main() — состав бетона NIST, плотность кирпича

// Удельные активности обычных ограждающих конструкций, Бк/кг (UNSCEAR, средние
// по миру; реальный разброс K 100..1600, Ra 10..100, Th 10..60). Используются
// только в режиме "all" (три серии сразу, для справочной картинки поля).
const double A_K40 = 400., A_RA226 = 40., A_TH232 = 30.;
const double RHO_CONCRETE = 2.30;  // г/см³, G4_CONCRETE

// РЕЖИМ ОДНОЙ СЕРИИ (argv[3] = K|Ra|Th|all). Причина: измеренный фон
// RC-103 показал излишек в полосе Tl-208 2614,5 кэВ против модели на
// усреднённом UNSCEAR-бетоне (см. план, раздел «Одним множителем не
// обойтись») — значит состав ограждающих конструкций ЭТОГО помещения
// отличается от справочного, и три активности надо подгонять по измерению,
// а не назначать. Подгонка линейна ПО АКТИВНОСТИ только если каждая серия
// посчитана отдельно с ЕДИНИЧНОЙ активностью (1 Бк/кг) — тогда результат
// этого прогона умножается на подобранную амплитуду напрямую, без повторного
// счёта Geant4. gSeries = -1 (all, три сразу, стандартные UNSCEAR) | 0 (K-40,
// 1 Бк/кг) | 1 (Ra-226 ряд, 1 Бк/кг) | 2 (Th-232 ряд, 1 Бк/кг).
int gSeries = -1;

// --- НУКЛИДНЫЙ РЕЖИМ (20.08.2026, канон geant4-spectrum-pipeline «метод 1»):
// шаблон на КАЖДОЕ ЗВЕНО цепочки, а не на цепочку целиком. Позволяет проверить
// РАВНОВЕСИЕ (Pb-214 против Bi-214 — радон; Ac-228 против Tl-208 — торий с
// ветвлением 35,94 %) вместо допущения о нём. Ключ nuc=<имя>, активность 1 Бк/кг
// родителя звена. Принадлежность линий — ENSDF/LNHB, явными списками (не по
// диапазонам энергий: диапазонное правило молча переклеивает линии при правке
// таблицы). Ra-224 240,0 слит с Pb-212: линии 238,6/240,0 на CsI НЕРАЗЛИЧИМЫ
// (FWHM ~30 кэВ), раздельные амплитуды были бы вырождены.
const double NUC_RA226[] = {186.2};
const double NUC_PB214[] = {241.9, 295.2, 351.9};
const double NUC_PB212[] = {238.6, 240.0, 300.1};
const double NUC_BI212[] = {727.3, 1620.5};
const double NUC_TL208[] = {510.7, 583.2, 860.6, 2614.5};
// series 1 без Ra-226/Pb-214 -> Bi-214;  series 2 без Pb-212/Bi-212/Tl-208 -> Ac-228

bool InList(double e, const double* arr, int n) {
  for (int i = 0; i < n; ++i) if (std::abs(e - arr[i]) < 0.5) return true;
  return false;
}

// -1 = режим выключен; иначе индекс в NUC_NAMES
int gNuc = -1;
const char* NUC_NAMES[] = {"K40", "Ra226", "Pb214", "Bi214",
                           "Pb212", "Ac228", "Bi212", "Tl208"};
const int N_NUC = 8;

// --- ИОННЫЙ РЕЖИМ (D-001, 21.08.2026) ---------------------------------------
// Ключ ion=<имя>: источником становится САМО ЯДРО, а не список его линий.
// Geant4 RDM берёт схему распада из ENSDF и даёт всё сразу: гаммы со своими
// интенсивностями БЕЗ порога отбора, бета-спектры, конверсионные электроны,
// характеристический рентген и ТОРМОЗНОЕ от бет. Таблица LINES не покрывает
// ничего из этого ниже 186,2 кэВ — там у неё нет ни одной линии вообще.
// Z и A по порядку NUC_NAMES; nucleusLimits ограничивает распад одним звеном,
// иначе цепочка пойдёт вниз и шаблон перестанет быть шаблоном звена.
const int NUC_Z[] = {19, 88, 82, 83, 82, 89, 83, 81};
const int NUC_A[] = {40, 226, 214, 214, 212, 228, 212, 208};
int gIonZ = 0, gIonA = 0;   // 0 = ионный режим выключен

int NucOfLine(double e_keV, int series) {
  if (series == 0) return 0;                                  // K-40
  if (series == 1) {
    if (InList(e_keV, NUC_RA226, 1)) return 1;                // Ra-226
    if (InList(e_keV, NUC_PB214, 3)) return 2;                // Pb-214
    return 3;                                                 // Bi-214
  }
  if (InList(e_keV, NUC_PB212, 3)) return 4;                  // Pb-212 (+Ra-224)
  if (InList(e_keV, NUC_BI212, 2)) return 6;                  // Bi-212
  if (InList(e_keV, NUC_TL208, 4)) return 7;                  // Tl-208
  return 5;                                                   // Ac-228
}

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
    // Слабые линии Bi-214 ВЫШЕ 2447,9 кэВ (добавлены 13.08.2026, задача №16).
    // До этого набор обрывался на 2447,9, а весь ряд Th — на 2614,5, из-за
    // чего модель НЕ МОГЛА дать ни одного отсчёта выше 2614,5 кэВ, и окно
    // подгонки мюонной амплитуды (2802-2840 кэВ) получало нулевое
    // гамма-предсказание ПО ПОСТРОЕНИЮ — весь измеренный там счёт
    // автоматически приписывался мюонам.
    // Источник: IAEA NDS / ENSDF, выгрузка 09.08.2026, файл репозитория
    // detectors/Gamma-1S/web-th232/data/ensdf_ra226_chain_lines.csv (✅).
    // Отобраны линии с интенсивностью >= 0,002 % (более слабые дают вклад
    // ниже статистической значимости прогона). Проценты ENSDF переведены в
    // доли: 0,0300 % -> 3,00e-4.
    // ⚠️ Оценка вклада (масштабированием отклика линии 1764,5 на отношение
    // выходов и эффективностей) даёт ~1e-6 cps против измеренных 1,63e-4
    // cps в окне — то есть эти линии окно НЕ объясняют и мюонную компоненту
    // не отменяют; они убирают нефизичный обрыв спектра, не более.
    {2472.9, 2.30e-5, 1},  {2505.5, 5.60e-5, 1},  {2694.7, 3.00e-4, 1},
    {2699.2, 2.77e-5, 1},  {2769.9, 2.45e-4, 1},  {2785.9, 5.50e-5, 1},
    {2827.0, 2.40e-5, 1},  {2880.4, 1.00e-4, 1},  {2893.6, 5.90e-5, 1},
    {2922.0, 1.36e-4, 1},  {2940.0, 3.60e-5, 1},  {2978.9, 1.36e-4, 1},
    {3000.0, 8.60e-5, 1},  {3053.9, 2.09e-4, 1},  {3081.8, 5.90e-5, 1},
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
  // gSeries>=0: только линии своей серии, активность = 1 Бк/кг — единичный
  // отклик для последующей подгонки амплитуды. Остальные линии получают вес
  // 0, а не выбрасываются из массива: розыгрыш просто никогда их не выберет,
  // и таблица LINES остаётся ЕДИНСТВЕННЫМ источником данных о линиях
  // (не дублируется под каждый режим).
  const double a[3] = {A_K40, A_RA226, A_TH232};
  double s = 0;
  for (int i = 0; i < NLINES; ++i) {
    double ai;
    if (gNuc >= 0) ai = (NucOfLine(LINES[i].E_keV, LINES[i].series) == gNuc) ? 1.0 : 0.0;
    else if (gSeries < 0) ai = a[LINES[i].series];
    else ai = (LINES[i].series == gSeries) ? 1.0 : 0.0;
    s += ai * LINES[i].yield * RHO_CONCRETE / 1000.0;  // 1/(см³·с)
    gCum[i] = s;
  }
  if (s <= 0) {
    G4Exception("wallfield::BuildWeights", "EMPTY_SERIES", FatalException,
               "выбранная серия не даёт ни одной линии с ненулевым весом");
  }
  gSvTotal = s;
  for (int i = 0; i < NLINES; ++i) gCum[i] /= s;

  if (gIonZ > 0) {
    // ИОННЫЙ РЕЖИМ: одно событие = ОДИН РАСПАД родителя, а не один фотон.
    // Поэтому объёмная скорость при 1 Бк/кг — это просто плотность распадов,
    // без всяких yield: они уже сидят внутри схемы распада ENSDF.
    // Оставить прежнюю нормировку значило бы поделить на сумму выходов
    // ещё раз и занизить шаблон во столько раз, сколько фотонов на распад.
    gSvTotal = RHO_CONCRETE / 1000.0;   // распад/(см^3 * с) при 1 Бк/кг
  }
}

G4LogicalVolume* gCavLV = nullptr;

}  // namespace

// --- геометрия ---------------------------------------------------------------

// Суммарные (наружная граница материала) и внутренние (граница воздушной
// полости до стены) полуразмеры + смещение центра box вдоль оси — используются
// и геометрией, и розыгрышем источника (#SHIELD-26 corner-режим).
struct BoxSpec { double halfCm, offsetCm; };
BoxSpec OuterSpec(const AxisSide& neg, const AxisSide& pos) {
  const double outNeg = neg.innerCm + neg.wallCm + neg.extraCm;
  const double outPos = pos.innerCm + pos.wallCm + pos.extraCm;
  return {0.5 * (outNeg + outPos), 0.5 * (outPos - outNeg)};
}
BoxSpec InnerSpec(const AxisSide& neg, const AxisSide& pos) {
  const double inNeg = neg.innerCm, inPos = pos.innerCm;
  return {0.5 * (inNeg + inPos), 0.5 * (inPos - inNeg)};
}
// Глобальные, нужны и геометрии, и генератору источника (симметричны по смыслу
// с R_WALL/R_CAV в сферическом режиме).
BoxSpec gOuterX, gOuterY, gOuterZ, gInnerX, gInnerY, gInnerZ;

class WallGeom : public G4VUserDetectorConstruction {
public:
  G4VPhysicalVolume* Construct() override {
    auto* nist = G4NistManager::Instance();
    if (!gCorner) {
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

    // --- corner-режим: асимметричный box, #SHIELD-26 -----------------------
    gOuterX = OuterSpec(gX_neg, gX_pos);
    gOuterY = OuterSpec(gY_neg, gY_pos);
    gOuterZ = OuterSpec(gZ_neg, gZ_pos);
    gInnerX = InnerSpec(gX_neg, gX_pos);
    gInnerY = InnerSpec(gY_neg, gY_pos);
    gInnerZ = InnerSpec(gZ_neg, gZ_pos);

    const double worldHalf = 1.1 * std::max({gOuterX.halfCm + std::abs(gOuterX.offsetCm),
                                              gOuterY.halfCm + std::abs(gOuterY.offsetCm),
                                              gOuterZ.halfCm + std::abs(gOuterZ.offsetCm)}) * cm;
    auto* world = new G4LogicalVolume(new G4Box("world", worldHalf, worldHalf, worldHalf),
                                      nist->FindOrBuildMaterial("G4_AIR"), "world");
    auto* pv = new G4PVPlacement(nullptr, {}, world, "world", nullptr, false, 0, true);

    auto* outerSolid = new G4Box("wall_outer", gOuterX.halfCm * cm, gOuterY.halfCm * cm,
                                 gOuterZ.halfCm * cm);
    auto* innerSolid = new G4Box("wall_inner", gInnerX.halfCm * cm, gInnerY.halfCm * cm,
                                 gInnerZ.halfCm * cm);
    // innerBox расположен внутри outerBox со смещением = разница их центров
    const G4ThreeVector relOffset(
        (gInnerX.offsetCm - gOuterX.offsetCm) * cm,
        (gInnerY.offsetCm - gOuterY.offsetCm) * cm,
        (gInnerZ.offsetCm - gOuterZ.offsetCm) * cm);
    auto* wallSolid = new G4SubtractionSolid("wall", outerSolid, innerSolid, nullptr, relOffset);

    auto* wall = new G4LogicalVolume(wallSolid, gWallMat, "wall");
    new G4PVPlacement(nullptr,
                       G4ThreeVector(gOuterX.offsetCm * cm, gOuterY.offsetCm * cm,
                                     gOuterZ.offsetCm * cm),
                       wall, "wall", world, false, 0, true);

    // Регистрирующая полость — маленький шар в ИСТИННОМ центре (0,0,0), как в
    // сферическом режиме; она заведомо внутри innerBox (R_CAV <= innerSpec.halfCm
    // с обеих сторон каждой оси по построению gX_neg/pos.innerCm >= R_CAV/cm).
    gCavLV = new G4LogicalVolume(new G4Orb("cav", R_CAV), nist->FindOrBuildMaterial("G4_AIR"), "cav");
    new G4PVPlacement(nullptr, {}, gCavLV, "cav", world, false, 0, true);
    return pv;
  }
};

// --- источник ----------------------------------------------------------------
class WallSource : public G4VUserPrimaryGeneratorAction {
  G4ParticleGun fGun{1};
  G4ParticleDefinition* fIon = nullptr;   // ионный режим, определяется лениво
public:
  WallSource() {
    fGun.SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("gamma"));
  }
  void GeneratePrimaries(G4Event* e) override {
    G4ThreeVector p;
    if (!gCorner) {
      // равномерно по бетону: розыгрыш в шаре с отбрасыванием полости
      do {
        p.set((2 * G4UniformRand() - 1) * R_WALL, (2 * G4UniformRand() - 1) * R_WALL,
              (2 * G4UniformRand() - 1) * R_WALL);
      } while (p.mag() > R_WALL || p.mag() < R_CAV);
    } else {
      // corner-режим: равномерно в outer box (со смещением центра), отбросить
      // если точка попала в innerBox (тоже со своим смещением) — #SHIELD-26.
      do {
        p.set(gOuterX.offsetCm * cm + (2 * G4UniformRand() - 1) * gOuterX.halfCm * cm,
              gOuterY.offsetCm * cm + (2 * G4UniformRand() - 1) * gOuterY.halfCm * cm,
              gOuterZ.offsetCm * cm + (2 * G4UniformRand() - 1) * gOuterZ.halfCm * cm);
      } while (std::abs((p.x() / cm - gInnerX.offsetCm)) < gInnerX.halfCm &&
               std::abs((p.y() / cm - gInnerY.offsetCm)) < gInnerY.halfCm &&
               std::abs((p.z() / cm - gInnerZ.offsetCm)) < gInnerZ.halfCm);
    }
    fGun.SetParticlePosition(p);

    if (gIonZ > 0) {
      // ИОННЫЙ РЕЖИМ (D-001): испускаем САМО ЯДРО в покое, остальное делает
      // RDM. Определение иона берётся лениво — таблица ионов существует
      // только после инициализации, в конструкторе её ещё нет.
      if (!fIon) fIon = G4IonTable::GetIonTable()->GetIon(gIonZ, gIonA, 0.0);
      fGun.SetParticleDefinition(fIon);
      fGun.SetParticleEnergy(0.0);
      fGun.SetParticleMomentumDirection({0, 0, 1});   // покоится, направление неважно
      fGun.GeneratePrimaryVertex(e);
      return;
    }

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
  // 2 кэВ на канал (было 10): ионный режим даёт ХРИ дочерних — Kα Po 77 и
  // Kβ 90 кэВ на 10-кэВной сетке сливались в два широких бина, форма
  // комплекса терялась ещё в источнике, до всякого детектора.
  static constexpr int kBins = 1500;
  static constexpr double kBinKeV = 2.0;
  std::vector<double> fLen{std::vector<double>(kBins + 1, 0.0)};
  G4String fOut = "wallfield.csv";

  void BeginOfRunAction(const G4Run*) override {
    std::fill(fLen.begin(), fLen.end(), 0.0);
  }
  void Add(double eKeV, double lenCm) {
    int b = static_cast<int>(eKeV / kBinKeV);
    if (b > kBins) b = kBins;  // последний канал — переполнение (E>=3000 кэВ,
                                // catch-all; см. WONTFIX D-004 в DECISIONS.md —
                                // make_m1_macros_v3.py трактует его как обычный
                                // узкий бин, эффект <0,005% полного потока)
    if (b >= 0) fLen[b] += lenCm;
  }
  void EndOfRunAction(const G4Run* run) override {
    const long N = run->GetNumberOfEvent();
    if (!N) return;
    double vWall;
    if (!gCorner) {
      vWall = 4. / 3. * pi * (std::pow(R_WALL / cm, 3) - std::pow(R_CAV / cm, 3));
    } else {
      const double vOuter = 8. * gOuterX.halfCm * gOuterY.halfCm * gOuterZ.halfCm;
      const double vInner = 8. * gInnerX.halfCm * gInnerY.halfCm * gInnerZ.halfCm;
      vWall = vOuter - vInner;  // см3, генератор уже отбрасывает inner-объём
    }
    const double vCav = 4. / 3. * pi * std::pow(R_CAV / cm, 3);
    // Реальная скорость испускания в объёме бетона, фотон/с
    const double rate = gSvTotal * vWall;
    // Флюенс = (средняя длина пробега на фотон) * (скорость) / объём полости
    const double norm = rate / (double(N) * vCav);

    double tot = 0;
    for (int i = 0; i <= kBins; ++i) tot += fLen[i] * norm;

    FILE* f = std::fopen(fOut.c_str(), "w");
    std::fprintf(f, "# поле ЕРН в помещении: флюенс в воздушной полости\n");
    std::fprintf(f, "# src_sha1 = %s\n", RCWF_SRC_SHA1);
    std::fprintf(f, "# git_describe = %s\n", RCWF_GIT_DESCRIBE);
    std::fprintf(f, "# build = %s %s\n", __DATE__, __TIME__);
    if (gNuc >= 0)
      std::fprintf(f, "# nuclide = %s, 1 Bq/kg parent, rho %.2f (nuclide template, pipeline method 1)\n",
                   NUC_NAMES[gNuc], RHO_CONCRETE);
    else if (gSeries < 0)
      std::fprintf(f, "# concrete: K-40 %.0f, Ra-226 %.0f, Th-232 %.0f Bq/kg, rho %.2f\n",
                   A_K40, A_RA226, A_TH232, RHO_CONCRETE);
    else
      std::fprintf(f, "# series = %s, активность 1 Bq/kg, rho %.2f (единичный "
                      "отклик для подгонки амплитуды, см. план раздел 4)\n",
                   gSeries == 0 ? "K-40" : gSeries == 1 ? "Ra-226" : "Th-232",
                   RHO_CONCRETE);
    if (!gCorner) {
      std::fprintf(f, "# R_wall_cm = %.1f  R_cav_cm = %.1f\n", R_WALL / cm, R_CAV / cm);
    } else {
      std::fprintf(f, "# CORNER #SHIELD-26: X- 1m/12cm brick +80cm(assumed more house), "
                      "X+ 80cm(=old model, no data), Y- 1.5m/50cm brick +0(loggia+street), "
                      "Y+ 80cm(=old model, no data), Z +-80cm(=old model, no data)\n");
      std::fprintf(f, "# outer_half_cm X=%.1f off=%.1f Y=%.1f off=%.1f Z=%.1f off=%.1f\n",
                   gOuterX.halfCm, gOuterX.offsetCm, gOuterY.halfCm, gOuterY.offsetCm,
                   gOuterZ.halfCm, gOuterZ.offsetCm);
      std::fprintf(f, "# inner_half_cm X=%.1f off=%.1f Y=%.1f off=%.1f Z=%.1f off=%.1f\n",
                   gInnerX.halfCm, gInnerX.offsetCm, gInnerY.halfCm, gInnerY.offsetCm,
                   gInnerZ.halfCm, gInnerZ.offsetCm);
      std::fprintf(f, "# wall_material = G4_CONCRETE composition, density %.2f g/cm3 (brick, "
                      "NOT measured)\n", RHO_BRICK);
    }
    std::fprintf(f, "# N = %ld\n", N);
    std::fprintf(f, "# Sv_total_per_cm3_s = %.6e\n", gSvTotal);
    std::fprintf(f, "# fluence_total_cm2_s = %.6e\n", tot);
    std::fprintf(f, "# bin_keV = %.3f  (последний канал = переполнение)\n", kBinKeV);
    std::fprintf(f, "E_keV,fluence_cm2_s\n");
    for (int i = 0; i <= kBins; ++i)
      if (fLen[i] > 0)
        std::fprintf(f, "%.1f,%.6e\n", (i + 0.5) * kBinKeV, fLen[i] * norm);
    std::fclose(f);
    G4cout << "RESULT N= " << N << " fluence_total= " << tot << " cm-2 s-1  file= "
           << fOut << G4endl;
    // СТОРОЖ: нулевой флюенс — это ПРОВАЛ прогона, а не «такой результат».
    // Прогон с нулём завершался кодом 0 и выглядел успешным (21.08, K-40 и
    // Ra-226 без порога thresholdForVeryLongDecayTime). Проверяем НАЛИЧИЕ
    // РЕЗУЛЬТАТА, а не наличие выхода.
    if (!(tot > 0)) {
      G4cerr << "FLUENCE_ZERO: прогон дал нулевой флюенс. Для долгоживущих ядер "
                "проверь /process/had/rdm/thresholdForVeryLongDecayTime." << G4endl;
      std::exit(3);
    }
    // СТОРОЖ ОТСЕЧЕНИЯ ЦЕПОЧКИ (21.08): печатаем фактический диапазон
    // nucleusLimits в шапку CSV, чтобы «шаблон звена» нельзя было принять за
    // таковой без доказательства. Без строки limits шаблон вбирает потомков,
    // и это НЕ ВИДНО ни по коду возврата, ни по флюенсу — только по составу
    // линий, который никто не смотрит.
    if (gIonZ > 0) {
      G4cout << "CHAIN_CUT: nucleusLimits A=" << gIonA << " Z=" << gIonZ
             << " (метод 1: распад ТОЛЬКО этого звена)" << G4endl;
    }
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
    if (gIonZ > 0) {
      // ИОННЫЙ РЕЖИМ (D-001). Без Decay+RDM ион просто стоит на месте.
      // Эталон: detectors/Gamma-1S/geometry/main.cc:82-84.
      RegisterPhysics(new G4DecayPhysics());
      RegisterPhysics(new G4RadioactiveDecayPhysics());
      // Порог 1 мм в бетоне отсекает электроны ниже ~350 кэВ: они не
      // трекаются, энергия кладётся локально — и ТОРМОЗНОГО ОТ НИХ НЕ БУДЕТ,
      // то есть ровно та компонента, ради которой режим и вводится, пропала
      // бы молча, а прогон при этом выглядел бы успешным.
      SetDefaultCutValue(0.05 * mm);
    } else {
      SetDefaultCutValue(1.0 * mm);  // линейный режим: важен транспорт гамма
    }
  }
};

int main(int argc, char** argv) {
  const long n = (argc > 1) ? std::atol(argv[1]) : 2000000;
  if (argc > 3) {
    const std::string s = argv[3];
    if (s == "K") gSeries = 0;
    else if (s == "Ra") gSeries = 1;
    else if (s == "Th") gSeries = 2;
    else if (s == "all") gSeries = -1;
    else {
      std::fprintf(stderr, "серия: ожидается K|Ra|Th|all, получено %s\n", s.c_str());
      return 2;
    }
  }
  // rwall=<см> — перекрытие радиуса бетонной сферы, только для проверки
  // насыщения (#SHIELD-16); corner=1 — асимметричная угловая геометрия
  // (#SHIELD-26). Аргумент argv[4], опционален.
  if (argc > 4) {
    const std::string a = argv[4];
    const std::string rwKey = "rwall=";
    const std::string cnKey = "corner=";
    if (a.rfind(rwKey, 0) == 0) R_WALL = std::atof(a.c_str() + rwKey.size()) * cm;
    else if (a.rfind(cnKey, 0) == 0) gCorner = std::atoi(a.c_str() + cnKey.size()) != 0;
    else if (a.rfind("nuc=", 0) == 0) {
      const std::string want = a.substr(4);
      for (int i = 0; i < N_NUC; ++i)
        if (want == NUC_NAMES[i]) { gNuc = i; break; }
      if (gNuc < 0) {
        std::fprintf(stderr, "nuc=: unknown nuclide\n");
        return 2;
      }
    }
    // ion=<имя> — ИОННЫЙ РЕЖИМ (D-001): полный распад ядра вместо таблицы
    // линий. Прежний nuc= сознательно оставлен рабочим: без него нельзя было
    // бы сравнить старый шаблон с новым на одних и тех же прогонах.
    else if (a.rfind("ion=", 0) == 0) {
      const std::string want = a.substr(4);
      for (int i = 0; i < N_NUC; ++i)
        if (want == NUC_NAMES[i]) { gNuc = i; gIonZ = NUC_Z[i]; gIonA = NUC_A[i]; break; }
      if (gIonZ <= 0) {
        std::fprintf(stderr, "ion=: unknown nuclide\n");
        return 2;
      }
    }
  }
  BuildWeights();
  if (gCorner) {
    gWallMat = G4NistManager::Instance()->BuildMaterialWithNewDensity(
        "G4_BRICK_approx", "G4_CONCRETE", RHO_BRICK * g / cm3);
  }

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
  if (gIonZ > 0) {
    // БЕЗ ЭТОЙ СТРОКИ ДОЛГОЖИВУЩИЕ ЯДРА НЕ РАСПАДАЮТСЯ ВООБЩЕ, а прогон
    // при этом завершается кодом 0 и пишет файл — с нулевым флюенсом.
    // Так и случилось 21.08: K-40 (1,25e9 лет) и Ra-226 (1600 лет) дали
    // fluence_total = 0 при EXITCODE=0. Порог RDM по умолчанию отсекает
    // распады с очень большим временем жизни.
    // Строка взята из эталона: Gamma-1S/macros/decay_th232_isotopes.mac:15.
    ui->ApplyCommand("/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+30 ns");
    // ОТСЕЧЕНИЕ ЦЕПОЧКИ — вторая обязательная строка метода 1. Без неё распад
    // идёт ВНИЗ ПО ВСЕЙ ВЕТВИ, и шаблон звена вбирает всех своих потомков:
    // у Ra-226 линия дочернего Bi-214 609 кэВ (6,3e-03) оказывалась СИЛЬНЕЕ
    // его собственной 186 кэВ (4,9e-03), а полный флюенс Ra-226 (0,145)
    // превышал флюенс Bi-214 (0,105) — при том что у радия одна слабая линия
    // с выходом 3,6 %. Отсюда корреляция шаблонов 0,993 и зануление амплитуд
    // в NNLS: шаблоны были не звеньями, а перекрывающимися кусками цепочек.
    // Метод 1 БЕЗ этой строки методом 1 не является (канон:
    // skills/geant4-spectrum-pipeline/SKILL.md:61-62).
    char lim[128];
    std::snprintf(lim, sizeof(lim),
                  "/process/had/rdm/nucleusLimits %d %d %d %d",
                  gIonA, gIonA, gIonZ, gIonZ);
    ui->ApplyCommand(lim);
    // Мутационная проверка сторожа проведена 21.08: со снятой строкой выше
    // K-40 даёт fluence_total = 0, сторож печатает FLUENCE_ZERO и завершает
    // прогон кодом 3 (лог wallion_mutation_K40.log). То есть сторож умеет
    // краснеть, а не только молчать при зелёном.
    // Проверка результата, а не кода возврата: нулевой флюенс — это провал,
    // даже когда всё «прошло успешно». Сторож ниже, после прогона.
  }
  ui->ApplyCommand("/run/beamOn " + std::to_string(n));
  delete rm;
  return 0;
}
