// Геометрия RadiaCode 101/102/103 + авторский сосуд Маринелли 200 мл.
//
// ОТКУДА ЧИСЛА
// ------------
// ПРИБОР, чертежи (drawings/rc101-103_case_*.png, разрезы Fusion):
//   корпус 123.00 x 34.00 x 17.50, кристалл CsI(Tl) 10x10x10 (1 см³, паспорт).
//   Размер 12.00 — от наружной плоскости носа до ЦЕНТРА кристалла: выносные
//   линии упираются в осевую линию кристалла, и пиксельный промер подсвеченного
//   кристалла в обеих проекциях даёт центр на 12.2 и 11.7 мм от торца.
//   Пара 8.20/9.30 — тоже до ЦЕНТРА от двух больших граней: их сумма ровно
//   17.50, промер даёт 8.09 и 9.42. Значит кристалл сдвинут на 0.55 мм к одной
//   большой грани, а по ширине стоит строго по центру (17.00 из 34.00).
//
// ПРИБОР, фото разборки (drawings/rc_teardown_*.jpg):
//   кристалл сидит в белой отражающей чашке с окном под фотоприёмник; чашка — в
//   отдельном чёрном модуле у носа, соединённом с платой шлейфом, то есть
//   основная плата под кристалл НЕ заходит. SiPM — одиночный кристалл ~6x6 мм
//   на круглой платке. Аккумулятор — Li-Po с маркировкой «602560 3.7V 1000mAh»,
//   то есть 6.0 x 25 x 60 мм; лежит в задней крышке, кристаллом не заслоняется.
//   Дисплей и плата — со стороны той грани, к которой ближе кристалл (8.20).
//
// СОСУД (STL v.2 Can/Cap, dnpro, thingiverse 6562353, CC BY; см. RCDetector.hh)
//   — промер лучевым сканированием:
//   наружный радиус 36.05, внутренний 33.24 (стенка 2.81), высота 69.30,
//   торцевая стенка со щелью 2.80, диск крышки 2.80.
//   Колодец — не цилиндр Ø30 (как на упрощённом эскизе объёма пробы), а тело по
//   форме прибора: полость 34.70 x 18.14, наружная поверхность гильзы
//   37.20 x 20.70, то есть стенка ровно 1.25 мм. Полость уходит на 47.81 мм от
//   наружной плоскости. Сечение — скруглённый прямоугольник: угловой профиль
//   радиуса даёт максимум 19.18 мм около 18°, тогда как прямой угол дал бы
//   sqrt(18.60²+10.35²) = 21.29, а полное скругление («стадион») — 18.60 при 0°.
//   Отсюда радиус углов сечения 6.75 мм.
//
// САМОПРОВЕРКИ
//   1. Полость стакана, залитая до среза горловины, даёт 200.15 см³ при высоте
//      пробы 66.50 мм — сходится и с надписью «200 ml» на сосуде, и с эскизом
//      объёма пробы Ø66.5 x 66.5. Построенное тело пробы должно дать столько же.
//   2. Центр кристалла при посадке прибора в колодец попадает в 0.2..0.4 мм от
//      центра объёма пробы, т.е. «кристалл в центре маринельки» выполняется.
//   3. Толщины стенок сошлись с ответом автора моделей в обсуждении thing:6562353
//      (комментарий от 11.05.2024): «Внешние стенки 2,8 мм ... Внутренняя стенка
//      1,2 мм». Промер дал 2.81 и 1.25 мм. Источник независимый: автор называет
//      проектные значения, промер идёт по выданной сетке треугольников.
//      Оттуда же — существенное для модели: внутреннюю стенку автор печатает со
//      100 % заполнением, то есть плотность стенки КОЛОДЦА (единственной, что
//      стоит между пробой и кристаллом) правомерно брать как у монолита.
//      Наружная стенка печатается с заполнением меньше 100 %, и сплошной пластик
//      там слегка завышает ослабление; путь через неё для отклика второстепенен.
//
// ЧТО НЕ МОДЕЛИРУЕТСЯ
//   Резьба и рёбра жёсткости сосуда (R до 40.25, ~6 см³ пластика при R>36.05) —
//   далеко от кристалла. Чёрный модуль детектора вокруг отражающей чашки: в
//   носовой части для него нет места при принятой толщине стенки, вклад ~1 мм
//   низкоатомного пластика. Заходная фаска колодца у щели.

#include "RCDetector.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4Polycone.hh"
#include "G4PhysicalConstants.hh"
#include "G4RotationMatrix.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"

#include <cmath>
#include <string>
#include <vector>

namespace {
G4Material* Mat(const G4String& n) {
  return G4NistManager::Instance()->FindOrBuildMaterial(n);
}

G4LogicalVolume* Put(G4VSolid* s, G4Material* m, const G4String& nm,
                     G4LogicalVolume* mother, const G4ThreeVector& pos,
                     const G4Colour& col, G4RotationMatrix* rot = nullptr) {
  auto* lv = new G4LogicalVolume(s, m, nm);
  auto* va = new G4VisAttributes(col);
  va->SetForceSolid(true);
  lv->SetVisAttributes(va);
  new G4PVPlacement(rot, pos, lv, nm, mother, false, 0, true);
  return lv;
}

// Пластина по мировым границам Z, толщиной t с центром по Y в yc.
G4Box* Slab(const G4String& nm, double x, double t, double z0, double z1) {
  return new G4Box(nm, 0.5 * x, 0.5 * t, 0.5 * (z1 - z0));
}
}  // namespace

RCDetector::RCDetector(bool withVessel) : fWithVessel(withVessel) {}

// ---------------------------------------------------------------------------
// Пресеты сосудов. ПРОМЕР STL: у обоих объём полости сошёлся с номиналом
// (200.15 и 498.9 см³), что и служит проверкой промера.
VesselGeom VesselGeom::Preset(const G4String& name) {
  VesselGeom v;
  if (name == "m200") return v;          // значения по умолчанию

  if (name != "m500") {
    G4Exception("VesselGeom::Preset", "rc001", FatalException,
                ("неизвестный сосуд: " + name).c_str());
  }
  v.name = "m500";
  v.outerR = -1;                          // обвод задан профилем
  // медианный по окружности наружный радиус, z_can от плоскости со щелью
  v.profZcan = {0.0, 9.5, 21.5, 33.5, 45.5, 57.5, 63.5, 66.5,
                78.5, 84.5, 87.5, 97.5};
  v.profR    = {47.1, 49.0, 50.4, 52.0, 53.8, 55.6, 54.7, 51.2,
                51.3, 51.1, 45.3, 46.1};
  v.innerR = 43.34;
  v.barrelH = 97.50;
  v.endWall = 3.30;
  // Колодец: наружная поверхность гильзы 40.50 x 23.60 (радиус углов 4.9),
  // полость 36.68 x 19.78 => стенка 1.91 мм — толще, чем 1.25 у m200.
  v.wellOutX = 20.25; v.wellOutY = 11.80; v.wellOutR = 4.90;
  v.wellInX  = 18.34; v.wellInY  =  9.89; v.wellInR  = 3.00;
  v.wellTip = 65.60; v.wellTipOut = 67.41;
  // как и у m200, зазор выбран из совпадения центров носовых куполов
  v.seatGap = v.wellInY - 8.75;           // 1.14 мм
  // Крышка (промер STL: 45.04 см³ пластика, габарит 110.9 x 110.4 x 14.3)
  // моделируется одним диском той же массы: юбка охватывает расширяющуюся книзу
  // бочку, толку от неё для отклика нет, а геометрию усложняет.
  v.capSkirtR = 55.45;
  v.capT = 45.04e3 / (pi * 55.45 * 55.45);   // 4.66 мм
  v.capH = v.capT;                            // юбки нет
  return v;
}

// ---------------------------------------------------------------------------
void RCDetector::DefineMaterials() {
  auto* nist = G4NistManager::Instance();
  auto* H = nist->FindOrBuildElement("H");
  auto* C = nist->FindOrBuildElement("C");
  auto* N = nist->FindOrBuildElement("N");
  auto* O = nist->FindOrBuildElement("O");
  auto* F = nist->FindOrBuildElement("F");
  auto* Al = nist->FindOrBuildElement("Al");
  auto* Si = nist->FindOrBuildElement("Si");
  auto* Ti = nist->FindOrBuildElement("Ti");
  auto* Br = nist->FindOrBuildElement("Br");
  auto* Co = nist->FindOrBuildElement("Co");
  auto* Cu = nist->FindOrBuildElement("Cu");
  auto* Li = nist->FindOrBuildElement("Li");

  auto* abs = new G4Material("ABS", 1.05 * g / cm3, 3);   // корпус прибора
  abs->AddElement(C, 0.851);
  abs->AddElement(H, 0.081);
  abs->AddElement(N, 0.068);

  // Белая отражающая чашка: полимер с наполнителем TiO2 (ДОПУЩЕНИЕ по виду)
  auto* refl = new G4Material("Reflector", 1.45 * g / cm3, 4);
  refl->AddElement(C, 0.55);
  refl->AddElement(H, 0.06);
  refl->AddElement(O, 0.29);
  refl->AddElement(Ti, 0.10);

  auto* pla = new G4Material("PLA", 1.24 * g / cm3, 3);   // C3H4O2
  pla->AddElement(C, 3);
  pla->AddElement(H, 4);
  pla->AddElement(O, 2);

  auto* petg = new G4Material("PETG", 1.27 * g / cm3, 3);  // C10H8O4
  petg->AddElement(C, 10);
  petg->AddElement(H, 8);
  petg->AddElement(O, 4);

  auto* fr4 = new G4Material("FR4", 1.85 * g / cm3, 5);   // типовой стеклотекстолит
  fr4->AddElement(Si, 0.2818);
  fr4->AddElement(O, 0.3937);
  fr4->AddElement(C, 0.2264);
  fr4->AddElement(H, 0.0281);
  fr4->AddElement(Br, 0.0700);

  auto* batt = new G4Material("LiPo", 2.00 * g / cm3, 7);  // суррогат Li-Po
  batt->AddElement(Al, 0.20);
  batt->AddElement(Cu, 0.10);
  batt->AddElement(C, 0.25);
  batt->AddElement(O, 0.20);
  batt->AddElement(Co, 0.15);
  batt->AddElement(Li, 0.02);
  batt->AddElement(F, 0.08);

  MakeMatrix(fVes.sampleMatrix, fVes.sampleDensity);

  (void)abs; (void)refl; (void)pla; (void)petg; (void)fr4; (void)batt;
}

// ---------------------------------------------------------------------------
// Матрицы пробы. Состав фиксирован, плотность — параметр: самопоглощение
// определяется произведением (плотность x массовый коэффициент ослабления),
// поэтому такой пары хватает, чтобы накрыть реальные пробы.
G4Material* RCDetector::MakeMatrix(const G4String& name, double rho,
                                   const G4String& g4name) {
  if (auto* have = G4Material::GetMaterial(g4name, false)) return have;
  auto* nist = G4NistManager::Instance();
  const double d = rho * g / cm3;

  if (name == "air") {
    auto* air = nist->FindOrBuildMaterial("G4_AIR");
    auto* m = new G4Material(g4name, air->GetDensity(), 1);
    m->AddMaterial(air, 1.0);
    return m;
  }
  if (name == "water") {
    auto* m = new G4Material(g4name, d, 2);
    m->AddElement(nist->FindOrBuildElement("H"), 2);
    m->AddElement(nist->FindOrBuildElement("O"), 1);
    return m;
  }
  if (name == "organic") {           // целлюлоза C6H10O5: зерно, сено, биомасса
    auto* m = new G4Material(g4name, d, 3);
    m->AddElement(nist->FindOrBuildElement("C"), 6);
    m->AddElement(nist->FindOrBuildElement("H"), 10);
    m->AddElement(nist->FindOrBuildElement("O"), 5);
    return m;
  }
  // Минеральные матрицы по массовым долям.
  struct Frac { const char* el; double w; };
  static const Frac soil[] = {{"O", 0.500}, {"Si", 0.280}, {"Al", 0.070},
                              {"Fe", 0.050}, {"Ca", 0.040}, {"K", 0.020},
                              {"Mg", 0.010}, {"Na", 0.010}, {"C", 0.020}};
  static const Frac ash[] = {{"O", 0.440}, {"Ca", 0.190}, {"Si", 0.130},
                             {"K", 0.090}, {"Mg", 0.040}, {"Al", 0.030},
                             {"Fe", 0.030}, {"P", 0.030}, {"C", 0.020},
                             {"S", 0.020}};
  const Frac* f = (name == "ash") ? ash : soil;
  const int n = (name == "ash") ? 10 : 9;
  auto* m = new G4Material(g4name, d, n);
  for (int i = 0; i < n; ++i)
    m->AddElement(nist->FindOrBuildElement(f[i].el), f[i].w);
  return m;
}

// ---------------------------------------------------------------------------
// Тело собирается ИЗ ПРИМИТИВОВ, а не выдавливанием: G4ExtrudedSolid в Geant4 —
// наследник тесселированного тела, и четыре таких тела в навигации давали сотни
// граней на каждом шаге (замер: 1.7 тыс. событий/с против 12 тыс. на примитивах).
//
// Сечение — скруглённый прямоугольник: два бруса крест-накрест плюс четыре
// цилиндра в углах. Нос — полуцилиндр радиусом hy с осью вдоль X плюс две сферы
// того же радиуса на его концах: в плоскости YZ это точная полуокружность (так
// нос и выглядит на продольном разрезе), в плане углы скруглены радиусом hy.
//
// У составного тела основание купола (стадион) уже, чем сечение прямой части
// (скруглённый прямоугольник), то есть в углах остаётся ступенька. Сама по себе
// она безобидна, но если у корпуса и у полости колодца ступеньки стоят в разных
// Z, корпус вылезает наружу. Поэтому посадочный зазор выбран так, чтобы центры
// куполов совпали: тогда полость объемлет корпус при любом Z (см. seatGap).
G4VSolid* RCDetector::Capsule(const G4String& nm, double hx, double hy, double r,
                              double len) {
  const double bz = 0.5 * (len - hy);        // полудлина прямой части
  const double cx = hx - r, cy = hy - r;

  G4VSolid* s = new G4Box(nm + "_bx", cx, hy, bz);
  s = new G4UnionSolid(nm + "_b2", s, new G4Box(nm + "_by", hx, cy, bz));
  auto* corner = new G4Tubs(nm + "_cr", 0., r, bz, 0., twopi);
  const double sx[4] = {1, -1, -1, 1}, sy[4] = {1, 1, -1, -1};
  for (int k = 0; k < 4; ++k)
    s = new G4UnionSolid(nm + "_c" + std::to_string(k), s, corner, nullptr,
                         G4ThreeVector(sx[k] * cx, sy[k] * cy, 0));

  const double cyl = hx - hy;                // полудлина носового полуцилиндра
  auto* rot = new G4RotationMatrix();
  rot->rotateY(90. * deg);
  s = new G4UnionSolid(nm + "_d", s, new G4Tubs(nm + "_dm", 0., hy, cyl, 0., twopi),
                       rot, G4ThreeVector(0, 0, -bz));
  auto* orb = new G4Orb(nm + "_ob", hy);
  s = new G4UnionSolid(nm + "_o1", s, orb, nullptr, G4ThreeVector(cyl, 0, -bz));
  return new G4UnionSolid(nm, s, orb, nullptr, G4ThreeVector(-cyl, 0, -bz));
}

// ---------------------------------------------------------------------------
void RCDetector::BuildDevice(G4LogicalVolume* world) {
  const DeviceGeom& d = fDev;
  const double hx = 0.5 * d.caseX, hy = 0.5 * d.caseY;
  const double noseTip = -d.crystalZ0;             // -12.00, кристалл в нуле
  const double cryY = -(hy - d.crystalToFace);     // -0.55

  // корпус
  auto* caseS = Capsule("case", hx, hy, d.caseEdgeR, d.caseZ);
  const double caseOZ = noseTip + CapOriginFromTip(hy, d.caseZ);
  auto* caseLV = Put(caseS, Mat("ABS"), "case", world, G4ThreeVector(0, 0, caseOZ),
                     G4Colour(0.25, 0.25, 0.28, 0.30));

  // внутренняя полость
  const double ihx = hx - d.wallSide, ihy = hy - d.wallFace;
  const double airTip = noseTip + d.wallNose;
  const double airTop = noseTip + d.caseZ - d.wallTail;
  const double airLen = airTop - airTip;
  const double airOZ = airTip + CapOriginFromTip(ihy, airLen);
  auto* airLV = Put(Capsule("caseAir", ihx, ihy, d.airEdgeR, airLen), Mat("G4_AIR"),
                    "caseAir", caseLV, G4ThreeVector(0, 0, airOZ - caseOZ),
                    G4Colour(0.8, 0.9, 1.0, 0.08));

  // мировые координаты -> система полости
  const auto inAir = [&](double y, double z) { return G4ThreeVector(0, y, z - airOZ); };

  // Белая отражающая чашка: закрывает кристалл со всех сторон, кроме -Y, где
  // оставлено окно под фотоприёмник (моделируется плёнкой 0.05 мм).
  const double win = 0.05;
  const double hc = 0.5 * d.crystal;
  const double rx = hc + d.reflector;              // по X и Z
  const double ryLo = hc + win, ryHi = hc + d.reflector;   // по -Y и +Y
  const double reflYc = cryY + 0.5 * (ryHi - ryLo);
  auto* reflLV = Put(new G4Box("reflector", rx, 0.5 * (ryLo + ryHi), rx),
                     Mat("Reflector"), "reflector", airLV, inAir(reflYc, 0.0),
                     G4Colour(0.97, 0.97, 0.95, 0.55));
  fCrystalLV = Put(new G4Box("crystal", hc, hc, hc), Mat("G4_CESIUM_IODIDE"), "crystal",
                   reflLV, G4ThreeVector(0, cryY - reflYc, 0),
                   G4Colour(0.35, 0.55, 0.95, 0.85));

  // фотоприёмник и круглая платка под окном (сторона -Y)
  const double yWin = cryY - ryLo;                 // наружная плоскость окна
  auto* rotY = new G4RotationMatrix();
  rotY->rotateX(90. * deg);                        // ось диска -> вдоль Y
  Put(new G4Box("sipm", 0.5 * d.sipmSide, 0.5 * d.sipm, 0.5 * d.sipmSide), Mat("G4_Si"),
      "sipm", airLV, inAir(yWin - 0.5 * d.sipm, 0.0), G4Colour(0.9, 0.7, 0.2));
  Put(new G4Tubs("sipmPcb", 0., 0.5 * d.subDia, 0.5 * d.subT, 0., twopi), Mat("FR4"),
      "sipmPcb", airLV, inAir(yWin - d.sipm - 0.5 * d.subT, 0.0),
      G4Colour(0.95, 0.95, 0.9), rotY);

  // Со стороны -Y к хвосту: дисплей у стенки, за ним основная плата.
  // Ширины дисплея и платы ограничены скруглением полости корпуса.
  const double yDisp = -ihy + 0.25 + 0.5 * d.dispT;
  Put(Slab("display", d.dispX, d.dispT, d.dispZ0, d.dispZ1), Mat("G4_GLASS_PLATE"),
      "display", airLV, inAir(yDisp, 0.5 * (d.dispZ0 + d.dispZ1)),
      G4Colour(0.15, 0.15, 0.18, 0.7));
  Put(Slab("pcb", d.pcbX, d.pcbT, d.pcbZ0, d.pcbZ1), Mat("FR4"), "pcb", airLV,
      inAir(yDisp + 0.5 * d.dispT + 0.5 * d.pcbT, 0.5 * (d.pcbZ0 + d.pcbZ1)),
      G4Colour(0.1, 0.45, 0.2));

  // Аккумулятор — в задней крышке, со стороны +Y.
  Put(Slab("batt", d.battX, d.battT, d.battZ0, d.battZ1), Mat("LiPo"), "batt", airLV,
      inAir(ihy - 0.50 - 0.5 * d.battT, 0.5 * (d.battZ0 + d.battZ1)),
      G4Colour(0.75, 0.7, 0.35));
}

// ---------------------------------------------------------------------------
void RCDetector::BuildVessel(G4LogicalVolume* world) {
  const VesselGeom& v = fVes;
  // Нос прибора в z = -crystalZ0. Вершина полости колодца — на seatGap ниже.
  const double zSlot = -fDev.crystalZ0 + v.wellTip - v.seatGap;    // 35.61
  const double zRim = zSlot - v.barrelH;                           // -33.69
  const double zSmpTop = zSlot - v.endWall;                        // 32.81
  const double smpH = zSmpTop - zRim;                              // 66.50

  const double zOut = 0.5 * (zSlot + zRim);
  const double zIn = 0.5 * (zSmpTop + zRim);

  // наружный обвод: цилиндр (m200) либо тело вращения по профилю (m500)
  G4VSolid* tOut = nullptr;
  if (v.outerR > 0) {
    tOut = new G4Tubs("tOut", 0., v.outerR, 0.5 * v.barrelH, 0., twopi);
  } else {
    const size_t n = v.profZcan.size();
    std::vector<double> zp(n), ri(n, 0.), ro(n);
    for (size_t i = 0; i < n; ++i) {
      // профиль задан от плоскости со щелью вниз; в локальной системе тела
      // ось Z направлена как в мире, поэтому порядок инвертируется
      const size_t k = n - 1 - i;
      zp[i] = (zSlot - v.profZcan[k]) - zOut;
      ro[i] = v.profR[k];
    }
    tOut = new G4Polycone("tOut", 0., twopi, static_cast<G4int>(n),
                          zp.data(), ri.data(), ro.data());
  }
  auto* tIn = new G4Tubs("tIn", 0., v.innerR, 0.5 * smpH, 0., twopi);

  const double cOutTip = zSlot - v.wellTipOut;
  const double cOutLen = zSlot - cOutTip;
  auto* cOut = Capsule("cOut", v.wellOutX, v.wellOutY, v.wellOutR, cOutLen);
  const double zCOut = cOutTip + CapOriginFromTip(v.wellOutY, cOutLen);

  const double cInTip = zSlot - v.wellTip;
  const double cInLen = (zSlot + 10.0) - cInTip;   // с запасом, чтобы прорезать торец
  auto* cIn = Capsule("cIn", v.wellInX, v.wellInY, v.wellInR, cInLen);
  const double zCIn = cInTip + CapOriginFromTip(v.wellInY, cInLen);

  // пластик стакана = ((наружный цилиндр - полость) + гильза) - колодец
  auto* s1 = new G4SubtractionSolid("can1", tOut, tIn, nullptr,
                                    G4ThreeVector(0, 0, zIn - zOut));
  auto* s2 = new G4UnionSolid("can2", s1, cOut, nullptr,
                              G4ThreeVector(0, 0, zCOut - zOut));
  auto* s3 = new G4SubtractionSolid("can", s2, cIn, nullptr,
                                    G4ThreeVector(0, 0, zCIn - zOut));
  Put(s3, Mat(v.plasticMat), "can", world, G4ThreeVector(0, 0, zOut),
      G4Colour(0.75, 0.75, 0.72, 0.25));

  // проба = полость стакана минус гильза колодца
  auto* smp = new G4SubtractionSolid("sampleS", tIn, cOut, nullptr,
                                     G4ThreeVector(0, 0, zCOut - zIn));
  fSampleLV = Put(smp, Mat("Sample"), "sample", world, G4ThreeVector(0, 0, zIn),
                  G4Colour(0.85, 0.65, 0.35, 0.25));
  fSampleVolumeCm3 = smp->GetCubicVolume() / cm3;

  // крышка: диск под срезом горловины + юбка снаружи стакана.
  // Радиус у горловины: у m200 это цилиндр, у m500 — последняя точка профиля.
  const double skirtH = v.capH - v.capT;
  const double rimR = (v.outerR > 0) ? v.outerR : v.profR.back();
  const double discR = (skirtH > 0.1) ? rimR : v.capSkirtR;
  Put(new G4Tubs("capDisc", 0., discR, 0.5 * v.capT, 0., twopi), Mat(v.plasticMat),
      "capDisc", world, G4ThreeVector(0, 0, zRim - 0.5 * v.capT),
      G4Colour(0.7, 0.7, 0.68, 0.4));
  if (skirtH > 0.1) {
    Put(new G4Tubs("capSkirt", rimR + 0.01, v.capSkirtR, 0.5 * skirtH, 0., twopi),
        Mat(v.plasticMat), "capSkirt", world,
        G4ThreeVector(0, 0, zRim + 0.5 * skirtH), G4Colour(0.7, 0.7, 0.68, 0.4));
  }
}

// ---------------------------------------------------------------------------
G4VPhysicalVolume* RCDetector::Construct() {
  DefineMaterials();

  auto* worldS = new G4Box("world", 150. * mm, 150. * mm, 200. * mm);
  auto* worldLV = new G4LogicalVolume(worldS, Mat("G4_AIR"), "world");
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto* worldPV =
      new G4PVPlacement(nullptr, {}, worldLV, "world", nullptr, false, 0, true);

  BuildDevice(worldLV);
  if (fWithVessel) BuildVessel(worldLV);

  G4cout << "\n=== RadiaCode: геометрия ===" << G4endl;
  G4cout << "  CsI(Tl) " << fDev.crystal << " мм куб, центр в начале координат"
         << G4endl;
  if (fWithVessel) {
    const double zSlot = -fDev.crystalZ0 + fVes.wellTip - fVes.seatGap;
    const double zc = 0.5 * ((zSlot - fVes.endWall) + (zSlot - fVes.barrelH));
    G4cout << "  сосуд " << fVes.name << " (" << fVes.plasticMat << "), проба "
           << fVes.sampleMatrix << " " << fVes.sampleDensity << " г/см³"
           << G4endl;
    G4cout << "  объём пробы: " << fSampleVolumeCm3 << " см³ (промер STL: "
           << (fVes.name == "m200" ? "200.15" : "498.9") << ")" << G4endl;
    G4cout << "  центр пробы: z = " << zc << " мм от центра кристалла" << G4endl;
  } else {
    G4cout << "  сосуд не построен: прибор в воздухе" << G4endl;
  }
  G4cout << G4endl;
  return worldPV;
}

