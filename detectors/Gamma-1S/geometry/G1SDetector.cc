// Геометрия спектрометра ГАММА-1С: УДС-ГЦ-63х63-USB в экране-защите «Экран-1СГ».
//
// ОТКУДА ЧИСЛА
// ------------
// ДЕТЕКТОР — чертёж «Чертеж 63х63.pdf» (БДС-Г, УДС-Г-63х63), подписи слоёв:
//   резина 1 мм / Al / резина 2 мм / Al / MgO, габариты Ø78, 74, 315, Ø71,
//   MgO 3,65 радиально и 6 по торцу, воздух 0,5.
//   Промер растра (1200 dpi, доля тёмных пикселей в строке) подтверждает, что
//   чертёж масштабный: диаметр кристалла 90,42 pt при осевой линии посередине,
//   то есть 1,435 pt/мм, и границы слоёв ложатся на подписанные толщины.
//   Две проставленные суммы сходятся:
//     Ø: 63 + 2*(3,65 + 0,5 + 2 + 1,5) = 78,3  против Ø78;
//     L:  63 + (6 + 0,5 + 2 + 1,5 + 0,5 + 1)  = 74,5  против 74.
//   Отсюда однозначное чтение: радиально резины 1 мм и воздуха НЕТ, они только
//   на входном торце; наружный корпус Al 1,5 мм.
//   Ø71 на чертеже относится к баллону ФЭУ и совпадает с наружным диаметром
//   банки кристалла (63 + 2*(3,65+0,5) = 71,3) — банка и ФЭУ состыкованы.
//
// ЭКРАН-ЗАЩИТА — состав задан оператором: сталь 3 мм, свинец 50 мм, кадмий,
//   медь. Толщины кадмия и меди не заданы; приняты по 1 мм. Порядок слоёв
//   от полости наружу: Cu -> Cd -> Pb -> сталь. Это стандартная градуированная
//   защита: ХРИ свинца 72–88 кэВ гасится кадмием, ХРИ кадмия 23 кэВ — медью,
//   ХРИ меди 8 кэВ уже ниже рабочего диапазона (паспорт: от 50 кэВ).
//
// РАЗМЕР ПОЛОСТИ — НЕ задан ни оператором, ни доступными документами
//   (руководство отсылает к отдельному РЭ на «Экран-1СГ» ДЦКИ.305179.038,
//   которого нет). Подобран по таблице 2.2 паспорта ДЦКИ.412131.001 ПС
//   «Содержание цветных металлов, кг, не менее»: свинец 165, медь 1,6,
//   кадмий 1,2, алюминий 1,7. Полость Ø200 x 190 мм даёт свинца 167,1 кг
//   и меди ровно 1,60 кг при вкладыше 1 мм — оба сходятся с паспортом.
//   Это подгонка под массу, а не чертёж: см. ПРОВЕРКА ПО МАССАМ в конце файла.
//
// ЧТО НЕ МОДЕЛИРУЕТСЯ
//   Делитель ФЭУ, плата USB-АЦП и разъёмы (заменены однородным блоком за
//   баллоном). Магнитный экран из пермаллоя — на чертеже его нет. Тележка,
//   петли и механизм крышки экрана. Кабельный ввод в свинце.

#include "G1SDetector.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4Polycone.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"

#include <cstdio>
#include <vector>

namespace {

G4Material* Mat(const G4String& n) {
  return G4NistManager::Instance()->FindOrBuildMaterial(n);
}

// Кольцо/диск по мировым границам z. rin=0 даёт сплошной диск.
G4LogicalVolume* Ring(const G4String& nm, double rin, double rout,
                      double z0, double z1, G4Material* m,
                      G4LogicalVolume* mother, const G4Colour& col,
                      G4LogicalVolume* reuse = nullptr) {
  auto* s = new G4Tubs(nm, rin * mm, rout * mm, 0.5 * (z1 - z0) * mm,
                       0.0, CLHEP::twopi);
  G4LogicalVolume* lv = reuse;
  if (!lv) {
    lv = new G4LogicalVolume(s, m, nm);
    auto* va = new G4VisAttributes(col);
    va->SetForceSolid(true);
    lv->SetVisAttributes(va);
  } else {
    // Отдельное тело того же материала: логический объём один на материал
    // не сделать, поэтому создаём новый, но с теми же атрибутами.
    lv = new G4LogicalVolume(s, m, nm);
    lv->SetVisAttributes(reuse->GetVisAttributes());
  }
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0.5 * (z0 + z1) * mm),
                    lv, nm, mother, false, 0, true);
  return lv;
}

double CylCm3(double rin, double rout, double h) {   // мм -> см³
  return CLHEP::pi * (rout * rout - rin * rin) * h / 1000.0;
}

}  // namespace

// ---------------------------------------------------------------------------
void G1SDetector::DefineMaterials() {
  auto* nist = G4NistManager::Instance();
  nist->FindOrBuildMaterial("G4_AIR");
  nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");
  nist->FindOrBuildMaterial("G4_Al");
  nist->FindOrBuildMaterial("G4_Pb");
  nist->FindOrBuildMaterial("G4_Cd");
  nist->FindOrBuildMaterial("G4_Cu");
  nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  nist->FindOrBuildMaterial("G4_Pyrex_Glass");
  nist->FindOrBuildMaterial("G4_RUBBER_NATURAL");
  nist->FindOrBuildMaterial("G4_Galactic");

  // Отражатель — НАСЫПНОЙ порошок MgO, а не спечённая керамика (ρ=3,58).
  // Плотность — параметр fHead.mgoDensity, см. пояснение в G1SDetector.hh.
  auto* mgo = nist->BuildMaterialWithNewDensity(
      "MgO_powder", "G4_MAGNESIUM_OXIDE", fHead.mgoDensity * g / cm3);
  (void)mgo;

  // Блок за баллоном ФЭУ: делитель, плата, разъёмы. ДОПУЩЕНИЕ: однородная
  // смесь текстолита и меди эффективной плотности 0,5 г/см³.
  auto* el = new G4Material("Electronics", 0.5 * g / cm3, 3);
  el->AddMaterial(nist->FindOrBuildMaterial("G4_BAKELITE"), 0.80);
  el->AddMaterial(nist->FindOrBuildMaterial("G4_Cu"), 0.15);
  el->AddMaterial(nist->FindOrBuildMaterial("G4_Al"), 0.05);

  // Имя материала — по роли («Sample»), а не по составу: матрица переключаема.
  MakeMatrix(fVessel.sampleMatrix, fVessel.sampleDensity, "Sample");
  nist->FindOrBuildMaterial("G4_POLYPROPYLENE");
}

// ---------------------------------------------------------------------------
// Матрица пробы ОИСН-16 — состав по массе из .efa ЛСРМ (Material=…Compound):
// H 0,022  C 0,206  N 0,009  O 0,049  Fe 0,714. Это имитатор насыпного образца
// на железной основе: 71 % железа по массе, НЕ грунт и НЕ органика.
// Отдельный статический метод нужен, чтобы mucalc считал mu по ТОЙ ЖЕ смеси.
G4Material* G1SDetector::MakeMatrix(const G4String& matrix, double rho,
                                    const G4String& g4name) {
  if (auto* have = G4Material::GetMaterial(g4name, false)) return have;
  auto* nist = G4NistManager::Instance();

  // Вода — вторая матрица, не для сверки с ЛСРМ, а для МИА: паспорт задаёт её
  // для сосуда Маринелли с ДИСТИЛЛИРОВАННОЙ ВОДОЙ. Считать прямо в воде лучше,
  // чем пересчитывать эффективность с ОИСН-16 формулой: у ОИСН-16 71 % железа,
  // и перенос между столь разными матрицами был бы экстраполяцией.
  if (matrix == "water") {
    auto* w = new G4Material(g4name, rho * g / cm3, 2);
    w->AddElement(nist->FindOrBuildElement("H"), 2);
    w->AddElement(nist->FindOrBuildElement("O"), 1);
    return w;
  }

  // Матрица смесевого источника РИСН-379 — состав по массе из ОРИГИНАЛЬНОГО
  // .spe (поле MATERIAL, «Поверка 2016»): лёгкая органо-минеральная основа
  // с 20 % кальция. НЕ вода: на 59,5 кэВ (Am-241) кальций поднимает
  // поглощение заметно.
  if (matrix == "risn379") {
    auto* m = new G4Material(g4name, rho * g / cm3, 7);
    m->AddElement(nist->FindOrBuildElement("H"), 0.0430);
    m->AddElement(nist->FindOrBuildElement("C"), 0.3303);
    m->AddElement(nist->FindOrBuildElement("N"), 0.0120);
    m->AddElement(nist->FindOrBuildElement("O"), 0.3484);
    m->AddElement(nist->FindOrBuildElement("Na"), 0.0410);
    m->AddElement(nist->FindOrBuildElement("Mg"), 0.0220);
    m->AddElement(nist->FindOrBuildElement("Ca"), 0.2033);
    return m;
  }
  if (matrix != "OISN16") {
    G4Exception("G1SDetector::MakeMatrix", "g1s002", FatalException,
                ("неизвестная матрица: " + matrix).c_str());
  }
  auto* m = new G4Material(g4name, rho * g / cm3, 5);
  m->AddElement(nist->FindOrBuildElement("H"), 0.022);
  m->AddElement(nist->FindOrBuildElement("C"), 0.206);
  m->AddElement(nist->FindOrBuildElement("N"), 0.009);
  m->AddElement(nist->FindOrBuildElement("O"), 0.049);
  m->AddElement(nist->FindOrBuildElement("Fe"), 0.714);
  return m;
}

// ---------------------------------------------------------------------------
G4VPhysicalVolume* G1SDetector::Construct() {
  DefineMaterials();

  auto* worldS = new G4Box("World", 400 * mm, 400 * mm, 700 * mm);
  auto* worldLV = new G4LogicalVolume(worldS, Mat("G4_AIR"), "World");
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World",
                                    nullptr, false, 0, true);

  if (fWithShield) BuildShield(worldLV);
  BuildHead(worldLV);
  if (fWithVessel) {
    // Маринелли — с колодцем, «Дента» и Петри — плоские кюветы на торце
    if (fVessel.wellInR > 0) BuildVessel(worldLV);
    else BuildCup(worldLV);
  }
  return worldPV;
}

// ---------------------------------------------------------------------------
// Сосуд Маринелли 1 л, надет на головку сверху. Габариты — таблица кювет ЛСРМ
// («Прецизионные измерения», с. 11): внешние Ø150, H = 110. Сосуд выше
// номинального объёма, поэтому уровень засыпки вычисляется из целевых 1000 мл;
// над пробой до крышки — воздух. Колодец Ø80 садится на головку Ø78,3; дно
// колодца лежит на резиновом протекторе торца (z = 43).
//
// ПРОВЕРКА 1 (объём): вычисленный уровень должен лечь НИЖЕ крышки — иначе
//   целевой объём в сосуд не помещается (печатается в ReportMasses).
// ПРОВЕРКА 2 (эффективная толщина): подгонка f(мю*ро*d) к отношению
//   эффективностей двух плотностей должна дать d = 26..31 мм (табличное
//   значение ЛСРМ и Thick из .efa этого экземпляра).
void G1SDetector::BuildVessel(G4LogicalVolume* w) {
  const VesselGeom& v = fVessel;
  const double zFace = 43.00;                 // наружная плоскость головки
  const double zWellFloor = zFace + v.wall;   // дно колодца (проба выше)
  const double zBot = zFace - v.wellDepth + v.wall;  // низ юбки сосуда
  const double zTop = zBot + v.height;               // верх сосуда
  const double rWellIn = v.wellInR;           // 40,00 — полость колодца
  const double rWellOut = rWellIn + v.wall;   // 42,00 — колодец снаружи
  const double rIn = v.outerR - v.wall;       // 73,00 — полость стакана
  const double rOut = v.outerR;               // 75,00

  // Уровень засыпки из целевого объёма: кольцо вокруг колодца заполняется
  // целиком, остаток ложится сплошным слоем над дном колодца.
  const double ringCm3 = CylCm3(rWellOut, rIn, zWellFloor - (zBot + v.wall));
  const double topH = (v.sampleCm3 - ringCm3) * 1000.0
                    / (CLHEP::pi * rIn * rIn);            // мм
  const double zFill = zWellFloor + topH;
  fSampleVolumeCm3 = ringCm3 + CylCm3(0, rIn, topH);
  fSampleFits = (zFill <= zTop - v.wall);

  const G4Colour cPP(0.9, 0.9, 0.6), cSm(0.55, 0.35, 0.15);
  auto* pp = Mat("G4_POLYPROPYLENE");

  // Пластик: обечайка, дно кольцом (вокруг колодца), крышка, гильза колодца,
  // донце колодца.
  Ring("V_side", rIn, rOut, zBot, zTop, pp, w, cPP);
  Ring("V_bottom", rWellOut, rIn, zBot, zBot + v.wall, pp, w, cPP);
  Ring("V_lid", 0, rIn, zTop - v.wall, zTop, pp, w, cPP);
  Ring("V_well", rWellIn, rWellOut, zBot, zFace + v.wall, pp, w, cPP);
  Ring("V_wellfloor", 0, rWellIn, zFace, zWellFloor, pp, w, cPP);

  // Проба — ОДНО тело (полигон вращения): кольцо вокруг колодца, переходящее
  // в сплошной слой до вычисленного уровня засыпки. Один физический объём
  // нужен GPS-розыгрышу (/gps/pos/confine принимает одно имя).
  auto* sm = Mat("Sample");
  const double zs[4] = {(zBot + v.wall) * mm, zWellFloor * mm,
                        zWellFloor * mm, zFill * mm};
  const double ri[4] = {rWellOut * mm, rWellOut * mm, 0, 0};
  const double ro[4] = {rIn * mm, rIn * mm, rIn * mm, rIn * mm};
  auto* smSolid = new G4Polycone("Sample", 0, CLHEP::twopi, 4, zs, ri, ro);
  fSampleLV = new G4LogicalVolume(smSolid, sm, "Sample");
  auto* va = new G4VisAttributes(cSm);
  va->SetForceSolid(true);
  fSampleLV->SetVisAttributes(va);
  new G4PVPlacement(nullptr, {}, fSampleLV, "Sample", w, false, 0, true);
}

// ---------------------------------------------------------------------------
// Пресеты сосудов комплекта поверки.
//
// Габариты — таблица измерительных кювет ЛСРМ («Прецизионные измерения»,
// с. 11): Маринелли 1,0 л — Ø150 H=110; «Дента» 0,12 л — Ø75 H=35;
// Петри 0,075 л — Ø88 H=14. Объёмы засыпки и матрица — из .efa комплекта:
// все три геометрии калиброваны на ОИСН-16 при ρ = 1,6, то есть на 1600, 192
// и 96 г соответственно — именно такие массы стоят у источника Th-232
// в описи комплекта, что и подтверждает прочтение.
//
// Дента и Петри ставятся ПРЯМО НА торец детектора (Distance = 0 в .efa),
// колодца у них нет.
VesselGeom VesselGeom::Preset(const G4String& n) {
  VesselGeom v;
  v.name = n;
  if (n == "marinelli") return v;          // значения по умолчанию

  if (n == "denta") {
    v.outerR = 37.50;      // ЛСРМ табл.: Ø75
    v.height = 35.00;      // ЛСРМ табл.: H = 35
    v.wall = 1.50;         // ДОПУЩЕНИЕ: стенка пластмассовой кюветы
    v.wellInR = -1;        // колодца нет
    v.sampleCm3 = 120.0;   // .efa: Volume = 120 мл
    return v;
  }
  if (n == "petri") {
    v.outerR = 44.00;      // ЛСРМ табл.: Ø88
    v.height = 14.00;      // ЛСРМ табл.: H = 14
    v.wall = 1.50;         // ДОПУЩЕНИЕ
    v.wellInR = -1;
    v.sampleCm3 = 60.0;    // .efa: Volume = 60 мл
    return v;
  }
  G4Exception("VesselGeom::Preset", "g1s003", FatalException,
              ("неизвестный сосуд: " + n).c_str());
  return v;
}

// ---------------------------------------------------------------------------
// Плоская кювета без колодца («Дента», Петри): стакан стоит дном на торце
// детектора, засыпка — от внутреннего дна до уровня, заданного объёмом.
void G1SDetector::BuildCup(G4LogicalVolume* w) {
  const VesselGeom& v = fVessel;
  const double zFace = 43.00;                 // наружная плоскость головки
  const double zBot = zFace;                  // сосуд стоит на торце
  const double zTop = zBot + v.height;
  const double rIn = v.outerR - v.wall;
  const double zIn = zBot + v.wall;           // внутреннее дно

  const double h = v.sampleCm3 * 1000.0 / (CLHEP::pi * rIn * rIn);  // мм
  const double zFill = zIn + h;
  fSampleVolumeCm3 = CylCm3(0, rIn, h);
  fSampleFits = (zFill <= zTop - v.wall);

  const G4Colour cPP(0.9, 0.9, 0.6), cSm(0.55, 0.35, 0.15);
  auto* pp = Mat("G4_POLYPROPYLENE");
  Ring("V_side", rIn, v.outerR, zBot, zTop, pp, w, cPP);
  Ring("V_bottom", 0, v.outerR, zBot, zIn, pp, w, cPP);
  Ring("V_lid", 0, rIn, zTop - v.wall, zTop, pp, w, cPP);

  auto* s = new G4Tubs("Sample", 0, rIn * mm, 0.5 * h * mm, 0, CLHEP::twopi);
  fSampleLV = new G4LogicalVolume(s, Mat("Sample"), "Sample");
  auto* va = new G4VisAttributes(cSm);
  va->SetForceSolid(true);
  fSampleLV->SetVisAttributes(va);
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0.5 * (zIn + zFill) * mm),
                    fSampleLV, "Sample", w, false, 0, true);
}

// ---------------------------------------------------------------------------
// Головка. Оси координат: центр кристалла — начало, +Z к входному торцу.
void G1SDetector::BuildHead(G4LogicalVolume* w) {
  const HeadGeom& h = fHead;

  const double rCry = 0.5 * h.cryDia;              // 31,50
  const double zCry = 0.5 * h.cryLen;              // 31,50
  const double rMgo = rCry + h.mgoRad;             // 35,15
  const double rCan = rMgo + h.alCan;              // 35,65  (Ø71,3 = Ø71 ФЭУ)
  const double rRub = rCan + h.rubber2;            // 37,65
  const double rCase = rRub + h.alCase;            // 39,15  (Ø78,3 = Ø78)

  // Границы по оси, от кристалла к входному торцу
  const double zMgoTop  = zCry + h.mgoFace;        // 37,50
  const double zCanTop  = zMgoTop + h.alCan;       // 38,00
  const double zRubTop  = zCanTop + h.rubber2;     // 40,00
  const double zCaseTop = zRubTop + h.alCase;      // 41,50
  const double zAirTop  = zCaseTop + h.airGap;     // 42,00
  const double zFace    = zAirTop + h.rubber1;     // 43,00 — наружная плоскость
  const double zWinBot  = -zCry - h.window;        // -36,50 — низ световода
  const double zPmtBot  = zWinBot - h.pmtLen;      // -156,50
  const double zTail    = zFace - h.unitLen;       // -272,00

  const G4Colour cCry(0.2, 0.9, 0.3), cMgo(1, 1, 1), cAl(0.7, 0.7, 0.75),
      cRub(0.15, 0.15, 0.15), cGlass(0.6, 0.8, 1.0), cEl(0.5, 0.3, 0.1);

  // Кристалл
  fCrystalLV = Ring("NaI", 0, rCry, -zCry, zCry, Mat("G4_SODIUM_IODIDE"), w, cCry);

  // Отражатель MgO: кольцо вдоль кристалла + шайба на входном торце
  Ring("MgO_side", rCry, rMgo, -zCry, zCry, Mat("MgO_powder"), w, cMgo);
  Ring("MgO_face", 0, rMgo, zCry, zMgoTop, Mat("MgO_powder"), w, cMgo);

  // Герметичная банка Al: обечайка и торцевая шайба; снизу банку закрывает
  // стеклянный световод, состыкованный с ФЭУ.
  Ring("AlCan_side", rMgo, rCan, zWinBot, zMgoTop, Mat("G4_Al"), w, cAl);
  Ring("AlCan_face", 0, rCan, zMgoTop, zCanTop, Mat("G4_Al"), w, cAl);
  Ring("Window", 0, rMgo, zWinBot, -zCry, Mat("G4_Pyrex_Glass"), w, cGlass);

  // Амортизатор «резина 2 мм»
  Ring("Rubber2_side", rCan, rRub, zPmtBot, zCanTop, Mat("G4_RUBBER_NATURAL"),
       w, cRub);
  Ring("Rubber2_face", 0, rRub, zCanTop, zRubTop, Mat("G4_RUBBER_NATURAL"),
       w, cRub);

  // Наружный корпус Al 1,5 мм — по всей длине устройства
  Ring("AlCase_side", rRub, rCase, zTail, zRubTop, Mat("G4_Al"), w, cAl);
  Ring("AlCase_face", 0, rCase, zRubTop, zCaseTop, Mat("G4_Al"), w, cAl);

  // Входной торец: воздух 0,5 и наружный протектор «резина 1 мм»
  Ring("FaceAir", 0, rCase, zCaseTop, zAirTop, Mat("G4_AIR"), w, cAl);
  Ring("Rubber1_face", 0, rCase, zAirTop, zFace, Mat("G4_RUBBER_NATURAL"),
       w, cRub);

  // ФЭУ: баллон Ø71 со стенкой 1,5 и вакуумом внутри
  const double rPmt = 0.5 * h.pmtDia;              // 35,50
  Ring("PMT_glass", rPmt - h.pmtGlass, rPmt, zPmtBot, zWinBot,
       Mat("G4_Pyrex_Glass"), w, cGlass);
  Ring("PMT_vac", 0, rPmt - h.pmtGlass, zPmtBot, zWinBot, Mat("G4_Galactic"),
       w, cGlass);

  // Делитель и плата
  Ring("Electronics", 0, rRub, zTail, zPmtBot, Mat("Electronics"), w, cEl);
}

// ---------------------------------------------------------------------------
// Экран-защита. Слои от полости наружу: Cu -> Cd -> Pb -> сталь.
void G1SDetector::BuildShield(G4LogicalVolume* w) {
  const ShieldGeom& s = fShield;

  const double rCav = 0.5 * s.cavityDia;                  // 100
  const double rCu = rCav + s.cu;                         // 101
  const double rCd = rCu + s.cd;                          // 102
  const double rPb = rCd + s.pb;                          // 152
  const double rSt = rPb + s.steel;                       // 155
  const double rBore = 0.5 * s.boreDia;                   // 41

  const double z0 = s.floorFromCryCentre;                 // -45  днище полости
  const double z1 = z0 + s.cavityH;                       // 145  потолок полости
  const double zCu0 = z0 - s.cu,       zCu1 = z1 + s.cu;
  const double zCd0 = zCu0 - s.cd,     zCd1 = zCu1 + s.cd;
  const double zPb0 = zCd0 - s.pb,     zPb1 = zCd1 + s.pb;
  const double zSt0 = zPb0 - s.steel,  zSt1 = zPb1 + s.steel;

  const G4Colour cCu(0.8, 0.5, 0.2), cCd(0.6, 0.6, 0.5), cPb(0.35, 0.35, 0.4),
      cSt(0.5, 0.55, 0.6);

  // Медь
  fCuLV = Ring("Cu_side", rCav, rCu, z0, z1, Mat("G4_Cu"), w, cCu);
  Ring("Cu_bottom", rBore, rCu, zCu0, z0, Mat("G4_Cu"), w, cCu);
  // Кадмий
  fCdLV = Ring("Cd_side", rCu, rCd, zCu0, zCu1, Mat("G4_Cd"), w, cCd);
  Ring("Cd_bottom", rBore, rCd, zCd0, zCu0, Mat("G4_Cd"), w, cCd);
  // Свинец
  fPbLV = Ring("Pb_side", rCd, rPb, zCd0, zCd1, Mat("G4_Pb"), w, cPb);
  Ring("Pb_bottom", rBore, rPb, zPb0, zCd0, Mat("G4_Pb"), w, cPb);
  // Сталь
  fSteelLV = Ring("St_side", rPb, rSt, zPb0, zPb1, Mat("G4_STAINLESS-STEEL"),
                  w, cSt);
  Ring("St_bottom", rBore, rSt, zSt0, zPb0, Mat("G4_STAINLESS-STEEL"), w, cSt);

  // Крышка. Для точечной геометрии 25 см она открыта — так и записано в
  // фоновом файле набора ЛСРМ (background_bg_2016_open_lid_point25cm).
  if (s.lidClosed) {
    Ring("Cu_top", 0, rCu, z1, zCu1, Mat("G4_Cu"), w, cCu);
    Ring("Cd_top", 0, rCd, zCu1, zCd1, Mat("G4_Cd"), w, cCd);
    Ring("Pb_top", 0, rPb, zCd1, zPb1, Mat("G4_Pb"), w, cPb);
    Ring("St_top", 0, rSt, zPb1, zSt1, Mat("G4_STAINLESS-STEEL"), w, cSt);
  }
}

// ---------------------------------------------------------------------------
// ПРОВЕРКА ПО МАССАМ. Паспорт ДЦКИ.412131.001 ПС, таблица 2.2, «не менее»:
//   свинец 165 кг, медь 1,6 кг, кадмий 1,2 кг, алюминий 1,7 кг.
// Массы считаются по формулам, а не по построенным телам, чтобы отчёт был
// доступен и без инициализации ядра.
void G1SDetector::ReportMasses() const {
  const ShieldGeom& s = fShield;
  const HeadGeom& h = fHead;
  const double rCav = 0.5 * s.cavityDia, rCu = rCav + s.cu, rCd = rCu + s.cd;
  const double rPb = rCd + s.pb, rSt = rPb + s.steel, rBore = 0.5 * s.boreDia;
  const double z0 = s.floorFromCryCentre, z1 = z0 + s.cavityH;
  const double hCu = z1 - z0, hCd = hCu + 2 * s.cu, hPb = hCd + 2 * s.cd;
  const double hSt = hPb + 2 * s.pb;

  const double vCu = CylCm3(rCav, rCu, hCu) + CylCm3(rBore, rCu, s.cu)
                   + (s.lidClosed ? CylCm3(0, rCu, s.cu) : 0.0);
  const double vCd = CylCm3(rCu, rCd, hCd) + CylCm3(rBore, rCd, s.cd)
                   + (s.lidClosed ? CylCm3(0, rCd, s.cd) : 0.0);
  const double vPb = CylCm3(rCd, rPb, hPb) + CylCm3(rBore, rPb, s.pb)
                   + (s.lidClosed ? CylCm3(0, rPb, s.pb) : 0.0);
  const double vSt = CylCm3(rPb, rSt, hSt) + CylCm3(rBore, rSt, s.steel)
                   + (s.lidClosed ? CylCm3(0, rSt, s.steel) : 0.0);
  const double vNaI = CylCm3(0, 0.5 * h.cryDia, h.cryLen);

  std::printf("\n--- массы построенных тел, кг (паспорт, «не менее») ---\n");
  std::printf("  свинец   %8.1f   (165)\n", vPb * 11.34 / 1000);
  std::printf("  медь     %8.2f   (1,6)\n", vCu * 8.96 / 1000);
  std::printf("  кадмий   %8.2f   (1,2)\n", vCd * 8.65 / 1000);
  std::printf("  сталь    %8.1f   (-)\n", vSt * 8.0 / 1000);
  std::printf("  NaI(Tl)  %8.3f   объём %.1f см³\n", vNaI * 3.667 / 1000, vNaI);
  std::printf("  полость  Ø%.0f x %.0f мм, крышка %s\n",
              s.cavityDia, s.cavityH, s.lidClosed ? "закрыта" : "открыта");
  if (fSampleVolumeCm3 > 0)
    std::printf("  проба    %8.1f см³  (цель %.0f, уровень %s)\n",
                fSampleVolumeCm3, fVessel.sampleCm3,
                fSampleFits ? "ниже крышки — ок" : "ВЫШЕ КРЫШКИ!");
  // Параметры, вынесенные для целевых прогонов, печатать обязательно: иначе
  // потом не восстановить, при каких значениях получен результат.
  std::printf("  MgO %.2f г/см³ ; колодец маринелли %.1f мм ; матрица %s\n",
              h.mgoDensity, fVessel.wellDepth, fVessel.sampleMatrix.c_str());
  std::printf("\n");
}
