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
  // Матрица ОИСН-06 — лёгкая (насыпная ро около 0,6), состав по массе из
  // поля MATERIAL оригинальных .spe поверки 2024 (партия Th-232 420-17031):
  // в отличие от ОИСН-16 железа всего 15 %, основа углеродная. Нужна для
  // прямого прогона распада без пересчёта f(мю*ро*d) между матрицами.
  if (matrix == "OISN06") {
    auto* m = new G4Material(g4name, rho * g / cm3, 5);
    m->AddElement(nist->FindOrBuildElement("H"), 0.064);
    m->AddElement(nist->FindOrBuildElement("C"), 0.612);
    m->AddElement(nist->FindOrBuildElement("N"), 0.021);
    m->AddElement(nist->FindOrBuildElement("O"), 0.151);
    m->AddElement(nist->FindOrBuildElement("Fe"), 0.151);
    return m;
  }
  // Матрица источника Ra-226 420-7-18 (поверка 2016) — ТРЕТЬЯ, отдельная
  // матрица под тем же насыпным классом "ОИСН-06", НЕ идентичная блоку
  // выше: поле MATERIAL этого файла (Ro=0,6) даёт H 8,13% C 60,91%
  // N 7,14% O 23,81%, БЕЗ железа — а блок выше (партия Th-232 420-17031,
  // поверка 2024) даёт 15,1% Fe. Тот же паттерн взаимно нерешаемых
  // деклараций MATERIAL, что у ОИСН-16 (см. коммент выше про OISN16 /
  // OISN16_2016) — оператор 10.08.2026 распорядился НИ ОДНОЙ из двух
  // деклараций не верить, и матрицу строить по плотности (Ro=0,6 —
  // совпадает в обоих файлах, не спорна) + типовому составу отверждённой
  // эпоксидной смолы (гипотеза оператора: насыпка из чешуек/крошки
  // эпоксидки, плотность монолита ~1,1-1,2 г/см³ против заявленных 0,6 —
  // отношение ~0,5-0,55 physически согласуется с рыхлой засыпкой).
  //
  // Состав — НЕ измерение конкретного источника (паспорт состава не
  // публикуется ни одним изготовителем, проверено sci-search 10.08.2026),
  // а расчёт брутто-формулы отверждения DGEBA + ПЭПА-аналог (триэтилен-
  // тетрамин, TETA) по стехиометрии раскрытия эпоксидного кольца амином:
  // C 70,98% H 7,77% N 4,80% O 16,45% (masses, PubChem CID 2286/5565,
  // ✅ прямое чтение). Альтернатива — измеренный EPOTEK-301-1 (PDG, ✅):
  // C 68,96% H 6,99% N 0,89% O 23,15% — расходится с расчётом по N/O
  // (тип отвердителя неизвестен), но Z/A обоих вариантов совпадает с
  // точностью 0,7% (0,538 против 0,534) — для комптон-доминированного
  // переноса в этом диапазоне энергий выбор варианта на результат почти
  // не влияет (см. ra226-remarks.md §14/15, verified-facts.jsonl).
  if (matrix == "OISN06_epoxy") {
    auto* m = new G4Material(g4name, rho * g / cm3, 4);
    m->AddElement(nist->FindOrBuildElement("H"), 0.0777);
    m->AddElement(nist->FindOrBuildElement("C"), 0.7098);
    m->AddElement(nist->FindOrBuildElement("N"), 0.0480);
    m->AddElement(nist->FindOrBuildElement("O"), 0.1645);
    return m;
  }
  // Матрица ОИСН-16, ВТОРАЯ объявленная версия состава. В комплекте под одним
  // именем «ОИСН-16» и одной плотностью Ro=1,6 лежат два разных состава:
  //   поверка 2024: H 0,022  C 0,206  N 0,009  O 0,049  Fe 0,714  (OISN16)
  //   поверка 2016: H 0,024826 C 0,218471 N 0,022840 O 0,078451 Fe 0,655412
  // Источник запаян 17.09.2007 и не вскрывался, поэтому состав между
  // поверками измениться не мог — верна не более чем одна запись.
  // Запись 2016 года — ТОЧНОЕ целочисленное отношение 25 : 220 : 23 : 79 : 660
  // при сумме 1007 (сходится до 5e-16), то есть нормировка рецептуры с целыми
  // частями. Запись 2024 года — круглые тысячные 22 : 206 : 9 : 49 : 714 при
  // сумме 1000, то есть набрано руками. Одно из другого не получается ни
  // округлением, ни пересчётом атомных долей в массовые.
  // Ни та, ни другая первичным документом НЕ подтверждена: паспорт эталона
  // (certificates/, имя изготовителя не публикуется) называет матрицу только
  // именем, а поле MATERIAL в .spe заполняется выбором из справочника ПО и,
  // по указанию владельца комплекта, содержать ошибки может. Достоверны там
  // масса, объём и заметки, а не состав.
  // Поэтому оба состава заведены как отдельные матрицы: разница между ними —
  // не выбор «правильного», а слагаемое бюджета систематики, и меряется
  // прогоном, а не оценкой через Z_эфф.
  if (matrix == "OISN16_2016") {
    auto* m = new G4Material(g4name, rho * g / cm3, 5);
    m->AddElement(nist->FindOrBuildElement("H"), 25.0 / 1007.0);
    m->AddElement(nist->FindOrBuildElement("C"), 220.0 / 1007.0);
    m->AddElement(nist->FindOrBuildElement("N"), 23.0 / 1007.0);
    m->AddElement(nist->FindOrBuildElement("O"), 79.0 / 1007.0);
    m->AddElement(nist->FindOrBuildElement("Fe"), 660.0 / 1007.0);
    return m;
  }

  if (matrix != "OISN16") {
    G4Exception("G1SDetector::MakeMatrix", "g1s002", FatalException,
                ("неизвестная матрица: " + matrix).c_str());
  }
  // Состав по объявлению поверки 2024 (см. разбор выше).
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
// Сосуд Маринелли 1 л, надет на головку сверху. Габариты — чертёж изготовителя
// источников (см. VesselGeom в заголовке): наружный Ø154, H = 112, колодец Ø97
// глубиной 65. Сосуд выше номинального объёма, поэтому уровень засыпки
// вычисляется из целевых 1000 мл; над пробой до крышки — воздух. Дно колодца
// лежит на крышке Al торца (z = 41, стек торца от оператора).
//
// ПРОВЕРКА 1 (объём): вычисленный уровень должен лечь НИЖЕ крышки — иначе
//   целевой объём в сосуд не помещается (печатается в ReportMasses).
// ПРОВЕРКА 2 (эффективная толщина) — ЧТО ОНА НА САМОМ ДЕЛЕ ЗНАЧИТ.
//   Подгонка f(мю*ро*d) к отношению эффективностей двух плотностей даёт
//   d_эфф нашей модели. Сравнивать его с табличным ЛСРМ (26(2) мм) и с
//   Thick из .efa (31(2) мм) МОЖНО, но это сверка двух ПОДГОНОЧНЫХ величин,
//   а не проверка геометрии: у ЛСРМ d_эфф — толщина воображаемой пластинки
//   в модели «точечный источник за поглотителем», и она «подбирается так,
//   чтобы кривые эффективности при пересчёте на материал совпали» (§8.5.2).
//   Что это параметр, а не размер, видно из их же таблиц: один и тот же
//   сосуд при одних габаритах имеет d_эфф 17(2) мм в табл. 8-1 и 26(2) мм
//   в табл. 8-2 — поменялась формула, не сосуд. Совпадение по d_эфф
//   геометрию НЕ подтверждает; геометрию подтверждает чертёж.
void G1SDetector::BuildVessel(G4LogicalVolume* w) {
  const VesselGeom& v = fVessel;
  // 41,0 = крышка Al 2 мм: стек торца Al2/воздух1/банка0,5/MgO6 (оператор)
  const double zFace = 41.00;                 // наружная плоскость головки
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

  // Прежние габариты Маринелли, снятые с таблицы ЛСРМ «Прецизионные
  // измерения», с. 11. Оставлены ОТДЕЛЬНЫМ пресетом, чтобы разницу с
  // чертежом можно было измерить прогоном, а не рассуждением: на этих
  // числах посчитано всё, что сделано в проекте до 08.08.2026.
  // Оператор: «Мы брали размеры из статьи ЛСРМ Прецизионные измерения и
  // там не понятно что.» Достоверность ниже чертежа изготовителя.
  if (n == "marinelli_lsrm") {
    v.outerR   = 75.00;    // таблица ЛСРМ: внешний Ø150
    v.height   = 110.00;   // таблица ЛСРМ: H = 110
    v.wall     = 2.00;     // ДОПУЩЕНИЕ
    v.wellInR  = 40.00;    // Ø80 — ДОПУЩЕНИЕ «колодец по головке Ø78,3»
    v.wellDepth = 74.00;   // ДОПУЩЕНИЕ «колодец садится на всю головку»
    v.sampleCm3 = 1000.0;
    return v;
  }

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
  // 41,0 = крышка Al (см. BuildVessel)
  const double zFace = 41.00;                 // наружная плоскость головки
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
  // Стенка начинается ОТ ВЕРХА ДНА (zIn), а не от zBot: дно — сплошной диск
  // на всё сечение, и стенка от zBot перекрывала бы его в кольце
  // [rIn, outerR] на всю толщину дна. Поймано прогоном /geometry/test/run:
  // перекрытие 0,74 мм, оба тела полипропилен — на ослабление не влияло,
  // но геометрия была формально невалидной.
  Ring("V_side", rIn, v.outerR, zIn, zTop, pp, w, cPP);
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
  // Торец по оператору: банка 0,5 -> воздух 1,0 -> крышка Al 2,0. Резины на
  // торце нет; амортизатор 2 мм — только радиально.
  const double zAirTop  = zCanTop + h.faceAir;     // 39,00 — верх зазора
  const double zFace    = zAirTop + h.alCaseFace;  // 41,00 — наружная плоскость
  const double zRubTop  = zCanTop;                 // радиальная резина до банки
  const double zWinBot  = -zCry - h.window;        // -32,00 — низ слоя геля
  const double zPmtBot  = zWinBot - h.pmtLen;      // -152,00
  const double zTail    = zFace - h.unitLen;       // -274,00

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
  // Стык кристалл->ФЭУ: тонкий оптический гель (силикон; по составу и
  // плотности близок к резине, берём готовый материал). Толщина h.window.
  Ring("Window", 0, rMgo, zWinBot, -zCry, Mat("G4_RUBBER_NATURAL"), w, cGlass);

  // Амортизатор «резина 2 мм» — ТОЛЬКО радиально (оператор: на торце резины
  // нет вовсе)
  Ring("Rubber2_side", rCan, rRub, zPmtBot, zRubTop, Mat("G4_RUBBER_NATURAL"),
       w, cRub);

  // Наружный корпус: бок Al 1,5 по всей длине, крышка торца Al 2,0.
  // Между банкой и крышкой — воздушный зазор 1,0 (чертёж, подтверждено
  // оператором).
  Ring("AlCase_side", rRub, rCase, zTail, zAirTop, Mat("G4_Al"), w, cAl);
  // Торцевой зазор банка<->крышка: в центре воздух, по периферии — кольцевая
  // прокладка «резина 1 мм» (оператор, 29.07.2026: фиксирует банку от
  // осевого смещения; прежнее чтение ошибочно вешало её обёрткой на бок).
  // Разница радиусов кольца — 5-6 мм по оператору (принято 5,5), наружный
  // край у стенки корпуса. Центр торца перед кристаллом остаётся воздухом —
  // на фронтальный тракт по оси кольцо не влияет, оно прикрывает только
  // периферию торца (~1 мм резины на косых путях).
  const double rSeal = rRub - h.faceSealW;
  Ring("FaceAir", 0, rSeal, zCanTop, zAirTop, Mat("G4_AIR"), w, cAl);
  if (h.faceSeal > 0)
    Ring("FaceSeal_ring", rSeal, rRub, zCanTop, zAirTop,
         Mat("G4_RUBBER_NATURAL"), w, cRub);
  else
    Ring("FaceAir_rim", rSeal, rRub, zCanTop, zAirTop, Mat("G4_AIR"), w, cAl);
  Ring("AlCase_face", 0, rCase, zAirTop, zFace, Mat("G4_Al"), w, cAl);

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
// Профили тел вращения защиты. Вынесены отдельно, потому что по ним И строится
// геометрия, И считается масса для сверки с паспортом: держать два независимых
// описания одной формы — верный способ разойтись после первой же правки.
namespace {

// Свинцовое тело без крышки, снизу вверх. lin — суммарная толщина вкладышей.
int PbBodyProfile(const ShieldGeom& s, double* z, double* ri, double* ro) {
  const double lin = s.cu + s.cd;
  int n = 0;
  auto add = [&](double zz, double a, double b) {
    z[n] = zz; ri[n] = a; ro[n] = b; ++n;
  };
  add(s.zBottom,      s.rBore, s.rNeck2);
  add(s.zLedge2,      s.rBore, s.rNeck2);
  add(s.zLedge2,      s.rBore, s.rNeck1);
  add(s.zLedge1,      s.rBore, s.rNeck1);
  add(s.zLedge1,      s.rBore, s.rPbOut);
  add(s.zFloor - lin, s.rBore, s.rPbOut);
  add(s.zFloor - lin, s.rCav + lin, s.rPbOut);
  add(s.zBlockTop,    s.rCav + lin, s.rPbOut);
  return n;
}

// Крышка: пробка радиусом по полости плюс диск на всю ширину блока.
int LidProfile(const ShieldGeom& s, double* z, double* ri, double* ro) {
  const double lin = s.cu + s.cd;
  int n = 0;
  auto add = [&](double zz, double a, double b) {
    z[n] = zz; ri[n] = a; ro[n] = b; ++n;
  };
  add(s.zCeil + lin, 0, s.rCav + lin);
  add(s.zBlockTop,   0, s.rCav + lin);
  add(s.zBlockTop,   0, s.rPbOut);
  add(s.zLidTop,     0, s.rPbOut);
  return n;
}

// Объём тела вращения, заданного профилем (мм³ -> см³).
double ProfileCm3(int n, const double* z, const double* ri, const double* ro) {
  double v = 0;
  for (int i = 0; i + 1 < n; ++i) {
    const double h = z[i + 1] - z[i];
    if (h <= 0) continue;               // вертикальная ступень объёма не даёт
    // Внутри сегмента радиусы постоянны (профиль ступенчатый).
    v += CLHEP::pi * (ro[i] * ro[i] - ri[i] * ri[i]) * h;
  }
  return v / 1000.0;
}

G4LogicalVolume* Body(const G4String& nm, int n, const double* z,
                      const double* ri, const double* ro, G4Material* m,
                      G4LogicalVolume* mother, const G4Colour& col) {
  std::vector<double> zz(n), a(n), b(n);
  for (int i = 0; i < n; ++i) { zz[i] = z[i] * mm; a[i] = ri[i] * mm; b[i] = ro[i] * mm; }
  auto* s = new G4Polycone(nm, 0, CLHEP::twopi, n, zz.data(), a.data(), b.data());
  auto* lv = new G4LogicalVolume(s, m, nm);
  auto* va = new G4VisAttributes(col);
  va->SetForceSolid(true);
  lv->SetVisAttributes(va);
  new G4PVPlacement(nullptr, {}, lv, nm, mother, false, 0, true);
  return lv;
}

}  // namespace

// ---------------------------------------------------------------------------
// Экран-защита по разрезу рисунка 1.1 РЭ. Слои от полости наружу:
// Cu -> Cd -> Pb -> стальной кожух. Форма ступенчатая: широкий блок с полостью
// под сосуд сверху, свинцовая шейка с каналом под головку снизу.
void G1SDetector::BuildShield(G4LogicalVolume* w) {
  const ShieldGeom& s = fShield;
  const double lin = s.cu + s.cd;
  const double rCu = s.rCav + s.cu;      // наружная граница меди
  const double rCd = s.rCav + lin;       // наружная граница кадмия

  const G4Colour cCu(0.8, 0.5, 0.2), cCd(0.6, 0.6, 0.5), cPb(0.35, 0.35, 0.4),
      cSt(0.5, 0.55, 0.6);

  // Свинец: тело и крышка отдельными объёмами — крышка снимается.
  double z[12], ri[12], ro[12];
  int n = PbBodyProfile(s, z, ri, ro);
  fPbLV = Body("Pb_body", n, z, ri, ro, Mat("G4_Pb"), w, cPb);

  // Вкладыши полости: обечайка и дно. Стоят внутри свинца, а не съедают
  // полость: радиус полости rCav снят с разреза по видимой стенке.
  fCuLV = Ring("Cu_side", s.rCav, rCu, s.zFloor, s.zCeil, Mat("G4_Cu"), w, cCu);
  fCdLV = Ring("Cd_side", rCu, rCd, s.zFloor, s.zCeil, Mat("G4_Cd"), w, cCd);
  Ring("Cu_floor", s.rBore, rCd, s.zFloor - s.cu, s.zFloor, Mat("G4_Cu"), w, cCu);
  Ring("Cd_floor", s.rBore, rCd, s.zFloor - lin, s.zFloor - s.cu, Mat("G4_Cd"),
       w, cCd);

  // Свинцовая пробка канала. Закрывает единственный путь в полость, не
  // перекрытый свинцом: снизу канал открыт наружу, и без пробки фон входил бы
  // прямо к хвосту устройства. Указана оператором, на разрезе видна как
  // заштрихованный блок в канале у палубы.
  Ring("Pb_plug", 0, s.rBore, s.zPlugTop - s.plugThick, s.zPlugTop,
       Mat("G4_Pb"), w, cPb);

  // Кожух — сталь по всей высоте, включая юбку ниже блока (на разрезе она
  // пустая внутри: между шейкой и кожухом воздух).
  fSteelLV = Ring("St_skin", s.rPbOut, s.rPbOut + s.steel, s.zBottom, s.zLidTop,
                  Mat("G4_STAINLESS-STEEL"), w, cSt);

  // Крышка. Для точечной геометрии 25 см она открыта — так и записано в
  // фоновом файле набора ЛСРМ (background_bg_2016_open_lid_point25cm).
  if (s.lidClosed) {
    n = LidProfile(s, z, ri, ro);
    Body("Pb_lid", n, z, ri, ro, Mat("G4_Pb"), w, cPb);
    Ring("Cu_top", 0, rCd, s.zCeil, s.zCeil + s.cu, Mat("G4_Cu"), w, cCu);
    Ring("Cd_top", 0, rCd, s.zCeil + s.cu, s.zCeil + lin, Mat("G4_Cd"), w, cCd);
    Ring("St_top", 0, s.rPbOut + s.steel, s.zLidTop, s.zLidTop + s.steel,
         Mat("G4_STAINLESS-STEEL"), w, cSt);
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
  const double lin = s.cu + s.cd;
  const double rCu = s.rCav + s.cu, rCd = s.rCav + lin;

  double z[12], ri[12], ro[12];
  int n = PbBodyProfile(s, z, ri, ro);
  double vPb = ProfileCm3(n, z, ri, ro) + CylCm3(0, s.rBore, s.plugThick);
  if (s.lidClosed) {
    n = LidProfile(s, z, ri, ro);
    vPb += ProfileCm3(n, z, ri, ro);
  }

  const double hCav = s.zCeil - s.zFloor;
  const double vCu = CylCm3(s.rCav, rCu, hCav) + CylCm3(s.rBore, rCd, s.cu)
                   + (s.lidClosed ? CylCm3(0, rCd, s.cu) : 0.0);
  const double vCd = CylCm3(rCu, rCd, hCav) + CylCm3(s.rBore, rCd, s.cd)
                   + (s.lidClosed ? CylCm3(0, rCd, s.cd) : 0.0);
  const double vSt = CylCm3(s.rPbOut, s.rPbOut + s.steel, s.zLidTop - s.zBottom)
                   + (s.lidClosed ? CylCm3(0, s.rPbOut + s.steel, s.steel) : 0.0);
  const double vNaI = CylCm3(0, 0.5 * h.cryDia, h.cryLen);

  std::printf("\n--- массы построенных тел, кг (паспорт, «не менее») ---\n");
  std::printf("  свинец   %8.1f   (165)\n", vPb * 11.34 / 1000);
  std::printf("  медь     %8.2f   (1,6)\n", vCu * 8.96 / 1000);
  std::printf("  кадмий   %8.2f   (1,2)\n", vCd * 8.65 / 1000);
  std::printf("  сталь    %8.1f   (-)\n", vSt * 8.0 / 1000);
  std::printf("  NaI(Tl)  %8.3f   объём %.1f см³\n", vNaI * 3.667 / 1000, vNaI);
  std::printf("  полость  Ø%.1f x %.1f мм, канал Ø%.1f, крышка %s\n",
              2 * s.rCav, hCav, 2 * s.rBore,
              s.lidClosed ? "закрыта" : "открыта");
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
