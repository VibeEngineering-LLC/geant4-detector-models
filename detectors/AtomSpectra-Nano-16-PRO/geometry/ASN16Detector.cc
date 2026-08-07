// Геометрия AtomSpectra Nano 16 PRO.
//
// ОТКУДА ЧИСЛА
// ------------
// Опубликованного чертежа прибора НЕТ. Размеры собраны из трёх источников,
// разобранных в reference/geometry-source.md:
//   1) слова оператора (29.07.2026, уточнения 05.08.2026) и его записка
//      «размеры.txt»: кристалл 18 x 15 x 60 (до 06.08.2026 стояло 57 —
//      размер соседней модели Nano 15), корпус 86 x 42 x 25, обёртка
//      ПТФЭ 1,0 + Al 0,1, кристалл упёрт в переднюю стенку;
//   2) чертёж профиля экструзии линейки Nano (Nano 5 PRO): рабочая стенка
//      1,20, дно 2,05, боковая (39,50 − 35,60)/2 = 1,95;
//   3) фотография открытого торца: плата по дну во всю длину, кристалл
//      опирается на неё; кристалл с обёрткой примерно по центру полости.
//
// ПОВЕРКА ПО ФОТО. Масштаб взят по внутренней ширине полости (38,10 мм).
// Кристалл с обёрткой: по фото ~20,5 x 17,3 мм, по принятой геометрии
// 20,20 x 17,20. Два независимых размера сошлись — довод в пользу того, что и
// обёртка 1,0 + 0,1, и внутренняя ширина полости приняты верно. Замер по
// растру без масштабной линейки, считать проверкой порядка, не чертежом.
//
// ЧТО НЕ МОДЕЛИРУЕТСЯ
//   Рёбра и винтовые бобышки профиля экструзии (в модели стенки гладкие —
//   это ЗАНИЖАЕТ массу алюминия и, значит, ослабление сбоку). Разъём и
//   элементы на плате. Светодиод. Тонкая проводка. Оптический контакт
//   кристалл—SiPM (моделируется прямым касанием без слоя геля).
//
// ЧЕГО НЕТ В ИСТОЧНИКАХ (перечислено в Nano16Geom со звёздочкой): сплав
// корпуса, толщина крышек и платы, состав крышек, платы и сборки SiPM,
// положение кристалла по X и Z внутри полости.

#include "ASN16Detector.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

#include <cstdio>

namespace {

G4Material* Mat(const G4String& n) {
  return G4NistManager::Instance()->FindOrBuildMaterial(n);
}

// Параллелепипед по мировым границам (x0..x1, y0..y1, z0..z1), помещённый в
// мать, у которой центр стоит в мировой точке `mother0`. Порядок границ
// произволен. Возвращает логический объём.
G4LogicalVolume* BoxAt(const G4String& nm,
                       double x0, double x1, double y0, double y1,
                       double z0, double z1,
                       G4Material* m, G4LogicalVolume* mother,
                       const G4ThreeVector& mother0, const G4Colour& col,
                       bool solidVis = true) {
  const double dx = std::abs(x1 - x0), dy = std::abs(y1 - y0),
               dz = std::abs(z1 - z0);
  auto* s = new G4Box(nm, 0.5 * dx * mm, 0.5 * dy * mm, 0.5 * dz * mm);
  auto* lv = new G4LogicalVolume(s, m, nm);
  auto* va = new G4VisAttributes(col);
  va->SetForceSolid(solidVis);
  lv->SetVisAttributes(va);
  const G4ThreeVector c(0.5 * (x0 + x1) * mm, 0.5 * (y0 + y1) * mm,
                        0.5 * (z0 + z1) * mm);
  new G4PVPlacement(nullptr, c - mother0, lv, nm, mother, false, 0, true);
  return lv;
}

double BoxCm3(double dx, double dy, double dz) {   // мм -> см³
  return dx * dy * dz / 1000.0;
}

}  // namespace

// ---------------------------------------------------------------------------
void ASN16Detector::DefineMaterials() {
  auto* nist = G4NistManager::Instance();
  nist->FindOrBuildMaterial("G4_AIR");
  nist->FindOrBuildMaterial("G4_CESIUM_IODIDE");
  nist->FindOrBuildMaterial("G4_TEFLON");
  nist->FindOrBuildMaterial(fGeom.matBody);
  nist->FindOrBuildMaterial(fGeom.matCap);
  nist->FindOrBuildMaterial(fGeom.matPcb);
  nist->FindOrBuildMaterial(fGeom.matSipm);
}

// --- опорные плоскости -------------------------------------------------------
// Начало координат — центр кристалла, поэтому все выводится от его габаритов.
double ASN16Detector::CrystalFrontZ() const { return +0.5 * fGeom.cryZ; }
double ASN16Detector::CrystalTopY() const { return +0.5 * fGeom.cryY; }

double ASN16Detector::BottomFaceY() const {
  return WorkFaceY() - fGeom.bodyY;
}

double ASN16Detector::FrontFaceZ() const {
  // обёртка (ПТФЭ + фольга) упёрта в крышку, крышка — в торец корпуса
  return CrystalFrontZ() + fGeom.alFoil + fGeom.ptfe + fGeom.wCap;
}

double ASN16Detector::WorkFaceY() const {
  return CrystalTopY() + fGeom.alFoil + fGeom.ptfe + fGeom.wFront;
}

// ---------------------------------------------------------------------------
G4VPhysicalVolume* ASN16Detector::Construct() {
  DefineMaterials();
  const Nano16Geom& g = fGeom;

  // --- границы по Y: снизу вверх, от дна корпуса к рабочей стенке ----------
  // ПОРЯДОК СЛОЁВ ОБЁРТКИ (ОПЕРАТОР, 05.08.2026): ПТФЭ лежит НА КРИСТАЛЛЕ,
  // алюминиевая фольга — НА ПТФЭ. До этой правки вложенность была обратной
  // (фольга на кристалле, ПТФЭ снаружи); суммарная поверхностная плотность
  // стека от порядка не зависит, поэтому пропускание было верным, но
  // рассеяние и флуоресценция фольги считались не с того места. Ошибку нашёл
  // независимый аудит, разрешил оператор.
  const double yCryT = +0.5 * g.cryY,            yCryB = -0.5 * g.cryY;
  const double yPtfeT = yCryT + g.ptfe,          yPtfeB = yCryB - g.ptfe;
  const double yFoilT = yPtfeT + g.alFoil,       yFoilB = yPtfeB - g.alFoil;
  const double yBodyT = yFoilT + g.wFront;       // наружная рабочая поверхность
  const double yBodyB = yBodyT - g.bodyY;
  const double yCavB  = yBodyB + g.wBot;         // внутренняя поверхность дна
  const double yPcbT  = yFoilB, yPcbB = yFoilB - g.pcbT;   // плата под обёрткой

  // --- границы по X --------------------------------------------------------
  const double xCry = 0.5 * g.cryX;
  const double xPtfe = xCry + g.ptfe, xFoil = xPtfe + g.alFoil;
  const double xBody = 0.5 * g.bodyX, xCav = xBody - g.wSide;

  // --- границы по Z --------------------------------------------------------
  // +Z — К ПЕРЕДНЕМУ ТОРЦУ И К ИСТОЧНИКУ (так же ориентирован Гамма-1С).
  // Это не косметика: конус GPS «/gps/ang/maxtheta» разыгрывается вокруг −Z,
  // поэтому источник обязан стоять на +Z, иначе конус светит в пустоту, а
  // прогон при этом честно отработает и запишет почти пустой спектр.
  const double zCryF = +0.5 * g.cryZ, zCryB = -0.5 * g.cryZ;
  const double zPtfeF = zCryF + g.ptfe;
  const double zFoilF = zPtfeF + g.alFoil;
  const double zCapFi = zFoilF;                  // внутренняя грань крышки
  const double zBodyF = zCapFi + g.wCap;         // наружный торец корпуса
  const double zBodyB = zBodyF - g.bodyZ;
  const double zCapBi = zBodyB + g.wCap;
  const double zSipmB = zCryB - g.sipmT;

  // --- мир -----------------------------------------------------------------
  // Запас должен вмещать точечный источник на 10 см от торца плюс воздух
  // вокруг него: рассеяние в воздухе перед прибором — часть измерения.
  // Мир расширяется только вместе с комнатой: пол на 700 мм ниже в коробку
  // 300 x 300 x 500 не помещается, а держать большой мир постоянно незачем.
  const double wHX = g.room ? g.roomHalf : 150.0;
  const double wHZ = g.room ? g.roomHalf : 250.0;
  auto* worldS = new G4Box("World", wHX * mm, wHX * mm, wHZ * mm);
  auto* worldLV = new G4LogicalVolume(worldS, Mat("G4_AIR"), "World");
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto* world = new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr,
                                  false, 0, true);
  const G4ThreeVector w0(0, 0, 0);

  // --- корпус: сплошная коробка, внутрь неё вставлена воздушная полость ----
  // Рёбра и бобышки профиля не воспроизводятся: стенки гладкие.
  // Корпус и крышки можно заменить вакуумом командой /asn16/caseOff on —
  // ТОЛЬКО ДИАГНОСТИКА. Разность двух прогонов одной ревизии показывает
  // полный вклад алюминиевого корпуса: и ослабление в пучке (дно 1,50 мм),
  // и рассеяние с боковых стенок и крышек, лежащих вне пучка. Обёртка
  // (AlFoil) при этом остаётся — она часть сборки кристалла, не корпуса.
  const G4String matCaseNow = g.caseOff ? G4String("G4_Galactic") : g.matBody;
  const G4String matCapNow = g.caseOff ? G4String("G4_Galactic") : g.matCap;
  auto* bodyLV = BoxAt("Body", -xBody, xBody, yBodyB, yBodyT, zBodyF, zBodyB,
                       Mat(matCaseNow), worldLV, w0,
                       G4Colour(0.60, 0.65, 0.68), false);
  const G4ThreeVector bodyC(0, 0.5 * (yBodyB + yBodyT) * mm,
                            0.5 * (zBodyF + zBodyB) * mm);

  // Полость по всей длине экструзии; закрывают её крышки, а не сам профиль.
  auto* cavLV = BoxAt("Cavity", -xCav, xCav, yCavB, yFoilT, zBodyF, zBodyB,
                      Mat("G4_AIR"), bodyLV, bodyC,
                      G4Colour(0.85, 0.90, 0.95), false);
  const G4ThreeVector cavC(0, 0.5 * (yCavB + yFoilT) * mm,
                           0.5 * (zBodyF + zBodyB) * mm);

  // --- торцевые крышки (АЛЮМИНИЙ 1,50 мм, ОПЕРАТОР 06.08.2026) ------------
  auto* capFrontLV = BoxAt("CapFront", -xCav, xCav, yCavB, yFoilT,
                           zBodyF, zCapFi, Mat(matCapNow), cavLV, cavC,
                           G4Colour(0.55, 0.60, 0.64));
  // ВХОДНОЕ ОКНО: фрезеровка передней крышки до 0,60 мм напротив кристалла
  // (ОПЕРАТОР, 06.08.2026). Сделана не булевой операцией, а вставкой
  // ВОЗДУШНОГО КАРМАНА в тело крышки: результат тот же — в пучке остаётся
  // g.wCapWin алюминия, — а тесселированных и булевых тел этот контур
  // избегает намеренно.
  // Карман примыкает к ВНУТРЕННЕЙ грани крышки; окно шире кристалла на
  // g.capWinPad с каждой стороны и УЖЕ обёртки (20,00 x 17,00 против
  // 20,20 x 17,20), поэтому обёртка упирается в ободок, а не проваливается
  // в выборку.
  if (g.capWindow) {
  const double xWin = 0.5 * g.cryX + g.capWinPad;
  const double yWinT = +0.5 * g.cryY + g.capWinPad;
  const double yWinB = -0.5 * g.cryY - g.capWinPad;
  const G4ThreeVector capFrontC(0, 0.5 * (yCavB + yFoilT) * mm,
                                0.5 * (zBodyF + zCapFi) * mm);
  BoxAt("CapWindow", -xWin, xWin, yWinB, yWinT,
        zCapFi, zCapFi + (g.wCap - g.wCapWin),
        Mat("G4_AIR"), capFrontLV, capFrontC, G4Colour(0.90, 0.95, 1.00),
        false);
  }
  BoxAt("CapBack", -xCav, xCav, yCavB, yFoilT, zCapBi, zBodyB,
        Mat(matCapNow), cavLV, cavC, G4Colour(0.30, 0.32, 0.34));

  // --- плата: по дну полости во всю длину между крышками -------------------
  // Наполнение платы (медь дорожек, припой, корпуса компонентов) вводится
  // ЭФФЕКТИВНЫМИ СПЛОШНЫМИ СЛОЯМИ поверх диэлектрика, со стороны кристалла:
  // реальная разводка неизвестна, поэтому моделируется не рисунок, а
  // эквивалент по массе. Порядок снизу вверх — диэлектрик, медь, припой,
  // корпуса компонентов; последние ближе всего к кристаллу, как на плате.
  // Стек ВСТРАИВАЕТСЯ В ТОЛЩИНУ платы, а не надстраивается над ней: иначе
  // кристалл сдвинулся бы вверх, и вместе с ним поехала бы вся геометрия
  // замера. Диэлектрику остаётся pcbT минус сумма трёх слоёв.
  const double tCu = g.pcbCuT, tSol = g.pcbSnPbT, tCmp = g.pcbCompT;
  const double tDie = g.pcbT - (tCu + tSol + tCmp);
  if (tDie <= 0.0) {
    G4Exception("ASN16Detector::Construct", "ASN16_PCB", FatalException,
                "слои платы толще самой платы");
  }
  const double yDieT = yPcbB + tDie;
  const double yCuT = yDieT + tCu, ySolT = yCuT + tSol;
  BoxAt("PCB", -xCav, xCav, yPcbB, yDieT, zCapFi, zCapBi,
        Mat(g.matPcb), cavLV, cavC, G4Colour(0.18, 0.48, 0.29));
  BoxAt("PcbCu", -xCav, xCav, yDieT, yCuT, zCapFi, zCapBi,
        Mat("G4_Cu"), cavLV, cavC, G4Colour(0.72, 0.45, 0.20));
  BoxAt("PcbSolder", -xCav, xCav, yCuT, ySolT, zCapFi, zCapBi,
        MakeSolder(), cavLV, cavC, G4Colour(0.62, 0.64, 0.67));
  BoxAt("PcbComp", -xCav, xCav, ySolT, yPcbT, zCapFi, zCapBi,
        Mat("G4_ALUMINUM_OXIDE"), cavLV, cavC, G4Colour(0.85, 0.85, 0.80));
  std::printf("  ПЛАТА: %.2f мм всего = диэлектрик %.3f %s + медь %.3f + "
              "припой %.3f %s + корпуса %.3f Al2O3 (эффективные сплошные "
              "слои, ДОПУЩЕНИЕ)\n",
              g.pcbT, tDie, g.matPcb.c_str(), tCu, tSol,
              g.pcbSolderPb ? "Sn63Pb37" : "SAC305", tCmp);

  // --- обёртка и кристалл ---------------------------------------------------
  // Вложение фольга -> ПТФЭ -> кристалл: ПТФЭ прилегает к кристаллу, фольга
  // лежит на ПТФЭ (ОПЕРАТОР). Задняя грань всех трёх в одной плоскости
  // z = zCryB, поэтому сзади обёртки нет — там SiPM.
  auto* foilLV = BoxAt("AlFoil", -xFoil, xFoil, yFoilB, yFoilT, zFoilF, zCryB,
                       Mat(g.matBody), cavLV, cavC,
                       G4Colour(0.78, 0.80, 0.82), false);
  const G4ThreeVector foilC(0, 0.5 * (yFoilB + yFoilT) * mm,
                            0.5 * (zFoilF + zCryB) * mm);

  auto* ptfeLV = BoxAt("PTFE", -xPtfe, xPtfe, yPtfeB, yPtfeT, zPtfeF, zCryB,
                       Mat("G4_TEFLON"), foilLV, foilC,
                       G4Colour(0.95, 0.95, 0.93), false);
  const G4ThreeVector ptfeC(0, 0.5 * (yPtfeB + yPtfeT) * mm,
                            0.5 * (zPtfeF + zCryB) * mm);

  fCrystalLV = BoxAt("Crystal", -xCry, xCry, yCryB, yCryT, zCryF, zCryB,
                     Mat("G4_CESIUM_IODIDE"), ptfeLV, ptfeC,
                     G4Colour(0.85, 0.68, 0.24));

  // --- SiPM на задней грани кристалла --------------------------------------
  BoxAt("SiPM", -xCry, xCry, yCryB, yCryT, zCryB, zSipmB,
        Mat(g.matSipm), cavLV, cavC, G4Colour(0.29, 0.44, 0.65));

  // --- ИСТОЧНИК: пачка электродов WT-20 под прибором -----------------------
  // Строится только по команде /asn16/wt20 on. Прибор лежит на пенале ДНОМ
  // (сторона платы), поэтому пачка стоит по −Y, а не по +Z, как точечный
  // источник опорного замера. Ослабление на пути: стенка пенала, дно корпуса
  // 1,50 мм Al, плата 1,60 мм — и только затем кристалл. Окна фрезеровки на
  // этом пути НЕТ: оно в передней ТОРЦЕВОЙ крышке.
  fPackW = 0.0;
  if (g.wt20 && g.wt20N > 0) {
    G4Material* wt20 = MakeWT20();
    const double r = 0.5 * g.wt20D;
    const double halfSpan = 0.5 * (g.wt20N - 1) * g.wt20Pitch;   // до центров
    // Прибор стоит на ЧЕТЫРЁХ НОЖКАХ, а не лежит дном: между дном корпуса и
    // крышкой пенала остаётся воздух высотой feetH (ОПЕРАТОР, 06.08.2026).
    const double yTopFace = BottomFaceY() - g.feetH;
    const double yInT = yTopFace - g.wt20Wall;      // внутренняя грань крышки
    const double yRod = yInT - g.wt20Gap - r;       // ось стержней
    const double yInB = yRod - r - g.wt20Gap;
    const double yOutB = yInB - g.wt20Wall;
    // Пачка центрирована под КОРПУСОМ («вдоль электродов, по центру»), а не
    // под кристаллом: кристалл упёрт в передний торец, и центр корпуса от него
    // смещён. Опорой служит то, на чём прибор лежит.
    const double zMid = 0.5 * (zBodyF + zBodyB);
    const double xHalfIn = halfSpan + r + g.wt20Gap;
    const double zHalfIn = 0.5 * g.wt20L + g.wt20Gap;

    auto* caseLV = BoxAt("Case", -(xHalfIn + g.wt20Wall), xHalfIn + g.wt20Wall,
                         yOutB, yTopFace,
                         zMid - (zHalfIn + g.wt20Wall),
                         zMid + (zHalfIn + g.wt20Wall),
                         Mat(g.matCase), worldLV, w0,
                         G4Colour(0.35, 0.55, 0.85), false);
    const G4ThreeVector caseC(0, 0.5 * (yOutB + yTopFace) * mm, zMid * mm);
    auto* airLV = BoxAt("CaseAir", -xHalfIn, xHalfIn, yInB, yInT,
                        zMid - zHalfIn, zMid + zHalfIn,
                        Mat("G4_AIR"), caseLV, caseC,
                        G4Colour(0.9, 0.9, 0.9), false);
    const G4ThreeVector airC(0, 0.5 * (yInB + yInT) * mm, zMid * mm);

    // Все стержни носят ОДНО имя «Rod»: розыгрыш источника задаётся командой
    // /gps/pos/confine Rod, и она сравнивает ИМЯ объёма в точке — значит одна
    // команда покрывает всю пачку. Разные имена потребовали бы десяти
    // источников GPS с ручными весами.
    auto* rodS = new G4Tubs("Rod", 0, r * mm, 0.5 * g.wt20L * mm, 0, 360 * deg);
    auto* rodLV = new G4LogicalVolume(rodS, wt20, "Rod");
    auto* rva = new G4VisAttributes(G4Colour(0.75, 0.20, 0.20));
    rva->SetForceSolid(true);
    rodLV->SetVisAttributes(rva);
    for (int i = 0; i < g.wt20N; ++i) {
      const double x = -halfSpan + i * g.wt20Pitch;
      new G4PVPlacement(nullptr,
                        G4ThreeVector(x * mm, yRod * mm, zMid * mm) - airC,
                        rodLV, "Rod", airLV, false, i, true);
    }
    // Плёнка осаждения дочерних торона: воздушная оболочка 25 мкм поверх
    // каждого стержня, имя RodSkin — под /gps/pos/confine RodSkin
    // (поверхностный источник, задача №9: горб ХРИ дочерних 75-95 кэВ
    // объёмным источником не воспроизводится — самопоглощение в вольфраме).
    // Материал — тот же воздух полости: на перенос слой не влияет и
    // существует только как область розыгрыша.
    auto* skinS = new G4Tubs("RodSkin", r * mm, (r + 0.025) * mm,
                             0.5 * g.wt20L * mm, 0, 360 * deg);
    auto* skinLV = new G4LogicalVolume(skinS, Mat("G4_AIR"), "RodSkin");
    skinLV->SetVisAttributes(G4VisAttributes::GetInvisible());
    for (int i = 0; i < g.wt20N; ++i) {
      const double x = -halfSpan + i * g.wt20Pitch;
      new G4PVPlacement(nullptr,
                        G4ThreeVector(x * mm, yRod * mm, zMid * mm) - airC,
                        skinLV, "RodSkin", airLV, false, i, true);
    }
    // Единица массы здесь пишется как `gram`, а не `g`: имя `g` в этой функции
    // занято ссылкой на геометрию, и `g/cm3` молча означало бы не то.
    const double vRod = 3.14159265358979 * r * r * g.wt20L / 1000.0;   // см³
    const double rhoRod = wt20->GetDensity() / (gram / cm3);
    fPackW = g.wt20N * vRod * rhoRod;
    fPackYRod = yRod;
    std::printf("--- ПАЧКА WT-20 ---\n");
    std::printf("  %d стержней %.2f x %.1f мм, шаг %.2f мм\n",
                g.wt20N, g.wt20D, g.wt20L, g.wt20Pitch);
    std::printf("  сплав %.3f г/см3 (W + %.1f %% масс. ThO2), объём пачки "
                "%.4f см3, масса %.2f г\n",
                rhoRod, g.wt20ThO2, g.wt20N * vRod, fPackW);
    // Печатаются ТРИ плоскости, а не две. Пока прибор «лежал», дно корпуса и
    // верх пенала совпадали, и одной строки хватало; с ножками это разные
    // уровни, и печатать верх пенала под именем дна корпуса значит подсунуть
    // в журнал верное число под чужим именем.
    std::printf("  ось стержней y = %+.2f, верх пенала y = %+.2f, дно корпуса "
                "y = %+.2f (ножки %.2f), просвет пенал-стержень %.2f мм\n",
                yRod, yTopFace, BottomFaceY(), g.feetH, yTopFace - yRod - r);
    std::printf("  на пути от стержня до кристалла: пенал %.2f мм %s, дно "
                "корпуса %.2f мм %s, плата %.2f мм %s\n",
                g.wt20Wall, g.matCase.c_str(), g.wBot, g.matBody.c_str(),
                g.pcbT, g.matPcb.c_str());

    // --- ножки прибора ------------------------------------------------------
    // Четыре пуговички по углам корпуса. На ослабление в пучке они не влияют
    // (стоят по углам, кристалл над серединой), но ПОДНИМАЮТ прибор — а это
    // прямо входит в телесный угол.
    {
      G4Material* rub = Mat(g.matFeet);
      if (!rub) rub = Mat(g.matCase);          // если марки нет в базе NIST
      const double fx = 0.5 * g.bodyX - g.feetInset - 0.5 * g.feetXY;
      const double fz = zMid + 0.5 * g.bodyZ - g.feetInset - 0.5 * g.feetXY;
      const double fz2 = zMid - 0.5 * g.bodyZ + g.feetInset + 0.5 * g.feetXY;
      for (int sx = -1; sx <= 1; sx += 2)
        for (int iz = 0; iz < 2; ++iz) {
          const double zc = iz ? fz : fz2;
          BoxAt("Foot", sx * fx - 0.5 * g.feetXY, sx * fx + 0.5 * g.feetXY,
                yTopFace, yTopFace + g.feetH,
                zc - 0.5 * g.feetXY, zc + 0.5 * g.feetXY,
                rub, worldLV, w0, G4Colour(0.25, 0.25, 0.28));
        }
      std::printf("  НОЖКИ: 4 x %.1f x %.1f мм, высота %.2f мм (%s) — прибор "
                  "поднят над пеналом\n", g.feetXY, g.feetXY, g.feetH,
                  rub ? rub->GetName().c_str() : "?");
    }

    // --- столешница под пеналом ------------------------------------------
    if (g.table) {
      G4Material* wood = G4Material::GetMaterial("Table", false);
      if (!wood)
        wood = G4NistManager::Instance()->BuildMaterialWithNewDensity(
            "Table", "G4_CELLULOSE_BUTYRATE", g.tableRho * gram / cm3);
      const double yT = yOutB - g.tableGap;
      BoxAt("Table", -0.5 * g.tableXY, 0.5 * g.tableXY, yT - g.tableT, yT,
            zMid - 0.5 * g.tableXY, zMid + 0.5 * g.tableXY,
            wood, worldLV, w0, G4Colour(0.55, 0.40, 0.25), false);
      std::printf("  СТОЛ: %.0f мм %s, плита %.0f x %.0f мм, верх y = %+.2f\n",
                  g.tableT, "целлюлоза (замена дерева)", g.tableXY, g.tableXY,
                  yT);
    }
  }

  // --- КОМНАТА: бетонный пол ------------------------------------------------
  // Ставится независимо от пачки: рассеиватель нужен и в опорных постановках.
  // Плита кроет мир по X и Z целиком — приближение бесконечной плоскости, у
  // которой возвращающийся поток от расстояния почти не зависит. Стен и
  // потолка НЕТ: пол — ближайшая и самая массивная поверхность, и вводить
  // остальное имеет смысл только если его вклада не хватит.
  if (g.room) {
    const double yTop = -g.floorDrop;
    if (yTop - g.floorT < -g.roomHalf) {
      G4Exception("ASN16Detector::Construct", "ASN16_ROOM", FatalException,
                  "пол не помещается в мир: увеличить roomHalf");
    }
    BoxAt("Floor", -0.98 * g.roomHalf, 0.98 * g.roomHalf,
          yTop - g.floorT, yTop,
          -0.98 * g.roomHalf, 0.98 * g.roomHalf,
          Mat("G4_CONCRETE"), worldLV, w0,
          G4Colour(0.45, 0.45, 0.45), false);
    std::printf("  КОМНАТА: пол бетон %.0f мм, верх y = %+.0f мм под прибором, "
                "мир расширен до +-%.0f мм (ДОПУЩЕНИЕ: высота стола не "
                "измерена)\n", g.floorT, yTop, g.roomHalf);
  }

  return world;
}

// ---------------------------------------------------------------------------
// Припой платы. Состав по прибору НЕ установлен, поэтому даётся выбор из двух
// типовых, и оба собираются из элементов по массовым долям:
//   Sn63Pb37 — эвтектический оловянно-свинцовый, 8,40 г/см³. Свинец в нём и
//              есть причина, по которой припой вообще попал в модель: Pb Kα1
//              74,97 и Kβ 84,9 кэВ ложатся в измеренный горб 75-95 кэВ, а
//              плата стоит вплотную к кристаллу;
//   SAC305   — Sn96,5 Ag3,0 Cu0,5, 7,38 г/см³, бессвинцовый: у него в этой
//              области линий нет (Sn Kα 25,3), и разница двух прогонов прямо
//              меряет вклад свинца припоя.
// Плотности — паспортные значения припоев, не выведены из правила смеси.
G4Material* ASN16Detector::MakeSolder() {
  const G4String name = fGeom.pcbSolderPb ? "SolderSnPb" : "SolderSAC";
  if (auto* have = G4Material::GetMaterial(name, false)) return have;
  auto* nist = G4NistManager::Instance();
  G4Material* m;
  if (fGeom.pcbSolderPb) {
    m = new G4Material(name, 8.40 * gram / cm3, 2);
    m->AddElement(nist->FindOrBuildElement("Sn"), 0.63);
    m->AddElement(nist->FindOrBuildElement("Pb"), 0.37);
  } else {
    m = new G4Material(name, 7.38 * gram / cm3, 3);
    m->AddElement(nist->FindOrBuildElement("Sn"), 0.965);
    m->AddElement(nist->FindOrBuildElement("Ag"), 0.030);
    m->AddElement(nist->FindOrBuildElement("Cu"), 0.005);
  }
  return m;
}

// ---------------------------------------------------------------------------
// Сплав электрода WT-20: вольфрам с 2 % масс. диоксида тория. Плотность НЕ
// взята из таблицы, а посчитана по правилу смеси из плотностей компонентов
// базы NIST: 1/rho = w1/rho1 + w2/rho2. Реальный спечённый пруток может быть
// плотнее или рыхлее — это допущение, и оно печатается при каждом запуске.
G4Material* ASN16Detector::MakeWT20() {
  if (auto* have = G4Material::GetMaterial("WT20", false)) return have;
  auto* nist = G4NistManager::Instance();
  G4Material* w = nist->FindOrBuildMaterial("G4_W");
  if (!w) {
    G4Exception("ASN16Detector::MakeWT20", "ASN16_MAT", FatalException,
                "в базе NIST нет G4_W");
    return nullptr;
  }
  // Диоксид тория собирается ИЗ ЭЛЕМЕНТОВ: в базе NIST его нет (есть оксид
  // урана и диоксид плутония, тория — нет). Стехиометрия ThO2 задана числом
  // атомов, а не долями массы: доли пришлось бы считать самому, и ошибка в
  // них выглядела бы как верный материал.
  G4Material* tho2 = G4Material::GetMaterial("ThO2", false);
  if (!tho2) {
    tho2 = new G4Material("ThO2", fGeom.rhoThO2 * gram / cm3, 2);
    tho2->AddElement(nist->FindOrBuildElement("Th"), 1);
    tho2->AddElement(nist->FindOrBuildElement("O"), 2);
  }
  const double f = fGeom.wt20ThO2 / 100.0;
  const double rw = w->GetDensity() / (gram / cm3);
  const double rt = tho2->GetDensity() / (gram / cm3);
  const double rho = 1.0 / ((1.0 - f) / rw + f / rt);
  auto* m = new G4Material("WT20", rho * gram / cm3, 2);
  m->AddMaterial(w, 1.0 - f);
  m->AddMaterial(tho2, f);
  return m;
}

double ASN16Detector::PackMassG() const { return fPackW; }

// ---------------------------------------------------------------------------
void ASN16Detector::ReportPlanes() const {
  const Nano16Geom& gm = fGeom;
  std::printf("--- ОПОРНЫЕ ПЛОСКОСТИ (мм, начало — центр кристалла) ---\n");
  std::printf("  наружный торец корпуса   z = %+8.2f  <- от него 10 см в замере\n",
              FrontFaceZ());
  std::printf("  передняя грань кристалла z = %+8.2f  (разница %.2f мм)\n",
              CrystalFrontZ(), CrystalFrontZ() - FrontFaceZ());
  std::printf("  наружная рабочая стенка  y = %+8.2f\n", WorkFaceY());
  std::printf("  рабочая грань кристалла  y = %+8.2f  (разница %.2f мм)\n",
              CrystalTopY(), WorkFaceY() - CrystalTopY());
  std::printf("  площадь торца  %.2f см², рабочей грани %.2f см², отношение %.2f\n",
              gm.cryX * gm.cryY / 100.0, gm.cryX * gm.cryZ / 100.0,
              gm.cryZ / gm.cryY);
}

void ASN16Detector::ReportMasses() const {
  const Nano16Geom& gm = fGeom;
  auto rho = [](const G4String& n) {
    G4Material* m = Mat(n);
    return m ? m->GetDensity() / (g / cm3) : 0.0;
  };
  const double rhoCsI = rho("G4_CESIUM_IODIDE");
  const double vCry = BoxCm3(gm.cryX, gm.cryY, gm.cryZ);
  std::printf("--- МАССЫ (г) ---\n");
  std::printf("  CsI(Tl)  %7.2f см³ x %.3f = %7.1f г\n", vCry, rhoCsI,
              vCry * rhoCsI);
  const double stack = gm.wFront / 10.0 * rho(gm.matBody)
                     + gm.ptfe / 10.0 * rho("G4_TEFLON")
                     + gm.alFoil / 10.0 * rho(gm.matBody);
  std::printf("  лицевой стек рабочей грани  %.6f г/см²\n", stack);
  // Торцевой стек печатается наравне с лицевым: в настоящей постановке в ПУЧКЕ
  // стоит именно он, а до 06.08.2026 не печатался нигде. Порядок слоёв тот же
  // (ПТФЭ на кристалле, фольга на ПТФЭ), но снаружи вместо стенки 1,20 стоит
  // торцевая крышка.
  const double stackWrap = gm.ptfe / 10.0 * rho("G4_TEFLON")
                         + gm.alFoil / 10.0 * rho(gm.matBody);
  const double stackWin = stackWrap + gm.wCapWin / 10.0 * rho(gm.matCap);
  const double stackRim = stackWrap + gm.wCap / 10.0 * rho(gm.matCap);
  if (gm.capWindow) {
    std::printf("  ТОРЦЕВОЙ стек В ОКНЕ        %.6f г/см²  <- В ПУЧКЕ\n",
                stackWin);
    std::printf("  ТОРЦЕВОЙ стек вне окна      %.6f г/см²\n", stackRim);
    std::printf("  окно фрезеровки  %.2f x %.2f мм, остаток крышки %.2f мм\n",
                gm.cryX + 2 * gm.capWinPad, gm.cryY + 2 * gm.capWinPad,
                gm.wCapWin);
  } else {
    std::printf("  ТОРЦЕВОЙ стек, крышка СПЛОШНАЯ  %.6f г/см²  <- В ПУЧКЕ\n",
                stackRim);
    std::printf("  ФРЕЗЕРОВКА ОТКЛЮЧЕНА (/asn16/capWindow off)\n");
  }
  std::printf("--- ВЕЩЕСТВА-ЗАМЕНИТЕЛИ (не подтверждены источником) ---\n");
  std::printf("  крышки: %s, фрезеровка до %.2f мм напротив кристалла\n",
              gm.matCap.c_str(), gm.wCapWin);
  std::printf("  плата:  %s вместо FR4 (масса платы ЗАНИЖЕНА)\n",
              gm.matPcb.c_str());
  std::printf("  SiPM:   %s, голый кремний вместо сборки\n", gm.matSipm.c_str());
  std::printf("  корпус: %s при заявленном «алюминиевый сплав»\n",
              gm.matBody.c_str());
}
