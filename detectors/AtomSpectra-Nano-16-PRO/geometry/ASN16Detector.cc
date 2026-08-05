// Геометрия AtomSpectra Nano 16 PRO.
//
// ОТКУДА ЧИСЛА
// ------------
// Опубликованного чертежа прибора НЕТ. Размеры собраны из трёх источников,
// разобранных в reference/geometry-source.md:
//   1) слова оператора (29.07.2026, уточнения 05.08.2026) и его записка
//      «размеры.txt»: кристалл 18 x 15 x 57, корпус 86 x 42 x 25, обёртка
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
  auto* worldS = new G4Box("World", 150 * mm, 150 * mm, 250 * mm);
  auto* worldLV = new G4LogicalVolume(worldS, Mat("G4_AIR"), "World");
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto* world = new G4PVPlacement(nullptr, {}, worldLV, "World", nullptr,
                                  false, 0, true);
  const G4ThreeVector w0(0, 0, 0);

  // --- корпус: сплошная коробка, внутрь неё вставлена воздушная полость ----
  // Рёбра и бобышки профиля не воспроизводятся: стенки гладкие.
  auto* bodyLV = BoxAt("Body", -xBody, xBody, yBodyB, yBodyT, zBodyF, zBodyB,
                       Mat(g.matBody), worldLV, w0,
                       G4Colour(0.60, 0.65, 0.68), false);
  const G4ThreeVector bodyC(0, 0.5 * (yBodyB + yBodyT) * mm,
                            0.5 * (zBodyF + zBodyB) * mm);

  // Полость по всей длине экструзии; закрывают её крышки, а не сам профиль.
  auto* cavLV = BoxAt("Cavity", -xCav, xCav, yCavB, yFoilT, zBodyF, zBodyB,
                      Mat("G4_AIR"), bodyLV, bodyC,
                      G4Colour(0.85, 0.90, 0.95), false);
  const G4ThreeVector cavC(0, 0.5 * (yCavB + yFoilT) * mm,
                           0.5 * (zBodyF + zBodyB) * mm);

  // --- торцевые крышки (пластик) -------------------------------------------
  BoxAt("CapFront", -xCav, xCav, yCavB, yFoilT, zBodyF, zCapFi,
        Mat(g.matCap), cavLV, cavC, G4Colour(0.30, 0.32, 0.34));
  BoxAt("CapBack", -xCav, xCav, yCavB, yFoilT, zCapBi, zBodyB,
        Mat(g.matCap), cavLV, cavC, G4Colour(0.30, 0.32, 0.34));

  // --- плата: по дну полости во всю длину между крышками -------------------
  BoxAt("PCB", -xCav, xCav, yPcbB, yPcbT, zCapFi, zCapBi,
        Mat(g.matPcb), cavLV, cavC, G4Colour(0.18, 0.48, 0.29));

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

  return world;
}

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
  std::printf("  лицевой стек рабочей грани  %.3f г/см²\n", stack);
  std::printf("--- ВЕЩЕСТВА-ЗАМЕНИТЕЛИ (не подтверждены источником) ---\n");
  std::printf("  крышки: %s вместо ABS\n", gm.matCap.c_str());
  std::printf("  плата:  %s вместо FR4 (масса платы ЗАНИЖЕНА)\n",
              gm.matPcb.c_str());
  std::printf("  SiPM:   %s, голый кремний вместо сборки\n", gm.matSipm.c_str());
  std::printf("  корпус: %s при заявленном «алюминиевый сплав»\n",
              gm.matBody.c_str());
}
