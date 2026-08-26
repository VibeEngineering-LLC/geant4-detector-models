// Геометрия RadiaCode-110 — реализация. См. шапку RC110Detector.hh: плоское
// размещение, координаты уже в мировом кадре (источник — RC110_detector.gdml,
// финальная перепроверенная версия от 26.08.2026, не пересчитывать заново).
// Единицы величин в задаче — мм, здесь везде домножено на `mm` явно.
//
// Донор паттерна: RadiaCode-103/geometry/RCDetector.cc (хелпер Put(),
// стиль G4VisAttributes+SetForceSolid+G4PVPlacement).

#include "RC110Detector.hh"

#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

namespace {

G4Material* Mat(const G4String& n) {
  return G4NistManager::Instance()->FindOrBuildMaterial(n);
}

// Тот же хелпер, что в RadiaCode-103/geometry/RCDetector.cc: логический
// объём + непрозрачная закраска (SetForceSolid) + размещение в матери.
G4LogicalVolume* Put(G4VSolid* s, G4Material* m, const G4String& nm,
                     G4LogicalVolume* mother, const G4ThreeVector& pos,
                     const G4Colour& col) {
  auto* lv = new G4LogicalVolume(s, m, nm);
  auto* va = new G4VisAttributes(col);
  va->SetForceSolid(true);
  lv->SetVisAttributes(va);
  new G4PVPlacement(nullptr, pos, lv, nm, mother, false, 0, true);
  return lv;
}

}  // namespace

// ---------------------------------------------------------------------------
void RC110Detector::DefineMaterials() {
  auto* nist = G4NistManager::Instance();

  // PCB. NIST-имени G4_G10 в этой сборке нет (стандартный список NIST-
  // материалов Geant4 не включает трейд-марку G10/FR4 напрямую — тот же
  // вывод, к которому уже пришли в RadiaCode-103: там для платы построен
  // собственный материал "FR4"). Запасной вариант здесь — идентичный по
  // составу материал (тот же донор), а не второй источник истины.
  fMatPcb = nist->FindOrBuildMaterial("G4_G10");
  if (!fMatPcb) {
    auto* Si = nist->FindOrBuildElement("Si");
    auto* O  = nist->FindOrBuildElement("O");
    auto* C  = nist->FindOrBuildElement("C");
    auto* H  = nist->FindOrBuildElement("H");
    auto* Br = nist->FindOrBuildElement("Br");
    auto* fr4 = new G4Material("FR4_RC110", 1.85 * g / cm3, 5);
    fr4->AddElement(Si, 0.2818);
    fr4->AddElement(O,  0.3937);
    fr4->AddElement(C,  0.2264);
    fr4->AddElement(H,  0.0281);
    fr4->AddElement(Br, 0.0700);
    fMatPcb = fr4;
    G4cout << "[RC110Detector] G4_G10 не найден в NIST DB этой сборки -> "
              "PCB собран как FR4 (состав как в "
              "RadiaCode-103/geometry/RCDetector.cc::DefineMaterials)"
           << G4endl;
  }

  // Батарея — ПЛЕЙСХОЛДЕР для картинки (vis-only сборка, BeamOn не
  // вызывается, точный состав Li-Po не требуется здесь). Основной вариант —
  // G4_POLYVINYL_ACETATE, запасной — G4_POLYETHYLENE (заведомо есть в NIST DB,
  // используется также для ESR-плёнки ниже).
  fMatBattery = nist->FindOrBuildMaterial("G4_POLYVINYL_ACETATE");
  if (!fMatBattery) {
    G4cout << "[RC110Detector] G4_POLYVINYL_ACETATE не найден в NIST DB этой "
              "сборки -> батарея (плейсхолдер) собрана как G4_POLYETHYLENE"
           << G4endl;
    fMatBattery = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  }
}

// ---------------------------------------------------------------------------
// Все 15 компонентов — по спецификации из брифа задачи (см. GEANT4-MODEL.md
// / RC110_detector.gdml). Половинные размеры G4Box передаются в порядке
// X,Y,Z, координаты центра — уже в мировом кадре.
void RC110Detector::BuildDevice(G4LogicalVolume* world) {
  // --- Корпус: оболочка, стенка 1.5 мм по всем трём осям --------------------
  auto* caseOuter = new G4Box("caseOuter", 63.3 * mm, 17.05 * mm, 10.85 * mm);
  auto* caseInner = new G4Box("caseInner", 61.8 * mm, 15.55 * mm, 9.35 * mm);
  auto* caseShell = new G4SubtractionSolid("CaseShell", caseOuter, caseInner);
  Put(caseShell, Mat("G4_POLYCARBONATE"), "CaseShell", world,
      G4ThreeVector(0, 0, 0), G4Colour(0.55, 0.55, 0.58, 0.35));

  // --- Кристалл CsI(Tl), куб 14 мм ------------------------------------------
  Put(new G4Box("Crystal", 7 * mm, 7 * mm, 7 * mm), Mat("G4_CESIUM_IODIDE"),
      "Crystal", world, G4ThreeVector(-50.25 * mm, 0, 0),
      G4Colour(1.0, 0.84, 0.0));

  // --- SiPM ------------------------------------------------------------------
  Put(new G4Box("SiPM", 3 * mm, 3 * mm, 0.4 * mm), Mat("G4_Si"), "SiPM",
      world, G4ThreeVector(-50.25 * mm, 0, -7.4 * mm),
      G4Colour(0.0, 0.7, 0.9));

  // --- Капсула кристалла: оболочка, стенка 1.5 мм ----------------------------
  auto* capOuter = new G4Box("CapsuleOuter", 9 * mm, 9 * mm, 9 * mm);
  auto* capInner = new G4Box("CapsuleInner", 7.5 * mm, 7.5 * mm, 7.5 * mm);
  auto* capShell = new G4SubtractionSolid("CapsuleShell", capOuter, capInner);
  Put(capShell, Mat("G4_POLYCARBONATE"), "CapsuleShell", world,
      G4ThreeVector(-49.9 * mm, 0, -0.325 * mm),
      G4Colour(0.20, 0.20, 0.20, 0.5));

  // --- ESR-обёртка (3M ESR, 65 мкм = 0.0325 мм полутолщины на грань) --------
  // Пять сплошных граней + шестая (-Z, к SiPM) с прорезанным окном под него.
  auto* esrMat = Mat("G4_POLYETHYLENE");
  const G4Colour esrCol(0.90, 0.90, 0.92, 0.80);
  Put(new G4Box("ESR_px", 0.0325 * mm, 7.065 * mm, 7.065 * mm), esrMat,
      "ESR_px", world, G4ThreeVector(-43.2175 * mm, 0, 0), esrCol);
  Put(new G4Box("ESR_nx", 0.0325 * mm, 7.065 * mm, 7.065 * mm), esrMat,
      "ESR_nx", world, G4ThreeVector(-57.2825 * mm, 0, 0), esrCol);
  Put(new G4Box("ESR_py", 7.065 * mm, 0.0325 * mm, 7.065 * mm), esrMat,
      "ESR_py", world, G4ThreeVector(-50.25 * mm, 7.0325 * mm, 0), esrCol);
  Put(new G4Box("ESR_ny", 7.065 * mm, 0.0325 * mm, 7.065 * mm), esrMat,
      "ESR_ny", world, G4ThreeVector(-50.25 * mm, -7.0325 * mm, 0), esrCol);
  Put(new G4Box("ESR_pz", 7.065 * mm, 7.065 * mm, 0.0325 * mm), esrMat,
      "ESR_pz", world, G4ThreeVector(-50.25 * mm, 0, 7.0325 * mm), esrCol);

  auto* esrNzOuter =
      new G4Box("ESR_nz_outer", 7.065 * mm, 7.065 * mm, 0.0325 * mm);
  auto* esrNzWindow = new G4Box("ESR_nz_window", 3.2 * mm, 3.2 * mm, 0.5 * mm);
  auto* esrNzShell =
      new G4SubtractionSolid("ESR_nz", esrNzOuter, esrNzWindow);
  Put(esrNzShell, esrMat, "ESR_nz", world,
      G4ThreeVector(-50.25 * mm, 0, -7.0325 * mm), esrCol);

  // --- Основная плата ---------------------------------------------------------
  Put(new G4Box("PCB", 50.75 * mm, 14.5 * mm, 0.5 * mm), fMatPcb, "PCB",
      world, G4ThreeVector(9.85 * mm, 0, -5.45 * mm),
      G4Colour(0.1, 0.45, 0.2));

  // --- Батарея (плейсхолдер) ---------------------------------------------------
  Put(new G4Box("Battery", 32.5 * mm, 12.5 * mm, 4.8 * mm), fMatBattery,
      "Battery", world, G4ThreeVector(23.8 * mm, 0, 4.55 * mm),
      G4Colour(0.75, 0.75, 0.75));

  // --- Окно дисплея -------------------------------------------------------------
  Put(new G4Box("DisplayWindow", 18.25 * mm, 7 * mm, 0.5 * mm),
      Mat("G4_POLYCARBONATE"), "DisplayWindow", world,
      G4ThreeVector(-19.25 * mm, 0, -8.95 * mm),
      G4Colour(0.10, 0.10, 0.12, 0.40));

  // --- LCD-панель дисплея -------------------------------------------------------
  Put(new G4Box("DisplayLCD", 17 * mm, 6.5 * mm, 1.1 * mm), Mat("G4_Si"),
      "DisplayLCD", world, G4ThreeVector(-19.25 * mm, 0, -7.35 * mm),
      G4Colour(0.05, 0.05, 0.15));

  // --- USB-разъём -----------------------------------------------------------
  Put(new G4Box("USB", 3.75 * mm, 4.45 * mm, 1.6 * mm), Mat("G4_Cu"), "USB",
      world, G4ThreeVector(58.75 * mm, 0, -4.0 * mm),
      G4Colour(0.85, 0.55, 0.20));
}

// ---------------------------------------------------------------------------
G4VPhysicalVolume* RC110Detector::Construct() {
  DefineMaterials();

  auto* worldS = new G4Box("world", fWorldHalfX * mm, fWorldHalfY * mm,
                           fWorldHalfZ * mm);
  auto* worldLV = new G4LogicalVolume(worldS, Mat("G4_AIR"), "world");
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto* worldPV =
      new G4PVPlacement(nullptr, {}, worldLV, "world", nullptr, false, 0, true);

  BuildDevice(worldLV);

  G4cout << "\n=== RadiaCode-110: геометрия (vis-only, плоское размещение) ==="
         << G4endl;
  G4cout << "  15 компонентов, мир " << 2 * fWorldHalfX << " x "
         << 2 * fWorldHalfY << " x " << 2 * fWorldHalfZ << " мм" << G4endl;
  G4cout << G4endl;

  return worldPV;
}