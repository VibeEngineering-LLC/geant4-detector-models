#pragma once
#include "G4VUserDetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4Exception.hh"
#include "globals.hh"
#include <cmath>
#include <string>
#include <initializer_list>
struct CabinParams {           // all lengths in cm
  double halfLen   = 400.0;    // half length of the sector along Y
  double rOut      = 197.5;    // outer radius of the fuselage (A320 family: 3.95 m)
  double rIn       = 185.0;    // inner radius of the wall (cabin width 3.70 m)
  double frameSurf = 1.85;     // g/cm2 skin+frames+stringers, Al 2024 equivalent; ESTIMATE +-30 % (skin 1.6 mm = 0.44 g/cm2, structure ~1.6 g/cm2 by fuselage mass / area)
  double blanketT  = 0.05;     // PET insulation film, rho 1.38
  double trimT     = 0.5;      // interior trim, rho 0.96
  double floorTop  = -60.0;    // Z of the floor top
  double floorT    = 1.1;      // Al floor (with beams, equivalent, ~3.0 g/cm2, ESTIMATE)
  double paxHalfX  = 170.0;    // half width of the people slab
  double paxTopZ   = 60.0;     // top of the people slab (Z)
  double paxRho    = 0.11;     // g/cm3 average density of the people layer (153 passengers x 84 kg / 117 m3, load factor 0.85)
  double cargoTopZ = -90.0, cargoBotZ = -165.0, cargoHalfXTop = 150.0, cargoHalfXBot = 75.0, cargoRho = 0.078;
  double tubeX = -150.0, tubeY = 0.0, tubeZ = 15.0, tubeR = 14.0, tubeHalfL = 150.0;
  double airOutRho = 0.364e-3; // g/cm3 air at 11 km (fills the whole world, cabin included)
  double airInRho  = 0.90e-3;  // g/cm3 cabin air (used only inside the tube)
  double worldHalf = 700.0;    // half size of the world box (>= 1.414 x source disk radius)
  bool   fuel      = false;    // вариант: крылья с топливом и двигатели (оценка их вклада)
  bool   empty     = false;    // проверочный режим: только мир и трубка, без конструкции (замыкающий тест двухэтапной схемы)
};
class CabinDet : public G4VUserDetectorConstruction {
 public:
  explicit CabinDet(const CabinParams& p = CabinParams()) : p_(p) {}
  G4VPhysicalVolume* Construct() override;
  const G4LogicalVolume* TubeLV() const { return tubeLV_; }
  const CabinParams& Params() const { return p_; }
 private:
  CabinParams p_;
  G4LogicalVolume* tubeLV_ = nullptr;
};
inline G4VPhysicalVolume* CabinDet::Construct() {
  auto nist = G4NistManager::Instance();
  auto airOut = nist->BuildMaterialWithNewDensity("AirOut", "G4_AIR", p_.airOutRho * g/cm3);
  auto airIn = nist->BuildMaterialWithNewDensity("AirIn", "G4_AIR", p_.airInRho * g/cm3);
  double rhoFrame = p_.frameSurf / (p_.rOut - p_.rIn - p_.blanketT) * g/cm3;   // толщина оболочки каркаса = rOut - rBlanketOut, см
  // Элементы с данными теплового рассеяния S(alpha,beta) (G4NDL ThermalScattering): имя G4Element обязано быть штатным TS_*,
  // иначе нейтроны рассеиваются на свободном газе. H в тканях и воде — ядро воды; в пластике, ПЭТ и целлюлозе — ядро полиэтилена (приближение).
  auto tsElem = [&](const char* tsName, const char* sym, int Z, std::initializer_list<int> Ns) {
    auto* e = new G4Element(tsName, sym, (G4int) Ns.size());
    double tot = 0.0;
    for (int N : Ns) tot += nist->GetIsotopeAbundance(Z, N);
    for (int N : Ns) {
      auto* iso = new G4Isotope(std::string(tsName) + "_iso" + std::to_string(N), Z, N, (double) N * g / mole);
      e->AddIsotope(iso, nist->GetIsotopeAbundance(Z, N) / tot);
    }
    return e;
  };
  auto* hWater = tsElem("TS_H_of_Water", "H", 1, {1, 2});
  auto* hPoly = tsElem("TS_H_of_Polyethylene", "H", 1, {1, 2});
  auto* alMetal = tsElem("TS_Aluminium_Metal", "Al", 13, {27});
  auto fillAl2024 = [&](G4Material* m) {
    m->AddElement(alMetal, 0.935);
    m->AddElement(nist->FindOrBuildElement("Cu"), 0.044);
    m->AddElement(nist->FindOrBuildElement("Mg"), 0.015);
    m->AddElement(nist->FindOrBuildElement("Mn"), 0.006);
  };
  auto frame = new G4Material("FrameAl", rhoFrame, 4);
  fillAl2024(frame);
  auto floorMat = new G4Material("Al2024", 2.78 * g/cm3, 4);
  fillAl2024(floorMat);
  auto pet = new G4Material("Mylar", 1.38 * g/cm3, 3);            // C10H8O4
  pet->AddElement(hPoly, 0.041960); pet->AddElement(nist->FindOrBuildElement("C"), 0.625017); pet->AddElement(nist->FindOrBuildElement("O"), 0.333023);
  auto trimMat = new G4Material("TrimPlastic", 0.96 * g/cm3, 2);   // полиэтилен CH2
  trimMat->AddElement(hPoly, 0.143711); trimMat->AddElement(nist->FindOrBuildElement("C"), 0.856289);
  auto paxMat = new G4Material("PaxLayer", p_.paxRho * g/cm3, 9);  // мягкая ткань ICRU-44 (состав), плотность слоя — усреднённая
  paxMat->AddElement(hWater, 0.102);
  const char* paxEl[8] = {"C", "N", "O", "Na", "P", "S", "Cl", "K"};
  const double paxW[8] = {0.143, 0.034, 0.708, 0.002, 0.003, 0.003, 0.002, 0.003};
  for (int k = 0; k < 8; ++k) paxMat->AddElement(nist->FindOrBuildElement(paxEl[k]), paxW[k]);
  auto cargoMat = new G4Material("Cargo", p_.cargoRho * g/cm3, 3);   // целлюлоза C6H10O5
  cargoMat->AddElement(hPoly, 0.062162); cargoMat->AddElement(nist->FindOrBuildElement("C"), 0.444450); cargoMat->AddElement(nist->FindOrBuildElement("O"), 0.493388);
  double rTrimOut = p_.rIn;
  double rTrimIn = p_.rIn - p_.trimT;
  double rBlanketOut = p_.rIn + p_.blanketT;
  auto worldLV = new G4LogicalVolume(new G4Box("World", p_.worldHalf*cm, p_.worldHalf*cm, p_.worldHalf*cm), airOut, "World");
  auto worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true);
  if (p_.empty) {   // пустой мир + трубка в воздухе: конструкции нет
    auto rotE = new G4RotationMatrix;
    rotE->rotateX(90*deg);
    auto tubeE = new G4LogicalVolume(new G4Tubs("Tube", 0, p_.tubeR*cm, p_.tubeHalfL*cm, 0, 360*deg), airIn, "Tube");
    new G4PVPlacement(rotE, G4ThreeVector(p_.tubeX*cm, p_.tubeY*cm, p_.tubeZ*cm), tubeE, "Tube", worldLV, false, 0, true);
    tubeLV_ = tubeE;
    return worldPV;
  }
  auto rot = new G4RotationMatrix;
  rot->rotateX(90*deg);
  auto trimLV = new G4LogicalVolume(new G4Tubs("Trim", rTrimIn*cm, rTrimOut*cm, p_.halfLen*cm, 0, 360*deg), trimMat, "Trim");
  new G4PVPlacement(rot, G4ThreeVector(), trimLV, "Trim", worldLV, false, 0, true);
  auto blanketLV = new G4LogicalVolume(new G4Tubs("Blanket", rTrimOut*cm, rBlanketOut*cm, p_.halfLen*cm, 0, 360*deg), pet, "Blanket");
  new G4PVPlacement(rot, G4ThreeVector(), blanketLV, "Blanket", worldLV, false, 0, true);
  auto skinLV = new G4LogicalVolume(new G4Tubs("Skin", rBlanketOut*cm, p_.rOut*cm, p_.halfLen*cm, 0, 360*deg), frame, "Skin");
  new G4PVPlacement(rot, G4ThreeVector(), skinLV, "Skin", worldLV, false, 0, true);
  auto cabAirLV = new G4LogicalVolume(new G4Tubs("CabinAir", 0, rTrimIn*cm, p_.halfLen*cm, 0, 360*deg), airIn, "CabinAir");
  new G4PVPlacement(rot, G4ThreeVector(), cabAirLV, "CabinAir", worldLV, false, 0, true);
  auto inCab = [&](G4LogicalVolume* lv, const char* nm, G4ThreeVector gpos) {   // daughter axis-aligned in global frame, mother CabinAir is rotated by rot
    new G4PVPlacement(G4Transform3D(*rot, (*rot) * gpos), lv, nm, cabAirLV, false, 0, true);
  };
  double xf = std::sqrt(rTrimIn*rTrimIn - (p_.floorTop - p_.floorT)*(p_.floorTop - p_.floorT)) - 2.0;
  auto floorLV = new G4LogicalVolume(new G4Box("Floor", xf*cm, p_.halfLen*cm, p_.floorT/2*cm), floorMat, "Floor");
  inCab(floorLV, "Floor", G4ThreeVector(0, 0, (p_.floorTop - p_.floorT/2)*cm));
  auto paxLV = new G4LogicalVolume(new G4Box("Pax", p_.paxHalfX*cm, p_.halfLen*cm, (p_.paxTopZ - p_.floorTop)/2*cm), paxMat, "Pax");
  if (p_.paxHalfX*p_.paxHalfX + p_.paxTopZ*p_.paxTopZ >= rTrimIn*rTrimIn) {
    G4Exception("CabinDet::Construct", "Cabin01", FatalException, "Pax slab pokes into the wall");
  }
  inCab(paxLV, "Pax", G4ThreeVector(0, 0, (p_.paxTopZ + p_.floorTop)/2*cm));
  auto cargoLV = new G4LogicalVolume(new G4Trd("Cargo", p_.cargoHalfXBot*cm, p_.cargoHalfXTop*cm, p_.halfLen*cm, p_.halfLen*cm, (p_.cargoTopZ - p_.cargoBotZ)/2*cm), cargoMat, "Cargo");
  inCab(cargoLV, "Cargo", G4ThreeVector(0, 0, (p_.cargoTopZ + p_.cargoBotZ)/2*cm));
  if (p_.fuel) {   // вариант: крылья (керосин 2 × 3,1 т + обшивка) и двигатели 2 × 2,3 т; прибор — над крылом
    auto ker = new G4Material("Kerosene", 0.80 * g/cm3, 2);
    ker->AddElement(hPoly, 0.154); ker->AddElement(nist->FindOrBuildElement("C"), 0.846);
    auto eng = nist->BuildMaterialWithNewDensity("EngineEq", "G4_Ti", 0.36 * g/cm3);   // 2,3 т на цилиндр R 90 см, L 250 см
    auto fuelLV = new G4LogicalVolume(new G4Box("Fuel", 325*cm, 200*cm, 7.5*cm), ker, "Fuel");
    auto wingLV = new G4LogicalVolume(new G4Box("WingSkin", 325*cm, 200*cm, 0.35*cm), floorMat, "WingSkin");
    auto engLV = new G4LogicalVolume(new G4Tubs("Engine", 0, 90*cm, 125*cm, 0, 360*deg), eng, "Engine");
    for (int sgn : {-1, 1}) {
      new G4PVPlacement(nullptr, G4ThreeVector(sgn*525*cm, 0, -170*cm), fuelLV, "Fuel", worldLV, false, 0, true);
      new G4PVPlacement(nullptr, G4ThreeVector(sgn*525*cm, 0, -162.15*cm), wingLV, "WingSkin", worldLV, false, 0, true);
      new G4PVPlacement(rot, G4ThreeVector(sgn*575*cm, 150*cm, -270*cm), engLV, "Engine", worldLV, false, 0, true);
    }
  }
  auto rotT = new G4RotationMatrix;
  rotT->rotateX(90*deg);
  double paxCentreZ = (p_.paxTopZ + p_.floorTop)/2;
  if (std::fabs(p_.tubeX) + p_.tubeR >= p_.paxHalfX || p_.tubeHalfL + std::fabs(p_.tubeY) > p_.halfLen ||
      p_.tubeZ - p_.tubeR <= p_.floorTop || p_.tubeZ + p_.tubeR >= p_.paxTopZ) {
    G4Exception("CabinDet::Construct", "Cabin02", FatalException, "Tube is not fully inside Pax");
  }
  auto tubeLV = new G4LogicalVolume(new G4Tubs("Tube", 0, p_.tubeR*cm, p_.tubeHalfL*cm, 0, 360*deg), airIn, "Tube");
  new G4PVPlacement(rotT, G4ThreeVector(p_.tubeX*cm, p_.tubeY*cm, (p_.tubeZ - paxCentreZ)*cm), tubeLV, "Tube", paxLV, false, 0, true);
  tubeLV_ = tubeLV;
  auto invisible = G4VisAttributes::GetInvisible();
  worldLV->SetVisAttributes(invisible);
  trimLV->SetVisAttributes(invisible);
  blanketLV->SetVisAttributes(invisible);
  skinLV->SetVisAttributes(invisible);
  floorLV->SetVisAttributes(invisible);
  paxLV->SetVisAttributes(invisible);
  cargoLV->SetVisAttributes(invisible);
  tubeLV->SetVisAttributes(invisible);
  return worldPV;
}
