You are a C++17 engineer working with Geant4 11.4.2 (MSVC). Write ONE complete header file `cabin_geom.hh` (all inline, no .cc). Output ONLY the code (no markdown fences, no prose). Comments in Russian, short.

## Purpose
Detector construction of a simplified aircraft cabin sector (Boeing 777-300ER) for a cosmic-ray simulation. Coordinates: X across the cabin (window on -X), Y along the fuselage axis, Z up. Origin = fuselage axis at the sector centre. All lengths in the parameter struct are in cm.

## Interface (exactly)
```cpp
#pragma once
#include "G4VUserDetectorConstruction.hh"
class G4LogicalVolume; class G4VPhysicalVolume;
struct CabinParams {           // all lengths in cm
  double halfLen   = 400.0;    // half length of the sector along Y
  double rOut      = 310.0;    // outer radius of the fuselage
  double rIn       = 293.0;    // inner radius of the wall (cabin half width)
  double frameSurf = 4.45;     // g/cm2 of the smeared skin + frame (Al 2024 equivalent)
  double blanketT  = 0.05;     // PET insulation film, rho 1.38
  double trimT     = 0.5;      // interior trim, rho 0.96
  double floorTop  = -55.0;    // Z of the floor top
  double floorT    = 1.3;      // Al floor (with beams, equivalent)
  double paxHalfX  = 280.0;    // half width of the people slab
  double paxTopZ   = 65.0;     // top of the people slab (Z)
  double paxRho    = 0.077;    // g/cm3 average density of the people layer
  double cargoTopZ = -100.0, cargoBotZ = -260.0, cargoHalfXTop = 225.0, cargoHalfXBot = 125.0, cargoRho = 0.122;
  double tubeX = -235.0, tubeY = 0.0, tubeZ = 20.0, tubeR = 14.0, tubeHalfL = 150.0;
  double airOutRho = 0.364e-3; // g/cm3 air at 11 km (fills the whole world, cabin included)
  double airInRho  = 0.90e-3;  // g/cm3 cabin air (used only inside the tube)
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
```
Define `Construct` after the class as `inline G4VPhysicalVolume* CabinDet::Construct()`.

## Materials (G4NistManager `nist = G4NistManager::Instance()`)
* `airOut = nist->BuildMaterialWithNewDensity("AirOut", "G4_AIR", p_.airOutRho * g/cm3)`; `airIn` the same with name "AirIn" and `airInRho`.
* `frame`: a new `G4Material("FrameAl", rhoFrame, 4)` with `rhoFrame = p_.frameSurf / (p_.rOut - p_.rIn + 0.0) * g/cm3` (frameSurf is g/cm2, thickness of the frame shell = rOut - rIn - trimT - blanketT, use that thickness in cm) and elements Al 0.935, Cu 0.044, Mg 0.015, Mn 0.006 (mass fractions, `AddElement(nist->FindOrBuildElement("Al"), 0.935)` etc.).
* `floorMat`: a second material "Al2024" with density 2.78 g/cm3 and the same four elements (same fractions).
* `pet = nist->FindOrBuildMaterial("G4_MYLAR")`.
* `trimMat = nist->BuildMaterialWithNewDensity("TrimPlastic", "G4_POLYETHYLENE", 0.96*g/cm3)`.
* `paxMat = nist->BuildMaterialWithNewDensity("PaxLayer", "G4_TISSUE_SOFT_ICRU-4", p_.paxRho * g/cm3)`.
* `cargoMat = nist->BuildMaterialWithNewDensity("Cargo", "G4_CELLULOSE_CELLOPHANE", p_.cargoRho * g/cm3)`.

## Geometry — logical volume names are EXACT (they are decoded by name later)
Define exactly (cm): `rTrimOut = rIn`, `rTrimIn = rIn - trimT`, `rBlanketOut = rIn + blanketT`, frame shell = [rBlanketOut, rOut]. Every solid uses `* cm`.
1. `World`: `G4Box` half sizes 560 cm in X, Y and Z, material `airOut`, LV name "World", physical volume placed with `new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, true)`; it is the return value.
2. Shells: `G4Tubs` with axis along local Z; every shell is placed directly in World with ONE shared rotation `rot` = `new G4RotationMatrix; rot->rotateX(90*deg);` so that the shell axis becomes the world Y axis, translation zero. Half length = halfLen. `Trim`: radii [rTrimIn, rTrimOut], trimMat. `Blanket`: radii [rTrimOut, rBlanketOut], pet. `Skin`: radii [rBlanketOut, rOut], frame. LV and PV names equal to the volume name.
3. `Floor`: G4Box half sizes X = `xf`, Y = halfLen, Z = floorT/2; centre Z = floorTop - floorT/2; material Al2024; `xf = std::sqrt(rTrimIn*rTrimIn - (floorTop - floorT)*(floorTop - floorT)) - 2.0`. No rotation.
4. `Pax`: G4Box half sizes X = paxHalfX, Y = halfLen, Z = (paxTopZ - floorTop)/2; centre Z = (paxTopZ + floorTop)/2; material paxMat; no rotation. If `paxHalfX*paxHalfX + paxTopZ*paxTopZ >= rTrimIn*rTrimIn` (corner would poke into the trim) call `G4Exception("CabinDet::Construct","Cabin01",FatalException,"Pax slab pokes into the wall")`.
5. `Cargo`: `G4Trd("Cargo", dx1, dx2, dy1, dy2, dz)` with dx1 = cargoHalfXBot (at -z), dx2 = cargoHalfXTop (at +z), dy1 = dy2 = halfLen, dz = (cargoTopZ - cargoBotZ)/2; centre Z = (cargoTopZ + cargoBotZ)/2; material cargoMat; no rotation. (G4Trd axis is Z, which is world Z, so no rotation is needed.)
6. `Tube` (recorder volume): `G4Tubs("Tube", 0, tubeR, tubeHalfL, 0, 360*deg)`, material `airIn`, LV name "Tube", a daughter of the `Pax` logical volume placed with rotation `rotT` = `new G4RotationMatrix; rotT->rotateX(90*deg);` (axis along world Y) and translation relative to the centre of Pax: (tubeX, tubeY, tubeZ - paxCentreZ) where paxCentreZ = (paxTopZ + floorTop)/2. Store the logical volume in `tubeLV_`. Check with `G4Exception` (FatalException) that the tube is fully inside Pax: `std::fabs(tubeX) + tubeR < paxHalfX`, `tubeHalfL + std::fabs(tubeY) <= halfLen`, `tubeZ - tubeR > floorTop`, `tubeZ + tubeR < paxTopZ`.
7. Every PV placement uses overlap checking `true`. Set every logical volume's vis attributes to `G4VisAttributes::GetInvisible()`.

## Requirements
* Include every Geant4 header used: G4NistManager, G4Material, G4Element, G4Box, G4Tubs, G4Trd, G4LogicalVolume, G4PVPlacement, G4RotationMatrix, G4SystemOfUnits, G4ThreeVector, G4VisAttributes, G4Exception, globals.hh, <cmath>.
* Convert every cm parameter with `* cm` where a Geant4 length is required; keep the parameter struct in plain cm doubles.
* No static mutable state.
