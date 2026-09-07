#pragma once

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"

class NpsmBenchDetectorConstruction : public G4VUserDetectorConstruction {
public:
    static double gCrystalXMm;
    static double gCrystalYMm;
    static double gCrystalZMm;

    // Режим независимой сверки (линия Г): кристаллу задаются оптические
    // свойства, и свет считает ШТАТНЫЙ G4Scintillation по интегральной
    // кривой L(E), а не наша модель по шагам. Умолчание false — прежние
    // прогоны не затрагиваются.
    static bool gScintLight;

    static double SourceRadiusMm();
    static G4LogicalVolume* GetCrystalLogicalVolume();

    G4VPhysicalVolume* Construct() override;

private:
    static G4LogicalVolume* fgCrystalLV;
};
