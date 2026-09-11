#pragma once

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include <string>

class NpsmBenchDetectorConstruction : public G4VUserDetectorConstruction {
public:
    static double gCrystalXMm;
    static double gCrystalYMm;
    static double gCrystalZMm;

    // Форма и объём кристалла ДЛЯ ШАПКИ CSV. Стенд строит только G4Box, и
    // объём считался произведением XYZ — верно для него и НЕВЕРНО для
    // проектов с чужой геометрией: run_g1s_npsm ставит сюда 63x63x63, а
    // строит цилиндр Ø63x63, и шапка объявляла 250 см³ вместо 196,4
    // (найдено 10.09.2026 при сборке МЕТОДа). Постановка обязана лежать в
    // файле (#CFG-2), поэтому форма называется явно, а объём при
    // gCrystalVolumeCm3 > 0 берётся заданным, а не выводится из габаритов.
    static std::string gCrystalShape;   // "box" | "cylinder"
    static double gCrystalVolumeCm3;    // 0 = считать как X*Y*Z

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
