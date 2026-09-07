#include "NpsmBenchDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4Region.hh"
#include "NpsmScintTable.hh"
#include <cstdio>
#include <iostream>
#include <vector>

bool NpsmBenchDetectorConstruction::gScintLight = false;
double NpsmBenchDetectorConstruction::gCrystalXMm = 102.0;
double NpsmBenchDetectorConstruction::gCrystalYMm = 102.0;
double NpsmBenchDetectorConstruction::gCrystalZMm = 406.0;

G4LogicalVolume* NpsmBenchDetectorConstruction::fgCrystalLV = nullptr;

double NpsmBenchDetectorConstruction::SourceRadiusMm() {
    // Расстояние от центра до диагонали кристалла + запас 30 мм
    double diag = 0.5 * std::sqrt(gCrystalXMm * gCrystalXMm + gCrystalYMm * gCrystalYMm + gCrystalZMm * gCrystalZMm) + 30.0;
    return diag;
}

G4LogicalVolume* NpsmBenchDetectorConstruction::GetCrystalLogicalVolume() {
    return fgCrystalLV;
}

G4VPhysicalVolume* NpsmBenchDetectorConstruction::Construct() {
    G4NistManager* nist = G4NistManager::Instance();

    // Определение материалов
    G4Material* galactic = nist->FindOrBuildMaterial("G4_Galactic");
    G4Material* sodiumIodide = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");

    // Создание мира
    double worldSize = SourceRadiusMm() + 50.0;
    G4Box* worldSolid = new G4Box("World", worldSize * mm, worldSize * mm, worldSize * mm);
    G4LogicalVolume* worldLV = new G4LogicalVolume(worldSolid, galactic, "World");
    G4VPhysicalVolume* worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0);

    // Оптические свойства кристалла — только для режима light=scint (линия Г).
    // Наша модель считает свет весом на КАЖДОМ шаге по dE/dx; штатный
    // G4Scintillation берёт готовую интегральную кривую L(E) и разыгрывает
    // число фотонов сам. Это два независимых пути к одной величине: если они
    // сойдутся, значит совпадает не набор параметров, а реализация.
    //
    // RESOLUTIONSCALE = 0 — собственное разрешение выключено: сравниваются
    // средние, а не ширины, и лишний источник разброса тут только мешает.
    // Фотоны не транспортируются: их убивает NpsmBenchStackingAction, считая.
    if (gScintLight) {
        auto* mpt = new G4MaterialPropertiesTable();
        std::vector<G4double> eArr(npsm_scint::kEnergyMeV,
                                   npsm_scint::kEnergyMeV + npsm_scint::kN);
        std::vector<G4double> yArr(npsm_scint::kYieldPhotons,
                                   npsm_scint::kYieldPhotons + npsm_scint::kN);
        for (auto& v : eArr) v *= MeV;
        // ⚠ В режиме «выход по типу частицы» Geant4 требует кривую для
        // КАЖДОГО типа, который выделяет энергию, иначе роняет прогон
        // (`Scint01: no correct entry in MaterialPropertiesTable`). На
        // 662 кэВ хватало одних электронов, а на 1173 кэВ начинается
        // рождение пар — и прогон упал (07.09.2026). Тяжёлым частицам
        // задаётся та же кривая: физически для них она неверна, но в наших
        // постановках их вклад отсутствует, а падать прогон не должен.
        // Если появится задача с протонами или альфа — кривые нужны свои.
        for (const char* name : {"ELECTRONSCINTILLATIONYIELD",
                                 "PROTONSCINTILLATIONYIELD",
                                 "DEUTERONSCINTILLATIONYIELD",
                                 "TRITONSCINTILLATIONYIELD",
                                 "ALPHASCINTILLATIONYIELD",
                                 "IONSCINTILLATIONYIELD"}) {
            mpt->AddProperty(name, eArr, yArr, true, true);
        }
        // Один фотон на единицу — абсолют уже зашит в таблицу.
        mpt->AddConstProperty("SCINTILLATIONYIELD", 1.0 / MeV, true);
        mpt->AddConstProperty("RESOLUTIONSCALE", 0.0, true);
        mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 250.0 * ns, true);
        mpt->AddConstProperty("SCINTILLATIONYIELD1", 1.0, true);
        // Спектр испускания: одна линия 415 нм (максимум NaI(Tl)). На счёт
        // фотонов не влияет, но без него G4Scintillation не работает.
        const G4double photonE[2] = {2.95 * eV, 3.02 * eV};
        const G4double spec[2] = {1.0, 1.0};
        const G4double rindex[2] = {1.85, 1.85};
        mpt->AddProperty("SCINTILLATIONCOMPONENT1",
                         std::vector<G4double>(photonE, photonE + 2),
                         std::vector<G4double>(spec, spec + 2), true, true);
        mpt->AddProperty("RINDEX", std::vector<G4double>(photonE, photonE + 2),
                         std::vector<G4double>(rindex, rindex + 2), true, true);
        sodiumIodide->SetMaterialPropertiesTable(mpt);
        std::printf("ScintLight: ELECTRONSCINTILLATIONYIELD задан, узлов %d, "
                    "L(1 МэВ)=%.0f фотонов\n",
                    npsm_scint::kN, npsm_scint::kYieldPhotons[npsm_scint::kN / 2]);
    }

    // Создание кристалла
    double halfX = gCrystalXMm / 2.0;
    double halfY = gCrystalYMm / 2.0;
    double halfZ = gCrystalZMm / 2.0;
    G4Box* crystalSolid = new G4Box("Crystal", halfX * mm, halfY * mm, halfZ * mm);
    G4LogicalVolume* crystalLV = new G4LogicalVolume(crystalSolid, sodiumIodide, "Crystal");
    new G4PVPlacement(nullptr, G4ThreeVector(), crystalLV, "Crystal", worldLV, false, 0, true);

    fgCrystalLV = crystalLV;

    // Регион кристалла. Нужен для линии В: деэкситация настраивается на
    // G4Region, а не на объём (civanch, форум, тема 2056), и это позволяет
    // держать PIXE/Оже там, где считается сигнал, не платя за них в остальной
    // геометрии. Регион создаётся ВСЕГДА, но пока никто не назначил ему
    // настройки, поведение не меняется — прежние прогоны воспроизводятся.
    {
        auto* reg = new G4Region("CrystalRegion");
        reg->AddRootLogicalVolume(crystalLV);
        std::printf("Region: CrystalRegion <- Crystal\n");
    }

    // Вывод информации о геометрии
    double volumeCm3 = (gCrystalXMm * gCrystalYMm * gCrystalZMm) / 1000.0;
    double density = sodiumIodide->GetDensity() / (g/cm3);
    std::printf("Crystal: %gx%gx%g mm, volume=%.2f cm3, density=%.4f g/cm3, world half-size=%.1f mm\n",
                gCrystalXMm, gCrystalYMm, gCrystalZMm, volumeCm3, density, worldSize);

    return worldPV;
}
