#pragma once
// Фантом тканей вместо прибора (режим CABIN_PHANTOM): шар мягкой ткани ICRU-4 радиусом 120 мм в воздухе.
// Слои учёта: кожа (0,07 мм от поверхности), глубина 10 мм (слой 1 мм), глубина 30 мм (слой 1 мм), весь шар.
#include <G4VUserDetectorConstruction.hh>
#include <G4NistManager.hh>
#include <G4Box.hh>
#include <G4Orb.hh>
#include <G4Sphere.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4SystemOfUnits.hh>
#include <vector>

struct PhantomDet : public G4VUserDetectorConstruction {
    static constexpr double kR = 120.0;                 // радиус шара, мм
    std::vector<G4LogicalVolume*> lv;                   // тела шара: 0 кожа, 1 ткань до 10 мм, 2 слой 10 мм, 3 ткань 10..30 мм, 4 слой 30 мм, 5 ядро
    std::vector<int> layerOf;                           // слой учёта тела: 0 кожа, 1 глубина 10 мм, 2 глубина 30 мм, -1 отдельно не учитывается
    std::vector<double> massKg;                         // масса тела
    G4VPhysicalVolume* Construct() override;
};
inline G4VPhysicalVolume* PhantomDet::Construct() {
    auto* nist = G4NistManager::Instance();
    auto* air = nist->FindOrBuildMaterial("G4_AIR");
    auto* tis = nist->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRU-4");
    auto* worldLV = new G4LogicalVolume(new G4Box("World", 300, 300, 1700), air, "World");
    auto* world = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, false);
    const double r[] = {kR, kR - 0.07, kR - 9.5, kR - 10.5, kR - 29.5, kR - 30.5};   // внешние границы тел, мм
    const int lay[] = {0, -1, 1, -1, 2, -1};
    const char* nm[] = {"Skin", "T1", "D10", "T2", "D30", "Core"};
    for (int i = 0; i < 6; ++i) {
        G4VSolid* sol = (i == 5) ? (G4VSolid*) new G4Orb(nm[i], r[i] * mm)
                                 : (G4VSolid*) new G4Sphere(nm[i], r[i + 1] * mm, r[i] * mm, 0, 360 * deg, 0, 180 * deg);
        auto* l = new G4LogicalVolume(sol, tis, nm[i]);
        new G4PVPlacement(nullptr, G4ThreeVector(), l, nm[i], worldLV, false, 0, false);
        lv.push_back(l); layerOf.push_back(lay[i]);
        massKg.push_back(sol->GetCubicVolume() * tis->GetDensity() / kg);
    }
    return world;
}
