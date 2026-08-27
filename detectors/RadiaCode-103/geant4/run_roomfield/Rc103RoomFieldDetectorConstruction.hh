// Геометрия: прямоугольная комната (воздух) в окружении шести плит
// ограждающих конструкций — четыре кирпичные стены и два бетонных перекрытия.
// Размеры берутся из gRoom/gSlabs (Rc103RoomFieldGeometry.hh), в код не зашиты.
//
// МАТЕРИАЛ КИРПИЧА — ПРИБЛИЖЕНИЕ. NIST-состава керамического кирпича в Geant4
// нет, поэтому строится BuildMaterialWithNewDensity("G4_BRICK_approx",
// "G4_CONCRETE", rho_brick): элементный состав бетона при плотности кирпича.
// Обоснование: в диапазоне 100..3000 кэВ доминирует комптоновское рассеяние,
// сечение которого определяется электронной плотностью, а отношение Z/A у
// силикатов бетона и керамики практически совпадает; расхождение существенно
// только ниже ~100 кэВ, где работает фотоэффект и важен Z. Это приближение
// печатается и в stdout, и в шапку выходного CSV — молчать о нём нельзя.
//
// Скоринг: воздушный шар радиуса ballRMm с центром в точке наблюдения,
// вложен в объём воздуха комнаты (тот же материал, физику не меняет).
#pragma once

#include "G4VUserDetectorConstruction.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

class Rc103RoomFieldDetectorConstruction : public G4VUserDetectorConstruction {
 public:
  Rc103RoomFieldDetectorConstruction() = default;
  G4VPhysicalVolume* Construct() override;

  // Шар-скорер. Разрешается только после Construct(); main обязан проверить
  // ненулевой указатель ПОСЛЕ Initialize() и до BeamOn.
  static G4LogicalVolume* GetBallLogicalVolume() { return fgBallLV; }
  // Плотность бетона, взятая фактом у материала NIST (г/см3). 0 до Construct().
  static double GetConcreteDensityGCm3() { return fgRhoConcrete; }

 private:
  static G4LogicalVolume* fgBallLV;
  static double fgRhoConcrete;
};
