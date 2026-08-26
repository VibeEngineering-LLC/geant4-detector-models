// RC103 case from Rc-103.stl via CADMesh (runtime tessellated solid).
// Requires: https://github.com/christopherpoole/CADMesh
//
// CMake (after find_package or add_subdirectory CADMesh):
//   target_link_libraries(your_target PRIVATE CADMesh::CADMesh)
//
// ВАЖНОЕ ОТЛИЧИЕ ОТ RC-110: Rc-103.stl экспортирован из Blender уже в системе
// устройства и в миллиметрах — bbox [-61.5..61.5] x [-17..17] x [-8.75..8.75],
// центр в нуле. Поэтому remap осей НЕ нужен, в отличие от RC110CADMeshCase, где
// применялся сдвиг 122.5/122.5/63.3.
// Число треугольников ЗАВИСИТ от того, какой файл резолвится (см.
// stl_to_gdml.py::resolve_stl() — берёт Rc-103.stl первым, .bak fallback):
// Rc-103.stl = 102252 tri (5112684 B), Rc-103.stl.bak = 18516 tri (925884 B) —
// проверено по факту 26.08.2026 (заголовок STL совпал с размером файла).
// Прежнее «2792 треугольника» в этом комментарии было неверным числом.
//
// Система координат: X — длина (USB в +X), Y — ширина, Z — толщина (−Z = лицо).

#ifndef RC103_CADMESH_CASE_HH
#define RC103_CADMESH_CASE_HH

#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class RC103CADMeshCase
{
public:
  // stlPath: абсолютный путь к Rc-103.stl
  // caseABS: материал корпуса (например, из G4NistManager)
  static G4LogicalVolume* BuildLogicalVolume(const G4String& stlPath,
                                             G4Material* caseABS,
                                             const G4ThreeVector& extraOffset =
                                               G4ThreeVector(0, 0, 0));

  // Габарит STL по проектному чертежу, мм. Для sanity-check после загрузки.
  static constexpr G4double kCaseLengthMM = 123.0;
  static constexpr G4double kCaseWidthMM  = 34.0;
  static constexpr G4double kCaseThickMM  = 17.5;
};

#endif
