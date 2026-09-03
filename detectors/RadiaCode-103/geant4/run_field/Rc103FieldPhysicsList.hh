// Физ-лист прогонов по полю: только ЭМ-физика.
//
// G4EmStandardPhysics_option4 — самый точный EM-конструктор Geant4 для задач
// низко-/среднеэнергетической гамма-дозиметрии.
//
// SetDefaultCutValue(0.05*mm) — здесь важна точность отклика в кристалле 10 мм
// и в тонких слоях обвязки (ESR-плёнка, SiPM, окно дисплея). В run_room/ стоял
// 1 мм, но там задача была транспорт через 80 см бетона — другая.
//
// Декей-физика НЕ регистрируется: первичные частицы этой программы — фотоны с
// разыгранной по CSV энергией, никаких ионов и распадов здесь нет.
//
// 03.09.2026: порог продукции и режим атомной деэкситации вынесены в параметры
// (ключи emcut= и deex= в main.cc) — проверяется утверждение «Geant4 упрощает
// физику на краях спектра ради скорости». Дефолты повторяют прежнее поведение
// бит в бит: cut 0.05 мм, deex=std.
//
// ⚠ Флаги деэкситации ставятся ПОСЛЕ конструктора G4EmStandardPhysics_option4:
// тот зовёт G4EmParameters::SetDefaults() и сбросил бы их обратно
// (G4EmStandardPhysics_option4.cc:113 в исходнике Geant4 11.2.1).
#pragma once

#include "G4VModularPhysicsList.hh"

#include <string>

class Rc103FieldPhysicsList : public G4VModularPhysicsList {
 public:
  // deexMode: "std"  — как было: только флуоресценция, включённая option4;
  //           "deex" — плюс Оже и DeexcitationIgnoreCut;
  //           "max"  — плюс PIXE (характеристический рентген вдоль торможения).
  explicit Rc103FieldPhysicsList(double cutMm = 0.05,
                                 const std::string& deexMode = "std");
  ~Rc103FieldPhysicsList() override = default;

  // Для шапки CSV: постановка обязана лежать В ФАЙЛЕ, а не в его имени.
  static double gCutMm;
  static std::string gDeexMode;
};