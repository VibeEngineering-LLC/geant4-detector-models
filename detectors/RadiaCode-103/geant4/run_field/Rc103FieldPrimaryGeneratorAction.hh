// Розыгрыш ИЗОТРОПНОГО поля фотонов внутри виртуальной сферы радиуса R_SRC.
//
// Виртуальная сфера — НЕ физический объём: никакого солида для неё не
// создаётся, она существует только как правило розыгрыша стартовой точки.
//
// R_SRC = 70 мм. Обоснование: полудиагональ корпуса 123x34x17.5 мм равна
// 64.4 мм, 70 мм даёт запас 5.6 мм. Мир (полуразмер 150 мм) её вмещает.
//
// КЛЮЧЕВАЯ ФИЗИКА: точка разыгрывается равномерно по поверхности сферы, а
// направление — по КОСИНУСНОМУ закону относительно ВНУТРЕННЕЙ нормали (не
// изотропно по направлению!). Только такая пара даёт изотропное поле внутри
// сферы, и только для неё верно тождество нормировки Ф = 4N/S.
#pragma once

#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4Event;
class Rc103FieldFluxSpectrum;

class Rc103FieldPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  static constexpr double kRSrcMm = 300.0;
  static constexpr double kRSrcCm = 30.0;

  explicit Rc103FieldPrimaryGeneratorAction(const Rc103FieldFluxSpectrum* spectrum);
  ~Rc103FieldPrimaryGeneratorAction() override = default;

  void GeneratePrimaries(G4Event* event) override;

 private:
  const Rc103FieldFluxSpectrum* fSpectrum;
  G4ParticleGun fGun;
};
