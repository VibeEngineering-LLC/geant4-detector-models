// Розыгрыш космических мюонов (mu-) через ГОРИЗОНТАЛЬНЫЙ ДИСК радиуса R_DISK
// на высоте Z_DISK над прибором; мюоны летят вниз.
//
// ГЕОМЕТРИЯ ПРИБОРА В НАШЕЙ МОДЕЛИ: корпус Case_outer 123 x 34 x 17.5 мм с
// центром в (0,0,0) — длинная ось по X, по Z прибор всего +-8.75 мм. Отсюда
// Z_DISK = 100 мм (заведомо выше) и R_DISK по умолчанию 150 мм при
// полудиагонали проекции sqrt(61.5^2 + 17^2) = 63.8 мм.
//
// R_DISK переопределяется ключом CLI `rdisk=<мм>` — это требование спеки,
// без него невозможна обязательная проверка насыщения по радиусу.
//
// ⚠ УГЛОВОЕ РАСПРЕДЕЛЕНИЕ — критическая деталь, здесь уже была ошибка.
// Интенсивность I(θ) ~ cos²θ задана на площадку, ПЕРПЕНДИКУЛЯРНУЮ треку.
// Через ГОРИЗОНТАЛЬНУЮ площадку проходит dN ~ I(θ)·cosθ·dΩ ~ cos³θ·sinθ·dθ,
// то есть плотность по cosθ есть p(cosθ) ~ cosθ³, а обратная функция
// распределения — U^(1/4). Поэтому в коде стоит std::pow(U, 0.25) и НИКОГДА
// не std::cbrt (тот дал бы p ~ cos²θ, мюоны стартовали бы наклоннее реальных
// и континуум искажался бы). Ошибка с cbrt найдена и исправлена 15.08.2026 —
// не повторять.
//
// Энергия и угол разыгрываются НЕЗАВИСИМО. Упрощение заявлено явно: в самой
// формуле Гайссера E и θ связаны через E·cosθ; разъединение оправдано тем, что
// вход и так экстраполяция (см. Rc103MuonSpectrum.hh), и добавлять точную 2D-
// связь означало бы ложную точность поверх приближения.
#pragma once

#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4Event;
class Rc103MuonSpectrum;

class Rc103MuonPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  static constexpr double kZDiskMm = 100.0;
  static constexpr double kRDiskDefaultMm = 150.0;

  Rc103MuonPrimaryGeneratorAction(const Rc103MuonSpectrum* spectrum,
                                  double rDiskMm);
  ~Rc103MuonPrimaryGeneratorAction() override = default;

  void GeneratePrimaries(G4Event* event) override;

 private:
  const Rc103MuonSpectrum* fSpectrum;
  double fRDiskMm;
  G4ParticleGun fGun;
};
