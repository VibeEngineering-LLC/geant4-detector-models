#pragma once

class NpsmLightYield {
 public:
  NpsmLightYield();  // параметры по умолчанию — значения для NaI(Tl)
  // Вес светового выхода. S — тормозная способность, МэВ/см.
  double Weight(double S_MeV_cm) const;
  // Настройка параметров; каждый обязан быть строго положительным,
  // eta — в интервале (0, 1]. Нарушение — громкий отказ, см. ниже.
  void SetParameters(double eta, double sOns, double sTrap, double sBirks);
  // Включена ли модель. Когда выключена, Weight ВСЕГДА возвращает 1.0.
  void SetEnabled(bool on) { fEnabled = on; }
  bool IsEnabled() const { return fEnabled; }
  // Для шапки CSV — постановка обязана лежать В ФАЙЛЕ, а не в его имени.
  double Eta() const; double SOns() const; double STrap() const; double SBirks() const;
  long long BadSCount() const { return fNBadS; }
  long long OutOfRangeCount() const { return fNOutOfRange; }

 private:
  bool fEnabled = false;   // ВЫКЛЮЧЕНА по умолчанию
  double fEta, fSOns, fSTrap, fSBirks;
  mutable long long fNBadS = 0;
  mutable long long fNOutOfRange = 0;
};

// Самопроверка модели: 5 утверждений, 0 при успехе. Определена в .cc;
// объявление здесь, чтобы main.cc не заводил свой extern-прототип.
int NpsmLightYieldSelfTest();
