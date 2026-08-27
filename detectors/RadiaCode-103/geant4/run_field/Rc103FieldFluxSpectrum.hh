// Чтение спектра флюенса поля ЕРН из CSV шага 1 двухшаговой схемы
// (detectors/RadiaCode-103/results/wallion/wf_m1_<нуклид>.csv).
//
// Формат входа (проверен фактом по wf_m1_K40.csv, 27.08.2026):
//   # ... комментарии, среди них строка
//   # fluence_total_cm2_s = 9.510101e-03
//   # bin_keV = 2.000  (последний канал = переполнение)
//   E_keV,fluence_cm2_s
//   1.0,1.730823e-07
//   ...
// Строки данных РАЗРЕЖЕНЫ (только ненулевые бины). E_keV берётся как есть —
// это центр бина шириной 2 кэВ.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

class Rc103FieldFluxSpectrum {
 public:
  // Возвращает false, если файл не открылся, не нашлась строка
  // fluence_total_cm2_s или не оказалось ни одной строки данных.
  bool Load(const std::string& path);

  // Ф из ЗАГОЛОВКА файла — это и есть номинальный флюенс нормировки.
  double HeaderTotalCm2S() const { return fHeaderTotal; }
  // Прямая сумма колонки fluence_cm2_s — для самопроверки чтения входа.
  double SumColumnCm2S() const { return fSumColumn; }
  std::size_t NBins() const { return fEnergyKeV.size(); }

  // Розыгрыш энергии по кумулятивному распределению колонки fluence_cm2_s.
  // Требует уже выполненного Load().
  double SampleEnergyKeV() const;

 private:
  std::vector<double> fEnergyKeV;
  std::vector<double> fCumulative;  // нормированная на 1 кумулята
  double fHeaderTotal = 0.0;
  double fSumColumn = 0.0;
};
