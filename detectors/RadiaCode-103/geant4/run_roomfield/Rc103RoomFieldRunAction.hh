// Накопление track-length спектра флюенса в шаре и запись CSV в формате,
// СОВМЕСТИМОМ с results/wallion/wf_m1_*.csv (шапка '#', строка
// 'E_keV,fluence_cm2_s', пары; поле '# fluence_total_cm2_s = ...' читает шаг 2).
//
// Нормировка (спека, раздел «Нормировка»):
//   S_v_i = A_уд * rho_i,  A_уд = 1 Бк/кг = 1e-3 Бк/г   [расп/(см3 с)]
//   R     = sum_i (S_v_i * V_i) = 1e-3 * M_полная[г]     [расп/с]
//   T     = N_событий / R                                [с]
//   Ф(бин) = len_бина[см] / (V_шара[см3] * T)            [1/(см2 с)]
//
// Самопроверка (обязательна): сумма НАПЕЧАТАННОЙ колонки сверяется с
// НАПЕЧАТАННЫМ значением из шапки — сравниваются те самые числа, что попали в
// файл, а не внутренние double. Расхождение > 1e-6 относительной — прогон
// объявляется проваленным (код возврата 4), а не «почти сошлось».
#pragma once

#include "G4UserRunAction.hh"
#include "globals.hh"

#include <string>
#include <vector>

class G4Run;

class Rc103RoomFieldRunAction : public G4UserRunAction {
 public:
  // 2 кэВ на канал, 0..3000 кэВ; индекс kBins — канал переполнения (catch-all),
  // ровно как в geometry/wallfield.cc, чей формат читает шаг 2.
  static constexpr int kBins = 1500;
  static constexpr double kBinKeV = 2.0;

  Rc103RoomFieldRunAction(std::string outputCsv, std::string nuclide, int ionZ,
                          int ionA, double ratePerS, double tRunS,
                          double ballVolumeCm3, double specificActivityBqPerKg);

  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run* run) override;

  void AddTrackLengthCm(double eKeV, double lenCm);

 private:
  std::string fOutputCsv;
  std::string fNuclide;
  int fIonZ;
  int fIonA;
  double fRatePerS;
  double fTRunS;
  double fBallVolumeCm3;
  double fSpecificActivityBqPerKg;
  std::vector<double> fLenCm;
};
