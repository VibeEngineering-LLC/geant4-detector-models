// Накопление спектра энерговыделения в Crystal_log, нормировка в cps на
// 1 Бк/кг родителя и запись CSV.
//
// Нормировка (тождество Ф = 4N/S для выпуклой поверхности в изотропном поле):
//   S     = 4*pi*R_SRC^2                 [см^2], R_SRC = 7.0 см
//   Ф     = fluence_total_cm2_s          [1/(см^2*с)] - ИЗ ЗАГОЛОВКА входа
//   rate  = Ф * S / 4.0                  [1/с]
//   T_run = N_events / rate              [с]
//   cps(канал) = counts(канал) / T_run   [1/с на 1 Бк/кг родителя]
//
// Режим --check-norm дополнительно копит суммарную длину треков ФОТОНОВ в
// контрольном шаре и считает Ф_изм = sumLen / (V * T_run).
#pragma once

#include "G4UserRunAction.hh"
#include "globals.hh"

#include <array>
#include <string>

class G4Run;

class Rc103FieldRunAction : public G4UserRunAction {
 public:
  static constexpr int kNBins = 3000;        // 0..2999 кэВ
  static constexpr double kBinWidthKeV = 1.0;  // сетка мелкая и НЕсвёрнутая

  Rc103FieldRunAction(std::string outputCsv, double fluxTotalCm2S,
                      double surfaceCm2, double ratePerS, double tRunS,
                      bool checkNormMode, double checkVolumeCm3);

  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run* run) override;

  // category: 0 = мимо свинца, 1 = рассеян в свинце, 2 = рождён в свинце
  void RecordEvent(G4double edepMeV, int category);
  void AddGammaTrackLengthMm(G4double lenMm) { fSumTrackLenMm += lenMm; }

 private:
  std::string fOutputCsv;
  double fFluxTotalCm2S;
  double fSurfaceCm2;
  double fRatePerS;
  double fTRunS;
  bool fCheckNormMode;
  double fCheckVolumeCm3;

  long long fNEvents = 0;
  long long fNHits = 0;
  double fSumTrackLenMm = 0.0;
  std::array<long long, kNBins> fHistogram{};
  // Разбиение того же счёта по происхождению (29.08.2026). Сумма трёх
  // гистограмм обязана совпадать с fHistogram побитово.
  static constexpr int kNCat = 3;
  std::array<std::array<long long, kNBins>, kNCat> fHistOrigin{};
  std::array<long long, kNCat> fNHitsOrigin{};
};
