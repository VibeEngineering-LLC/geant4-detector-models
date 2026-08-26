// Гистограмма депонированной в Crystal_log энергии (1 кэВ/канал, 0..699
// кэВ) + счётчики (события/hits/фотопик), печать в G4cout + запись CSV в
// конце прогона. Простой smoke-test паттерн — БЕЗ разложения по процессам,
// БЕЗ Y-биннинга (в отличие от старой линии detectors/RadiaCode-103/
// geometry/, которую этот код НЕ копирует и не читал).
#pragma once

#include "G4UserRunAction.hh"
#include "globals.hh"

#include <array>
#include <string>

class G4Run;

class Rc103RunRunAction : public G4UserRunAction {
 public:
  // outputCsv — путь итогового CSV; solidAngleFractionCone — Omega_cone/4pi
  // конуса источника (для пересчёта эффективности «в конусе» в абсолютную).
  Rc103RunRunAction(std::string outputCsv, double solidAngleFractionCone);
  ~Rc103RunRunAction() override = default;

  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

  void RecordEvent(G4double edepMeV);

  static constexpr int kNBins = 700;             // 1 кэВ/канал, 0..699 кэВ
  static constexpr double kBinWidthKeV = 1.0;
  static constexpr double kLineKeV = 661.657;     // Cs-137
  static constexpr double kWindowHalfKeV = 10.0;  // +-10 кэВ окно фотопика

 private:
  std::string fOutputCsv;
  double fSolidAngleFractionCone;
  long long fNEvents = 0;
  long long fNHits = 0;
  long long fNPhotopeak = 0;
  std::array<long long, kNBins> fHistogram{};
};
