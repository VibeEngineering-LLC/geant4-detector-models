// Гистограмма депонированной в Crystal_log энергии (2 кэВ/канал, 0..3000 кэВ,
// та же сетка что wallfield.cc::Fluence — для прямой сравнимости бинов) +
// счётчики событий/hits, печать в G4cout + запись CSV. Смоук-тест паттерн:
// БЕЗ биасинга, статистика в кристалле честно может оказаться мала (задание,
// см. main.cc и предупреждение ниже) — так и печатается, не приукрашивается.
#pragma once

#include "G4UserRunAction.hh"
#include "globals.hh"

#include <array>
#include <string>

class G4Run;

class Rc103RoomRunAction : public G4UserRunAction {
 public:
  // outputCsv — путь итогового CSV; wallEmissionRatePerSec — полная объёмная
  // скорость испускания фотонов во ВСЕЙ толще бетона (фотон/с, из
  // Rc103RoomPrimaryGeneratorAction::GetWallEmissionRatePerSec()) — нужна
  // для пересчёта hits/N_events в АБСОЛЮТНУЮ прогнозную скорость счёта в
  // кристалле, не только относительную эффективность.
  Rc103RoomRunAction(std::string outputCsv, double wallEmissionRatePerSec);
  ~Rc103RoomRunAction() override = default;

  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

  void RecordEvent(G4double edepMeV);

  static constexpr int kNBins = 1500;         // 2 кэВ/канал, 0..2999 кэВ
  static constexpr double kBinWidthKeV = 2.0;

 private:
  std::string fOutputCsv;
  double fWallEmissionRatePerSec;
  long long fNEvents = 0;
  long long fNHits = 0;
  std::array<long long, kNBins> fHistogram{};
};
