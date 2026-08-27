// Накопление спектра энерговыделения в Crystal_log от космических мюонов и
// запись CSV.
//
// НОРМИРОВКА — НА ОДИН МЮОН, ПЕРЕСЁКШИЙ ДИСК. Абсолютный поток НЕ добивается:
// в выход идёт per_muon = counts / n_events, то есть вероятность отклика.
// Амплитуда (сколько таких мюонов в секунду) подбирается позже NNLS вместе с
// активностями нуклидов.
//
// Сетка гистограммы 1 кэВ/бин, 0..2999 кэВ — та же, что у нуклидных шаблонов
// run_field/, иначе их нельзя сложить в одну матрицу.
//
// ⚠ ПЕРЕПОЛНЕНИЕ. У мюона депозиты доходят до десятков МэВ (сквозной трек
// через ~1 см CsI даёт ~5-8 МэВ). Всё, что >= 3000 кэВ, идёт в отдельный
// счётчик fNOverflow и НЕ сваливается в последний бин и НЕ теряется молча.
#pragma once

#include "G4UserRunAction.hh"
#include "globals.hh"

#include <array>
#include <string>

class G4Run;

class Rc103MuonRunAction : public G4UserRunAction {
 public:
  static constexpr int kNBins = 3000;          // 0..2999 кэВ
  static constexpr double kBinWidthKeV = 1.0;  // сетка мелкая и НЕсвёрнутая

  Rc103MuonRunAction(std::string outputCsv, double rDiskMm, double zDiskMm,
                     double diskAreaCm2, double pdgExpectedPerS, double eLoGeV,
                     double eHiGeV);

  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run* run) override;

  void RecordEvent(G4double edepMeV);

 private:
  std::string fOutputCsv;
  double fRDiskMm;
  double fZDiskMm;
  double fDiskAreaCm2;
  double fPdgExpectedPerS;
  double fELoGeV;
  double fEHiGeV;

  long long fNEvents = 0;
  long long fNHits = 0;
  long long fNOverflow = 0;
  double fMaxEdepKeV = 0.0;
  std::array<long long, kNBins> fHistogram{};
};
