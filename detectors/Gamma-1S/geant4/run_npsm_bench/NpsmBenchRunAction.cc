#include "NpsmBenchRunAction.hh"
// Шапка CSV читает статические поля геометрии и физ-листа напрямую — без этих
// заголовков компилятор видел «<error type>» на SourceRadiusMm() (сборка 05.09).
#include "NpsmBenchDetectorConstruction.hh"
#include "Rc103FieldPhysicsList.hh"
#include "G4EmParameters.hh"
#include "G4SystemOfUnits.hh"
#include "NpsmLightYield.hh"
#include "NpsmBenchPrimaryGeneratorAction.hh"
#include "NpsmBenchSteppingAction.hh"
#include "G4RunManager.hh"
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <cstdlib>

std::string NpsmBenchRunAction::gPrimaryKind = "gamma";
int NpsmBenchRunAction::gIonZ = 0;
int NpsmBenchRunAction::gIonA = 0;
int NpsmBenchRunAction::gThreads = 0;

NpsmBenchRunAction::NpsmBenchRunAction(std::string outCsv, double energyKeV, long long nEventsRequested, long seed, const NpsmLightYield* lightYield)
    : fOutCsv(outCsv), fEnergyKeV(energyKeV), fNEventsRequested(nEventsRequested), fSeed(seed),
      fLightYield(lightYield) {
    if (!fLightYield) {
        std::fprintf(stderr, "NpsmBenchRunAction: FATAL lightYield == nullptr\n");
        std::abort();
    }
}

void NpsmBenchRunAction::BeginOfRunAction(const G4Run*) {
    // Ничего не делаем
}

// Объект прогона создаёт Geant4 — по одному на поток в MT и один в Serial.
// Возвращать здесь свой тип обязательно: иначе накопителей просто не будет,
// а слияние потоков окажется нечем выполнять.
G4Run* NpsmBenchRunAction::GenerateRun() {
    return new NpsmBenchRun();
}

void NpsmBenchRunAction::RecordEvent(int nCompt, int nRayl, bool phot, bool conv, bool escaped, double edepMeV, double edepLightMeV) {
    // Пишем в объект ТЕКУЩЕГО прогона (свой у каждого потока), а не в поля
    // этого действия. G4RunManager::GetNonConstCurrentRun даёт его и на
    // рабочем потоке. Нулевой указатель здесь означал бы, что события идут
    // мимо накопителя — молчать об этом нельзя.
    auto* raw = G4RunManager::GetRunManager()->GetNonConstCurrentRun();
    auto* run = dynamic_cast<NpsmBenchRun*>(raw);
    if (!run) {
        std::fprintf(stderr, "FATAL: текущий G4Run не NpsmBenchRun — "
                             "событие некуда записать, статистика ПОТЕРЯНА\n");
        std::fflush(stderr);
        std::abort();
    }
    run->RecordEvent(nCompt, nRayl, phot, conv, escaped, edepMeV, edepLightMeV,
                     fEnergyKeV);
}

void NpsmBenchRunAction::EndOfRunAction(const G4Run* aRun) {
    // На мастере aRun — уже СЛИТЫЙ объект (Geant4 вызывает Merge для каждого
    // рабочего потока до этого момента). На рабочем потоке файл писать нельзя:
    // потоки затирали бы результат друг друга, и в файле осталась бы доля
    // статистики при внешне успешном прогоне.
    if (!IsMaster()) {
        return;
    }
    const auto* runPtr = dynamic_cast<const NpsmBenchRun*>(aRun);
    if (!runPtr) {
        std::fprintf(stderr, "FATAL: EndOfRunAction получил не NpsmBenchRun — "
                             "писать нечего, результат ПОТЕРЯН\n");
        std::fflush(stderr);
        std::abort();
    }
    const NpsmBenchRun& run = *runPtr;

    WriteCSV(run);

    double mean = run.fNAbsorbed > 0 ? run.fSumComptAbsorbed / run.fNAbsorbed : 0.0;
    double sem = 0.0;
    if (run.fNAbsorbed > 0) {
        double mean2 = run.fSumCompt2Absorbed / run.fNAbsorbed;
        double variance = mean2 - mean * mean;
        if (variance >= 0 && run.fNAbsorbed > 1) {
            sem = std::sqrt(variance / (run.fNAbsorbed - 1));
        }
    }

    // Вывод в stdout
    std::printf("npsm_bench: E=%.1f keV  N_absorbed=%lld  mean_ncompt=%.4f +- %.4f  escaped=%lld\n",
                fEnergyKeV, run.fNAbsorbed, mean, sem, run.fNEscaped);
    
    double ratio = (run.fSumEdepMeV > 0.0) ? run.fSumLightMeV / run.fSumEdepMeV : 0.0;
    const double meanEdepKeV  = (run.fNWithEdep > 0) ? run.fSumEdepMeV  * 1000.0 / run.fNWithEdep : 0.0;
    const double meanLightKeV = (run.fNWithEdep > 0) ? run.fSumLightMeV * 1000.0 / run.fNWithEdep : 0.0;
    std::printf("npsm_bench: n_with_edep=%lld mean_edep_keV=%.3f mean_light_keV=%.3f ratio=%.6f\n",
                run.fNWithEdep, meanEdepKeV, meanLightKeV, ratio);
}

void NpsmBenchRunAction::WriteCSV(const NpsmBenchRun& run) {
    std::ofstream file(fOutCsv, std::ios::out);
    if (!file.is_open()) {
        // ТИХИЙ ОТКАЗ ЗДЕСЬ УБИТ 07.09.2026 (W-067). Прежде функция молча
        // возвращалась: прогон отрабатывал до конца, печатал статистику в
        // stdout, возвращал ноль — и не создавал файла. Так были потеряны
        // прогоны 07.09, запущенные с относительным путём `out=out/...`, когда
        // рабочим каталогом процесса оказался каталог без подпапки `out`.
        // Считать такой прогон успешным нельзя: результата нет.
        std::fprintf(stderr,
                     "FATAL: не удалось открыть для записи '%s'. "
                     "Проверьте, что каталог существует; относительный путь "
                     "считается от рабочего каталога процесса, а не от места "
                     "запуска .cmd. Результат прогона ПОТЕРЯН.\n",
                     fOutCsv.c_str());
        std::fflush(stderr);
        std::abort();
    }

    file << "# npsm_bench geometry/transport benchmark\n";
    file << "energy_keV," << fEnergyKeV << "\n";
    file << "n_events_requested," << fNEventsRequested << "\n";
    file << "n_events_processed," << run.fNEvents << "\n";
    file << "seed," << fSeed << "\n";
    // Режим счёта. Serial и MT — разные выборки одних и тех же событий,
    // поэтому сравнивать между собой можно только прогоны одного режима.
    // Без этой строки два таких файла внешне неразличимы (класс W-068).
    file << "threads," << gThreads << "\n";
    // Каким путём считался свет и сколько фотонов родила штатная
    // сцинтилляция (0, если режим не включён) — линия Г.
    file << "light," << (NpsmBenchDetectorConstruction::gScintLight ? "scint" : "npsm") << "\n";
    file << "n_scint_photons," << run.fNScintPhotons << "\n";
    file << "em_cut_mm," << Rc103FieldPhysicsList::gCutMm << "\n";
    file << "em_deex," << Rc103FieldPhysicsList::gDeexMode << "\n";
    file << "em_option," << Rc103FieldPhysicsList::gEmOption << "\n";
    file << "deex_region," << (Rc103FieldPhysicsList::gDeexRegion.empty()
                                   ? "-" : Rc103FieldPhysicsList::gDeexRegion)
         << "\n";
    // Ключи, принимавшиеся командной строкой, но НЕ попадавшие в результат
    // (W-068, скан класса): модели сечений PIXE меняют физику деэкситации,
    // а ветвь local/continuous и restricted dE/dx — саму формулу света.
    // Значения PIXE берутся ИЗ G4EmParameters (что применилось), а не из
    // наших переменных (W-052); при deex != max Geant4 их игнорирует, и в
    // файле будет видно именно это.
    file << "pixe_model," << G4EmParameters::Instance()->PIXECrossSectionModel() << "\n";
    file << "pixe_e_model," << G4EmParameters::Instance()->PIXEElectronCrossSectionModel() << "\n";
    file << "npsm_local_keV," << NpsmBenchSteppingAction::gLocalEnergyKeV << "\n";
    file << "restricted_dedx," << (NpsmBenchSteppingAction::gRestrictedDedx ? 1 : 0) << "\n";
    file << "lowest_electron_energy_keV," << G4EmParameters::Instance()->LowestElectronEnergy() / keV << "\n";
    file << "crystal_mm," << NpsmBenchDetectorConstruction::gCrystalXMm << "x" << NpsmBenchDetectorConstruction::gCrystalYMm << "x" << NpsmBenchDetectorConstruction::gCrystalZMm << "\n";
    file << "crystal_volume_cm3," << (NpsmBenchDetectorConstruction::gCrystalXMm * NpsmBenchDetectorConstruction::gCrystalYMm * NpsmBenchDetectorConstruction::gCrystalZMm) / 1000.0 << "\n";
    file << "r_src_mm," << NpsmBenchDetectorConstruction::SourceRadiusMm() << "\n";
    // Постановка первички. `particle` — поле генератора СТЕНДА и в проекте с
    // чужим генератором (run_g1s_npsm) остаётся на умолчании; поэтому рядом
    // печатается gPrimaryKind, который выставляет main, и — из самого физ-листа,
    // не из наших переменных (W-052) — распад и угловые корреляции. Без этих
    // строк два ионных прогона с corr_gamma=0/1 были неразличимы по шапке (07.09).
    file << "particle," << NpsmBenchPrimaryGeneratorAction::gParticle << "\n";
    file << "primary," << gPrimaryKind << "\n";
    if (gPrimaryKind == "ion") {
        file << "ion_z," << gIonZ << "\n";
        file << "ion_a," << gIonA << "\n";
    }
    file << "decay," << (Rc103FieldPhysicsList::gDecay ? 1 : 0) << "\n";
    file << "correlated_gamma," << (Rc103FieldPhysicsList::gCorrGamma ? 1 : 0) << "\n";
    file << "beam," << (NpsmBenchPrimaryGeneratorAction::gPencilBeam ? "pencil" : "iso") << "\n";
    file << "npsm_enabled," << (fLightYield->IsEnabled() ? 1 : 0) << "\n";
    file << "npsm_eta," << fLightYield->Eta() << "\n";
    file << "npsm_s_ons," << fLightYield->SOns() << "\n";
    file << "npsm_s_trap," << fLightYield->STrap() << "\n";
    file << "npsm_s_birks," << fLightYield->SBirks() << "\n";
    file << "npsm_bad_s_count," << fLightYield->BadSCount() << "\n";
    file << "npsm_out_of_range_count," << fLightYield->OutOfRangeCount() << "\n";
    file << "n_with_edep," << run.fNWithEdep << "\n";
    const std::streamsize prevPrec = file.precision();
    file << std::setprecision(17);
    file << "sum_edep_MeV," << run.fSumEdepMeV << "\n";
    file << "sum_light_MeV," << run.fSumLightMeV << "\n";
    file << "sum_light2_MeV2," << run.fSumLight2MeV2 << "\n";
    file << std::setprecision(prevPrec);
    file << "n_overflow_edep," << run.fNOverflowEdep << "\n";
    file << "n_overflow_light," << run.fNOverflowLight << "\n";
    file << "n_absorbed_phot," << run.fNAbsorbed << "\n";
    file << "n_full_edep_1keV," << run.fNFullEdep << "\n";
    file << std::setprecision(17);
    file << "sum_light_full_MeV," << run.fSumLightFullMeV << "\n";
    file << "sum_light2_full_MeV2," << run.fSumLight2FullMeV2 << "\n";
    file << std::setprecision(6);
    file << "n_conv," << run.fNConv << "\n";
    file << "n_escaped," << run.fNEscaped << "\n";
    file << "n_other," << run.fNOther << "\n";
    file << "mean_ncompt_absorbed," << (run.fNAbsorbed > 0 ? run.fSumComptAbsorbed / run.fNAbsorbed : 0.0) << "\n";
    file << "sem_ncompt_absorbed," << (run.fNAbsorbed > 0 ? std::sqrt((run.fSumCompt2Absorbed / run.fNAbsorbed - (run.fSumComptAbsorbed / run.fNAbsorbed) * (run.fSumComptAbsorbed / run.fNAbsorbed)) / run.fNAbsorbed) : 0.0) << "\n";
    file << "n_compt,count_absorbed,count_all\n";

    for (int i = 0; i <= NpsmBenchRun::kMaxCompt; ++i) {
        file << i << "," << run.fHistAbsorbed[i] << "," << run.fHistAll[i] << "\n";
    }

    // Вторая таблица — спектры
    file << "\nbin_keV,count_edep,count_light\n";
    for (int i = 0; i < NpsmBenchRun::kNBins; ++i) {
        if (run.fSpecEdep[i] != 0 || run.fSpecLight[i] != 0) {
            const double binCenterKeV = (i + 0.5) * NpsmBenchRun::kBinKeV;
            file << binCenterKeV << "," << run.fSpecEdep[i] << "," << run.fSpecLight[i] << "\n";
        }
    }

    file.close();
}
