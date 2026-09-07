#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstring>

#include "G4RunManagerFactory.hh"
#include "G4Threading.hh"
#include "NpsmBenchActionInitialization.hh"
#include "G4EmParameters.hh"

// Число рабочих потоков. **Умолчание -1 = все ядра машины** (решение
// оператора 07.09.2026: «проще включить и пересчитать что нужно заново»).
// 0 — однопоточный режим, которым считались прогоны до 07.09.2026;
// N>0 — ровно N потоков.
//
// Почему это безопасно: слияние потоков проверено (1/2/4/8 дают идентичный
// результат до последнего целого), согласие с однопоточным режимом
// подтверждено на десяти зёрнах (1,0σ по энергии). Serial и MT дают РАЗНЫЕ
// выборки одних событий, поэтому режим пишется в шапку CSV, а merge_runs
// отказывается складывать отрезки разных режимов.
static int gThreads = -1;
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "Randomize.hh"  // G4Random::setTheSeed объявлен здесь, а не в G4Random.hh (его нет в 11.2.1)
#include "G4SystemOfUnits.hh"

#include "NpsmBenchDetectorConstruction.hh"
#include "NpsmBenchPrimaryGeneratorAction.hh"
#include "NpsmBenchEventAction.hh"
#include "NpsmBenchRunAction.hh"
#include "NpsmBenchSteppingAction.hh"
#include "Rc103FieldPhysicsList.hh"
#include "NpsmLightYield.hh"
#include "G4Electron.hh"
#include "G4EmCalculator.hh"
#include <cfloat>

int main(int argc, char** argv) {
    if (argc >= 2 && std::string(argv[1]) == "selftest") {
        return NpsmLightYieldSelfTest();
    }

    if (argc < 4) {
        std::fprintf(stderr, "Usage: %s <energy_keV> <n_events> <out_csv> [seed=<N>] [emcut=<mm>] [deex=std|deex|max] [crystal=<X>x<Y>x<Z>] [particle=gamma|e-] [beam=iso|pencil] [npsm=on|off] [npsm_eta=<value>] [npsm_ons=<value>] [npsm_trap=<value>] [npsm_birks=<value>]\n", argv[0]);
        return 2;
    }

    double energyKeV = std::stod(argv[1]);
    if (energyKeV <= 0) {
        std::fprintf(stderr, "Ошибка: энергия должна быть положительной\n");
        return 2;
    }

    long long nEvents = std::stoll(argv[2]);
    if (nEvents <= 0) {
        std::fprintf(stderr, "Ошибка: количество событий должно быть положительным\n");
        return 2;
    }

    std::string outCsv = argv[3];

    long seed = 0;
    double emCutMm = 0.05;
    std::string deexMode = "max";
    std::string crystalStr = "102x102x406";
    std::string particle = "gamma";
    std::string beam = "iso";
    std::string npsm = "off";
    double npsmEta = 0.0;
    double npsmOns = 0.0;
    double npsmTrap = 0.0;
    double npsmBirks = 0.0;

    for (int i = 4; i < argc; ++i) {
        const char* arg = argv[i];
        if (std::strncmp(arg, "seed=", 5) == 0) {
            seed = std::stoll(arg + 5);
        } else if (std::strncmp(arg, "emcut=", 6) == 0) {
            emCutMm = std::stod(arg + 6);
            if (emCutMm <= 0) {
                std::fprintf(stderr, "Ошибка: emcut должен быть положительным\n");
                return 2;
            }
        } else if (std::strncmp(arg, "deex=", 5) == 0) {
            deexMode = arg + 5;
            if (deexMode != "std" && deexMode != "fluo" &&
                deexMode != "pixefluo" && deexMode != "deex" &&
                deexMode != "max") {
                std::fprintf(stderr, "Ошибка: неизвестный режим deex\n");
                return 2;
            }
        } else if (std::strncmp(arg, "light=", 6) == 0) {
            // npsm  — наша модель: вес на каждом шаге по dE/dx (умолчание);
            // scint — штатный G4Scintillation по интегральной кривой L(E).
            // Второй путь нужен для независимой сверки (линия Г): он считает
            // ту же величину другим способом, поэтому расхождение указывает
            // на реализацию, а не на параметры.
            const std::string v(arg + 6);
            if (v == "scint") {
                NpsmBenchDetectorConstruction::gScintLight = true;
                Rc103FieldPhysicsList::gScintLight = true;
            } else if (v != "npsm") {
                std::fprintf(stderr, "Ошибка: light= принимает npsm или scint\n");
                return 2;
            }
        } else if (std::strncmp(arg, "threads=", 8) == 0) {
            // 0 (умолчание) — Serial, как считались все прежние прогоны.
            // >0 — многопоточный режим с этим числом рабочих потоков.
            gThreads = std::atoi(arg + 8);
            if (gThreads < 0) {
                std::fprintf(stderr, "Ошибка: threads не может быть отрицательным\n");
                return 2;
            }
            // В шапку значение кладётся ниже, при создании менеджера: там
            // известно ФАКТИЧЕСКОЕ число потоков (умолчание -1 разворачивается
            // в число ядер), а не то, что было написано в командной строке.
        } else if (std::strncmp(arg, "deex_region=", 12) == 0) {
            // Имя региона для адресной деэкситации (линия В). Значение
            // попадает в шапку CSV; проверка существования региона — в
            // Geant4, он ругается на неизвестное имя.
            Rc103FieldPhysicsList::gDeexRegion = arg + 12;
        } else if (std::strncmp(arg, "em=", 3) == 0) {
            // opt4 (умолчание) либо opt0 — суррогат «без доплеровского
            // уширения» для линии Д. Значение попадает в шапку CSV.
            Rc103FieldPhysicsList::gEmOption = arg + 3;
            if (Rc103FieldPhysicsList::gEmOption != "opt4" &&
                Rc103FieldPhysicsList::gEmOption != "opt0") {
                std::fprintf(stderr, "Ошибка: em= принимает opt4 или opt0\n");
                return 2;
            }
        } else if (std::strncmp(arg, "crystal=", 8) == 0) {
            crystalStr = arg + 8;
        } else if (std::strncmp(arg, "npsm_local_keV=", 15) == 0) {
            // Граница ветвей continuous/local; умолчание 1 кэВ (Algorithm S1).
            NpsmBenchSteppingAction::gLocalEnergyKeV = std::stod(arg + 15);
            if (NpsmBenchSteppingAction::gLocalEnergyKeV < 0) {
                std::fprintf(stderr, "Ошибка: npsm_local_keV не может быть отрицательным\n");
                return 2;
            }
        } else if (std::strncmp(arg, "restricted_dedx=", 16) == 0) {
            NpsmBenchSteppingAction::gRestrictedDedx = (std::string(arg + 16) == "1");
        } else if (std::strncmp(arg, "dedx_table=", 11) == 0) {
            NpsmBenchSteppingAction::gPrintDedxTable = (std::string(arg + 11) == "1");
        } else if (std::strncmp(arg, "pixe_model=", 11) == 0) {
            // Ключи разбора W-055. Пустое значение = умолчание Geant4,
            // применяются только при deex=max.
            Rc103FieldPhysicsList::gPixeModel = arg + 11;
        } else if (std::strncmp(arg, "pixe_e_model=", 13) == 0) {
            Rc103FieldPhysicsList::gPixeElecModel = arg + 13;
        } else if (std::strncmp(arg, "particle=", 9) == 0) {
            particle = arg + 9;
            if (particle != "gamma" && particle != "e-") {
                std::fprintf(stderr, "Ошибка: неизвестный тип частицы\n");
                return 2;
            }
        } else if (std::strncmp(arg, "beam=", 5) == 0) {
            beam = arg + 5;
            if (beam != "iso" && beam != "pencil") {
                std::fprintf(stderr, "Ошибка: неизвестный тип пучка\n");
                return 2;
            }
        } else if (std::strncmp(arg, "npsm=", 5) == 0) {
            npsm = arg + 5;
            if (npsm != "on" && npsm != "off") {
                std::fprintf(stderr, "Ошибка: неизвестный режим NPSM\n");
                return 2;
            }
        } else if (std::strncmp(arg, "npsm_eta=", 9) == 0) {
            npsmEta = std::stod(arg + 9);
            if (npsmEta <= 0 || npsmEta > 1.0) {
                std::fprintf(stderr, "Ошибка: npsm_eta должен быть в интервале (0, 1]\n");
                return 2;
            }
        } else if (std::strncmp(arg, "npsm_ons=", 9) == 0) {
            npsmOns = std::stod(arg + 9);
            if (npsmOns <= 0) {
                std::fprintf(stderr, "Ошибка: npsm_ons должен быть положительным\n");
                return 2;
            }
        } else if (std::strncmp(arg, "npsm_trap=", 10) == 0) {
            npsmTrap = std::stod(arg + 10);
            if (npsmTrap <= 0) {
                std::fprintf(stderr, "Ошибка: npsm_trap должен быть положительным\n");
                return 2;
            }
        } else if (std::strncmp(arg, "npsm_birks=", 11) == 0) {
            npsmBirks = std::stod(arg + 11);
            if (npsmBirks <= 0) {
                std::fprintf(stderr, "Ошибка: npsm_birks должен быть положительным\n");
                return 2;
            }
        } else {
            std::fprintf(stderr, "Ошибка: неизвестный параметр %s\n", arg);
            return 2;
        }
    }

    // Парсинг размеров кристалла
    long long x, y, z;
    if (std::sscanf(crystalStr.c_str(), "%lldx%lldx%lld", &x, &y, &z) != 3) {
        std::fprintf(stderr, "Ошибка: неверный формат crystal=<X>x<Y>x<Z>\n");
        return 2;
    }

    NpsmBenchDetectorConstruction::gCrystalXMm = static_cast<double>(x);
    NpsmBenchDetectorConstruction::gCrystalYMm = static_cast<double>(y);
    NpsmBenchDetectorConstruction::gCrystalZMm = static_cast<double>(z);

    // Создание объекта NpsmLightYield и установка параметров
    NpsmLightYield lightYield;
    double eta = lightYield.Eta();
    double ons = lightYield.SOns();
    double trap = lightYield.STrap();
    double birks = lightYield.SBirks();

    if (npsmEta != 0.0) eta = npsmEta;
    if (npsmOns != 0.0) ons = npsmOns;
    if (npsmTrap != 0.0) trap = npsmTrap;
    if (npsmBirks != 0.0) birks = npsmBirks;

    lightYield.SetParameters(eta, ons, trap, birks);
    lightYield.SetEnabled(npsm == "on");

    // Инициализация Геанта. Умолчание — Serial: все прежние прогоны контура
    // считаны так, и их воспроизводимость менять нельзя. Многопоточность
    // включается только явным ключом threads=N (линия И, 07.09.2026).
    //
    // ⚠ Тип запрошен именно MT, а не Default: по умолчанию 11.x берёт режим
    // Tasking, а в нём (тема форума, разработчик) в 11.4.x есть регрессия —
    // параметр luxury генератора не инициализируется в рабочих потоках.
    G4RunManager* runManager = nullptr;
    if (gThreads != 0) {
        // -1 — все ядра машины; иначе ровно столько потоков, сколько задано.
        const int nThreads = (gThreads < 0)
                                 ? G4Threading::G4GetNumberOfCores()
                                 : gThreads;
        runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::MT);
        runManager->SetNumberOfThreads(nThreads);
        NpsmBenchRunAction::gThreads = nThreads;  // в шапку идёт ФАКТ, не «-1»
    } else {
        runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
        NpsmBenchRunAction::gThreads = 0;
    }
    runManager->SetVerboseLevel(0);

    NpsmBenchPrimaryGeneratorAction::gParticle = particle;
    NpsmBenchPrimaryGeneratorAction::gPencilBeam = (beam == "pencil");

    auto detectorConstruction = new NpsmBenchDetectorConstruction();
    runManager->SetUserInitialization(detectorConstruction);

    auto physicsList = new Rc103FieldPhysicsList(emCutMm, deexMode);
    runManager->SetUserInitialization(physicsList);

    // Действия создаёт инициализатор: в MT они нужны по экземпляру на поток,
    // а мастеру — только RunAction (он получает слитый результат и пишет
    // файл). Прямые SetUserAction для этого не годятся. В Serial поведение
    // то же самое: Build() зовётся один раз.
    runManager->SetUserInitialization(new NpsmBenchActionInitialization(
        outCsv, energyKeV, nEvents, seed, &lightYield, 0.0));

    runManager->Initialize();

    if (!NpsmBenchDetectorConstruction::GetCrystalLogicalVolume()) {
        std::fprintf(stderr, "Ошибка: кристалл не был создан\n");
        return 3;
    }

    // Флаги деэкситации ПОСЛЕ инициализации — прочитанные из G4EmParameters.
    // Строка физ-листа печатается в его конструкторе, то есть ДО того, как
    // отработали остальные конструкторы физики, и подтверждает лишь наше
    // намерение. В полной геометрии это уже поймано: G4RadioactiveDecayPhysics
    // безусловно включает Оже и игнорирование порогов (W-068 доп. 1). Здесь
    // распад не регистрируется, но строка нужна по тому же правилу:
    // постановка проверяется по состоянию Geant4 в момент счёта.
    {
        auto* em = G4EmParameters::Instance();
        std::printf("POST_INIT_EM: fluo=%d auger=%d pixe=%d ignore_cut=%d "
                    "apply_cuts=%d lowest_e_keV=%.3f\n",
                    static_cast<int>(em->Fluo()), static_cast<int>(em->Auger()),
                    static_cast<int>(em->Pixe()),
                    static_cast<int>(em->DeexcitationIgnoreCut()),
                    static_cast<int>(em->ApplyCuts()),
                    em->LowestElectronEnergy() / keV);
        std::fflush(stdout);
    }

    // Установка UI команд
    auto ui = G4UImanager::GetUIpointer();
    ui->ApplyCommand("/run/verbose 1");
    ui->ApplyCommand("/event/verbose 0");
    ui->ApplyCommand("/tracking/verbose 0");

    long long printProgress = nEvents / 10;
    if (printProgress < 1) printProgress = 1;
    std::ostringstream oss;
    oss << "/run/printProgress " << printProgress;
    ui->ApplyCommand(oss.str());

    // Вывод информации перед запуском
    double srcRadius = NpsmBenchDetectorConstruction::SourceRadiusMm();
    std::printf("npsm_bench: energy_keV=%.3f n_events=%lld out_csv=%s crystal_mm=%gx%gx%g r_src_mm=%g seed=%ld\n",
                energyKeV, nEvents, outCsv.c_str(),
                NpsmBenchDetectorConstruction::gCrystalXMm,
                NpsmBenchDetectorConstruction::gCrystalYMm,
                NpsmBenchDetectorConstruction::gCrystalZMm,
                srcRadius, seed);

    std::printf("npsm_bench: particle=%s beam=%s npsm=%s eta=%g s_ons=%g s_trap=%g s_birks=%g\n",
                particle.c_str(), beam.c_str(), npsm.c_str(),
                eta, ons, trap, birks);

    // Установка зерна случайных чисел
    if (seed != 0) {
        G4Random::setTheSeed(seed);
    }

    runManager->BeamOn(nEvents);

    delete runManager;

    std::printf("EXITCODE=0\n");
    return 0;
}
