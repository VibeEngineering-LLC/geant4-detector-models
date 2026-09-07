#include "G4RunManagerFactory.hh"
#include "G4Threading.hh"
#include "G1sNpsmActionInitialization.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"

#include "G1SDetector.hh"
#include "Rc103FieldPhysicsList.hh"
#include "NpsmLightYield.hh"
#include "NpsmBenchRunAction.hh"
#include "NpsmBenchEventAction.hh"
#include "NpsmBenchSteppingAction.hh"
#include "NpsmBenchDetectorConstruction.hh"
#include "G1sNpsmPrimaryGeneratorAction.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <string>
#include <cstdlib>

// Параметры по умолчанию
std::string gOutCsv = "g1s_npsm.csv";
long long gEvents = 1000000;
long gSeed = 12345;
double gEnergyKeV = 661.657;
double gSourceZmm = 91.0;
std::string gNpsm = "on";
double gEta = 0.453;
double gSOns = 36.4;
double gSTrap = 12.0;
double gSBirks = 185.3;
double gCutMm = 0.05;
std::string gDeexMode = "deex";
std::string gPixeModel = "";
std::string gPixeElecModel = "";
int gShield = 1;
int gMasses = 1;
// primary=gamma|ion; при ion первичкой служит ядро (ion_z, ion_a),
// которое распадается само — тогда в спектре есть каскадное суммирование.
std::string gPrimary = "gamma";
int gIonZ = 27;
int gIonA = 60;
int gCorrGamma = 0;
// Число рабочих потоков: -1 (умолчание) — все ядра, 0 — однопоточный режим,
// N>0 — ровно N. Ключ threads=N.
int gThreads = -1;

void ParseArgs(int argc, char** argv) {
    std::map<std::string, std::string> args;
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        size_t pos = arg.find('=');
        if (pos == std::string::npos) {
            std::cerr << "Неизвестный ключ: " << arg << std::endl;
            exit(2);
        }
        std::string key = arg.substr(0, pos);
        std::string value = arg.substr(pos + 1);
        args[key] = value;
    }

    // ⚠ Ключи лежат в словаре БЕЗ знака '=' (он отрезан выше). Искать их со
    // знаком — значит не найти ни одного и молча уйти на умолчания: прогон
    // 06.09 с events=10000 отработал 10⁶ событий и записал файл не туда.
    // Тихих умолчаний быть не должно, поэтому ниже ещё и проверка на
    // неизвестный ключ: опечатка обязана останавливать прогон, а не менять
    // постановку незаметно.
    const std::set<std::string> known = {
        "out", "events", "seed", "energy_keV", "src_z_mm", "npsm", "npsm_eta",
        "npsm_s_ons", "npsm_s_trap", "npsm_s_birks", "cut_mm", "deex",
        "pixe_model", "pixe_e_model", "shield", "masses",
        "primary", "ion_z", "ion_a", "corr_gamma", "deex_region", "em",
        "threads"};
    for (const auto& kv : args) {
        if (known.find(kv.first) == known.end()) {
            std::cerr << "Неизвестный ключ: " << kv.first << std::endl;
            exit(2);
        }
    }

    if (args.count("out")) gOutCsv = args["out"];
    if (args.count("events")) gEvents = std::stoll(args["events"]);
    if (args.count("seed")) gSeed = std::stoll(args["seed"]);
    if (args.count("energy_keV")) gEnergyKeV = std::stod(args["energy_keV"]);
    if (args.count("src_z_mm")) gSourceZmm = std::stod(args["src_z_mm"]);
    if (args.count("npsm")) gNpsm = args["npsm"];
    if (args.count("npsm_eta")) gEta = std::stod(args["npsm_eta"]);
    if (args.count("npsm_s_ons")) gSOns = std::stod(args["npsm_s_ons"]);
    if (args.count("npsm_s_trap")) gSTrap = std::stod(args["npsm_s_trap"]);
    if (args.count("npsm_s_birks")) gSBirks = std::stod(args["npsm_s_birks"]);
    if (args.count("cut_mm")) gCutMm = std::stod(args["cut_mm"]);
    if (args.count("deex")) gDeexMode = args["deex"];
    if (args.count("pixe_model")) gPixeModel = args["pixe_model"];
    if (args.count("pixe_e_model")) gPixeElecModel = args["pixe_e_model"];
    if (args.count("shield")) gShield = (args["shield"] == "1") ? 1 : 0;
    if (args.count("masses")) gMasses = (args["masses"] == "1") ? 1 : 0;
    if (args.count("primary")) gPrimary = args["primary"];
    if (args.count("ion_z")) gIonZ = std::stoi(args["ion_z"]);
    if (args.count("ion_a")) gIonA = std::stoi(args["ion_a"]);
    if (args.count("corr_gamma")) gCorrGamma = (args["corr_gamma"] == "1") ? 1 : 0;
    // Адресная деэкситация (линия В) и вариант ЭМ-конструктора (линия Д).
    // Оба ключа уже есть в стенде; здесь добавлены, чтобы постановка на
    // полной геометрии задавалась тем же способом и попадала в шапку CSV.
    if (args.count("threads")) gThreads = std::stoi(args["threads"]);
    if (args.count("deex_region"))
        Rc103FieldPhysicsList::gDeexRegion = args["deex_region"];
    if (args.count("em")) {
        Rc103FieldPhysicsList::gEmOption = args["em"];
        if (Rc103FieldPhysicsList::gEmOption != "opt4" &&
            Rc103FieldPhysicsList::gEmOption != "opt0") {
            std::cerr << "Некорректное значение em: "
                      << Rc103FieldPhysicsList::gEmOption << std::endl;
            exit(2);
        }
    }

    // Проверка корректности значений
    if (gPrimary != "gamma" && gPrimary != "ion") {
        std::cerr << "Некорректное значение primary: " << gPrimary << std::endl;
        exit(2);
    }
    if (gNpsm != "on" && gNpsm != "off") {
        std::cerr << "Некорректное значение npsm: " << gNpsm << std::endl;
        exit(2);
    }
}

int main(int argc, char** argv) {
    ParseArgs(argc, argv);

    // Установка зерна
    G4Random::setTheSeed(gSeed);

    // Установка статических полей физ-листа
    Rc103FieldPhysicsList::gCutMm = gCutMm;
    Rc103FieldPhysicsList::gDeexMode = gDeexMode;
    Rc103FieldPhysicsList::gPixeModel = gPixeModel;
    Rc103FieldPhysicsList::gPixeElecModel = gPixeElecModel;
    Rc103FieldPhysicsList::gDecay = (gPrimary == "ion");
    Rc103FieldPhysicsList::gCorrGamma = (gCorrGamma == 1);

    // Менеджер прогона создаётся ФАБРИКОЙ — так предписывает руководство
    // (Book For Application Developers, пример main()); прежний прямой
    // `new G4RunManager` был единственным местом, где этот проект расходился
    // со стендом, и всплывал подозреваемым при разборе W-066.
    //
    // 07.09.2026: многопоточность включена по умолчанию (решение оператора).
    // Прежний комментарий про «Serial — не отступление, а условие
    // правильности» относился к состоянию, когда NpsmBenchRunAction копил
    // суммы в своих полях без слияния по потокам. Это устранено: накопители
    // переехали в NpsmBenchRun с Merge(), слияние проверено (1/2/4/8/32
    // потока дают идентичные целочисленные счётчики), согласие с
    // однопоточным режимом подтверждено на десяти зёрнах.
    G4RunManager* runManager = nullptr;
    if (gThreads != 0) {
        const int nThreads = (gThreads < 0) ? G4Threading::G4GetNumberOfCores()
                                            : gThreads;
        runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::MT);
        runManager->SetNumberOfThreads(nThreads);
        NpsmBenchRunAction::gThreads = nThreads;
    } else {
        runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
        NpsmBenchRunAction::gThreads = 0;
    }

    // Шапку CSV пишет NpsmBenchRunAction, и поле crystal_mm он берёт из
    // статических полей геометрии СТЕНДА. Здесь геометрия своя (G1SDetector),
    // стендовая не строится, но её умолчания 102x102x406 попали бы в шапку и
    // соврали бы о постановке — #CFG-2 требует, чтобы постановка лежала в
    // файле. Приводим поля к реальному кристаллу Гамма-1С Ø63x63
    // (G1SDetector.hh:21-22, чертёж «Чертеж 63х63.pdf»).
    NpsmBenchDetectorConstruction::gCrystalXMm = 63.0;
    NpsmBenchDetectorConstruction::gCrystalYMm = 63.0;
    NpsmBenchDetectorConstruction::gCrystalZMm = 63.0;

    // Создание детектора
    G1SDetector* detector = new G1SDetector;
    detector->fWithShield = (gShield == 1);
    detector->fWithVessel = false;

    runManager->SetUserInitialization(detector);

    // Создание физического листа
    Rc103FieldPhysicsList* physicsList = new Rc103FieldPhysicsList(gCutMm, gDeexMode);
    runManager->SetUserInitialization(physicsList);

    // Инициализация
    // ⚠ Инициализатор действий регистрируется ДО Initialize(): в
    // многопоточном режиме рабочие потоки берут действия именно из
    // него, и если поставить его после, потоки стартуют без
    // генератора первички — G4Exception Run0032 «G4VUserPrimaryGeneratorAction
    // is not defined», по одному на каждый поток (поймано 07.09.2026).
    // В однопоточном режиме порядок был безразличен, поэтому дефект
    // проявился только при включении многопоточности.
    // Модель светового выхода. Объект создаётся ВСЕГДА, а выключение идёт
    // через SetEnabled — ровно как в стенде run_npsm_bench (его main.cc:158-170).
    // Нулевой указатель сюда передавать нельзя: NpsmBenchSteppingAction и
    // NpsmBenchRunAction обрывают работу на нулевой модели.
    static NpsmLightYield lightYield;
    lightYield.SetParameters(gEta, gSOns, gSTrap, gSBirks);
    lightYield.SetEnabled(gNpsm == "on");

    // Постановка первички — в шапку CSV (см. NpsmBenchRunAction.hh: без этого
    // ионные прогоны выглядели в файле как обычный гамма-прогон).
    NpsmBenchRunAction::gPrimaryKind = gPrimary;
    NpsmBenchRunAction::gIonZ = gIonZ;
    NpsmBenchRunAction::gIonA = gIonA;

    // Поля генератора статические: задаются ДО построения объекта — то есть
    // до того, как инициализатор создаст его на каждом рабочем потоке.
    G1sNpsmPrimaryGeneratorAction::gSourceZmm = gSourceZmm;
    G1sNpsmPrimaryGeneratorAction::gEnergyKeV = gEnergyKeV;
    G1sNpsmPrimaryGeneratorAction::gPrimary = gPrimary;
    G1sNpsmPrimaryGeneratorAction::gIonZ = gIonZ;
    G1sNpsmPrimaryGeneratorAction::gIonA = gIonA;

    // Действия создаёт инициализатор: в многопоточном режиме они нужны по
    // экземпляру на поток, а мастеру — только RunAction, который получает
    // слитый результат и пишет файл. Прямые SetUserAction для этого негодны.
    runManager->SetUserInitialization(new G1sNpsmActionInitialization(
        gOutCsv, gEnergyKeV, gEvents, gSeed, &lightYield));


    runManager->Initialize();

    // Флаги деэкситации ПОСЛЕ инициализации — из самого G4EmParameters.
    // Физ-лист печатает их в конструкторе, то есть ДО того, как другие
    // конструкторы физики успели их переставить. По сообщению разработчика на
    // форуме (тема о G4RadioactiveDecayPhysics в 11.4) распад принудительно
    // включает Оже и игнорирование порогов при инициализации; строка ниже
    // проверяет это фактом на нашей сборке, а не верит форуму на слово.
    {
        auto* em = G4EmParameters::Instance();
        // apply_cuts: по словам разработчика (форум, тема 9388 пакета D) флаг
        // действует ТОЛЬКО на гамма-процессы; для ионизации и тормозного
        // излучения работает лишь значение порога. Печатается, чтобы постановка
        // была видна, а не предполагалась.
        std::printf("POST_INIT_EM: fluo=%d auger=%d pixe=%d ignore_cut=%d "
                    "apply_cuts=%d lowest_e_keV=%.3f\n",
                    static_cast<int>(em->Fluo()), static_cast<int>(em->Auger()),
                    static_cast<int>(em->Pixe()),
                    static_cast<int>(em->DeexcitationIgnoreCut()),
                    static_cast<int>(em->ApplyCuts()),
                    em->LowestElectronEnergy() / keV);
    }

    // Получение логического объема кристалла и установка в stepping action
    G4LogicalVolume* crystalLV = detector->fCrystalLV;
    if (!crystalLV) {
        std::cerr << "Ошибка: fCrystalLV не установлен" << std::endl;
        return 3;
    }
    NpsmBenchSteppingAction::SetCrystalLogicalVolume(crystalLV);

    // Порог «очень долгого распада». Geant4 11.x МОЛЧА считает стабильными
    // нуклиды с периодом полураспада длиннее порога (по умолчанию около года),
    // и Co-60 (5,27 года) просто не распадается: прогон отрабатывает успешно и
    // пишет спектр из нулей. Ровно это и произошло на первом дымовом прогоне
    // 07.09.2026. Приём взят у Am6er (BecqMoni, ветка pie, tools/g4cf/g4cf.cc:760).
    // Команду можно подавать только ПОСЛЕ инициализации ядра.
    if (gPrimary == "ion") {
        G4UImanager::GetUIpointer()->ApplyCommand(
            "/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year");
    }

    // Печать масс
    if (gMasses) {
        detector->ReportMasses();
    }

    // Печать постановки
    std::cout << "G1S_NPSM_SETUP: events=" << gEvents
              << " seed=" << gSeed
              << " energy_keV=" << gEnergyKeV
              << " src_z_mm=" << gSourceZmm
              << " npsm=" << gNpsm
              << " eta=" << gEta
              << " ons=" << gSOns
              << " trap=" << gSTrap
              << " birks=" << gSBirks
              << " cut_mm=" << gCutMm
              << " primary=" << gPrimary
              << " ion=" << gIonZ << "/" << gIonA
              << " corr_gamma=" << gCorrGamma
              << " deex=" << gDeexMode
              << " shield=" << gShield
              << std::endl;

    // Запуск моделирования
    runManager->BeamOn(gEvents);

    delete runManager;
    return 0;
}
