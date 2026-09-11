#include "G4RunManagerFactory.hh"
#include "G4Threading.hh"
#include "G1sNpsmActionInitialization.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

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
#include <cstdio>
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
// Кювета с пробой. "none" (умолчание) — прежнее поведение, голая головка и
// точечный источник. Иначе — пресет сосуда по VesselGeom::Preset и первичка
// из объёма пробы (src=sample). Матрица и плотность задаются явно: у каждого
// комплекта поверки они свои и в шапку CSV обязаны попасть (#CFG-2).
std::string gVessel = "none";
std::string gMatrix = "";      // пусто — умолчание пресета (OISN16)
double gRho = -1.0;            // <0 — умолчание пресета
double gSampleCm3 = -1.0;      // <0 — умолчание пресета
std::string gSrcMode = "point";
// chain=Amin:Amax:Zmin:Zmax — диапазон распадающихся ядер (nucleusLimits).
// Пусто — ограничения нет (прежнее поведение).
std::string gChain = "";
double gMgoRho = -1.0;   // <0 — умолчание геометрии (0,80 г/см³)
// Зазор проба↔стенка колодца, мм. <0 — умолчание пресета (0, вплотную).
double gWellGap = -1.0;
// Область розыгрыша внутри пробы (проверка неоднородности): доля высоты снизу
// и доля радиального размаха снаружи. 1,0 — вся проба, прежнее поведение.
double gSrcZFrac = 1.0;
double gSrcRFrac = 1.0;

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
        "threads",
        // Проба в кювете (11.09.2026, линия смеси AmTiCsEu). До этого проект
        // строил только «голую» головку: источник точечный, сосуда нет. Для
        // поверочного комплекта это упрощение недопустимо — активность
        // засыпана в кювету, и самопоглощение в матрице входит в отклик.
        "vessel", "matrix", "rho", "sample_cm3", "src",
        // Ограничение цепочки распада (11.09.2026). Без него порог
        // thresholdForVeryLongDecayTime=1e60 year, поставленный ради Co-60,
        // распускает ВСЮ цепочку: у Am-241 немедленно распадается дочерний
        // Np-237 (T½ = 2,14 млн лет), которого в реальном источнике за 24 года
        // накопилось ничтожно мало. Шаблон нуклида обязан содержать только его
        // собственный вклад — тот же приём nucleusLimits, что в макросе
        // decay_amticseu_isotopes.mac.
        "chain",
        // Плотность порошкового отражателя MgO. Ключ нужен НЕ для подгонки:
        // калибровать её по мягким линиям оператор запретил намеренно
        // (G1SDetector.hh:71-74 — подогнанная плотность впитает всё прочее
        // неверное в торце и станет эффективной массовой толщиной под чужим
        // именем). Ключ нужен, чтобы ИЗМЕРИТЬ систематику от физического
        // разброса 0,5…1,2 г/см³, а не оценивать её по формуле.
        "mgo_rho", "well_gap", "src_z_frac", "src_r_frac"};
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
    if (args.count("vessel")) gVessel = args["vessel"];
    if (args.count("matrix")) gMatrix = args["matrix"];
    if (args.count("rho")) gRho = std::stod(args["rho"]);
    if (args.count("sample_cm3")) gSampleCm3 = std::stod(args["sample_cm3"]);
    if (args.count("src")) gSrcMode = args["src"];
    if (args.count("chain")) gChain = args["chain"];
    if (args.count("mgo_rho")) gMgoRho = std::stod(args["mgo_rho"]);
    if (args.count("well_gap")) gWellGap = std::stod(args["well_gap"]);
    if (args.count("src_z_frac")) gSrcZFrac = std::stod(args["src_z_frac"]);
    if (args.count("src_r_frac")) gSrcRFrac = std::stod(args["src_r_frac"]);
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
    // Пресеты сосудов — те же имена, что понимает VesselGeom::Preset. Список
    // закрытый: опечатка в имени иначе молча дала бы пресет по умолчанию
    // (Preset возвращает marinelli на неизвестное имя), то есть НЕ ту
    // геометрию при внешне успешном прогоне.
    if (gVessel != "none" && gVessel != "marinelli" &&
        gVessel != "marinelli_lsrm" && gVessel != "denta" && gVessel != "petri") {
        std::cerr << "Некорректное значение vessel: " << gVessel << std::endl;
        exit(2);
    }
    if (gSrcMode != "point" && gSrcMode != "sample") {
        std::cerr << "Некорректное значение src: " << gSrcMode << std::endl;
        exit(2);
    }
    // Розыгрыш по объёму пробы без построенной кюветы невозможен — отказ
    // здесь, до прогона, а не G4Exception на первом событии.
    if (gSrcMode == "sample" && gVessel == "none") {
        std::cerr << "src=sample требует vessel≠none" << std::endl;
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
    // ⚠ 10.09.2026: габариты сами по себе шапку не спасали. Поля названы
    // X/Y/Z, и объём считался их произведением — 250 см³ вместо реальных
    // 196,4 у цилиндра Ø63x63 (завышение 27 %). Форма и объём задаются явно.
    NpsmBenchDetectorConstruction::gCrystalXMm = 63.0;
    NpsmBenchDetectorConstruction::gCrystalYMm = 63.0;
    NpsmBenchDetectorConstruction::gCrystalZMm = 63.0;
    NpsmBenchDetectorConstruction::gCrystalShape = "cylinder_D63xH63";
    {
        const double rMm = 0.5 * 63.0, hMm = 63.0;
        NpsmBenchDetectorConstruction::gCrystalVolumeCm3 =
            CLHEP::pi * rMm * rMm * hMm / 1000.0;   // 196,35 см³
    }

    // Создание детектора
    G1SDetector* detector = new G1SDetector;
    detector->fWithShield = (gShield == 1);
    detector->fWithVessel = (gVessel != "none");
    if (gMgoRho > 0) detector->fHead.mgoDensity = gMgoRho;
    if (detector->fWithVessel) {
        detector->fVessel = VesselGeom::Preset(gVessel);
        if (!gMatrix.empty()) detector->fVessel.sampleMatrix = gMatrix;
        if (gRho > 0) detector->fVessel.sampleDensity = gRho;
        if (gSampleCm3 > 0) detector->fVessel.sampleCm3 = gSampleCm3;
        if (gWellGap >= 0) detector->fVessel.sampleGap = gWellGap;
    }

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
    // Постановка пробы — в шапку CSV. Значения берутся из ПОСТРОЕННОГО
    // пресета, а не из ключей: пресет мог подставить умолчание там, где ключ
    // не задан, и в файле должно лежать применённое, а не запрошенное (W-052).
    NpsmBenchRunAction::gVessel = gVessel;
    NpsmBenchRunAction::gSrcMode = gSrcMode;
    NpsmBenchRunAction::gChainLimits = gChain;
    if (detector->fWithVessel) {
        NpsmBenchRunAction::gMatrix = detector->fVessel.sampleMatrix;
        NpsmBenchRunAction::gRhoSample = detector->fVessel.sampleDensity;
        NpsmBenchRunAction::gSampleCm3 = detector->fVessel.sampleCm3;
    }

    // Поля генератора статические: задаются ДО построения объекта — то есть
    // до того, как инициализатор создаст его на каждом рабочем потоке.
    G1sNpsmPrimaryGeneratorAction::gSourceZmm = gSourceZmm;
    G1sNpsmPrimaryGeneratorAction::gEnergyKeV = gEnergyKeV;
    G1sNpsmPrimaryGeneratorAction::gPrimary = gPrimary;
    G1sNpsmPrimaryGeneratorAction::gIonZ = gIonZ;
    G1sNpsmPrimaryGeneratorAction::gIonA = gIonA;
    G1sNpsmPrimaryGeneratorAction::gSourceMode = gSrcMode;
    G1sNpsmPrimaryGeneratorAction::gSrcZFrac = gSrcZFrac;
    G1sNpsmPrimaryGeneratorAction::gSrcRFrac = gSrcRFrac;

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
        // Ограничение цепочки — ПОСЛЕ снятия порога и до BeamOn. Порядок
        // важен: порог распускает всю цепочку, nucleusLimits её обратно
        // сужает до заданного нуклида.
        if (!gChain.empty()) {
            int aMin = 0, aMax = 0, zMin = 0, zMax = 0;
            if (std::sscanf(gChain.c_str(), "%d:%d:%d:%d",
                            &aMin, &aMax, &zMin, &zMax) != 4) {
                std::cerr << "Некорректный chain (нужно Amin:Amax:Zmin:Zmax): "
                          << gChain << std::endl;
                return 5;
            }
            char cmd[128];
            std::snprintf(cmd, sizeof(cmd),
                          "/process/had/rdm/nucleusLimits %d %d %d %d",
                          aMin, aMax, zMin, zMax);
            G4UImanager::GetUIpointer()->ApplyCommand(cmd);
            std::printf("CHAIN_LIMITS: %s\n", cmd);
        }
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
              << " mgo_rho=" << detector->fHead.mgoDensity
              << " vessel=" << gVessel
              << " src=" << gSrcMode
              << " src_z_frac=" << gSrcZFrac << " src_r_frac=" << gSrcRFrac;
    if (detector->fWithVessel) {
        std::cout << " matrix=" << detector->fVessel.sampleMatrix
                  << " rho=" << detector->fVessel.sampleDensity
                  << " sample_cm3=" << detector->fVessel.sampleCm3
                  << " sample_built_cm3=" << detector->fSampleVolumeCm3
                  << " well_gap=" << detector->fVessel.sampleGap
                  << " sample_fits=" << (detector->fSampleFits ? 1 : 0);
    }
    std::cout << std::endl;
    // Уровень засыпки выше крышки — геометрия построена, но физически
    // невозможна. Отказ, а не предупреждение: прогон дал бы правдоподобный
    // спектр не той пробы (#SA-6 — «построилось» не значит «верно»).
    if (detector->fWithVessel && !detector->fSampleFits) {
        std::cerr << "ОТКАЗ: уровень засыпки выше крышки сосуда" << std::endl;
        return 4;
    }

    // Запуск моделирования
    runManager->BeamOn(gEvents);

    delete runManager;
    return 0;
}
