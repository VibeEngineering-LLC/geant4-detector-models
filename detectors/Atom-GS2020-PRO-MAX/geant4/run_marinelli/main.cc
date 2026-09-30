#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <atomic>
#include <memory>

#include "G4RunManagerFactory.hh"
#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmParameters.hh"
#include "G4ThreeVector.hh"
#include "G4IonTable.hh"
#include "G4UImanager.hh"
#include "G4UserStackingAction.hh"
#include "G4Track.hh"
#include "G4Gamma.hh"
#include "G4VProcess.hh"
#include "G4Exception.hh"
#include "decay_tally.hh"  // #GS-49 (30.09): GS2020_DECAY_TALLY=<файл> — учёт состава фотонов/e⁻ конверсии по распадам, без переноса
#include "ext_source.hh"   // #GS-44 (29.09): внешнее изотропное поле со сферы, env GS2020_EXT_R_MM / _CZ_MM / _HEMI

// GS2020_GDML (оператор 26.09, #CFG-1: проверка d_эфф против ЛСРМ при разной плотности матрицы Epoxy_crumb) —
// путь к GDML переопределяется, если задан; иначе дефолт как раньше.
std::string ResolveGdmlPath() {
    const char* env = std::getenv("GS2020_GDML");
    return (env && *env) ? std::string(env) : std::string("C:/g4work/gs2020/GS2020_marinelli_th232.gdml");
}
const std::string kGdmlPath = ResolveGdmlPath();
static G4VSolid* gSourceSolid = nullptr;
static G4LogicalVolume* gSourceLV = nullptr;
static G4LogicalVolume* gCrystalLV = nullptr;
// GS2020_PRIMARY_TABLE (29.09.2026): первичный источник — фотон, энергия которого
// разыгрывается по табличному спектру CSV "k_keV,dNdk_per_decay_per_keV" (формат и
// проверки — донор run_g1s_npsm/G1sNpsmPrimaryGeneratorAction.{hh,cc}, LoadSpectrumTable/
// SampleTableEnergyKeV/PrintSpectrumDrawCounts, перенесено дословно). Позиция и
// направление — код PrimaryGeneratorAction::GeneratePrimaries (не дублируются).
// Если переменная окружения не задана, gUseTable=false, ниже ничего не исполняется —
// поведение и расход RNG побитово прежние.
namespace GammaTable {
std::vector<double> gK, gCdf;
double gStep = 0.0;
std::unique_ptr<std::atomic<long long>[]> gDrawn;

void Fail(const std::string& code, const std::string& msg) {
    G4Exception("GammaTable::Load", code.c_str(), FatalException, msg.c_str());
}
std::string Trim(const std::string& s) {
    const size_t b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    const size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}
bool ParseDouble(const std::string& s, double& out) {
    const std::string t = Trim(s);
    if (t.empty()) return false;
    char* end = nullptr;
    out = std::strtod(t.c_str(), &end);
    return end == t.c_str() + t.size() && std::isfinite(out);
}

void Load(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) Fail("TableNoFile", "GS2020_PRIMARY_TABLE: файл не открывается: " + path);
    std::vector<double> k, dens;
    std::string line; bool header = false; int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        const std::string t = Trim(line);
        if (t.empty() || t[0] == '#') continue;
        if (!header) {
            if (t != "k_keV,dNdk_per_decay_per_keV")
                Fail("TableBadHeader", "GS2020_PRIMARY_TABLE: строка " + std::to_string(lineNo) + ": ожидался заголовок k_keV,dNdk_per_decay_per_keV, получено: " + t);
            header = true; continue;
        }
        const size_t c = t.find(',');
        double kv = 0.0, dv = 0.0;
        if (c == std::string::npos || !ParseDouble(t.substr(0, c), kv) || !ParseDouble(t.substr(c + 1), dv))
            Fail("TableBadRow", "GS2020_PRIMARY_TABLE: строка " + std::to_string(lineNo) + " не разбирается как два числа: " + t);
        if (dv < 0.0)
            Fail("TableNegative", "GS2020_PRIMARY_TABLE: строка " + std::to_string(lineNo) + ": отрицательная плотность " + t);
        if (!k.empty() && !(kv > k.back()))
            Fail("TableNotIncreasing", "GS2020_PRIMARY_TABLE: строка " + std::to_string(lineNo) + ": энергия не возрастает " + t);
        if (kv <= 0.0)
            Fail("TableBadEnergy", "GS2020_PRIMARY_TABLE: строка " + std::to_string(lineNo) + ": энергия не положительна " + t);
        k.push_back(kv); dens.push_back(dv);
    }
    if (!header) Fail("TableBadHeader", "GS2020_PRIMARY_TABLE: нет заголовка: " + path);
    if (k.size() < 2) Fail("TableTooShort", "GS2020_PRIMARY_TABLE: строк данных меньше 2 (шаг сетки не определить): " + path);
    double step = k[1] - k[0];
    for (size_t i = 2; i < k.size(); ++i) step = std::min(step, k[i] - k[i - 1]);
    if (k.front() - 0.5 * step < 0.0)
        Fail("TableBadEnergy", "GS2020_PRIMARY_TABLE: левый край первого бина отрицателен: " + path);
    std::vector<double> cdf(k.size());
    double acc = 0.0;
    for (size_t i = 0; i < k.size(); ++i) { acc += dens[i] * step; cdf[i] = acc; }
    if (!(acc > 0.0) || !std::isfinite(acc))
        Fail("TableZeroIntegral", "GS2020_PRIMARY_TABLE: нулевой интеграл плотности: " + path);
    gK = std::move(k); gCdf = std::move(cdf); gStep = step;
    gDrawn.reset(new std::atomic<long long>[gK.size()]);
    for (size_t i = 0; i < gK.size(); ++i) gDrawn[i] = 0;
    std::printf("GAMMA_TABLE_LOADED: rows=%d kmin_keV=%.6g kmax_keV=%.6g step_keV=%.6g integral_per_decay=%.9g\n",
                (int)gK.size(), gK.front(), gK.back(), gStep, acc);
}

double Sample() {
    const double u = G4UniformRand() * gCdf.back();
    size_t i = static_cast<size_t>(std::upper_bound(gCdf.begin(), gCdf.end(), u) - gCdf.begin());
    if (i >= gCdf.size()) i = gCdf.size() - 1;
    gDrawn[i].fetch_add(1, std::memory_order_relaxed);
    return gK[i] + (G4UniformRand() - 0.5) * gStep;
}
long long DrawnTotal() {
    if (gK.empty()) return 0;
    long long total = 0;
    for (size_t i = 0; i < gK.size(); ++i) total += gDrawn[i].load();
    return total;
}
void PrintDrawCounts() {
    if (gK.empty()) return;
    std::printf("gamma_table: drawn_total=%lld\n", DrawnTotal());
    if (gK.size() > 20) return;
    for (size_t i = 0; i < gK.size(); ++i)
        std::printf("gamma_table: bin=%zu k_keV=%.6g n=%lld\n", i, gK[i], gDrawn[i].load());
}
}  // namespace GammaTable

static std::string ResolvePrimaryTablePath() {
    const char* env = std::getenv("GS2020_PRIMARY_TABLE");
    return (env && *env) ? std::string(env) : std::string();
}
const std::string gPrimaryTablePath = ResolvePrimaryTablePath();
const bool gUseTable = !gPrimaryTablePath.empty();

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    G4VPhysicalVolume* Construct() override {
        G4GDMLParser parser;
        parser.SetOverlapCheck(false);
        parser.Read(kGdmlPath, false);

        G4VPhysicalVolume* world = parser.GetWorldVolume();

        G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();
        gSourceLV = store->GetVolume("LV_SourceMatrix");
        gCrystalLV = store->GetVolume("LV_NaI_crystal");

        if (!gSourceLV || !gCrystalLV) {
            std::fprintf(stderr, "FATAL: volume not found\n");
            std::exit(3);
        }

        gSourceSolid = gSourceLV->GetSolid();
        ExtSource::CheckEnclosure(world);   // в режиме внешнего поля: сфера внутри мира и охватывает все объёмы, иначе FATAL

        return world;
    }
};

class GS2020PhysicsList : public G4VModularPhysicsList {
public:
    GS2020PhysicsList() {
        RegisterPhysics(new G4EmStandardPhysics_option4());
        RegisterPhysics(new G4DecayPhysics());
        RegisterPhysics(new G4RadioactiveDecayPhysics());

        SetDefaultCutValue(0.05*mm);
        G4EmParameters* p = G4EmParameters::Instance();
        p->SetFluo(true);
        p->SetAuger(true);
        p->SetPixe(false);
        p->SetDeexcitationIgnoreCut(true);
    }

    void SetCuts() override {
        G4VUserPhysicsList::SetCuts();
    }
};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
private:
    G4ParticleGun* fGun;
    double fEnergyKeV;
    int fIonZ;
    int fIonA;
    G4ParticleDefinition* fIon;

public:
    PrimaryGeneratorAction(double energyKeV, int ionZ, int ionA) : fEnergyKeV(energyKeV), fIonZ(ionZ), fIonA(ionA), fIon(nullptr) {
        fGun = new G4ParticleGun(1);
        if (fIonZ > 0 && !gUseTable) {
            // ион берётся лениво в GeneratePrimaries: до Initialize() таблица ионов не готова
        } else {
            G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
            G4String particleName = "gamma";
            G4ParticleDefinition* particleDef = particleTable->FindParticle(particleName);
            fGun->SetParticleDefinition(particleDef);
        }
    }

    void GeneratePrimaries(G4Event* event) override {
        if (ExtSource::Get().on) {   // #GS-44: фотон с поверхности сферы внутрь по косинусному закону (ext_source.hh)
            G4ThreeVector p, d;
            ExtSource::Sample(p, d);
            fGun->SetParticlePosition(p);
            fGun->SetParticleMomentumDirection(d);
            fGun->SetParticleEnergy(fEnergyKeV * keV);
            fGun->GeneratePrimaryVertex(event);
            return;
        }
        const double R = 70.5 * mm;
        const double HALF_Z = 52.5 * mm;

        G4ThreeVector position;
        for (;;) {
            double r = R * std::sqrt(G4UniformRand());
            double phi = CLHEP::twopi * G4UniformRand();
            double x = r * std::cos(phi);
            double y = r * std::sin(phi);
            double z = (2.0 * G4UniformRand() - 1.0) * HALF_Z;
            G4ThreeVector pLocal(x, y, z);

            if (gSourceSolid->Inside(pLocal) == kInside) {
                position = pLocal;
                break;
            }
        }

        // Translate to world coordinates
        G4ThreeVector worldPos = position + G4ThreeVector(0, 0, 40.5 * mm);

        double ct = 2.0 * G4UniformRand() - 1.0;
        double st = std::sqrt(1.0 - ct * ct);
        double ph2 = CLHEP::twopi * G4UniformRand();
        G4ThreeVector dir(st * std::cos(ph2), st * std::sin(ph2), ct);

        fGun->SetParticlePosition(worldPos);
        fGun->SetParticleMomentumDirection(dir);
        if (fIonZ > 0 && !gUseTable) {
            if (!fIon) {
                fIon = G4IonTable::GetIonTable()->GetIon(fIonZ, fIonA, 0.0);
                if (!fIon) { std::fprintf(stderr, "FATAL: ion Z=%d A=%d not found\n", fIonZ, fIonA); std::abort(); }
                fGun->SetParticleDefinition(fIon);
            }
            fGun->SetParticleEnergy(0.0 * keV);
        } else if (gUseTable) {
            // GS2020_PRIMARY_TABLE: энергия из таблицы (донор: SampleTableEnergyKeV);
            // позиция/направление — код выше, общий с обычным режимом.
            fGun->SetParticleEnergy(GammaTable::Sample() * keV);
        } else {
            fGun->SetParticleEnergy(fEnergyKeV * keV);
        }
        fGun->GeneratePrimaryVertex(event);
    }
};

class SteppingAction : public G4UserSteppingAction {
public:
    static double gEdepThisEvent;

    void UserSteppingAction(const G4Step* step) override {
        G4LogicalVolume* lv = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
        if (lv == gCrystalLV) {
            gEdepThisEvent += step->GetTotalEnergyDeposit() / keV;
        }
    }
};

double SteppingAction::gEdepThisEvent = 0.0;

// #GS-32 (29.09): GS2020_BETA_ONLY=1 — β/e⁻-компонента для метода 2: фотоны, рождённые самим распадом (γ и рентген —
// они уже в библиотеке линий), убиваются при рождении; e⁻/e⁺ переносятся, тормозное и возбуждённый рентген рождаются.
static const bool gBetaOnly = std::getenv("GS2020_BETA_ONLY") && std::string(std::getenv("GS2020_BETA_ONLY")) == "1";
static long long gKilledDecayPhotons = 0;
class BetaOnlyStackingAction : public G4UserStackingAction {
public:
    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* t) override {
        if (t->GetDefinition() == G4Gamma::Definition() && t->GetCreatorProcess()
            && t->GetCreatorProcess()->GetProcessName() == "Radioactivation") { ++gKilledDecayPhotons; return fKill; }
        return fUrgent;
    }
};

class RunAction : public G4UserRunAction {
private:
    std::vector<long long> fHist;
    std::string fOutCsv;
    double fEnergyKeV;
    long long fSeed;
    long long fNEvents;
    int fIonZ;
    int fIonA;
    long long fNWithEdep;
    std::string fPrimaryTable;

public:
    RunAction(const std::string& outCsv, double energyKeV, long long seed, long long nEvents, int ionZ, int ionA, const std::string& primaryTable)
        : fOutCsv(outCsv), fEnergyKeV(energyKeV), fSeed(seed), fNEvents(nEvents), fIonZ(ionZ), fIonA(ionA), fPrimaryTable(primaryTable) {
        fHist.resize(3700, 0);
        fNWithEdep = 0;
    }

    void Fill(double edepKeV) {
        if (edepKeV > 0.5) {
            fNWithEdep++;
            int ch = (int)std::floor(edepKeV);
            if (ch >= 0 && ch <= 3699) {
                fHist[ch]++;
            }
        }
    }

    void EndOfRunAction(const G4Run* run) override {
        {
            std::ofstream f(fOutCsv);
            f << "# gs2020_marinelli volume source in Marinelli, GDML " << kGdmlPath << "\n";
            f << "energy_keV," << fEnergyKeV << "\n";
            f << "n_events_requested," << fNEvents << "\n";
            f << "n_events_processed," << run->GetNumberOfEvent() << "\n";
            f << "seed," << fSeed << "\n";
            f << "em_cut_mm,0.05\n";
            f << "em_deex,deex\n";
            f << "em_option,opt4\n";
            f << "particle," << ((fIonZ > 0 && !gUseTable) ? "ion" : "gamma") << "\n";
            f << "ion_z," << fIonZ << "\n";
            f << "ion_a," << fIonA << "\n";
            f << "decay," << ((fIonZ > 0 && !gUseTable) ? 1 : 0) << "\n";
            f << "beta_only," << (gBetaOnly ? 1 : 0) << "\n";
            f << "killed_decay_photons," << gKilledDecayPhotons << "\n";
            f << "primary_table," << fPrimaryTable << "\n";
            f << "table_drawn," << GammaTable::DrawnTotal() << "\n";
            f << "vessel,marinelli\n";
            if (ExtSource::Get().on) {   // #GS-44: шапка режима внешнего поля; в режиме пробы строки прежние (побитовая регрессия)
                f << "src_mode,ext_sphere\n";
                f << "ext_R_mm," << ExtSource::Get().R / mm << "\n";
                f << "ext_cz_mm," << ExtSource::Get().cz / mm << "\n";
                f << "ext_hemi," << ExtSource::Get().hemi << "\n";
            } else {
                f << "src_mode,sample\n";
            }
            f << "sample_matrix," << gSourceLV->GetMaterial()->GetName() << "\n";   // из GDML (было зашито Epoxy_crumb)
            f << "sample_rho_g_cm3," << (gSourceLV->GetMaterial()->GetDensity() / (g / cm3)) << "\n";  // читаем из GDML, не хардкод (#CFG-1, 26.09)
            f << "npsm_enabled,0\n";
            f << "n_with_edep," << fNWithEdep << "\n";
            f << "bin_keV,count_edep,count_light\n";

            for (int i = 0; i < 3700; ++i) {
                f << i + 0.5 << "," << fHist[i] << "," << fHist[i] << "\n";
            }
            f.close();

            long long total = 0;
            for (long long x : fHist) total += x;
            std::cout << "[gs2020_marinelli] written " << fOutCsv << " total_hits=" << total << "\n";
        }
    }
};

class EventAction : public G4UserEventAction {
private:
    RunAction* fRunAction;

public:
    EventAction(RunAction* runAction) : fRunAction(runAction) {}

    void BeginOfEventAction(const G4Event*) override {
        SteppingAction::gEdepThisEvent = 0.0;
        if (DecayTally::Enabled()) DecayTally::BeginEvent();
    }

    void EndOfEventAction(const G4Event*) override {
        fRunAction->Fill(SteppingAction::gEdepThisEvent);
        if (DecayTally::Enabled()) DecayTally::EndEvent();
    }
};

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "usage: gs2020_marinelli.exe <energy_keV|ion:Z:A> <n_events> <out_csv> [seed]\n");
        return 2;
    }

    double energyKeV = 0.0;
    int ionZ = 0;
    int ionA = 0;

    std::string arg1(argv[1]);
    if (arg1.rfind("ion:", 0) == 0) {
        size_t pos1 = arg1.find(':', 4);
        if (pos1 != std::string::npos) {
            ionZ = std::stoi(arg1.substr(4, pos1 - 4));
            ionA = std::stoi(arg1.substr(pos1 + 1));
        } else {
            std::fprintf(stderr, "FATAL: invalid ion format\n");
            return 3;
        }
    } else {
        energyKeV = std::stod(argv[1]);
    }

    long long nEvents = std::stoll(argv[2]);
    std::string outCsv = std::string(argv[3]);
    unsigned long seed = (argc > 4) ? std::stoul(argv[4]) : (unsigned long)std::time(nullptr);

    G4Random::setTheSeed(seed);
    if (ExtSource::Get().on && (ionZ > 0 || gUseTable)) {   // #GS-44: внешнее поле — только монолиния фотонов
        std::fprintf(stderr, "FATAL: GS2020_EXT_R_MM with ion or primary table\n");
        return 6;
    }

    // GS2020_PRIMARY_TABLE: таблица читается один раз, до BeamOn (донор грузит её
    // на мастере до создания рабочих потоков; здесь run manager Serial — тем более
    // до всякого использования GammaTable::Sample() в GeneratePrimaries).
    if (gUseTable) {
        GammaTable::Load(gPrimaryTablePath);
    }

    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
    runManager->SetUserInitialization(new DetectorConstruction());
    runManager->SetUserInitialization(new GS2020PhysicsList());

    std::cout << "cut_mm=0.0500 deex=deex fluo=1 auger=1 pixe=0 ignore_cut=1\n";

    runManager->SetUserAction(new PrimaryGeneratorAction(energyKeV, ionZ, ionA));
    RunAction* runAction = new RunAction(outCsv, energyKeV, seed, nEvents, ionZ, ionA, gPrimaryTablePath);
    runManager->SetUserAction(runAction);
    runManager->SetUserAction(new EventAction(runAction));
    runManager->SetUserAction(new SteppingAction());
    if (gBetaOnly) runManager->SetUserAction(new BetaOnlyStackingAction());
    if (DecayTally::Enabled()) {   // #GS-49: режим учёта распада, перенос выключен (убиваются все не-ионы при рождении)
        if (gBetaOnly || ionZ <= 0 || gUseTable) { std::fprintf(stderr, "REFUSED: DECAY_TALLY only for ion:Z:A without BETA_ONLY/table\n"); return 4; }
        runManager->SetUserAction(new DecayTally::StackingAction());
    }
    std::cout << "beta_only=" << (gBetaOnly ? 1 : 0) << "\n";
    std::cout << "primary_table=" << gPrimaryTablePath << "\n";

    runManager->Initialize();

    if (ionZ > 0 && !gUseTable) {
        G4UImanager* ui = G4UImanager::GetUIpointer();
        ui->ApplyCommand("/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+30 ns");
        ui->ApplyCommand("/process/had/rdm/nucleusLimits " + std::to_string(ionA) + " " + std::to_string(ionA) + " " + std::to_string(ionZ) + " " + std::to_string(ionZ));
        std::cout << "CHAIN_CUT: nucleusLimits A=" << ionA << " Z=" << ionZ << "\n";
    }

    runManager->BeamOn(nEvents);
    if (DecayTally::Enabled()) DecayTally::Write(seed, ionZ, ionA);
    if (gUseTable) GammaTable::PrintDrawCounts();
    ExtSource::Report();   // #GS-44: доля траекторий с прицельным параметром < R/4 против a²/R² (косинусный закон)

    delete runManager;
    return 0;
}
