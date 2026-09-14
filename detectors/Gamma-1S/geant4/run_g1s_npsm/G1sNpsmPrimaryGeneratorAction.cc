#include "G1sNpsmPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4IonTable.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Exception.hh"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <algorithm>
#ifdef _WIN32
// Последним: макросы windows.h не должны задевать заголовки Geant4 выше.
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

double G1sNpsmPrimaryGeneratorAction::gSourceZmm = 91.0;
double G1sNpsmPrimaryGeneratorAction::gEnergyKeV = 661.657;
std::string G1sNpsmPrimaryGeneratorAction::gPrimary = "gamma";
int G1sNpsmPrimaryGeneratorAction::gIonZ = 27;   // Co
int G1sNpsmPrimaryGeneratorAction::gIonA = 60;   // Co-60
std::string G1sNpsmPrimaryGeneratorAction::gSourceMode = "point";
double G1sNpsmPrimaryGeneratorAction::gSrcZFrac = 1.0;
double G1sNpsmPrimaryGeneratorAction::gSrcRFrac = 1.0;
std::string G1sNpsmPrimaryGeneratorAction::gSpectrumCsv = "";
std::vector<double> G1sNpsmPrimaryGeneratorAction::gTabK;
std::vector<double> G1sNpsmPrimaryGeneratorAction::gTabCdf;
double G1sNpsmPrimaryGeneratorAction::gTabStep = 0.0;
std::unique_ptr<std::atomic<long long>[]> G1sNpsmPrimaryGeneratorAction::gTabDrawn;

namespace {
void TableFail(const std::string& code, const std::string& msg) {
    // Громкий отказ: таблица, прочитанная «как-нибудь», дала бы правдоподобный
    // спектр не той формы при внешне успешном прогоне.
    G4Exception("G1sNpsmPrimaryGeneratorAction::LoadSpectrumTable",
                code.c_str(), FatalException, msg.c_str());
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
}  // namespace

void G1sNpsmPrimaryGeneratorAction::LoadSpectrumTable(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open())
        TableFail("TableNoFile", "spectrum_csv: файл не открывается: " + path);
    std::vector<double> k, dens;
    std::string line;
    bool header = false;
    int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        const std::string t = Trim(line);
        if (t.empty() || t[0] == '#') continue;
        if (!header) {
            if (t != "k_keV,dNdk_per_decay_per_keV")
                TableFail("TableBadHeader", "spectrum_csv: строка " + std::to_string(lineNo)
                          + ": ожидался заголовок k_keV,dNdk_per_decay_per_keV, получено: " + t);
            header = true;
            continue;
        }
        const size_t c = t.find(',');
        double kv = 0.0, dv = 0.0;
        if (c == std::string::npos || !ParseDouble(t.substr(0, c), kv)
            || !ParseDouble(t.substr(c + 1), dv))
            TableFail("TableBadRow", "spectrum_csv: строка " + std::to_string(lineNo)
                      + " не разбирается как два числа: " + t);
        if (dv < 0.0)
            TableFail("TableNegative", "spectrum_csv: строка " + std::to_string(lineNo)
                      + ": отрицательная плотность " + t);
        if (!k.empty() && !(kv > k.back()))
            TableFail("TableNotIncreasing", "spectrum_csv: строка " + std::to_string(lineNo)
                      + ": энергия не возрастает " + t);
        if (kv <= 0.0)
            TableFail("TableBadEnergy", "spectrum_csv: строка " + std::to_string(lineNo)
                      + ": энергия не положительна " + t);
        k.push_back(kv);
        dens.push_back(dv);
    }
    if (!header) TableFail("TableBadHeader", "spectrum_csv: нет заголовка: " + path);
    if (k.size() < 2)
        TableFail("TableTooShort", "spectrum_csv: строк данных меньше 2 (шаг сетки не определить): " + path);
    // Шаг сетки — наименьшее расстояние между соседними центрами.
    double step = k[1] - k[0];
    for (size_t i = 2; i < k.size(); ++i) step = std::min(step, k[i] - k[i - 1]);
    if (k.front() - 0.5 * step < 0.0)
        TableFail("TableBadEnergy", "spectrum_csv: левый край первого бина отрицателен: " + path);
    std::vector<double> cdf(k.size());
    double acc = 0.0;
    for (size_t i = 0; i < k.size(); ++i) {
        acc += dens[i] * step;   // вес бина = плотность × ширина
        cdf[i] = acc;
    }
    if (!(acc > 0.0) || !std::isfinite(acc))
        TableFail("TableZeroIntegral", "spectrum_csv: нулевой интеграл плотности: " + path);
    gTabK = std::move(k);
    gTabCdf = std::move(cdf);
    gTabStep = step;
    gTabDrawn.reset(new std::atomic<long long>[gTabK.size()]);
    for (size_t i = 0; i < gTabK.size(); ++i) gTabDrawn[i] = 0;
    std::printf("GAMMA_TABLE_LOADED: rows=%d kmin_keV=%.6g kmax_keV=%.6g step_keV=%.6g "
                "integral_per_decay=%.9g\n",
                SpectrumRows(), SpectrumKminKeV(), SpectrumKmaxKeV(), gTabStep, acc);
}

double G1sNpsmPrimaryGeneratorAction::SampleTableEnergyKeV() {
    // Обратная функция распределения: равномерное число на [0, полный вес),
    // первый бин, у которого накопленный вес его превышает.
    const double u = G4UniformRand() * gTabCdf.back();
    size_t i = static_cast<size_t>(
        std::upper_bound(gTabCdf.begin(), gTabCdf.end(), u) - gTabCdf.begin());    if (i >= gTabCdf.size()) i = gTabCdf.size() - 1;   // u == полный вес
    gTabDrawn[i].fetch_add(1, std::memory_order_relaxed);
    return gTabK[i] + (G4UniformRand() - 0.5) * gTabStep;
}

std::string G1sNpsmPrimaryGeneratorAction::SpectrumCsvUtf8() {
#ifdef _WIN32
    const std::string& s = gSpectrumCsv;
    if (s.empty()) return s;
    int wn = MultiByteToWideChar(CP_ACP, 0, s.c_str(), -1, nullptr, 0);
    if (wn <= 0) return s;
    std::wstring w(static_cast<size_t>(wn), L'\0');
    MultiByteToWideChar(CP_ACP, 0, s.c_str(), -1, &w[0], wn);
    int un = WideCharToMultiByte(CP_UTF8, 0, w.c_str(), -1, nullptr, 0, nullptr, nullptr);
    if (un <= 0) return s;
    std::string u(static_cast<size_t>(un), '\0');
    WideCharToMultiByte(CP_UTF8, 0, w.c_str(), -1, &u[0], un, nullptr, nullptr);
    u.resize(static_cast<size_t>(un - 1));
    return u;
#else
    return gSpectrumCsv;
#endif
}

void G1sNpsmPrimaryGeneratorAction::PrintSpectrumDrawCounts() {
    if (gTabK.empty()) return;
    long long total = 0;
    for (size_t i = 0; i < gTabK.size(); ++i) total += gTabDrawn[i].load();
    std::printf("gamma_table: drawn_total=%lld\n", total);
    if (gTabK.size() > 20) return;
    for (size_t i = 0; i < gTabK.size(); ++i)
        std::printf("gamma_table: bin=%zu k_keV=%.6g n=%lld\n",
                    i, gTabK[i], gTabDrawn[i].load());
}

G4ThreeVector G1sNpsmPrimaryGeneratorAction::SamplePointInSample() {
    if (!fSampleSolid) {
        G4LogicalVolume* lv =
            G4LogicalVolumeStore::GetInstance()->GetVolume("Sample", false);
        if (!lv) {
            // Громкий отказ, а не тихий возврат нуля: молча вылетающая из
            // центра мира первичка дала бы правдоподобный, но неверный спектр.
            G4Exception("G1sNpsmPrimaryGeneratorAction::SamplePointInSample",
                        "NoSample", FatalException,
                        "src=sample, но логического объёма Sample нет: "
                        "сосуд не построен (vessel=none?)");
        }
        fSampleSolid = lv->GetSolid();
        fSampleSolid->BoundingLimits(fSampleMin, fSampleMax);
    }
    // Отбор с отклонением: тело пробы — G4Polycone (кольцо вокруг колодца плюс
    // слой над ним), аналитического равномерного розыгрыша для него нет, а
    // Inside() у солида точен. Счётчик срыва защищает от бесконечного цикла.
    for (int i = 0; i < 10000; ++i) {
        const G4ThreeVector p(
            fSampleMin.x() + (fSampleMax.x() - fSampleMin.x()) * G4UniformRand(),
            fSampleMin.y() + (fSampleMax.y() - fSampleMin.y()) * G4UniformRand(),
            fSampleMin.z() + (fSampleMax.z() - fSampleMin.z()) * G4UniformRand());
        if (fSampleSolid->Inside(p) != kInside) continue;
        // Ограничение области: снизу по высоте и снаружи по радиусу.
        if (gSrcZFrac < 1.0) {
            const double zTop = fSampleMin.z()
                + (fSampleMax.z() - fSampleMin.z()) * gSrcZFrac;
            if (p.z() > zTop) continue;
        }
        if (gSrcRFrac < 1.0) {
            const double rMax = std::max(std::abs(fSampleMax.x()),
                                         std::abs(fSampleMin.x()));
            const double rMin = rMax * (1.0 - gSrcRFrac);
            if (p.perp() < rMin) continue;
        }
        return p;
    }
    G4Exception("G1sNpsmPrimaryGeneratorAction::SamplePointInSample",
                "RejectionFailed", FatalException,
                "10000 попыток отбора не дали точки внутри пробы");
    return G4ThreeVector();
}

G1sNpsmPrimaryGeneratorAction::G1sNpsmPrimaryGeneratorAction()
    : fGun(1) {
    // В режиме "ion" частицу здесь НЕ задаём: таблица ионов доступна только
    // после инициализации ядра, и обращение к ней из конструктора действия
    // вернёт нуль. Ион создаётся при первом событии, см. GeneratePrimaries.
    if (gPrimary != "ion") {
        G4ParticleDefinition* gamma = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
        fGun.SetParticleDefinition(gamma);
    }
}

void G1sNpsmPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    // Позиция источника: точка на оси либо равномерно по объёму пробы.
    const G4ThreeVector position = (gSourceMode == "sample")
        ? SamplePointInSample()
        : G4ThreeVector(0., 0., gSourceZmm * mm);
    fGun.SetParticlePosition(position);

    if (gPrimary == "eplus" || gPrimary == "eplus_gamma") {
        // β⁺ (Sc-44): позитрон покоится в точке рождения и аннигилирует сам —
        // две встречные 511 кэВ; при eplus_gamma из той же точки в том же
        // событии — изотропный квант gEnergyKeV. Пробег позитрона и аннигиляция
        // на лету не учитываются — приближение, названное в отчёте.
        static G4ParticleDefinition* ep = G4ParticleTable::GetParticleTable()->FindParticle("e+");
        static G4ParticleDefinition* gm = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
        fGun.SetParticleDefinition(ep);
        fGun.SetParticleEnergy(0.0);
        fGun.SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
        fGun.GeneratePrimaryVertex(anEvent);
        if (gPrimary == "eplus") return;
        fGun.SetParticleDefinition(gm);
        fGun.SetParticleEnergy(gEnergyKeV * keV);
        const G4double ct = 1.0 - 2.0 * G4UniformRand(), ph = 2.0 * M_PI * G4UniformRand();
        const G4double st = std::sqrt(1.0 - ct * ct);
        fGun.SetParticleMomentumDirection(G4ThreeVector(st * std::cos(ph), st * std::sin(ph), ct));
        fGun.GeneratePrimaryVertex(anEvent);
        return;
    }

    if (gPrimary == "ion") {
        // Ядро создаётся один раз, при первом событии: таблица ионов готова
        // только после Initialize. Энергия нулевая — ядро покоится, распад
        // разыгрывает G4RadioactiveDecay, он же испускает весь каскад.
        if (!fIon) {
            fIon = G4IonTable::GetIonTable()->GetIon(gIonZ, gIonA, 0.0);
            if (!fIon) {
                std::fprintf(stderr, "Не найден ион Z=%d A=%d\n", gIonZ, gIonA);
                std::abort();
            }
            fGun.SetParticleDefinition(fIon);
            fGun.SetParticleCharge(0.0);
        }
        fGun.SetParticleEnergy(0.0);
    } else if (gPrimary == "gamma_table") {
        // Фотон с энергией из таблицы; частица gamma задана в конструкторе,
        // позиция и изотропное направление — общие с режимом gamma (ниже).
        fGun.SetParticleEnergy(SampleTableEnergyKeV() * keV);
    } else {
        // Задаем энергию
        fGun.SetParticleEnergy(gEnergyKeV * keV);
    }

    // Генерируем изотропное направление
    // Разыгрываем косинус полярного угла для равномерного распределения по телесному углу
    G4double cosTheta = 1.0 - 2.0 * G4UniformRand();
    G4double phi = 2.0 * M_PI * G4UniformRand();
    G4double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);

    G4ThreeVector direction(sinTheta * std::cos(phi),
                            sinTheta * std::sin(phi),
                            cosTheta);
    fGun.SetParticleMomentumDirection(direction);

    // Генерируем первичную вершину
    fGun.GeneratePrimaryVertex(anEvent);
}
