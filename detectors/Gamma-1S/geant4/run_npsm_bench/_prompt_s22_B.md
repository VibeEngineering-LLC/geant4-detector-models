# Задание B: накопление второй величины и два спектра на выходе

Ты — генератор кода C++ для проекта Geant4 11.2.1 (Windows, MSVC, C++17).

## Что вернуть

**ЧЕТЫРЕ ПОЛНЫХ ФАЙЛА** — целиком, от первой строки до последней, готовые к
компиляции. НЕ фрагменты, НЕ диффы, НЕ «...остальное без изменений».

Формат вывода — строго такой, без единого слова вне блоков:

```
=== FILE: NpsmBenchEventAction.hh ===
<полный текст>
=== FILE: NpsmBenchEventAction.cc ===
<полный текст>
=== FILE: NpsmBenchRunAction.hh ===
<полный текст>
=== FILE: NpsmBenchRunAction.cc ===
<полный текст>
```

Текущий текст этих файлов приведён в конце задания под заголовком «БАЗА».
Бери его за основу: всё, что не названо изменением, обязано сохраниться
дословно, включая комментарии на русском языке. Комментарии пиши по-русски.

---

## Файл 1-2: NpsmBenchEventAction

### Изменение B1. Вторая накапливаемая величина

Поле `double fEdepLightMeV;` с ТЕМ ЖЕ жизненным циклом, что `fEdepMeV`:
инициализация нулём в конструкторе, обнуление в `BeginOfEventAction`,
передача в `RecordEvent` дополнительным (последним) аргументом.

Мутатор: `inline void AddEdepLight(double v) { fEdepLightMeV += v; }`.

Всё остальное — дословно как есть, включая громкий `abort()` в конструкторе
и комментарий о том, почему проверки на `nullptr` нет в `EndOfEventAction`.

---

## Файл 3-4: NpsmBenchRunAction

### Изменение B2. Зависимость от модели

Конструктор получает ПЯТЫЙ аргумент — `const NpsmLightYield* lightYield`,
сохраняется в поле, проверяется на `nullptr` с печатью в `stderr` и
`std::abort()`. Нужен для шапки CSV: постановка обязана лежать В ФАЙЛЕ, а не
в его имени.

### Изменение B3. Сигнатура RecordEvent

```cpp
void RecordEvent(int nCompt, int nRayl, bool phot, bool conv, bool escaped,
                 double edepMeV, double edepLightMeV);
```

Вся существующая логика (взаимоисключающие исходы первичного фотона,
гистограммы рассеяний, проверка полного поглощения) сохраняется дословно.

### Изменение B4. Два спектра

```cpp
static constexpr int kNBins = 3200;
static constexpr double kBinKeV = 1.0;
std::array<long long, kNBins> fSpecEdep{};
std::array<long long, kNBins> fSpecLight{};
long long fNOverflowEdep = 0;
long long fNOverflowLight = 0;
```

В `RecordEvent` каждая из двух величин переводится в кэВ
(`value_keV = value_MeV * 1000.0`) и попадает в бин
`static_cast<int>(value_keV / kBinKeV)`. Если бин вне `[0, kNBins)` —
величина в спектр НЕ идёт, инкрементируется соответствующий счётчик
переполнения. Две величины обрабатываются независимо друг от друга.

### Изменение B5. Суммы для кривой отклика

Это ГЛАВНАЯ выходная величина этапа — по ней строится кривая
«свет / энергия», сравниваемая с независимым расчётом. Поля:

```cpp
double fSumEdepMeV = 0.0;
double fSumLightMeV = 0.0;
double fSumLight2MeV2 = 0.0;   // для оценки погрешности среднего
long long fNWithEdep = 0;      // события с edepMeV > 0
```

Накапливаются в `RecordEvent` ТОЛЬКО для событий с `edepMeV > 0`.
Порядок сложения для двух величин обязан быть одинаковым: при выключенной
модели вес ровно 1.0, и суммы обязаны совпасть ПОБИТОВО — это контрольная
проверка приёмки, а не украшение.

### Изменение B6. Шапка CSV

После существующих строк шапки, ДО строки `n_compt,count_absorbed,count_all`,
добавляются (значения — из `fLightYield` и из статических полей генератора
`NpsmBenchPrimaryGeneratorAction::gParticle` и `::gPencilBeam`):

```
particle,<gamma|e->
beam,<iso|pencil>
npsm_enabled,<0|1>
npsm_eta,<...>
npsm_s_ons,<...>
npsm_s_trap,<...>
npsm_s_birks,<...>
npsm_bad_s_count,<...>
npsm_out_of_range_count,<...>
n_with_edep,<...>
sum_edep_MeV,<...>
sum_light_MeV,<...>
sum_light2_MeV2,<...>
n_overflow_edep,<...>
n_overflow_light,<...>
```

Числовые значения сумм печатать с полной точностью:
`file << std::setprecision(17)` перед ними (побитовое сравнение двух сумм
невозможно при печати с точностью по умолчанию — на этом проекте уже был
ложный вывод о расхождении из-за недостаточной точности печати).

### Изменение B7. Вторая таблица

После существующей таблицы гистограммы рассеяний (она заканчивается циклом
по `kMaxCompt`) добавляется вторая таблица:

```
bin_keV,count_edep,count_light
```

Строка печатается только если ХОТЯ БЫ ОДИН из двух счётчиков бина ненулевой
(иначе файл распухнет на 3200 строк при десятке заполненных). В первой
колонке — центр бина: `(i + 0.5) * kBinKeV`.

### Изменение B8. Строка в stdout

К существующей печати в `EndOfRunAction` добавить вторую строку:

```
npsm_bench: n_with_edep=<...> mean_edep_keV=<...> mean_light_keV=<...> ratio=<...>
```

где `ratio = fSumLightMeV / fSumEdepMeV` при ненулевом знаменателе, иначе 0.

---

## Ограничения

- C++17, Geant4 11.2.1. Заголовки — только реально существующие в 11.2.1.
- `NpsmLightYield` — свой класс проекта, заголовок `NpsmLightYield.hh`,
  от Geant4 не зависит. Его интерфейс приведён в БАЗЕ.
- Не добавлять ничего, о чём не сказано в задании.

---

## БАЗА — текущий текст файлов

### NpsmBenchEventAction.hh

```cpp
#pragma once

#include "G4UserEventAction.hh"
#include "G4Event.hh"

class NpsmBenchRunAction;

class NpsmBenchEventAction : public G4UserEventAction {
public:
    // Указатель на RunAction ОБЯЗАТЕЛЕН. В первой версии конструктор был без
    // аргументов, fRunAction оставался nullptr, а EndOfEventAction молча
    // пропускал каждое событие через `if (fRunAction)` — прогон 10 000 событий
    // отработал успешно и записал CSV со всеми нулями (05.09).
    explicit NpsmBenchEventAction(NpsmBenchRunAction* runAction);
    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

    // Мутаторы
    inline void AddCompt() { fNCompt++; }
    inline void AddRayl() { fNRayl++; }
    inline void SetPhotAbsorbed() { fPhotAbsorbed = true; }
    inline void SetConv() { fConv = true; }
    inline void SetEscaped() { fEscaped = true; }
    inline void AddEdep(double edepMeV) { fEdepMeV += edepMeV; }

private:
    int fNCompt;
    int fNRayl;
    bool fPhotAbsorbed;
    bool fConv;
    bool fEscaped;
    double fEdepMeV;

    NpsmBenchRunAction* fRunAction;
};
```

### NpsmBenchEventAction.cc

```cpp
#include "NpsmBenchEventAction.hh"
#include "NpsmBenchRunAction.hh"

#include <cstdio>
#include <cstdlib>

NpsmBenchEventAction::NpsmBenchEventAction(NpsmBenchRunAction* runAction)
    : fNCompt(0), fNRayl(0), fPhotAbsorbed(false), fConv(false), fEscaped(false), fEdepMeV(0.0),
      fRunAction(runAction) {
    // Отказ громкий, а не тихий: без RunAction прогон не имеет смысла, и лучше
    // не стартовать вовсе, чем отработать 10 000 событий с нулевой статистикой.
    if (!fRunAction) {
        std::fprintf(stderr, "NpsmBenchEventAction: FATAL runAction == nullptr\n");
        std::abort();
    }
}

void NpsmBenchEventAction::BeginOfEventAction(const G4Event*) {
    fNCompt = 0;
    fNRayl = 0;
    fPhotAbsorbed = false;
    fConv = false;
    fEscaped = false;
    fEdepMeV = 0.0;
}

void NpsmBenchEventAction::EndOfEventAction(const G4Event*) {
    // Проверки на nullptr здесь намеренно нет: он исключён в конструкторе, а
    // тихая проверка в этом месте и была причиной нулевой статистики.
    fRunAction->RecordEvent(fNCompt, fNRayl, fPhotAbsorbed, fConv, fEscaped, fEdepMeV);
}
```

### NpsmBenchRunAction.hh

```cpp
#pragma once

#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include <array>
#include <string>

class NpsmBenchRunAction : public G4UserRunAction {
public:
    explicit NpsmBenchRunAction(std::string outCsv, double energyKeV, long long nEventsRequested, long seed);
    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

    void RecordEvent(int nCompt, int nRayl, bool phot, bool conv, bool escaped, double edepMeV);

private:
    std::string fOutCsv;
    double fEnergyKeV;
    long long fNEventsRequested;
    long fSeed;

    static constexpr int kMaxCompt = 32;
    std::array<long long, kMaxCompt+1> fHistAbsorbed{};
    std::array<long long, kMaxCompt+1> fHistAll{};
    long long fNEvents;
    long long fNAbsorbed;
    long long fNConv;
    long long fNEscaped;
    long long fNOther;
    long long fNFullEdep;
    double fSumComptAbsorbed;
    double fSumCompt2Absorbed;

    void WriteCSV();
};
```

### NpsmBenchRunAction.cc

```cpp
#include "NpsmBenchRunAction.hh"
// Шапка CSV читает статические поля геометрии и физ-листа напрямую — без этих
// заголовков компилятор видел «<error type>» на SourceRadiusMm() (сборка 05.09).
#include "NpsmBenchDetectorConstruction.hh"
#include "Rc103FieldPhysicsList.hh"
#include "G4EmParameters.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>
#include <iomanip>
#include <cmath>

NpsmBenchRunAction::NpsmBenchRunAction(std::string outCsv, double energyKeV, long long nEventsRequested, long seed)
    : fOutCsv(outCsv), fEnergyKeV(energyKeV), fNEventsRequested(nEventsRequested), fSeed(seed),
      fNEvents(0), fNAbsorbed(0), fNConv(0), fNEscaped(0), fNOther(0), fNFullEdep(0),
      fSumComptAbsorbed(0.0), fSumCompt2Absorbed(0.0) {}

void NpsmBenchRunAction::BeginOfRunAction(const G4Run*) {
    // Ничего не делаем
}

void NpsmBenchRunAction::RecordEvent(int nCompt, int nRayl, bool phot, bool conv, bool escaped, double edepMeV) {
    ++fNEvents;

    // Исходы первичного фотона — ВЗАИМОИСКЛЮЧАЮЩИЕ категории. Сумма
    // n_absorbed_phot + n_conv + n_escaped + n_other обязана равняться
    // n_events_processed — это проверяет анализатор. В сгенерированной версии
    // «other» инкрементировался и для поглощённых событий (ветка else после
    // проверки conv/escaped внутри блока phot) — счётчик дублировал n_absorbed.
    if (phot) {
        ++fNAbsorbed;
        fSumComptAbsorbed += nCompt;
        fSumCompt2Absorbed += static_cast<double>(nCompt) * nCompt;
        const int bin = (nCompt > kMaxCompt) ? kMaxCompt : nCompt;
        ++fHistAbsorbed[static_cast<std::size_t>(bin)];
    } else if (conv) {
        ++fNConv;
    } else if (escaped) {
        ++fNEscaped;
    } else {
        ++fNOther;  // ни одного из трёх флагов: иной конец трека
    }

    // Все события
    int bin = (nCompt > kMaxCompt) ? kMaxCompt : nCompt;
    ++fHistAll[bin];

    // Проверка полного поглощения по энергии
    double diff = std::abs(edepMeV - fEnergyKeV / 1000.0);
    if (diff < 0.001) { // 1 keV
        ++fNFullEdep;
    }
}

void NpsmBenchRunAction::EndOfRunAction(const G4Run*) {
    WriteCSV();

    double mean = fNAbsorbed > 0 ? fSumComptAbsorbed / fNAbsorbed : 0.0;
    double sem = 0.0;
    if (fNAbsorbed > 0) {
        double mean2 = fSumCompt2Absorbed / fNAbsorbed;
        double variance = mean2 - mean * mean;
        if (variance >= 0 && fNAbsorbed > 1) {
            sem = std::sqrt(variance / (fNAbsorbed - 1));
        }
    }

    // Вывод в stdout
    std::printf("npsm_bench: E=%.1f keV  N_absorbed=%lld  mean_ncompt=%.4f +- %.4f  escaped=%lld\n",
                fEnergyKeV, fNAbsorbed, mean, sem, fNEscaped);
}

void NpsmBenchRunAction::WriteCSV() {
    std::ofstream file(fOutCsv, std::ios::out);
    if (!file.is_open()) {
        return;
    }

    file << "# npsm_bench geometry/transport benchmark\n";
    file << "energy_keV," << fEnergyKeV << "\n";
    file << "n_events_requested," << fNEventsRequested << "\n";
    file << "n_events_processed," << fNEvents << "\n";
    file << "seed," << fSeed << "\n";
    file << "em_cut_mm," << Rc103FieldPhysicsList::gCutMm << "\n";
    file << "em_deex," << Rc103FieldPhysicsList::gDeexMode << "\n";
    file << "lowest_electron_energy_keV," << G4EmParameters::Instance()->LowestElectronEnergy() / keV << "\n";
    file << "crystal_mm," << NpsmBenchDetectorConstruction::gCrystalXMm << "x" << NpsmBenchDetectorConstruction::gCrystalYMm << "x" << NpsmBenchDetectorConstruction::gCrystalZMm << "\n";
    file << "crystal_volume_cm3," << (NpsmBenchDetectorConstruction::gCrystalXMm * NpsmBenchDetectorConstruction::gCrystalYMm * NpsmBenchDetectorConstruction::gCrystalZMm) / 1000.0 << "\n";
    file << "r_src_mm," << NpsmBenchDetectorConstruction::SourceRadiusMm() << "\n";
    file << "n_absorbed_phot," << fNAbsorbed << "\n";
    file << "n_full_edep_1keV," << fNFullEdep << "\n";
    file << "n_conv," << fNConv << "\n";
    file << "n_escaped," << fNEscaped << "\n";
    file << "n_other," << fNOther << "\n";
    file << "mean_ncompt_absorbed," << (fNAbsorbed > 0 ? fSumComptAbsorbed / fNAbsorbed : 0.0) << "\n";
    file << "sem_ncompt_absorbed," << (fNAbsorbed > 0 ? std::sqrt((fSumCompt2Absorbed / fNAbsorbed - (fSumComptAbsorbed / fNAbsorbed) * (fSumComptAbsorbed / fNAbsorbed)) / fNAbsorbed) : 0.0) << "\n";
    file << "n_compt,count_absorbed,count_all\n";

    for (int i = 0; i <= kMaxCompt; ++i) {
        file << i << "," << fHistAbsorbed[i] << "," << fHistAll[i] << "\n";
    }

    file.close();
}
```

### NpsmLightYield.hh

```cpp
#pragma once

class NpsmLightYield {
 public:
  NpsmLightYield();  // параметры по умолчанию — значения для NaI(Tl)
  // Вес светового выхода. S — тормозная способность, МэВ/см.
  double Weight(double S_MeV_cm) const;
  // Настройка параметров; каждый обязан быть строго положительным,
  // eta — в интервале (0, 1]. Нарушение — громкий отказ, см. ниже.
  void SetParameters(double eta, double sOns, double sTrap, double sBirks);
  // Включена ли модель. Когда выключена, Weight ВСЕГДА возвращает 1.0.
  void SetEnabled(bool on) { fEnabled = on; }
  bool IsEnabled() const { return fEnabled; }
  // Для шапки CSV — постановка обязана лежать В ФАЙЛЕ, а не в его имени.
  double Eta() const; double SOns() const; double STrap() const; double SBirks() const;
  long long BadSCount() const { return fNBadS; }
  long long OutOfRangeCount() const { return fNOutOfRange; }

 private:
  bool fEnabled = false;   // ВЫКЛЮЧЕНА по умолчанию
  double fEta, fSOns, fSTrap, fSBirks;
  mutable long long fNBadS = 0;
  mutable long long fNOutOfRange = 0;
};

// Самопроверка модели: 5 утверждений, 0 при успехе. Определена в .cc;
// объявление здесь, чтобы main.cc не заводил свой extern-прототип.
int NpsmLightYieldSelfTest();
```

### NpsmBenchPrimaryGeneratorAction.hh

```cpp
#pragma once

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include <string>

class NpsmBenchPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    // Постановка задаётся ДО построения объекта (main разбирает ключи раньше),
    // поэтому поля статические: "gamma" либо "e-"; iso-поток либо пучок из центра.
    static std::string gParticle;
    static bool gPencilBeam;

    explicit NpsmBenchPrimaryGeneratorAction(double energyKeV);
    void GeneratePrimaries(G4Event* anEvent) override;

private:
    G4ParticleGun fGun;
    double fEnergyKeV;
};
```
