# Задание C: ключи командной строки и порядок сборки объектов в main.cc

Ты — генератор кода C++ для проекта Geant4 11.2.1 (Windows, MSVC, C++17).

## Что вернуть

**ОДИН ПОЛНЫЙ ФАЙЛ** — целиком, от первой строки до последней, готовый к
компиляции. НЕ фрагмент, НЕ дифф, НЕ «...остальное без изменений».

Формат вывода — строго такой, без единого слова вне блока:

```
=== FILE: main.cc ===
<полный текст>
```

Текущий текст файла приведён в конце задания под заголовком «БАЗА». Бери его
за основу: всё, что не названо изменением, обязано сохраниться дословно,
включая комментарии на русском языке. Комментарии пиши по-русски.

---

## Изменение C1. Режим самопроверки — ПЕРВЫМ ДЕЛОМ

В самом начале `main`, ДО проверки `argc < 4`:

```cpp
if (argc >= 2 && std::string(argv[1]) == "selftest") {
    return NpsmLightYieldSelfTest();
}
```

Функция объявлена в `NpsmLightYield.hh`, возвращает 0 при успехе. Порядок
важен: самопроверка не должна требовать ни энергии, ни числа событий, ни
имени выходного файла — она не запускает перенос вовсе.

## Изменение C2. Новые ключи

К существующим `seed=`, `emcut=`, `deex=`, `crystal=` добавляются:

| Ключ | Значения | Умолчание |
|---|---|---|
| `particle=` | `gamma` либо `e-` | `gamma` |
| `beam=` | `iso` либо `pencil` | `iso` |
| `npsm=` | `on` либо `off` | `off` |
| `npsm_eta=` | число из (0, 1] | из класса |
| `npsm_ons=` | число > 0 | из класса |
| `npsm_trap=` | число > 0 | из класса |
| `npsm_birks=` | число > 0 | из класса |

Неизвестный ключ или недопустимое значение — сообщение в `stderr` и
возврат 2, как это уже сделано для существующих ключей. Строку `Usage:`
дополнить новыми ключами.

Умолчания четырёх числовых параметров НЕ зашивать числами в `main.cc`:
создать объект `NpsmLightYield`, прочитать его умолчания методами `Eta()`,
`SOns()`, `STrap()`, `SBirks()` в локальные переменные, поверх наложить
значения ключей, если они заданы, и одним вызовом `SetParameters` записать
обратно. Иначе умолчания разъедутся между двумя файлами.

## Изменение C3. Порядок построения объектов

СТРОГО такой (нарушение порядка даст либо `nullptr` в чужом конструкторе,
либо шапку CSV с умолчаниями вместо заданных значений):

```
1. Разбор всех ключей.
2. Статические поля геометрии — как сейчас.
3. NpsmBenchPrimaryGeneratorAction::gParticle   = <particle>;
   NpsmBenchPrimaryGeneratorAction::gPencilBeam = (<beam> == "pencil");
4. NpsmLightYield lightYield;              // на стеке main, живёт до конца
   lightYield.SetParameters(eta, ons, trap, birks);
   lightYield.SetEnabled(npsm == "on");
5. RunManager, DetectorConstruction, PhysicsList — как сейчас.
6. runAction  = new NpsmBenchRunAction(outCsv, energyKeV, nEvents, seed, &lightYield);
7. eventAction = new NpsmBenchEventAction(runAction);
8. primaryGenerator = new NpsmBenchPrimaryGeneratorAction(energyKeV);
9. steppingAction = new NpsmBenchSteppingAction(eventAction, &lightYield);
10. runManager->Initialize();  далее как сейчас.
```

## Изменение C4. Печать постановки перед прогоном

К существующей строке `npsm_bench: energy_keV=... ` добавить вторую строку
ДО `BeamOn`:

```
npsm_bench: particle=<...> beam=<...> npsm=<on|off> eta=<...> s_ons=<...> s_trap=<...> s_birks=<...>
```

Это проверяемый артефакт: по логу прогона обязано быть видно, с какими
параметрами он шёл, не заглядывая в CSV.

## Изменение C5. Заголовки

Добавить `#include "NpsmLightYield.hh"` и `#include <cstring>` (последний
нужен для `std::strncmp`, который уже используется, но подключался
транзитивно — на другом компиляторе это сломается).

---

## Ограничения

- C++17, Geant4 11.2.1.
- Не добавлять ничего, о чём не сказано в задании.
- Сравнение строк ключей — как в существующем коде, через `std::strncmp`.

---

## БАЗА — текущий текст файла

### main.cc

```cpp
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iomanip>

#include "G4RunManagerFactory.hh"
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

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "Usage: %s <energy_keV> <n_events> <out_csv> [seed=<N>] [emcut=<mm>] [deex=std|deex|max] [crystal=<X>x<Y>x<Z>]\n", argv[0]);
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
            if (deexMode != "std" && deexMode != "deex" && deexMode != "max") {
                std::fprintf(stderr, "Ошибка: неизвестный режим deex\n");
                return 2;
            }
        } else if (std::strncmp(arg, "crystal=", 8) == 0) {
            crystalStr = arg + 8;
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

    // Инициализация Геанта
    auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
    runManager->SetVerboseLevel(0);

    auto detectorConstruction = new NpsmBenchDetectorConstruction();
    runManager->SetUserInitialization(detectorConstruction);

    auto physicsList = new Rc103FieldPhysicsList(emCutMm, deexMode);
    runManager->SetUserInitialization(physicsList);

    auto runAction = new NpsmBenchRunAction(outCsv, energyKeV, nEvents, seed);
    runManager->SetUserAction(runAction);

    auto eventAction = new NpsmBenchEventAction(runAction);
    runManager->SetUserAction(eventAction);

    auto primaryGenerator = new NpsmBenchPrimaryGeneratorAction(energyKeV);
    runManager->SetUserAction(primaryGenerator);

    auto steppingAction = new NpsmBenchSteppingAction(eventAction);
    runManager->SetUserAction(steppingAction);

    runManager->Initialize();

    if (!NpsmBenchDetectorConstruction::GetCrystalLogicalVolume()) {
        std::fprintf(stderr, "Ошибка: кристалл не был создан\n");
        return 3;
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

    // Установка зерна случайных чисел
    if (seed != 0) {
        G4Random::setTheSeed(seed);
    }

    runManager->BeamOn(nEvents);

    delete runManager;

    std::printf("EXITCODE=0\n");
    return 0;
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

### NpsmBenchRunAction.hh

```cpp
#pragma once

#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include <array>
#include <string>

class NpsmLightYield;

class NpsmBenchRunAction : public G4UserRunAction {
public:
    explicit NpsmBenchRunAction(std::string outCsv, double energyKeV, long long nEventsRequested, long seed, const NpsmLightYield* lightYield);
    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

    void RecordEvent(int nCompt, int nRayl, bool phot, bool conv, bool escaped, double edepMeV, double edepLightMeV);

private:
    std::string fOutCsv;
    double fEnergyKeV;
    long long fNEventsRequested;
    long fSeed;
    const NpsmLightYield* fLightYield;

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

    static constexpr int kNBins = 3200;
    static constexpr double kBinKeV = 1.0;
    std::array<long long, kNBins> fSpecEdep{};
    std::array<long long, kNBins> fSpecLight{};
    long long fNOverflowEdep = 0;
    long long fNOverflowLight = 0;

    double fSumEdepMeV = 0.0;
    double fSumLightMeV = 0.0;
    double fSumLight2MeV2 = 0.0;   // для оценки погрешности среднего
    long long fNWithEdep = 0;      // события с edepMeV > 0

    void WriteCSV();
};
```

### NpsmBenchSteppingAction.hh

```cpp
#pragma once

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "G4EmCalculator.hh"

class NpsmBenchEventAction;
class NpsmLightYield;

class NpsmBenchSteppingAction : public G4UserSteppingAction {
public:
    NpsmBenchSteppingAction(NpsmBenchEventAction* eventAction,
                            const NpsmLightYield* lightYield);
    void UserSteppingAction(const G4Step*) override;

private:
    NpsmBenchEventAction* fEventAction;
    const NpsmLightYield* fLightYield;
    G4EmCalculator fEmCalc;
    static constexpr double kMinStepMm = 1e-4;
};
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
