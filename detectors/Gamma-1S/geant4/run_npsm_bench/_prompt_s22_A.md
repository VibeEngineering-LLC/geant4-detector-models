# Задание A: генератор первичных частиц и взвешивание в SteppingAction

Ты — генератор кода C++ для проекта Geant4 11.2.1 (Windows, MSVC, C++17).

## Что вернуть

**ЧЕТЫРЕ ПОЛНЫХ ФАЙЛА** — целиком, от первой строки до последней, готовые к
компиляции. НЕ фрагменты, НЕ диффы, НЕ «...остальное без изменений».

Формат вывода — строго такой, без единого слова вне блоков:

```
=== FILE: NpsmBenchPrimaryGeneratorAction.hh ===
<полный текст>
=== FILE: NpsmBenchPrimaryGeneratorAction.cc ===
<полный текст>
=== FILE: NpsmBenchSteppingAction.hh ===
<полный текст>
=== FILE: NpsmBenchSteppingAction.cc ===
<полный текст>
```

Текущий текст этих файлов приведён в конце задания под заголовком
«БАЗА». Бери его за основу: всё, что не названо изменением, обязано
сохраниться дословно, включая комментарии на русском языке. Комментарии
пиши по-русски.

---

## Файл 1-2: NpsmBenchPrimaryGeneratorAction

### Изменение A1. Выбор частицы

Два новых статических поля класса (объявить в `.hh`, определить в `.cc`):

```cpp
static std::string gParticle;   // "gamma" по умолчанию, либо "e-"
static bool gPencilBeam;        // false по умолчанию
```

В конструкторе частица ищется по имени `gParticle`:
`G4ParticleTable::GetParticleTable()->FindParticle(gParticle)`. Если
вернулся `nullptr` — печать в `stderr` и `std::abort()`. Отказ громкий, а не
тихий: молча подставленная гамма испортила бы весь прогон (инцидент W-053 в
этом же проекте — тихая защита дала прогон с нулевой статистикой).

### Изменение A2. Режим карандашного пучка

В `GeneratePrimaries`: если `gPencilBeam == true`, то позиция `(0,0,0)`
(центр кристалла), направление `(0,0,1)`. Никакого розыгрыша по сфере в этом
режиме не делается вовсе.

Иначе — ровно существующий изотропный розыгрыш со сферы, без изменений.

Комментарием обосновать: для проверки кривой электронного отклика нужен
электрон, тормозящийся ЦЕЛИКОМ внутри кристалла; пробег электрона 3 МэВ в
NaI около 4 мм при полутолщине кристалла 51 мм — запас более чем
десятикратный.

---

## Файл 3-4: NpsmBenchSteppingAction

### Изменение A3. Зависимости через конструктор

Конструктор принимает ДВА указателя:

```cpp
NpsmBenchSteppingAction(NpsmBenchEventAction* eventAction,
                        const NpsmLightYield* lightYield);
```

ОБА проверяются на `nullptr` в конструкторе: печать в `stderr` и
`std::abort()`. Проверок на `nullptr` в `UserSteppingAction` НЕ ставить —
именно тихая проверка в горячем методе была причиной инцидента W-053.

Член класса `G4EmCalculator fEmCalc;` — создаётся ОДИН раз как поле, а не на
каждом шаге (создание на шаге дорого). Заголовок `G4EmCalculator.hh`.

### Изменение A4. Расчёт веса светового выхода

Внутри уже существующей ветви накопления энергии
(`edep > 0 && preVol->GetLogicalVolume() == crystalLV`) вместо одной строки
`fEventAction->AddEdep(edep / MeV)` делается следующее, СТРОГО в этом порядке:

```
1. edepMeV = edep / MeV
2. w = 1.0
3. Если fLightYield->IsEnabled() И частица шага есть e- либо e+:
      w = fLightYield->Weight(S), где S вычисляется по п.4
   Иначе w остаётся ровно 1.0, и S НЕ вычисляется вовсе
   (иначе счётчики модели в контрольном прогоне npsm=off перестали бы
    быть нулевыми, а этот ноль — критерий приёмки).
4. Вычисление S (удельные потери, МэВ/см), две ветви:
   ветвь «шаг конечной длины»: если stepLength > kMinStepMm * mm
        S = edepMeV / (stepLength / cm)
   ветвь «локальный депозит»: иначе
        dedx = fEmCalc.ComputeElectronicDEDX(Ek, particleDef, material, DBL_MAX)
        S    = dedx / (MeV / cm)
   где Ek = step->GetPreStepPoint()->GetKineticEnergy(),
       material = step->GetPreStepPoint()->GetMaterial(),
       particleDef = step->GetTrack()->GetDefinition().
5. fEventAction->AddEdep(edepMeV);
   fEventAction->AddEdepLight(edepMeV * w);
```

Константа: `static constexpr double kMinStepMm = 1e-4;` — порог, ниже
которого шаг считается локальным депозитом.

Тип частицы сравнивать с `G4Electron::Definition()` и
`G4Positron::Definition()` (заголовки `G4Electron.hh`, `G4Positron.hh`),
НЕ по строковому имени.

### Почему ComptonElectronicDEDX, а не GetDEDX — воспроизвести это в комментарии

`G4EmCalculator::GetDEDX(...)` возвращает ОГРАНИЧЕННУЮ тормозную способность:
потери выше порога продукции отнесены к дельта-электронам и в неё не входят,
поэтому значение зависит от `emcut`. Модель непропорциональности описывает
плотность возбуждений вдоль трека — величину физическую, от вычислительного
порога не зависящую. Поэтому берётся ПОЛНАЯ электронная тормозная
способность: `ComputeElectronicDEDX(Ek, particle, material, DBL_MAX)`.

### Изменение A5. Остальное без изменений

Блок счёта процессов первичного фотона (compt / Rayl / phot / conv /
escaped) сохраняется дословно как есть.

---

## Ограничения

- C++17, Geant4 11.2.1. Заголовки — только реально существующие в 11.2.1.
- `DBL_MAX` требует `<cfloat>`.
- Единицы Geant4 внутренние: `MeV = 1`, `mm = 1`, `cm = 10`. Деление на
  `MeV`/`cm` обязательно там, где величина уходит наружу в физических
  единицах.
- Никаких `std::cout` в горячем методе.
- Не добавлять ничего, о чём не сказано в задании.

---

## БАЗА — текущий текст файлов

### NpsmBenchPrimaryGeneratorAction.hh

```cpp
#pragma once

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

class NpsmBenchPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    explicit NpsmBenchPrimaryGeneratorAction(double energyKeV);
    void GeneratePrimaries(G4Event* anEvent) override;

private:
    G4ParticleGun fGun;
    double fEnergyKeV;
};
```

### NpsmBenchPrimaryGeneratorAction.cc

```cpp
#include "NpsmBenchPrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4Event.hh"

#include "NpsmBenchDetectorConstruction.hh"

NpsmBenchPrimaryGeneratorAction::NpsmBenchPrimaryGeneratorAction(double energyKeV)
    : fGun(1), fEnergyKeV(energyKeV) {
    G4ParticleDefinition* gamma = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
    fGun.SetParticleDefinition(gamma);
}

void NpsmBenchPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    const double R_SRC = NpsmBenchDetectorConstruction::SourceRadiusMm() * mm;

    // 1) точка равномерно по сфере радиуса R_SRC
    const double cosT = 2.0 * G4UniformRand() - 1.0;
    const double sinT = std::sqrt(1.0 - cosT * cosT);
    const double phi  = twopi * G4UniformRand();
    const G4ThreeVector n(sinT * std::cos(phi), sinT * std::sin(phi), cosT); // наружная нормаль
    const G4ThreeVector pos = R_SRC * n;

    // 2) направление ВНУТРЬ, косинусное относительно -n
    const G4ThreeVector e3 = -n;
    G4ThreeVector e1 = e3.orthogonal().unit();
    const G4ThreeVector e2 = e3.cross(e1).unit();
    const double ct  = std::sqrt(G4UniformRand());
    const double st  = std::sqrt(1.0 - ct * ct);
    const double psi = twopi * G4UniformRand();
    const G4ThreeVector dir = st * std::cos(psi) * e1 + st * std::sin(psi) * e2 + ct * e3;

    fGun.SetParticlePosition(pos);
    fGun.SetParticleMomentumDirection(dir);
    fGun.SetParticleEnergy(fEnergyKeV * keV);

    fGun.GeneratePrimaryVertex(anEvent);
}
```

### NpsmBenchSteppingAction.hh

```cpp
#pragma once

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"

class NpsmBenchEventAction;

class NpsmBenchSteppingAction : public G4UserSteppingAction {
public:
    explicit NpsmBenchSteppingAction(NpsmBenchEventAction* eventAction);
    void UserSteppingAction(const G4Step*) override;

private:
    NpsmBenchEventAction* fEventAction;
};
```

### NpsmBenchSteppingAction.cc

```cpp
#include "NpsmBenchSteppingAction.hh"
#include "NpsmBenchEventAction.hh"
#include "NpsmBenchDetectorConstruction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHandle.hh"
#include "G4SystemOfUnits.hh"

NpsmBenchSteppingAction::NpsmBenchSteppingAction(NpsmBenchEventAction* eventAction)
    : fEventAction(eventAction) {}

void NpsmBenchSteppingAction::UserSteppingAction(const G4Step* step) {
    // Получаем предшествующий объем
    const G4VPhysicalVolume* preVol = step->GetPreStepPoint()->GetPhysicalVolume();
    if (!preVol) return;

    // Получаем логический объем кристалла
    G4LogicalVolume* crystalLV = NpsmBenchDetectorConstruction::GetCrystalLogicalVolume();
    if (!crystalLV) return;

    // Энергия, депонированная в кристалле
    double edep = step->GetTotalEnergyDeposit();
    // У G4StepPoint нет GetLogicalVolume(); логический объём берётся через
    // физический объём пре-шага (preVol получен выше).
    if (edep > 0 && preVol->GetLogicalVolume() == crystalLV) {
        fEventAction->AddEdep(edep / MeV);
    }

    // Только для первичного трека (parent ID == 0)
    if (step->GetTrack()->GetParentID() != 0) return;

    // Получаем процесс, определивший шаг
    const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
    std::string name = "";
    if (process) {
        name = process->GetProcessName();
    }

    // Проверяем, находится ли постшаг внутри кристалла
    const G4VPhysicalVolume* postVol = step->GetPostStepPoint()->GetPhysicalVolume();
    bool inCrystal = (postVol && postVol->GetLogicalVolume() == crystalLV);

    if (inCrystal) {
        // Если шаг завершен внутри кристалла, обновляем статистику
        if (name == "compt") {
            fEventAction->AddCompt();
        } else if (name == "Rayl") {
            fEventAction->AddRayl();  // Почему Rayleigh считается отдельно: это упругий процесс, меняющий только направление, а не количество Compton-взаимодействий
        } else if (name == "phot") {
            fEventAction->SetPhotAbsorbed();
        } else if (name == "conv") {
            fEventAction->SetConv();
        }
    } else if (!postVol) {
        // Если трек покинул мир, помечаем как ушедший
        fEventAction->SetEscaped();
    }
}
```

### NpsmBenchDetectorConstruction.hh

```cpp
#pragma once

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"

class NpsmBenchDetectorConstruction : public G4VUserDetectorConstruction {
public:
    static double gCrystalXMm;
    static double gCrystalYMm;
    static double gCrystalZMm;

    static double SourceRadiusMm();
    static G4LogicalVolume* GetCrystalLogicalVolume();

    G4VPhysicalVolume* Construct() override;

private:
    static G4LogicalVolume* fgCrystalLV;
};
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
