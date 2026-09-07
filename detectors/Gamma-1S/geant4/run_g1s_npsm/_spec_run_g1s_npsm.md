# Спека: проект `run_g1s_npsm` — Гамма-1С с моделью непропорциональности

Ты — генератор кода C++17 (Geant4 11.2.1) и CMake. Верни ТРИ ПОЛНЫХ файла,
без пояснений вне кода. Комментарии — по-русски. Формат строго такой:

```
=== FILE: G1sNpsmPrimaryGeneratorAction.hh ===
<полный текст>
=== FILE: G1sNpsmPrimaryGeneratorAction.cc ===
<полный текст>
=== FILE: main.cc ===
<полный текст>
```

## Назначение

Этап 5: расчёт отклика спектрометра Гамма-1С на аттестованный точечный
источник ¹³⁷Cs с включённой моделью непропорциональности сцинтилляции, для
сверки с измеренным спектром.

## Что УЖЕ существует и должно быть использовано БЕЗ ИЗМЕНЕНИЙ

Эти классы собираются в тот же бинарник (пути даны в CMake, который писать НЕ
нужно — он уже готов). Их заголовки включать, реализацию не дублировать:

| класс | заголовок | назначение |
|---|---|---|
| `G1SDetector` | `G1SDetector.hh` | полная геометрия: кристалл NaI 63×63, банка, отражатель, ФЭУ, свинцовый экран |
| `Rc103FieldPhysicsList` | `Rc103FieldPhysicsList.hh` | физ-лист option4, конструктор `(double cutMm, const std::string& deexMode)` |
| `NpsmLightYield` | `NpsmLightYield.hh` | модель Пейна |
| `NpsmBenchRunAction` | `NpsmBenchRunAction.hh` | накопление спектра и запись CSV |
| `NpsmBenchEventAction` | `NpsmBenchEventAction.hh` | сбор события |
| `NpsmBenchSteppingAction` | `NpsmBenchSteppingAction.hh` | взвешивание светового выхода |

**Точные сигнатуры, которые нельзя менять:**

```cpp
NpsmBenchRunAction(std::string outCsv, double energyKeV, long long nEventsRequested,
                   long seed, const NpsmLightYield* lightYield);
NpsmBenchEventAction(NpsmBenchRunAction* runAction);          // аргумент ОБЯЗАТЕЛЕН
NpsmBenchSteppingAction(NpsmBenchEventAction* eventAction,
                        const NpsmLightYield* lightYield);     // оба ОБЯЗАТЕЛЬНЫ
Rc103FieldPhysicsList(double cutMm = 0.05, const std::string& deexMode = "std");
```

🔴 Конструкторов по умолчанию у этих трёх классов НЕТ. Вызов
`new NpsmBenchEventAction()` или `new NpsmBenchSteppingAction()` без
аргументов — ошибка компиляции; так уже было сделано в первой генерации.

Статические поля физ-листа: `Rc103FieldPhysicsList::gCutMm`,
`::gDeexMode`, `::gPixeModel`, `::gPixeElecModel`.

Класс `G1SDetector` — наследник `G4VUserDetectorConstruction`, поля:
`fWithShield` (bool), `fWithVessel` (bool), `fCrystalLV` (`G4LogicalVolume*`,
заполняется в `Construct()`), метод `ReportMasses() const`.

## Файл 1: `G1sNpsmPrimaryGeneratorAction.hh` / Файл 2: `.cc`

Класс `G1sNpsmPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction`.

**Точечный ИЗОТРОПНЫЙ источник** в точке `(0, 0, gSourceZmm)` мм мировой
системы (начало координат — центр кристалла, +Z к входному торцу).

- статические поля: `static double gSourceZmm;` (умолчание 91.0),
  `static double gEnergyKeV;` (умолчание 661.657);
- конструктор без аргументов, внутри `G4ParticleGun fGun{1};` частица `gamma`;
- в `GeneratePrimaries(G4Event*)`: задать позицию
  `G4ThreeVector(0, 0, gSourceZmm * mm)`, энергию `gEnergyKeV * keV` и
  **изотропное направление**, разыгранное правильно:
  `cosTheta = 1 - 2*G4UniformRand()`, `phi = 2*pi*G4UniformRand()`,
  `sinTheta = sqrt(1 - cosTheta*cosTheta)`, направление
  `(sinTheta*cos(phi), sinTheta*sin(phi), cosTheta)`.
  ⚠ Разыгрывать надо КОСИНУС полярного угла, а не сам угол: `theta` из
  равномерного распределения даёт сгущение к полюсам и неверный телесный угол.
- затем `fGun.GeneratePrimaryVertex(anEvent)`.

## Файл 3: `main.cc`

Точка входа. Разбор ключей вида `имя=значение` из `argv` (порядок любой,
неизвестный ключ — сообщение и возврат 2):

| ключ | умолчание | смысл |
|---|---|---|
| `out=` | `g1s_npsm.csv` | путь выходного CSV |
| `events=` | `1000000` | число историй |
| `seed=` | `12345` | зерно |
| `energy_keV=` | `661.657` | энергия фотона |
| `src_z_mm=` | `91.0` | положение источника по Z |
| `npsm=` | `on` | `on`/`off` — модель непропорциональности |
| `npsm_eta=` | `0.453` | параметр η модели |
| `npsm_s_ons=` | `36.4` | S_Онз |
| `npsm_s_trap=` | `12.0` | S_захв |
| `npsm_s_birks=` | `185.3` | S_Биркс |
| `cut_mm=` | `0.05` | порог продукции |
| `deex=` | `deex` | режим деэкситации |
| `pixe_model=` | пусто | модель сечений PIXE |
| `pixe_e_model=` | пусто | электронная модель PIXE |
| `shield=` | `1` | `1`/`0` — свинцовый экран |
| `masses=` | `1` | `1`/`0` — печатать массы тел |

**Умолчания параметров модели — ТАБЛИЧНЫЕ PAYNE** (η=0,453; 36,4; 12,0; 185,3),
а НЕ значения Breitenmoser: этап 5 согласован именно на них.

Порядок в `main`:

1. Разобрать ключи. Значения `npsm_*` и `cut_mm`/`deex`/`pixe_*` положить в
   соответствующие статические поля до создания объектов.
2. Создать `G4RunManager`. Тип — последовательный (`G4RunManager`), не MT.
3. Создать `G1SDetector`, задать `fWithShield` по ключу `shield`,
   `fWithVessel = false` (точечный источник, сосуда нет).
   Передать в `SetUserInitialization`.
4. Создать и передать `Rc103FieldPhysicsList(cutMm, deexMode)`.
5. `Initialize()`.
6. **После `Initialize()`** взять `G4LogicalVolume*` кристалла из
   `G1SDetector::fCrystalLV` и передать его в
   `NpsmBenchSteppingAction::SetCrystalLogicalVolume(...)` — статический
   сеттер (он добавлен в донора отдельно, считай, что он есть).
   Если указатель нулевой — сообщение и возврат 3, молча не продолжать.
7. Если `masses=1` — вызвать `ReportMasses()`.
8. Создать объект `NpsmLightYield lightYield;` (именно ОБЪЕКТ, не указатель),
   задать `lightYield.SetParameters(eta, ons, trap, birks)` и
   `lightYield.SetEnabled(npsm == "on")`.
   🔴 **Нулевой указатель модели передавать ЗАПРЕЩЕНО:** и
   `NpsmBenchSteppingAction`, и `NpsmBenchRunAction` на нулевой модели
   аварийно завершаются. Выключение модели делается ТОЛЬКО через
   `SetEnabled(false)`, объект при этом существует.
9. Зарегистрировать `NpsmBenchRunAction`, `NpsmBenchEventAction`,
   `NpsmBenchSteppingAction`, `G1sNpsmPrimaryGeneratorAction`.
10. `BeamOn(events)`.
11. Вернуть 0.

**Печать постановки.** До `BeamOn` напечатать одной строкой:
`G1S_NPSM_SETUP: events=… seed=… energy_keV=… src_z_mm=… npsm=… eta=… ons=… trap=… birks=… cut_mm=… deex=… shield=…`
Это проверяемый артефакт #CFG-2: постановка обязана лежать в логе, а не в
имени файла.

## Ограничения

- Только заголовки Geant4 и стандартная библиотека.
- Никаких `using namespace std;`.
- Никаких тихих умолчаний при неразобранном ключе — сообщение и ненулевой код.
- Не изобретать своих классов геометрии, физ-листа, RunAction — они заданы выше.
