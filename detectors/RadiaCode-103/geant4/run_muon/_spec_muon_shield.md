Сгенерируй ОДИН файл C++ — `Rc103MuonDetectorConstruction.cc` ЦЕЛИКОМ, взамен существующего. Ответ: только код, без пояснений и без markdown-ограждения.

# Что меняется и зачем

В мюонный прогон добавляется ОПЦИОНАЛЬНЫЙ свинцовый домик — тот же, что уже
реализован в `run_field/Rc103FieldDetectorConstruction.cc` (эталон, строки 80–194).
Без него нельзя сравнить измеренный фон в домике с расчётом: свинец гамма-фон
ослабляет примерно в 20 раз, а космические мюоны — нет, поэтому в жёсткой части
спектра остаётся именно мюонная компонента, и её надо посчитать в той же защите.

Всё остальное поведение файла обязано остаться ДОСЛОВНО прежним: мир из `G4_AIR`,
чтение прибора из GDML, поиск объёмов с fallback, посадка прибора в (0,0,0),
`pSurfChk=true`, статический указатель на кристалл. При `shield=off` (дефолт)
новый код не должен исполняться вообще.

# Исходный файл (сохранить как есть)

Ниже — текущий файл целиком; правки вносятся ТОЧЕЧНО, остальное переносится без
изменений, включая все комментарии и тексты сообщений.

```cpp
#include "Rc103MuonDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"

#include <cctype>
#include <cstdio>
#include <cstdlib>

G4LogicalVolume* Rc103MuonDetectorConstruction::fgCrystalLV = nullptr;

Rc103MuonDetectorConstruction::Rc103MuonDetectorConstruction(
    const G4String& gdmlPath, double worldHalfMm)
    : fGdmlPath(gdmlPath), fWorldHalfMm(worldHalfMm) {}

namespace {
// Детерминированный поиск логического объёма по имени с фоллбэком на
// ЕДИНСТВЕННОЕ substring-совпадение — паттерн дословно из
// run_field/Rc103FieldDetectorConstruction.cc (эталон), не гадание.
G4LogicalVolume* FindVolumeOrFallback(const char* exactName,
                                      const char* substrLower) {
  auto* lvStore = G4LogicalVolumeStore::GetInstance();
  G4LogicalVolume* lv = lvStore->GetVolume(exactName, false);
  if (lv) return lv;

  std::fprintf(stderr,
               "Rc103MuonDetectorConstruction: '%s' NOT FOUND verbatim after "
               "GDML parse. Dumping ALL logical volumes in store (%zu total):\n",
               exactName, lvStore->size());
  for (auto* v : *lvStore) {
    std::fprintf(stderr, "  - '%s'\n", v->GetName().c_str());
  }

  G4LogicalVolume* candidate = nullptr;
  int nCandidates = 0;
  for (auto* v : *lvStore) {
    G4String lower = v->GetName();
    for (auto& c : lower) {
      c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    if (lower.find(substrLower) != G4String::npos) {
      candidate = v;
      ++nCandidates;
    }
  }
  if (nCandidates == 1) {
    std::fprintf(stderr,
                 "Rc103MuonDetectorConstruction: fallback matched EXACTLY one "
                 "candidate by substring '%s': '%s'\n",
                 substrLower, candidate->GetName().c_str());
    return candidate;
  }
  std::fprintf(stderr,
               "Rc103MuonDetectorConstruction: fallback found %d candidates "
               "(need exactly 1) - refusing to guess.\n",
               nCandidates);
  return nullptr;
}
}  // namespace

G4VPhysicalVolume* Rc103MuonDetectorConstruction::Construct() {
  auto* nist = G4NistManager::Instance();

  const double halfWorld = fWorldHalfMm * mm;
  auto* worldLV = new G4LogicalVolume(
      new G4Box("world_muon", halfWorld, halfWorld, halfWorld),
      nist->FindOrBuildMaterial("G4_AIR"), "world_muon");
  auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "world_muon", nullptr,
                                    false, 0, true);

  fParser.SetOverlapCheck(true);   // не глушить известные наложения
  fParser.Read(fGdmlPath, false);  // false = не валидировать против XSD

  // Берём "RC103_device_log" (сам прибор), НЕ "World" из парсера — тот
  // остаётся orphan-объёмом в сторе и из Construct() не возвращается.
  G4LogicalVolume* deviceLV =
      FindVolumeOrFallback("RC103_device_log", "device_log");
  if (!deviceLV) {
    std::fprintf(stderr,
                 "Rc103MuonDetectorConstruction: FATAL - could not resolve "
                 "RC103_device_log after GDML parse. Aborting.\n");
    std::abort();
  }

  G4LogicalVolume* crystalLV = FindVolumeOrFallback("Crystal_log", "crystal");
  if (!crystalLV) {
    std::fprintf(stderr,
                 "Rc103MuonDetectorConstruction: FATAL - could not resolve "
                 "Crystal_log after GDML parse. Aborting.\n");
    std::abort();
  }
  fgCrystalLV = crystalLV;

  // Прибор в (0,0,0) без поворота: корпус Case_outer 123 x 34 x 17.5 мм,
  // длинная ось по X, по Z всего +-8.75 мм. pSurfChk=true последним
  // аргументом — обязательная проверка наложений с миром.
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), deviceLV,
                    "pv_rc103_in_muon", worldLV, false, 0, true);

  std::fprintf(stdout,
               "Rc103MuonDetectorConstruction: world half=%.1f mm, device "
               "placed at origin, Crystal_log resolved.\n",
               fWorldHalfMm);
  return worldPV;
}

```

# Правка 1 — конструктор

Сигнатура становится:

```cpp
Rc103MuonDetectorConstruction::Rc103MuonDetectorConstruction(
    const G4String& gdmlPath, double worldHalfMm, bool shieldOn, double zDiskMm)
    : fGdmlPath(gdmlPath), fWorldHalfMm(worldHalfMm), fShieldOn(shieldOn),
      fZDiskMm(zDiskMm) {}
```

Поля `fShieldOn` и `fZDiskMm` объявлены в заголовке (его правит человек, не ты).

# Правка 2 — метод `BuildLeadShield`

Добавь приватный метод

```cpp
void Rc103MuonDetectorConstruction::BuildLeadShield(G4LogicalVolume* worldLV,
                                                    G4LogicalVolume* deviceLV)
```

Геометрия — ДОСЛОВНО как в эталоне, числа НЕ пересчитывать и не «улучшать»:

- `kShieldPbMm = 50.0` — стенки и дно;
- полость `kShieldCavityXMm = kShieldCavityYMm = 150.0`, `kShieldCavityZMm = 385.0`;
- наружный габарит `250 x 250 x 435` мм (`X = cavity + 2*Pb`, `Z = cavity + Pb`: дно есть, крышки нет);
- домик центрирован в (0,0,0), прибор тоже в (0,0,0);
- полуразмеры: `outHX = outHY = 125`, `outHZ = 217.5`, `cavHX = cavHY = 75`;
- `cavZMin = -outHZ + kShieldPbMm` (= −167.5, верх свинцового дна), `cavZMax = +outHZ` (= +217.5, срез открытого верха).

Эти константы объявляются в заголовке как `static constexpr` — в .cc их не дублировать,
обращаться по именам.

Тело метода:

1. `auto* nist = G4NistManager::Instance();`
2. Полуразмеры по формулам выше.
3. **Проверка вместимости.** Солид прибора привести `dynamic_cast<G4Box*>(deviceLV->GetSolid())`;
   если каст не удался — `std::fprintf(stderr, ...)` с именем солида и `std::abort()`
   (габарит нельзя проверить фактом, молча строить нельзя). Иначе взять
   `devHX/devHY/devHZ` через `GetXHalfLength()/mm` и проверить, что прибор с центром
   в (0,0,0) укладывается в полость по всем трём осям (по z — в интервал
   `[cavZMin, cavZMax]`). Напечатать габарит прибора, размеры полости и все четыре
   зазора; при несоответствии — `FATAL` в stderr и `std::abort()`.
4. **Проверка вложенности в мир.** Перебрать все 8 углов наружного габарита, взять
   максимум расстояния от (0,0,0), напечатать его. Если `outHX >= fWorldHalfMm` или
   `outHY >= fWorldHalfMm` или `outHZ >= fWorldHalfMm` — `FATAL` с обоими наборами
   чисел и `std::abort()`.
5. **Проверка, что диск старта мюонов ВЫШЕ домика** — это НОВАЯ проверка, в эталоне
   её нет и быть не могло (там источник — сфера). Если `fZDiskMm <= cavZMax` —
   `FATAL` в stderr с текстом о том, что диск старта оказался бы внутри полости и
   мюоны рождались бы за защитой, с печатью `fZDiskMm` и `cavZMax`, затем `std::abort()`.
   Запас не требовать, достаточно строгого неравенства.
6. **Солид.** Наружный `G4Box("shield_pb_outer", outHX*mm, outHY*mm, outHZ*mm)`
   минус вычитаемый `G4Box("shield_pb_cut", cavHX*mm, cavHY*mm, cutHZ*mm)`, где
   `cutHZ = 0.5*(cavZMax - cavZMin) + kShieldPbMm` и центр выреза
   `cutZ0 = cavZMin + cutHZ`. Вычитаемый бокс НАРОЧНО выступает за верхнюю грань —
   так верх остаётся открытым, а дно 50 мм сохраняется. Смещение задать
   `G4ThreeVector(0, 0, cutZ0*mm)` в `G4SubtractionSolid`.
7. **Объём и постановка.** `G4LogicalVolume(shieldSolid, nist->FindOrBuildMaterial("G4_Pb"), "shield_pb_log")`,
   затем `new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), shieldLV, "pv_shield_pb", worldLV, false, 0, true)`
   — последний аргумент `true` есть обязательная проверка наложений.
   Записать `fgShieldLV = shieldLV;`.
8. Итоговый `std::fprintf(stdout, ...)`: материал, толщина стенок и дна, размеры
   полости, «верх ОТКРЫТ», центр полости по z и то, что прибор в (0,0,0).

Все сообщения печатать с префиксом `Rc103MuonDetectorConstruction: SHIELD ...`.

# Правка 3 — вызов в `Construct()`

После постановки прибора (`pv_rc103_in_muon`) и ДО финального `std::fprintf`
добавить:

```cpp
  if (fShieldOn) {
    BuildLeadShield(worldLV, deviceLV);
  }
```

В финальном сообщении дополнительно печатать состояние домика — словом `ON` либо
`OFF`, и высоту диска старта `fZDiskMm`.

# Правка 4 — статический указатель

Рядом с `fgCrystalLV` определить `G4LogicalVolume* Rc103MuonDetectorConstruction::fgShieldLV = nullptr;`

# Заголовки

Добавить к существующим включениям: `G4SubtractionSolid.hh`, `G4NistManager.hh`
(если его ещё нет), `<algorithm>`, `<cmath>`. Ничего лишнего не подключать.

# Провенанс — обязательный комментарий

В шапке файла, после существующих комментариев, добавить абзац: геометрия домика
перенесена 29.08.2026 из `run_field/Rc103FieldDetectorConstruction.cc` (строки
80–194), числа сообщены оператором 15.08.2026; копия сделана потому, что классы
живут в разных сборочных проектах, а вынесение в общий модуль потребовало бы
пересборки `run_field` во время идущих прогонов. При правке одной копии обязана
правиться вторая. Файлы `geometry/RCDetector.*` источником не являются — они
запрещены оператором 27.08.2026.
