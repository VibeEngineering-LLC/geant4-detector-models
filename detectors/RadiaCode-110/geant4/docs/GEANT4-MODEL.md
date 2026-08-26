# Geant4-модель RadiaCode-110

Полный комплект геометрии RC-110 для Monte Carlo (Geant4): корпус, электроника, детектор.

**Источники:** `RC110_model.blend`, `Rc-110.stl`, рентген/фото, [спека RC-110](https://radiacode.com/docs/en/100-series/devices/100-series-introduction/technical-specification).

**Обновлено:** 2026-08-26 — закрыты 3 геометрических дефекта #AUD-6 (капсула/SiPM/смещение узла)
и расхождение ESR↔GDML (6-я пластина на −Z с вырезом 6.4×6.4 под SiPM). Ранее (2026-08-25):
добавлены плата, LiPo, дисплей, USB; файлы разложены по подпапкам.  
**Session state:** `SESSION-STATE.md` в референсной папке RC-110 (локальные
данные оператора, вне репозитория).

---

## Структура каталогов

```
geant4/
├── docs/
│   └── GEANT4-MODEL.md          ← этот файл
├── gdml/
│   ├── detector/
│   │   └── RC110_detector.gdml  ← SSOT: все внутренности + полый прямоугольный корпус
│   ├── full/
│   │   └── RC110_full.gdml      ← генерируется: STL-корпус + те же внутренности
│   └── case_mesh/
│       └── RC110_case_mesh.gdml ← генерируется: только STL корпус
├── cadmesh/
│   ├── RC110CADMeshCase.hh/.cc  ← загрузчик STL в рантайме (CADMesh)
│   ├── RC110HybridExample.cc
│   └── CMakeLists.example.txt
├── geometry/
│   ├── RC110Detector.hh/.cc     ← геометрия НАПРЯМУЮ в C++ (обход отсутствия GDML)
│   ├── vis_render.cc            ← headless-рендер через TOOLSSG_OFFSCREEN
│   └── CMakeLists.txt           ← таргет rc110_vis_render
├── scripts/
│   ├── stl_to_gdml.py           ← STL → case_mesh + full
│   └── verify_align.py         ← matplotlib, числовая проверка (НЕ генерирует repo-картинку)
└── verify/
    └── RC110_align_check.png    ← рендер Geant4 (rc110_vis_render), не генерируется автоматически
```

**В корне RC110** (не в `geant4/`):

| Файл | Назначение |
|---|---|
| `Rc-110.stl` | Исходный mesh корпуса (110 666 tri) |
| `RC110_model.blend` | Blender-модель с внутренностями |

---

## Система координат

| Ось | Направление |
|---|---|
| **X** | Длина прибора; **USB → +X** |
| **Y** | Ширина |
| **Z** | Толщина; **−Z = лицевая** (значок радиации) |

**Якорь:** центр корпуса `(0, 0, 0)`.

**Remap STL → device** (в `stl_to_gdml.py` и CADMesh):

```
device_x = stl_z − 63.3
device_y = stl_x − 122.5
device_z = stl_y − 122.5
```

---

## Три варианта корпуса

| Вариант | Файл / код | Корпус | Когда использовать |
|---|---|---|---|
| **A. Примитив** | `gdml/detector/RC110_detector.gdml` | Полая ABS-оболочка (бокс минус бокс), 126.6×34.1×21.7, стенка 1.5 мм | **Физика MC** — нет перекрытий, быстро |
| **B. STL в GDML** | `gdml/full/RC110_full.gdml` | 110 666 треугольников из `Rc-110.stl` | Vis / детальная форма; **STL сплошной** → перекрытие с внутренностями |
| **C. CADMesh** | `cadmesh/RC110CADMeshCase.*` | Тот же STL, загрузка в C++ | Гибрид: STL-vis + GDML-примитивы |

---

## Внутренние узлы (все три варианта B/C наследуют из detector GDML)

| Узел | Размер, мм | Центр (X, Y, Z) | Материал (GDML) |
|---|---|---|---|
| Плата FR4 | 101.5 × 29 × 1.0 | (9.85, 0, −5.45) | `PCB_FR4` |
| LiPo DTP962565 | 65 × 25 × 9.6 | (23.8, 0, 4.55) | `LiPo_DTP962565` |
| Окно дисплея | 36.5 × 14 × 1.0 | (−19.25, 0, −8.95) | `Display_window_PC` |
| LCD | 34 × 13 × 2.2 | (−19.25, 0, −7.35) | `LCD_stack` |
| USB Type-C | 7.5 × 8.9 × 3.2 | (58.75, 0, −4.0) | `USB_connector` |
| Капсула | 18³ (стенка 1.5) | модуль @ −50.25 | ABS |
| Кристалл CsI(Tl) | 14³ | @ −50.25 (в капсуле +0.35) | `G4_CESIUM_IODIDE` |
| SiPM | 6 × 6 × 0.8 | грань **−Z**, z=−7.4 локально | `G4_Si` |
| ESR 3M | 65 µm, 6 граней (−Z с вырезом 6.4×6.4 под SiPM) | вокруг кристалла | `ESR_foil` |

**Не моделируется:** отдельные SMD на плате (по решению оператора).

---

## Иерархия GDML

```
World
 └── RC110_device_log
      ├── Case_*          (shell или STL)
      └── Case_interior_log
           ├── PCB_log
           ├── Battery_log
           ├── Display_window_log / Display_LCD_log
           ├── USB_log
           └── DetectorModule_log @ (−50.25, 0, 0)
                └── капсула → кристалл + SiPM + ESR
```

---

## Сборка и проверка

```powershell
cd <референсная папка RC110>\geant4\scripts
python stl_to_gdml.py      # ~30 s, пишет gdml/full и gdml/case_mesh
```

**Требования:** Python 3.

⚠ **`stl_to_gdml.py` требует `Rc-110.stl`** в корне референсной папки
(`geant4/../Rc-110.stl`). STL в репозиторий не входит (5,5 МБ исходного
меша) — брать из референсной папки оператора.

⚠ **`geant4/scripts/verify_align.py` — чужой файл**, под этим именем и
путём лежит скрипт ДРУГОГО детектора (докстрока «Сверка геометрии
AtomSpectra PRO») — для RC-110 нерабочий, не удалять бездумно (возможно,
чужая рабочая копия, не моя зона §12). Настоящий, рабочий `verify_align.py`
для RC-110 лежит ОДНИМ УРОВНЕМ ВЫШЕ — в корне референсной папки
(`<референсная RC-110>\verify_align.py`, НЕ в `geant4/`). Он рисует
2D/3D-схему matplotlib (оси, легенда, цифры совмещения) — полезен для
ЧИСЛОВОЙ проверки, но сам по себе в репозиторий не входит и картинку сюда
больше не пишет (см. ниже).

**`verify/RC110_align_check.png` — теперь настоящий рендер Geant4** (с
26.08, заменил прежнюю matplotlib-версию по решению оператора). Построен
СВОИМ движком Geant4 через `geant4/geometry/RC110Detector.{hh,cc}` —
геометрия описана напрямую в C++ (обход отсутствия GDML-модуля в prebuilt-
сборке, см. раздел «Geant4 на этой машине» ниже) — и `vis_render.cc`,
рендер драйвером `TOOLSSG_OFFSCREEN` (НЕ RayTracer — тот крашится
`0xC0000005` на этой сборке всегда, факт проверен независимо от геометрии).
Стенка корпуса и капсулы полупрозрачны (alpha 0.35/0.5), поэтому видно
кристалл CsI (золотой), плату (зелёная), аккумулятор (серый), USB
(коричневый). Пересобрать: `build/RadiaCode-110/_build.cmd`, таргет
`rc110_vis_render` (CMake+Ninja, требует VS2022 BuildTools, см. `g4setup.ps1`).
Координаты компонентов в `RC110Detector.cc` — плоские (мировой кадр),
взяты из `RC110_detector.gdml` того же дня, независимо перепроверены.
⚠ SiPM/ESR (доли мм) на этом масштабе визуально не различимы — это
особенность рендера, не отсутствие в геометрии. Числовая проверка
совмещения (с осями/легендой) — по-прежнему через matplotlib-скрипт выше,
если нужны координаты как текст, не картинка.

---

## Geant4 на этой машине

- Установка: `C:\geant4` 11.2.1, vis OpenGL есть, **GDML OFF** в prebuilt
  (`Geant4_gdml_FOUND OFF` зашито в `lib\cmake\Geant4\Geant4Config.cmake`
  на этапе конфигурации — не runtime-переключатель).
- Рабочая папка MC: `C:\g4work` (пока RC-103).
- **Визуализация без GDML — сделано:** геометрия перенесена в C++
  (`geant4/geometry/RC110Detector.{hh,cc}`, см. раздел «Сборка и проверка»
  выше) — картинка Geant4 не требует GDML вовсе.
- **Пересборка Geant4 с `GEANT4_USE_GDML=ON`** (для загрузки самого
  `.gdml`-файла, не только для картинки) — отдельная задача, в отдельный
  префикс `C:\geant4-gdml`, текущую `C:\geant4` не трогает; статус — см.
  `SESSION-STATE.md` контура (локально), не входит в этот репозиторий.

### Рекомендации по MC

1. **Self-absorption / спектрометрия:** `gdml/detector/RC110_detector.gdml` (полый бокс + все внутренности).
2. **Визуализация формы корпуса:** `gdml/full/RC110_full.gdml` или CADMesh + detector GDML (без дублирования корпуса в MC).
3. **STL в full.gdml** — tessellated solid = сплошной ABS; для физики overlapping с детектором. Не использовать full для transport без доработки.

### CADMesh гибрид

```cpp
// cadmesh/RC110HybridExample.cc
// 1) RC110CADMeshCase::BuildLogicalVolume(stlPath) — корпус из STL
// 2) G4GDMLParser::Read(".../gdml/detector/RC110_detector.gdml")
//    → разместить RC110_device_log или только Case_interior_log
```

Пути по умолчанию:

- STL: `../../Rc-110.stl`
- GDML: `../gdml/detector/RC110_detector.gdml`

---

## Ограничения модели

- Корпус-примитив: без скруглений, USB-выреза, окна дисплея в ABS-оболочке.
- STL: замкнутый mesh → в GDML считается сплошным.
- Материалы LiPo/LCD/USB — **упрощённые** смеси (не datasheet-grade).
- ESR optical skin на −Z в C++ может потребовать отдельного thin skin (в GDML теперь 6 slab, включая −Z с вырезом 6.4×6.4 под SiPM — закрыто 2026-08-26, ранее грани −Z не было вовсе).

---

## Провенанс

| Утверждение | Источник |
|---|---|
| Корпус 126.6×34.1×21.7 | `Rc-110.stl` bbox |
| Кристалл −50.25 | крестик на задней +Z, SESSION-STATE |
| SiPM на −Z | правка оператора 2026-08-24 |
| PCB/Battery coords | `RC110_model.blend` rebuild (jsonl:186) |
| GDML detector | `geant4/gdml/detector/RC110_detector.gdml` |
