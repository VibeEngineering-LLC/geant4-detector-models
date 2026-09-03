Напиши ОДИН файл Python 3.11 `calib_nonlinear.py`. Выведи только код, без пояснений и без markdown-ограждений.

# Назначение
#CAL-1 (03.09.2026). Проверка, нужна ли нелинейная энергетическая шкала измеренному спектру RadiaCode-103: по реперным линиям фона измеряются центроиды пиков в КАНАЛАХ, затем сравниваются модели шкалы «энергия от канала»: полиномы степени 2, 3, 4 и модель с насыщением SiPM. Штатное ПО прибора ограничено степенью 2; в работе допустима любая до 4 (решение оператора 03.09.2026). Прогонов Geant4 нет.

# CLI
`python calib_nonlinear.py [--spectrum <путь.xml>] [--out <путь.json>] [--selftest]`
- `--spectrum` по умолчанию: `os.path.join(os.environ.get("G4MODELS_MEASURED", r"C:\g4work\measured\RadiaCode-103"), "Фон 7 дней без домика.xml")`.
- `--out` по умолчанию `<каталог скрипта>/out/calib_nonlinear.json` (каталог создать).
- `--selftest` — самопроверка без спектра (см. ниже), код возврата 0 при успехе, 1 при провале.

# Точные сигнатуры донора (переиспользование, НЕ переписывать)
```python
HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "analysis"))          # read_rcxml
DONOR = r"D:\GoogleDrive\Дозиметрия\ИИ\1 Скилы\0_Work\gamma-spectrum-analysis\scripts"
sys.path.insert(0, DONOR)
import read_rcxml
import gamma.calibration.fwhm_measure as FM
from gamma.calibration.energy_fit import polynomial_energy_fit
```
- `read_rcxml.read(path)` возвращает список объектов; брать `[0]`. У объекта: `.counts` (numpy-массив отсчётов по каналам, ПОСЛЕДНИЙ канал — переполнение, использовать `counts[:-1]`), `.coef` — кортеж `(c0, c1, c2)` коэффициентов прибора «энергия = c0 + c1·ch + c2·ch²», `.live` — живое время, с.
- `FM.measure_fwhm(counts, *, energy_keV=<кэВ>, energy_cal=<список коэффициентов от МЛАДШЕГО>, window_factor=1.25, methods=(FM.METHOD_HALF_MAX, FM.METHOD_GAUSSIAN), known_lines_keV=<список кэВ всех реперов>, reject_blended=False)` возвращает объект с полями `.passed` (bool), `.centroid_channel` (float или None), `.centroid_keV`, `.fwhm_keV`, `.fwhm_uncertainty_keV`, `.significance_sigma`. Если `.passed` ложно или `centroid_channel is None` — линия НЕ используется в подгонке, но заносится в отчёт с пометкой «не измерена».
- `polynomial_energy_fit(channels, energies, *, max_degree=d, min_degree=d, target_residual_keV=None)` возвращает объект `EnergyFitResult`; у него есть поля `coefficients` (список от МЛАДШЕГО к старшему, совместим с `numpy.polynomial.polynomial.polyval`), `residuals_keV` (список невязок E_ref − E_fit по точкам) и `degree`. Если каких-то из этих полей нет — обращаться через `getattr(res, имя, None)` и при отсутствии `residuals_keV` считать невязки самому через `numpy.polynomial.polynomial.polyval(channels, coefficients)`.

# Реперные линии (кэВ, имя)
```
(74.2,  "Pb K-alpha (72.80+74.97, смесь)"),
(85.5,  "Pb K-beta (84.94+87.3, смесь)"),
(186.0, "Ra-226 186.2 + U-235 185.7 (смесь)"),
(238.63,"Pb-212"), (295.22,"Pb-214"), (351.93,"Pb-214"),
(511.0, "аннигиляция"), (583.19,"Tl-208"), (609.31,"Bi-214"),
(911.20,"Ac-228"), (1120.29,"Bi-214"), (1460.82,"K-40"),
(1764.49,"Bi-214"), (2614.51,"Tl-208"),
```
Центроиды искать, используя ШКАЛУ ПРИБОРА `smp.coef` как `energy_cal` (она нужна только чтобы найти окно пика). Для подгонки использовать `centroid_channel`.

# Модели шкалы (подгонка «канал → энергия» по прошедшим линиям)
1. `poly2`, `poly3`, `poly4` — через `polynomial_energy_fit` с `min_degree=max_degree=d`.
2. `sat` — модель насыщения SiPM: `E = a0 + a1·g + a2·g²`, где `g = −ch0·ln(1 − ch/ch0)`, `ch0 > max(ch)`. Подгонять `scipy.optimize.least_squares` по параметрам `(a0, a1, a2, ch0)`, стартовые `a0,a1,a2` = коэффициенты прибора, `ch0 = 5·max(ch)`, границы: `ch0` в `[1.05·max(ch), 1e6]`.
Для каждой модели: коэффициенты, невязка по каждой линии `E_ref − E_model(ch)` в кэВ и в долях FWHM этой линии (если FWHM измерена), RMS невязки, максимум |невязки| отдельно для линий ниже 300 кэВ и выше 1500 кэВ, число параметров, χ² = Σ(невязка/σ)² с σ = `fwhm_keV/2.355/sqrt(значимость)` если доступны, иначе σ = 1 кэВ.
Дополнительно две «референсные» шкалы без подгонки: `device` (коэффициенты прибора) и `cal_room = [-3.711311, 2.444318, 0.000321]` (шкала, принятая в конвейере) — те же невязки по линиям.

# Вывод
- В stdout: таблица линий (имя, E_ref, канал центроида, FWHM, значимость, статус), затем таблица моделей (модель, число параметров, RMS, max|res|<300, max|res|>1500, χ²), затем для каждой модели строка коэффициентов.
- В JSON (`--out`): всё то же структурой `{"spectrum": ..., "device_coef": [...], "lines": [...], "models": {...}}`, `ensure_ascii=False`, отступ 2.

# Самопроверка `--selftest` (обязательна, должна уметь краснеть)
1. Позитив: сгенерировать 12 точек `ch` от 30 до 900, `E = 2.0 + 2.4·ch + 3e-4·ch² + 1e-7·ch³` (точные, без шума). `poly3` и `poly4` обязаны дать max|невязка| < 1e-6 кэВ; `poly2` обязана дать max|невязка| > 0.5 кэВ (кубический член существенный). Если хоть одно условие нарушено — напечатать «SELFTEST FAIL: ...» и вернуть 1.
2. Негатив-мутация: испортить одну точку на +15 кэВ и убедиться, что `poly3` теперь даёт max|невязка| > 1 кэВ (детектор чувствует порчу). Иначе «SELFTEST FAIL».
3. Модель `sat`: сгенерировать `E` точно по формуле насыщения с `ch0=2000, a0=0, a1=2.4, a2=0`, подогнать, требовать |ch0_fit − 2000| / 2000 < 0.05 и max|невязка| < 0.05 кэВ.
При успехе напечатать «SELFTEST OK» и вернуть 0.

# Технические требования
- `# -*- coding: utf-8 -*-`, docstring по-русски (назначение, CLI, модели, самопроверка), `sys.stdout.reconfigure(encoding="utf-8")` до первого print.
- numpy, scipy, json, argparse; matplotlib НЕ использовать.
- Функции: `measure_lines(counts, device_coef) -> list[dict]`, `fit_models(lines) -> dict`, `residual_table(model_name, coef_or_params, lines) -> dict`, `selftest() -> int`, `main(argv) -> int`, блок `if __name__ == "__main__": sys.exit(main(sys.argv[1:]))`.
- Импорты донора (read_rcxml, FM, polynomial_energy_fit) делать ВНУТРИ `main` после разбора аргументов, чтобы `--selftest` работал без донора и без спектра; в `selftest` для полиномов использовать `numpy.polynomial.polynomial.polyfit`/`polyval` напрямую.
- Все сообщения — по-русски.
