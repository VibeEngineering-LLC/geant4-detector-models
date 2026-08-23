Задача: написать Python-скрипт `regression_check.py` — объективный шлюз §34-петли
«регрессия физики на эталонных задачах» контура GEANT4. Верни ПОЛНЫЙ код файла
целиком, без markdown fences, без пояснений до/после — чистый Python.

## Контекст

Скрипт запускает уже собранный бинарник `rc_curves.exe` на двух коротких
эталонных макросах (моноэнергетический источник, кристалл CsI(Tl) детектора
RadiaCode-103) и сравнивает число зарегистрированных событий (`hits`) с
сохранённым эталоном. Если бинарник не пересобирался (код физики/геометрии не
менялся) — RNG-seed НЕ задаётся явно нигде в проекте, поэтому Geant4 использует
детерминированный дефолтный движок: результат при том же коде ПОЛНОСТЬЮ
воспроизводим (не просто статистически близок). Из-за этого допуск сравнения
(3σ по Пуассону от эталонного `hits`) в норме будет давать точное совпадение —
он существует "с запасом" на случай, если seed когда-нибудь станет
рандомизированным, а не потому что реальная дисперсия измерена и ненулевая.
Это должно быть явно написано в комментарии в шапке файла — не выдавать
допуск за измеренную величину (правило #AH-1 проекта: конфигурационное
ожидание — не то же самое, что измерение).

## Константы

```python
GEANT4_ENV = {
    "GEANT4_ROOT": r"C:\geant4",
    "G4NEUTRONHPDATA": r"C:\geant4\share\data\G4NDL4.7",
    "G4LEDATA": r"C:\geant4\share\data\G4EMLOW8.5",
    "G4LEVELGAMMADATA": r"C:\geant4\share\data\PhotonEvaporation5.7",
    "G4RADIOACTIVEDATA": r"C:\geant4\share\data\RadioactiveDecay5.6",
    "G4PARTICLEXSDATA": r"C:\geant4\share\data\G4PARTICLEXS4.0",
    "G4SAIDXSDATA": r"C:\geant4\share\data\G4SAIDDATA2.0",
    "G4ABLADATA": r"C:\geant4\share\data\G4ABLA3.3",
    "G4INCLDATA": r"C:\geant4\share\data\G4INCL1.2",
    "G4ENSDFSTATEDATA": r"C:\geant4\share\data\G4ENSDFSTATE2.3",
}
EXE = r"D:\Claude_files\repos\geant4-detector-models\build\RadiaCode-103\rc_curves.exe"
HERE = папка самого скрипта (os.path.dirname(os.path.abspath(__file__)))
BASELINE_PATH = HERE / "baseline.json"
STATE_PATH = HERE / "state.json"
REFS = [
    {"name": "cs137_662", "macro": "ref_cs137.mac"},
    {"name": "k40_1461", "macro": "ref_k40.mac"},
]
SIGMA_K = 3  # множитель допуска (3 сигма Пуассона)
```

## Функция `run_macro(macro_name)`

Собирает `env` копией `os.environ`, добавляет в PATH `C:\geant4\bin;` в начало,
устанавливает переменные из `GEANT4_ENV`. Запускает через
`subprocess.run([EXE, macro_path, "bare"], env=env, capture_output=True,
text=True, timeout=120)` (`macro_path` = `os.path.join(HERE, macro_name)`).
Ищет в `result.stdout` строку, начинающуюся с `RESULT` (например
`RESULT E_keV= 662 N= 30000 hits= 8116 eff_total= 0.270533 file= ...`),
парсит `hits=` числом после `hits=` до следующего пробела (regex
`r"hits=\s*(\d+)"`), возвращает `int`. Если процесс упал (returncode != 0)
или строка `RESULT` не найдена — `raise RuntimeError` с текстом, включающим
`result.stdout[-2000:]` и `result.stderr[-2000:]` (для диагностики).

## Функция `load_baseline()`

Читает `BASELINE_PATH` (JSON), если файла нет — возвращает `None`.
Формат файла: `{"cs137_662": {"hits": 8116, "date": "2026-08-24"}, "k40_1461": {...}}`.

## Функция `save_state(results)`

Пишет `STATE_PATH` (JSON, `ensure_ascii=False`, `indent=2`): список по каждому
эталону `{"name":..., "baseline_hits":..., "new_hits":..., "tolerance":...,
"status": "OK"|"FAIL"}`, плюс верхнеуровневое поле `"overall": "OK"|"FAIL"` и
`"checked_at"` (принять как аргумент функции строку, не вычислять `datetime.now()`
внутри — вызывающий код передаёт готовую строку времени).

## `main()`

1. Разобрать `sys.argv`: если есть флаг `--update-baseline` — режим обновления
   эталона (см. п.5), иначе обычная проверка.
2. `baseline = load_baseline()`. Если `baseline is None` и НЕ `--update-baseline` —
   вывести `"БАЗОВЫЙ ЭТАЛОН НЕ НАЙДЕН — запустите с --update-baseline один раз"`,
   `sys.exit(2)`.
3. Для каждого `ref in REFS`: `hits = run_macro(ref["macro"])`, вывести строку
   `"[%s] hits=%d" % (ref["name"], hits)`.
4. **Обычный режим:** для каждого эталона взять `base_hits = baseline[ref["name"]]["hits"]`,
   допуск `tol = max(1, round(SIGMA_K * (base_hits ** 0.5)))`, `status = "OK" if
   abs(hits - base_hits) <= tol else "FAIL"`. Печатать построчно:
   `"  %-12s baseline=%-8d new=%-8d допуск=±%-6d -> %s"`. Собрать список
   результатов, вызвать `save_state(results_with_status, время)`
   (`time.strftime("%Y-%m-%d %H:%M:%S")`). Если хотя бы один `FAIL` —
   напечатать `"РЕГРЕССИЯ: физика изменилась за пределы допуска"`,
   `sys.exit(1)`. Иначе `"OK: все эталоны в допуске"`, `sys.exit(0)`.
5. **`--update-baseline` режим:** записать текущие `hits` как новый эталон в
   `BASELINE_PATH` (JSON, тот же формат, добавить `"date"` = сегодняшняя дата
   строкой), напечатать `"Эталон обновлён: <path>"`, `sys.exit(0)`. Этот режим
   НЕ сравнивает ничего — просто фиксирует новое базовое значение (используется
   человеком осознанно после проверенного изменения физики, не автоматически).

## Требования к коду

- `sys.stdout.reconfigure(encoding="utf-8")` в начале файла.
- Только стандартная библиотека (`subprocess`, `json`, `os`, `sys`, `re`, `time`,
  `pathlib`) — без внешних зависимостей.
- Все печатаемые строки — по-русски (диагностика для оператора).
- Комментарий в шапке файла (docstring) — 5-8 строк: что скрипт делает, почему
  допуск "с запасом", а не измеренная величина (см. раздел "Контекст" выше),
  ссылка "часть §34-петли контура GEANT4, карточка audit/loops.md".
