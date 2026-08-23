# SPEC: analysis/make_m1_macros.py — версия 2, с автоматическими «стенами»

Сгенерируй ОДИН файл Python 3. Выведи ТОЛЬКО код, без пояснений и без markdown-
обрамления (без ```). Первая строка `# -*- coding: utf-8 -*-`. Комментарии — по-русски.
Используй `os`, `sys`.

## Контекст (для докстроки модуля, напечатать ДОСЛОВНО как есть)

```
"""Готовит пары макросов (спектр-источник + run) для METODA 1: единичный
(1 Bq/kg) polnyi raspad ODNOGO zvena, nucleusLimits otsekaet dochernie.

VERSIA 2 (21.08, posle /compact): dobavlena funktsia add_walls() - konnrolnyi
progon podtverdil, chto Geant4 GPS (/gps/hist/type arb + /gps/hist/inter Lin)
traktuet tochki (E,weight) kak uzly NEPRERYVNOI plotnosti (trapetsevidnaya
interpolyatsia MEZHDU sosednimi tochkami po ih energeticheskomu RASSTOYANIYU),
a NE kak "soderzhimoe bina", v otlichie ot wallfield.cc, kotoryi pishet dannye
imenno v formate "soderzhimoe 2-keV bina". Rezkaya liniya, okruzhennaya
razrezhennym kontinuumom (sosedняя nenulevaya tochka daleko), rasplyvaetsya
GPS v shirokii treugolnik i (dlya poslednei tochki gistogrammy) OBRYVAETSYA
sprava - proverено faktom na K-40 1460.8 keV: original (bez sten) daval
model_cps sootvetstvuyushchiy 96.4 Bq/kg, s dobavlennymi stenami vplotnuyu
(1.459/1.463 MeV) - 131.1 Bq/kg, protiv 122.1 u nezavisimogo monootklika
(analiz DECISIONS.md D-002/D-003). add_walls() vstavlyaet nulevye tochki
VPLOTNUYU k rodnoi setke (2 keV nizhe SPLIT_KEV, 10 keV vyshe) na kazhdom
razryve mezhdu sosednimi nenulevymi tochkami i na kraiah (posledn'aya tochka
vsegda poluchaet "stenu" sprava - inache GPS obryvaet raspredelenie).

Zapusk: python make_m1_macros.py
"""
```

## Импорты и sys.path (дословно, как в существующей версии)

```python
sys.stdout.reconfigure(encoding="utf-8")

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(_HERE, "..", "..", ".."))
BUILD = os.path.join(REPO, "build", "RadiaCode-103")
RESULTS = os.path.normpath(os.path.join(_HERE, "..", "results"))
WALLION = os.path.join(RESULTS, "wallion")

NUCS = ["K40", "Ra226", "Pb214", "Bi214", "Pb212", "Ac228", "Bi212", "Tl208"]
SRC_TAG = "m1"
```

## Константы (комментарий рядом — дословно, это ЖЁСТКИЙ ЛИМИТ этой сборки Geant4)

```python
# ЖЁСТКИЙ ЛИМИТ /gps/hist/point в этой сборке Geant4 (проверено фактом 21.08,
# бинарным перебором на живом прогоне: 1024 точки -> код 0, 1025 -> STATUS_
# STACK_BUFFER_OVERRUN). Похоже на фиксированный массив где-то в GPS Arb
# этой версии G4, не в нашем коде. MAX_HIST_POINTS - с запасом от границы.
MAX_HIST_POINTS = 900

# Родная сетка wallfield.cc: 2 кэВ бины по всему диапазону 0-3000 кэВ
# (kBinKeV в geometry/wallfield.cc). rebin_for_gps огрубляет её выше
# SPLIT_KEV до COARSE_KEV, чтобы уложиться в MAX_HIST_POINTS.
FINE_KEV = 2.0
SPLIT_KEV = 200.0
COARSE_KEV = 10.0
```

## Функция `rebin_for_gps` — БЕЗ ИЗМЕНЕНИЙ, переписать дословно как есть

```python
def rebin_for_gps(energies_kev, fluences, split_kev=SPLIT_KEV, coarse_kev=COARSE_KEV):
    """Переменное разрешение: точная 2-кэВ сетка НИЖЕ split_kev (там K-X
    комплекс, ради которого сетку мельчили), грубые окна ВЫШЕ. Сумма внутри
    окна - интеграл сохраняется точно; крупность там, где детектор с ПШПВ
    в десятки кэВ её всё равно не различит."""
    fine_e, fine_f, bucket = [], [], {}
    for e, f in zip(energies_kev, fluences):
        if f <= 0:
            continue
        if e < split_kev:
            fine_e.append(e)
            fine_f.append(f)
        else:
            w = int(e // coarse_kev)
            bucket[w] = bucket.get(w, 0.0) + f
    keys = sorted(bucket)
    return fine_e + [(w + 0.5) * coarse_kev for w in keys], fine_f + [bucket[w] for w in keys]
```

## НОВАЯ функция `add_walls(energies_kev, fluences, split_kev=SPLIT_KEV, fine_kev=FINE_KEV, coarse_kev=COARSE_KEV)`

Docstring (по-русски): «Вставляет нулевые точки ВПЛОТНУЮ к родной сетке на
каждом разрыве между соседними ненулевыми точками (родной шаг — fine_kev
ниже split_kev, coarse_kev выше) и на правом краю последней точки — иначе
GPS Arb+Lin трактует разрыв как широкий треугольник плотности, а последнюю
точку — как жёсткий обрыв распределения (см. докстрока модуля, D-002/D-003).
Список energies_kev предполагается уже отсортированным по возрастанию.»

Тело (реализуй ТОЧНО эту логику):

1. Если `energies_kev` пуст — вернуть `(list(energies_kev), list(fluences))` как есть.
2. Завести `out_e = []`, `out_f = []`.
3. Пройти по индексам `i` от `0` до `len(energies_kev)-1`:
   - `e = energies_kev[i]`, `f = fluences[i]`.
   - Определить родной шаг ДЛЯ ЭТОЙ точки: `step = fine_kev if e < split_kev else coarse_kev`.
   - **Левый край**: если `i == 0` — если `e - step > 0` (не у самого нуля шкалы),
     добавить в `out_e`/`out_f` стену `(e - step, 0.0)` ПЕРЕД текущей точкой.
     Если `i > 0` — если `e - energies_kev[i-1] > step + 1e-9` (разрыв больше
     родного шага; допуск `1e-9` — числа с плавающей точкой), добавить стену
     `(e - step, 0.0)` ПЕРЕД текущей точкой, но ТОЛЬКО если `e - step` строго
     больше `energies_kev[i-1]` (чтобы не вставить точку левее предыдущей и
     не нарушить возрастающий порядок; если условие не выполняется — стену
     не добавлять, разрыв слишком мал для отдельной стены).
   - Добавить саму точку: `out_e.append(e)`, `out_f.append(f)`.
   - **Правый край**: если `i == len(energies_kev) - 1` (это ПОСЛЕДНЯЯ точка
     во всём списке) — ВСЕГДА добавить стену `(e + step, 0.0)` ПОСЛЕ неё
     (последняя точка гистограммы GPS обрывает распределение — без стены
     энергия выше неё не сэмплируется никогда, см. D-002/D-003).
     Если `i < len(energies_kev) - 1` — если
     `energies_kev[i+1] - e > step + 1e-9`, добавить стену `(e + step, 0.0)`
     ПОСЛЕ текущей точки, ТОЛЬКО если `e + step` строго меньше
     `energies_kev[i+1]` (та же защита от нарушения порядка).
4. Вернуть `(out_e, out_f)`.

## Функция `main()` — как в существующей версии, ОДНО изменение: вызов `add_walls`

Структура ИДЕНТИЧНА существующей версии (чтение CSV, гейт на `MAX_HIST_POINTS`
до и после `rebin_for_gps`, запись `field_spectrum_m1_%s.mac`, копирование
run-макроса из архива, печать итоговой таблицы, `sys.exit`) — приведи её
полностью, но со следующей ОДНОЙ вставкой: **сразу после блока ребиннинга
(после `rebin_for_gps` и до записи макроса в файл) вызвать**:

```python
n_before_walls = len(energies)
energies, fluences = add_walls(energies, fluences)
n_after_walls = len(energies)
if n_after_walls > MAX_HIST_POINTS:
    print(u"ОШИБКА: add_walls добавила стены, и точек стало %d > %d для %s "
          u"(лимит /gps/hist/point). Нужно расширить бюджет ребиннинга."
          % (n_after_walls, MAX_HIST_POINTS, nuc))
    sys.exit(3)
if n_after_walls != n_before_walls:
    print(u"  %s: точек %d -> %d (add_walls: стены вокруг резких линий/разрывов)"
          % (nuc, n_before_walls, n_after_walls))
```

Всё остальное (чтение CSV, гейт `n_before > MAX_HIST_POINTS` на исходных
точках ДО ребиннинга, запись `/gps/hist/point`, копирование run-макроса,
итоговая таблица, коды возврата) — переписать ДОСЛОВНО так же, как в
базовой версии ниже (структура и тексты сообщений сохраняются один в один,
меняется только вставка `add_walls` в указанном месте):

```python
def main():
    written_pairs = 0
    summary = []

    bg_dir = os.path.join(RESULTS, "bare", "background")
    os.makedirs(bg_dir, exist_ok=True)

    for nuc in NUCS:
        csv_path = os.path.join(WALLION, "wf_%s_%s.csv" % (SRC_TAG, nuc))

        if not os.path.exists(csv_path):
            print(u"ПРЕДУПРЕЖДЕНИЕ: Отсутствует файл %s" % csv_path)
            continue

        energies = []
        fluences = []
        with open(csv_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split(",")
                if len(parts) != 2:
                    continue
                try:
                    e_kev = float(parts[0])
                    fluence = float(parts[1])
                    energies.append(e_kev)
                    fluences.append(fluence)
                except ValueError:
                    continue

        total_fluence = sum(fluences)
        if total_fluence <= 0:
            print(u"ОШИБКА: Общая флюенция для %s равна нулю или меньше" % nuc)
            sys.exit(3)

        n_before = sum(1 for f in fluences if f > 0)
        if n_before > MAX_HIST_POINTS:
            energies, fluences = rebin_for_gps(energies, fluences)
            drift = abs(sum(fluences) - total_fluence) / total_fluence
            if drift > 1e-9:
                print(u"ОШИБКА: ребиннинг исказил сумму потока для %s (дрейф %.2e)"
                      % (nuc, drift))
                sys.exit(3)
            print(u"  %s: точек %d -> %d (ребиннинг >200 кэВ до 10 кэВ, сумма потока сохранена)"
                  % (nuc, n_before, sum(1 for f in fluences if f > 0)))
        if sum(1 for f in fluences if f > 0) > MAX_HIST_POINTS:
            print(u"ОШИБКА: после ребиннинга всё ещё %d точек > %d для %s"
                  % (sum(1 for f in fluences if f > 0), MAX_HIST_POINTS, nuc))
            sys.exit(3)

        # (здесь вставка add_walls, см. блок выше)

        spectrum_path = os.path.join(RESULTS, "field_spectrum_m1_%s.mac" % nuc)
        with open(spectrum_path, "w", encoding="utf-8") as f:
            f.write(u"# Единичный (1 Бк/кг) отклик звена %s, wallfield.exe.\n" % nuc)
            f.write(u"# МЕТОД 1: полный распад ТОЛЬКО этого звена, nucleusLimits\n")
            f.write(u"# отсекает дочерние (канон geant4-spectrum-pipeline, D-001).\n")
            f.write(u"# ВЕРСИЯ 2 (21.08): add_walls вокруг резких линий/разрывов (D-002/D-003).\n")
            f.write(u"# FLUENCE_TOTAL_CM2_S = %.6f\n" % total_fluence)
            stamp = os.environ.get("G4MODELS_STAMP", "")
            if stamp:
                f.write(u"# Сгенерировано %s\n" % stamp)
            f.write(u"/gps/particle gamma\n")
            f.write(u"/gps/ene/type Arb\n")
            f.write(u"/gps/hist/type arb\n")
            for e_kev, fluence in zip(energies, fluences):
                e_meV = e_kev / 1000.0
                f.write(u"/gps/hist/point %.6f %.6e\n" % (e_meV, fluence))
            f.write(u"/gps/hist/inter Lin\n")

        run_source_path = os.path.join(BUILD, "_attic_table_method_20260821",
                                       "field_run_nucb_%s.mac" % nuc)
        if not os.path.exists(run_source_path):
            print(u"ОШИБКА: Не найден исходный файл для run-макроса: %s" % run_source_path)
            sys.exit(3)

        run_dest_path = os.path.join(BUILD, "field_run_m1b_%s.mac" % nuc)
        with open(run_source_path, "r", encoding="utf-8") as src_f:
            with open(run_dest_path, "w", encoding="utf-8") as dst_f:
                for line in src_f:
                    if line.strip().startswith("/control/execute"):
                        abs_spectrum = os.path.abspath(spectrum_path).replace("\\", "/")
                        dst_f.write(u"/control/execute %s\n" % abs_spectrum)
                    elif line.strip().startswith("/rc/outFile"):
                        out_file = os.path.join(RESULTS, "bare", "background",
                                                "bg_bare_field_m1_%s.csv" % nuc).replace("\\", "/")
                        dst_f.write(u"/rc/outFile %s\n" % out_file)
                    else:
                        dst_f.write(line)

        hist_points = len(energies)
        summary.append((nuc, total_fluence, hist_points, spectrum_path, run_dest_path))
        written_pairs += 1

    print(u"Сгенерировано пар макросов: %d из %d" % (written_pairs, len(NUCS)))
    if written_pairs == len(NUCS):
        print(u"Все макросы успешно созданы.")
        print(u"Итоговая таблица:")
        print(u"Нуклид\t\tОбщая флюенция\tТочек гистограммы\tПуть к спектру\t\t\t\t\t\tПуть к run-макросу")
        for nuc, fluence, points, spec_path, run_path in summary:
            print(u"%s\t\t%.6f\t\t%d\t\t%s\t\t%s" % (nuc, fluence, points, spec_path, run_path))
        sys.exit(0)
    else:
        missing = len(NUCS) - written_pairs
        print(u"ОШИБКА: Пропущено %d макросов." % missing)
        sys.exit(1)


if __name__ == "__main__":
    main()
```

ВАЖНО: точка вставки `add_walls` — СРАЗУ ПОСЛЕ блока проверки
`if sum(1 for f in fluences if f > 0) > MAX_HIST_POINTS: ... sys.exit(3)`
и ПЕРЕД `spectrum_path = os.path.join(...)`. Используй код вставки из
раздела «Функция main() — ... вызов add_walls» выше буквально на этом месте.
