# -*- coding: utf-8 -*-
"""Готовит пары макросов (спектр-источник + run) для МЕТОДА 1: единичный
(1 Бк/кг) полный распад ОДНОГО звена, nucleusLimits отсекает дочерние.

ВЕРСИЯ 2 (21.08, после /compact): добавлена add_walls() — контрольный прогон
подтвердил, что Geant4 GPS (/gps/hist/type arb + /gps/hist/inter Lin) трактует
точки (E,weight) как узлы НЕПРЕРЫВНОЙ плотности (трапецеидальная интерполяция
между соседними точками по их энергетическому РАССТОЯНИЮ), а не как «содержимое
бина» — резкая линия среди разреженного континуума размывалась в широкий
треугольник и (для последней точки гистограммы) обрывалась справа. Проверено
фактом на K-40 1460,8 кэВ: без стен модель давала 96,4 Бк/кг, со стенами
вплотную (1,459/1,463 МэВ) — 131,1, против 122,1 у независимого моноотклика
(DECISIONS.md D-002/D-003). add_walls() вставляет нулевые точки вплотную к
родной сетке на каждом разрыве и на правом краю последней точки.

Запуск: python make_m1_macros_v2.py
"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(_HERE, "..", "..", ".."))
BUILD = os.path.join(REPO, "build", "RadiaCode-103")
RESULTS = os.path.normpath(os.path.join(_HERE, "..", "results"))
WALLION = os.path.join(RESULTS, "wallion")

NUCS = ["K40", "Ra226", "Pb214", "Bi214", "Pb212", "Ac228", "Bi212", "Tl208"]
SRC_TAG = "m1"

# ЖЁСТКИЙ ЛИМИТ /gps/hist/point в этой сборке Geant4 (проверено фактом 21.08,
# бинарным перебором на живом прогоне: 1024 точки -> код 0, 1025 -> STATUS_
# STACK_BUFFER_OVERRUN). Похоже на фиксированный массив где-то в GPS Arb
# этой версии G4, не в нашем коде. MAX_HIST_POINTS - с запасом от границы.
MAX_HIST_POINTS = 950  # 21.08: 900 не хватило Ac228 после add_walls (926) — подняли
# запас; граница 1024/1025 по-прежнему проверена бинарным перебором, 950 оставляет
# ещё 74 точки запаса.

# Родная сетка wallfield.cc: 2 кэВ бины по всему диапазону 0-3000 кэВ
# (kBinKeV в geometry/wallfield.cc). rebin_for_gps огрубляет её выше
# SPLIT_KEV до COARSE_KEV, чтобы уложиться в MAX_HIST_POINTS.
FINE_KEV = 2.0
SPLIT_KEV = 200.0
COARSE_KEV = 10.0

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

def add_walls(energies_kev, fluences, rebinned, split_kev=SPLIT_KEV, fine_kev=FINE_KEV, coarse_kev=COARSE_KEV):
    """Вставляет нулевые точки ВПЛОТНУЮ к родной сетке на каждом разрыве между соседними
    ненулевыми точками (родной шаг — fine_kev ниже split_kev; выше — coarse_kev, но ТОЛЬКО
    если rebin_for_gps реально огрублял сетку (rebinned=True), иначе fine_kev везде — у K40
    (640 точек < MAX_HIST_POINTS) ребиннинг не срабатывал вовсе, и весь спектр, включая
    линию 1460,8, остался на родной 2-кэВ сетке; жёсткое предположение coarse_kev=10 там
    маскировало 6-кэВ разрыв как "недостаточный" для стены — найдено при вычитке 21.08) и на
    правом краю последней точки — иначе GPS Arb+Lin трактует разрыв как широкий треугольник
    плотности, а последнюю точку — как жёсткий обрыв распределения (D-002/D-003).
    Список energies_kev предполагается уже отсортированным по возрастанию."""
    if not energies_kev:
        return list(energies_kev), list(fluences)
    
    out_e = []
    out_f = []
    
    for i in range(len(energies_kev)):
        e = energies_kev[i]
        f = fluences[i]
        
        # Определить родной шаг ДЛЯ ЭТОЙ точки (coarse только если ребиннинг реально был)
        step = fine_kev if (e < split_kev or not rebinned) else coarse_kev
        
        # Левый край
        if i == 0:
            # Если не у самого нуля шкалы, добавить стену перед текущей точкой
            if e - step > 0:
                out_e.append(e - step)
                out_f.append(0.0)
        else:
            # Если разрыв больше родного шага, добавить стену перед текущей точкой
            if e - energies_kev[i-1] > step + 1e-9:
                if e - step > energies_kev[i-1]:
                    out_e.append(e - step)
                    out_f.append(0.0)
        
        # Добавить саму точку
        out_e.append(e)
        out_f.append(f)
        
        # Правый край
        if i == len(energies_kev) - 1:
            # Всегда добавить стену после последней точки
            out_e.append(e + step)
            out_f.append(0.0)
        else:
            # Если разрыв больше родного шага, добавить стену после текущей точки
            if energies_kev[i+1] - e > step + 1e-9:
                if e + step < energies_kev[i+1]:
                    out_e.append(e + step)
                    out_f.append(0.0)
    
    return out_e, out_f

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
        was_rebinned = n_before > MAX_HIST_POINTS
        if was_rebinned:
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

        # Вставка add_walls
        n_before_walls = len(energies)
        energies, fluences = add_walls(energies, fluences, was_rebinned)
        n_after_walls = len(energies)
        if n_after_walls > MAX_HIST_POINTS:
            print(u"ОШИБКА: add_walls добавила стены, и точек стало %d > %d для %s "
                  u"(лимит /gps/hist/point). Нужно расширить бюджет ребиннинга."
                  % (n_after_walls, MAX_HIST_POINTS, nuc))
            sys.exit(3)
        if n_after_walls != n_before_walls:
            print(u"  %s: точек %d -> %d (add_walls: стены вокруг резких линий/разрывов)"
                  % (nuc, n_before_walls, n_after_walls))

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
