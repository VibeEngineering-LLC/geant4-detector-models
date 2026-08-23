# -*- coding: utf-8 -*-
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(_HERE, "..", "..", ".."))
BUILD = os.path.join(REPO, "build", "RadiaCode-103")
RESULTS = os.path.normpath(os.path.join(_HERE, "..", "results"))
WALLION = os.path.join(RESULTS, "wallion")

NUCS = ["K40", "Ra226", "Pb214", "Bi214", "Pb212", "Ac228", "Bi212", "Tl208"]
SRC_TAG = "m1"     # wf_<SRC_TAG>_<nuc>.csv from the 2-keV ion runs

# ЖЁСТКИЙ ЛИМИТ /gps/hist/point в этой сборке Geant4 (проверено фактом 21.08,
# бинарным перебором на живом прогоне: 1024 точки -> код 0, 1025 -> STATUS_
# STACK_BUFFER_OVERRUN). Похоже на фиксированный массив где-то в GPS Arb
# этой версии G4, не в нашем коде. MAX_HIST_POINTS - с запасом от границы.
MAX_HIST_POINTS = 900

def rebin_for_gps(energies_kev, fluences, split_kev=200.0, coarse_kev=10.0):
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


def main():
    written_pairs = 0
    summary = []
    
    # Убедиться, что директория RESULTS/bare/background существует
    bg_dir = os.path.join(RESULTS, "bare", "background")
    os.makedirs(bg_dir, exist_ok=True)
    
    for nuc in NUCS:
        csv_path = os.path.join(WALLION, "wf_%s_%s.csv" % (SRC_TAG, nuc))
        
        if not os.path.exists(csv_path):
            print(u"ПРЕДУПРЕЖДЕНИЕ: Отсутствует файл %s" % csv_path)
            continue
            
        # Читаем данные из CSV
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

        # ГЕЙТ на 1024-лимит (см. MAX_HIST_POINTS выше): не пишем макрос,
        # который упадёт в Geant4 - сначала ужимаем переменным разрешением.
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

        # Записываем SOURCE SPECTRUM макрос
        spectrum_path = os.path.join(RESULTS, "field_spectrum_m1_%s.mac" % nuc)
        with open(spectrum_path, "w", encoding="utf-8") as f:
            f.write(u"# Единичный (1 Бк/кг) отклик звена %s, wallfield.exe.\n" % nuc)
            f.write(u"# МЕТОД 1: полный распад ТОЛЬКО этого звена, nucleusLimits\n")
            f.write(u"# отсекает дочерние (канон geant4-spectrum-pipeline, D-001).\n")
            f.write(u"# FLUENCE_TOTAL_CM2_S = %.6f\n" % total_fluence)
            stamp = os.environ.get("G4MODELS_STAMP", "")
            if stamp:
                f.write(u"# Сгенерировано %s\n" % stamp)
            f.write(u"/gps/particle gamma\n")
            f.write(u"/gps/ene/type Arb\n")
            f.write(u"/gps/hist/type arb\n")
            for i, (e_kev, fluence) in enumerate(zip(energies, fluences)):
                if fluence > 0:
                    e_meV = e_kev / 1000.0
                    f.write(u"/gps/hist/point %.6f %.6e\n" % (e_meV, fluence))
            f.write(u"/gps/hist/inter Lin\n")
            
        # Записываем RUN макрос
        # Шаблон берётся из архива старого метода: сами его прогоны отменены,
        # но ГЕОМЕТРИЯ источника (цилиндр полости) от метода не зависит и
        # обязана совпадать с той, на которой считался wallfield. Выдумывать
        # её здесь нельзя — молча разошлась бы с флюенсом.
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
                        # Путь с прямыми слэшами — как во всех существующих макросах
                        abs_spectrum = os.path.abspath(spectrum_path).replace("\\", "/")
                        dst_f.write(u"/control/execute %s\n" % abs_spectrum)
                    elif line.strip().startswith("/rc/outFile"):
                        out_file = os.path.join(RESULTS, "bare", "background",
                                                "bg_bare_field_m1_%s.csv" % nuc).replace("\\", "/")
                        dst_f.write(u"/rc/outFile %s\n" % out_file)
                    else:
                        dst_f.write(line)
                        
        # Собираем информацию для таблицы
        hist_points = sum(1 for f in fluences if f > 0)
        summary.append((nuc, total_fluence, hist_points, spectrum_path, run_dest_path))
        written_pairs += 1
        
    # Выводим итоговую таблицу
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
