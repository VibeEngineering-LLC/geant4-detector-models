# -*- coding: utf-8 -*-
"""Генерирует Geant4 GPS макросы для METHOD 1 (полный распад цепочки одного нуклида).
Переход с /gps/hist/type arb + /gps/hist/inter Lin на /gps/ene/type User + /gps/hist/type energy.
См. DECISIONS.md D-003 за подробностями о нестабильности Arb+Lin и необходимости ступенчатой гистограммы."""
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

MAX_HIST_POINTS = 950  # безопасный предел, оставлен для совместимости с v2
NATIVE_STEP_KEV = 2.0  # шаг бина в wallfield.cc


def read_wallfield_csv(csv_path):
    # header_fluence — сумма из собственного заголовка wallfield.cc
    # (# fluence_total_cm2_s = ...), раньше не читалась вообще, хотя уже
    # лежит в файле — дешёвая сверка целостности против усечённого
    # прогона (найдено циклом 4/5 стерильного аудита 22.08). None, если
    # заголовка нет (старый формат файла).
    energies = []
    fluences = []
    header_fluence = None
    with open(csv_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                if "fluence_total_cm2_s" in line and "=" in line:
                    try:
                        header_fluence = float(line.split("=", 1)[1].strip())
                    except ValueError:
                        pass
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
    return energies, fluences, header_fluence


def build_full_grid(energies_kev, fluences, step_kev):
    min_e = min(energies_kev)
    max_e = max(energies_kev)
    n_bins = int(round((max_e - min_e) / step_kev)) + 1
    grid_energies = [min_e + i * step_kev for i in range(n_bins)]
    grid_fluences = [0.0] * n_bins

    # Используем словарь для точного сопоставления по округленной энергии.
    # НАКОПЛЕНИЕ (+=), не перезапись — если у двух входных строк энергия
    # совпадёт после округления до 3 знаков, вторая раньше молча стирала
    # первую и теряла флюенс без единого предупреждения (найдено циклом
    # 3/5 стерильного аудита 22.08, воспроизведено синтетически; на
    # текущих 8 production-CSV дублей нет, риск был теоретическим).
    energy_dict = {}
    for e, f in zip(energies_kev, fluences):
        key = round(e, 3)
        energy_dict[key] = energy_dict.get(key, 0.0) + f

    for i, e in enumerate(grid_energies):
        key = round(e, 3)
        if key in energy_dict:
            grid_fluences[i] = energy_dict[key]

    return grid_energies, grid_fluences


def coarsen_grid(grid_energies_kev, grid_fluences, native_step_kev, max_points):
    # Возвращает ВЕРХНИЕ ГРАНИЦЫ бинов (не центр+общий шаг) — так граница
    # последнего, возможно НЕПОЛНОГО блока (n_before % k != 0) считается
    # верно сама по себе, без предположения о его ширине. Баг найден
    # стерильным аудитом 22.08: единый used_step_kev на все бины давал
    # разрыв/нахлёст ровно на последней паре точек Bi214/Tl208 (n_before=1501
    # не делится на k=2 нацело) — .../3.000000 -> 3.003000 вместо шага 4 кэВ.
    n_before = len(grid_energies_kev)
    if n_before <= max_points - 1:
        upper_edges = [e + native_step_kev / 2.0 for e in grid_energies_kev]
        return upper_edges, grid_fluences

    # Индексное блочное суммирование (НЕ сравнение float-границ — округление
    # там теряло крайние точки, найдено на живом прогоне Bi214 22.08): найти
    # минимальное k (число родных бинов на один укрупнённый), затем разбить
    # ОТСОРТИРОВАННЫЙ список на подряд идущие блоки по k точек и просуммировать
    # каждый блок. Каждая fine-точка попадает РОВНО в один блок — сумма
    # гарантированно сохраняется, последний блок может быть короче k.
    for k in range(2, 10000):
        if (n_before + k - 1) // k <= max_points - 1:
            break
    else:
        raise RuntimeError("Не удалось найти подходящий коэффициент ребиннинга")

    coarse_upper_edges = []
    coarse_fluences = []
    for start in range(0, n_before, k):
        block_e = grid_energies_kev[start:start + k]
        block_f = grid_fluences[start:start + k]
        # Верхняя граница = правый край ПОСЛЕДНЕЙ fine-точки блока — верно
        # для любого размера блока, включая укороченный последний.
        coarse_upper_edges.append(block_e[-1] + native_step_kev / 2.0)
        coarse_fluences.append(sum(block_f))

    # Проверка сохранения суммы. RuntimeError (не sys.exit!) — main() ловит
    # именно этот тип для всей функции; sys.exit(3) здесь раньше давал
    # необработанный SystemExit (except RuntimeError его не перехватывает),
    # обрывая весь batch вопреки заявленному "единому мягкому стилю" —
    # найдено циклом 3/5 стерильного аудита 22.08. sum_orig<=0 защищено
    # отдельно от ZeroDivisionError (после фикса накопления в
    # build_full_grid недостижимо в реальном потоке, но самостоятельная
    # функция не должна на это полагаться).
    sum_orig = sum(grid_fluences)
    sum_coarse = sum(coarse_fluences)
    if sum_orig <= 0:
        raise RuntimeError("сумма исходных флюенций <= 0, ребиннинг невозможен")
    if abs(sum_coarse - sum_orig) / sum_orig >= 1e-9:
        raise RuntimeError("сумма флюенций не сохранена при ребиннинге")

    return coarse_upper_edges, coarse_fluences


def write_energy_macro(spectrum_path, upper_edges_kev, coarse_fluences,
                        lower_edge_first_bin_kev, total_fluence, nuc):
    with open(spectrum_path, "w", encoding="utf-8") as f:
        f.write(u"# Единичный (1 Бк/кг) отклик звена %s, wallfield.exe.\n" % nuc)
        f.write(u"# МЕТОД 1: полный распад ТОЛЬКО этого звена, nucleusLimits\n")
        f.write(u"# отсекает дочерние (канон geant4-spectrum-pipeline, D-001).\n")
        f.write(u"# ВЕРСИЯ 3 (22.08): /gps/hist/type energy (ступенчатая гистограмма), НЕ arb+Lin —\n")
        f.write(u"# см. DECISIONS.md D-003 (нестабильность Arb+Lin, дефолт-откат/NaN на реальных прогонах).\n")
        f.write(u"# FLUENCE_TOTAL_CM2_S = %.6e\n" % total_fluence)
        f.write("/gps/particle gamma\n")
        f.write("/gps/ene/type User\n")
        f.write("/gps/hist/type energy\n")

        # Начальный нулевой бин (нижняя граница ПЕРВОЙ точки родной сетки —
        # первый блок при ребиннинге всегда полный, эта граница от него не зависит)
        lower_edge = max(lower_edge_first_bin_kev, 0.0)
        f.write("/gps/hist/point %.6f 0.0\n" % (lower_edge / 1000.0))

        # Остальные бины — верхние границы уже готовы (см. coarsen_grid)
        for upper_edge, weight in zip(upper_edges_kev, coarse_fluences):
            f.write("/gps/hist/point %.6f %.6e\n" % (upper_edge / 1000.0, weight))


def main():
    written_pairs = 0
    summary = []

    for nuc in NUCS:
        csv_path = os.path.join(WALLION, "wf_%s_%s.csv" % (SRC_TAG, nuc))
        if not os.path.exists(csv_path):
            print(u"ПРЕДУПРЕЖДЕНИЕ: Отсутствует файл %s" % csv_path)
            continue

        # Ниже — единый мягкий стиль для ВСЕХ проблем ОДНОГО нуклида (не
        # только отсутствия файла-шаблона): предупреждение + continue, а не
        # sys.exit(3), обрывающий необработанные нуклиды дальше по NUCS.
        # Найдено циклом 2/5 стерильного аудита 22.08 — асимметрия осталась
        # в трёх местах после первого фикса (только для шаблона).
        energies, fluences, header_fluence = read_wallfield_csv(csv_path)
        total_fluence = sum(fluences)
        if total_fluence <= 0:
            print(u"ПРЕДУПРЕЖДЕНИЕ: Общая флюенция для %s равна нулю или меньше" % nuc)
            continue

        # Сверка с собственным заголовком wallfield.cc — дешёвая защита от
        # усечённого (частично записанного) прогона; найдено циклом 4/5.
        # header_fluence<=0 (не только None) — ТОЖЕ подозрительно (битый
        # заголовок), а не молчаливый пропуск проверки: найдено циклом 5/5.
        if header_fluence is not None:
            if header_fluence <= 0:
                print(u"ПРЕДУПРЕЖДЕНИЕ: %s — заголовок CSV содержит "
                      u"некорректную сумму (%.6e), файл повреждён" %
                      (nuc, header_fluence))
                continue
            rel = abs(total_fluence - header_fluence) / header_fluence
            if rel >= 1e-6:
                print(u"ПРЕДУПРЕЖДЕНИЕ: %s — сумма из данных не совпадает "
                      u"с заголовком CSV, файл похож на усечённый" % nuc)
                continue

        grid_energies, grid_fluences = build_full_grid(energies, fluences, NATIVE_STEP_KEV)
        # Защита от молчаливой потери флюенса при энергии, не лежащей ровно
        # на решётке NATIVE_STEP_KEV (round(e,3)-ключ не найдёт совпадения
        # в grid_energies) — найдено циклом 4/5 стерильного аудита 22.08.
        grid_sum = sum(grid_fluences)
        if abs(grid_sum - total_fluence) / total_fluence >= 1e-9:
            print(u"ПРЕДУПРЕЖДЕНИЕ: %s — при укладке на сетку %.1f кэВ потерян "
                  u"флюенс — есть энергии вне решётки" % (nuc, NATIVE_STEP_KEV))
            continue
        try:
            coarse_upper_edges, coarse_fluences = coarsen_grid(
                grid_energies, grid_fluences, NATIVE_STEP_KEV, MAX_HIST_POINTS)
        except RuntimeError as exc:
            print(u"ПРЕДУПРЕЖДЕНИЕ: ребиннинг %s не удался (%s)" % (nuc, exc))
            continue

        n_before = len(grid_energies)
        n_after = len(coarse_upper_edges)
        if n_before > n_after:
            print(u"  %s: точек %d -> %d (ребиннинг, шаг ~%.1f кэВ, сумма потока сохранена)" %
                  (nuc, n_before, n_after, NATIVE_STEP_KEV * n_before / n_after))

        if len(coarse_upper_edges) > MAX_HIST_POINTS - 1:
            print(u"ПРЕДУПРЕЖДЕНИЕ: после ребиннинга для %s точек больше допустимого %d" %
                  (nuc, MAX_HIST_POINTS))
            continue

        spectrum_path = os.path.join(RESULTS, "field_spectrum_m1_%s.mac" % nuc)
        lower_edge_first_bin = grid_energies[0] - NATIVE_STEP_KEV / 2.0
        write_energy_macro(spectrum_path, coarse_upper_edges, coarse_fluences,
                           lower_edge_first_bin, total_fluence, nuc)

        # Генерация run-макроса. Отсутствие шаблона — как отсутствие CSV
        # (пропуск нуклида, не обрыв всего batch): раньше был sys.exit(3),
        # что молча теряло ещё не обработанные нуклиды в списке NUCS —
        # асимметрия с обработкой отсутствующего CSV, найдено аудитом 22.08.
        template_path = os.path.join(BUILD, "_attic_table_method_20260821", "field_run_nucb_%s.mac" % nuc)
        if not os.path.exists(template_path):
            print(u"ПРЕДУПРЕЖДЕНИЕ: Отсутствует шаблон %s" % template_path)
            continue

        run_path = os.path.join(BUILD, "field_run_m1b_%s.mac" % nuc)
        with open(template_path, "r", encoding="utf-8") as src_f:
            with open(run_path, "w", encoding="utf-8") as dst_f:
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

        summary.append((nuc, total_fluence, len(coarse_upper_edges), spectrum_path, run_path))
        written_pairs += 1

    print(u"Сгенерировано пар макросов: %d из %d" % (written_pairs, len(NUCS)))

    # Таблица печатается по накопленному summary ВСЕГДА (не только при
    # written_pairs==len(NUCS)) — иначе частичный сбой (1 из 8 нуклидов
    # пропущен) прятал диагностику по остальным 7 успешным. Найдено
    # циклом 2/5 стерильного аудита 22.08.
    if summary:
        print(u"  Звено     Флюенция      Точек   Путь к спектру")
        print(u"  ---------------------------------------------------")
        for nuc, fluence, points, spec_path, run_path in summary:
            print(u"  %-8s  %10.3e  %6d  %s" % (nuc, fluence, points, os.path.basename(spec_path)))

    if written_pairs == len(NUCS):
        sys.exit(0)
    else:
        missing = len(NUCS) - written_pairs
        print(u"ОШИБКА: Пропущено %d пар макросов" % missing)
        sys.exit(1)


if __name__ == "__main__":
    main()
