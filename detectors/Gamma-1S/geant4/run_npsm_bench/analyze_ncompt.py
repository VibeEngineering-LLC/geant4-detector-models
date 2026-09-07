import csv
import math
import glob
import argparse
import pathlib
import sys
import json

sys.stdout.reconfigure(encoding="utf-8")

def parse_csv(file_path):
    data = {}
    table = []
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                if "," in line:
                    key, value = line.strip().split(",", 1)
                    # Приводим к числу ОДИН РАЗ при чтении — все потребители
                    # ниже получают уже числа. Точечная конвертация по списку
                    # ключей оставляла mean/sem строками и роняла критерии K1/K2.
                    if key in ("em_deex", "crystal_mm"):
                        data[key] = value            # заведомо строковые поля
                    else:
                        try:
                            data[key] = int(value)
                        except ValueError:
                            try:
                                data[key] = float(value)
                            except ValueError:
                                data[key] = value
                elif line.strip() == "n_compt,count_absorbed,count_all":
                    continue
                elif line.strip().startswith("0,"):
                    break
            for line in f:
                if "," in line and not line.startswith("#"):
                    parts = line.strip().split(",")
                    if len(parts) >= 3:
                        table.append([int(parts[0]), int(parts[1]), int(parts[2])])
    except Exception as e:
        print(f"ПРЕДУПРЕЖДЕНИЕ: Не удалось прочитать файл {file_path}: {e}")
        return None
    if "energy_keV" not in data or "mean_ncompt_absorbed" not in data:
        print(f"ПРЕДУПРЕЖДЕНИЕ: В файле {file_path} отсутствуют необходимые поля")
        return None
    data["table"] = table
    return data

def weighted_mean(table):
    total_weighted = 0
    total_weight = 0
    for n_compt, count_absorbed, _ in table:
        total_weighted += n_compt * count_absorbed
        total_weight += count_absorbed
    if total_weight == 0:
        return float('nan')
    return total_weighted / total_weight

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--glob", default="out/ncompt_*.csv")
    parser.add_argument("--ref-lo", nargs=3, type=float, default=[500, 1.4, 0.15])
    parser.add_argument("--ref-hi", nargs=3, type=float, default=[2000, 2.3, 0.20])
    parser.add_argument("--json", type=str)
    args = parser.parse_args()

    files = glob.glob(args.glob)
    if not files:
        print("ПРЕДУПРЕЖДЕНИЕ: Не найдено файлов по шаблону")
        sys.exit(2)

    results = []
    for f in files:
        data = parse_csv(f)
        if data is None:
            continue
        results.append(data)

    if len(results) < 2:
        print("ПРЕДУПРЕЖДЕНИЕ: Загружено меньше двух файлов")
        sys.exit(2)

    # Ключи словаря — строки из CSV; сортировать обязательно по ЧИСЛУ, иначе
    # порядок лексикографический ("1100" < "500") и проверка монотонности ложна.
    results.sort(key=lambda x: float(x["energy_keV"]))

    print("Энергия, N событий, N поглощенных, N полной энергии, N ушедших, среднее ± SEM, отношение")
    print("-----------------------------------------------------------------------------------------------")

    all_means = []
    all_sems = []
    all_energies = []

    for data in results:
        # Все поля шапки приходят строками — приводим к числам здесь, один раз.
        energy = float(data["energy_keV"])
        n_events = int(data["n_events_processed"])
        n_absorbed = int(data["n_absorbed_phot"])
        n_full_edep = int(data["n_full_edep_1keV"])
        n_escaped = int(data["n_escaped"])
        mean = float(data["mean_ncompt_absorbed"])
        sem = float(data["sem_ncompt_absorbed"])
        ratio = n_absorbed / n_full_edep if n_full_edep != 0 else float('nan')

        print(f"{energy:6.0f}, {n_events:8d}, {n_absorbed:7d}, {n_full_edep:7d}, {n_escaped:7d}, {mean:6.3f} ± {sem:5.3f}, {ratio:6.3f}")

        all_means.append(mean)
        all_sems.append(sem)
        all_energies.append(energy)

        # Проверка согласованности среднего
        w_mean = weighted_mean(data["table"])
        if not math.isnan(w_mean) and abs(mean - w_mean) > 1e-3:
            print(f"ПРОТИВОРЕЧИЕ В ФАЙЛЕ {data['energy_keV']}")

    # Отсутствие нужной точки — НЕ несовпадение: на неполной серии ближайший
    # файл может отстоять на сотни кэВ. Такой случай = N/A, в итог не идёт.
    def check_ref(tag, ref):
        e0, m0, tol = ref
        c = min(results, key=lambda x: abs(x["energy_keV"] - e0))
        de = abs(c["energy_keV"] - e0)
        if de > 0.05 * e0:
            print(f"{tag}: N/A — нет точки вблизи {e0:.0f} кэВ "
                  f"(ближайшая {c['energy_keV']:.0f})")
            return None
        d = abs(c["mean_ncompt_absorbed"] - m0)
        print(f"{tag}: {'PASS' if d <= tol else 'FAIL'} | {c['energy_keV']:.0f} кэВ: "
              f"{c['mean_ncompt_absorbed']:.4f} ± {c['sem_ncompt_absorbed']:.4f} "
              f"против {m0:.4f}, допуск {tol:.4f}, расхождение {d:.4f}")
        return d <= tol

    k1_pass = check_ref("K1 низкая энергия", args.ref_lo)
    k2_pass = check_ref("K2 высокая энергия", args.ref_hi)

    # Критерий K3
    k3_pass = True
    for i in range(1, len(all_means)):
        diff = all_means[i] - all_means[i-1]
        sem_sum = all_sems[i] + all_sems[i-1]
        if diff < -2 * sem_sum:
            print(f"K3 монотонность: FAIL")
            print(f"  Пара {all_energies[i-1]} и {all_energies[i]}: разница = {diff:.4f}, "
                  f"сумма SEM = {sem_sum:.4f}")
            k3_pass = False
            break

    if k3_pass:
        print("K3 монотонность: PASS")

    # Критерий K4
    k4_pass = True
    worst_ratio = 0.0
    worst_file = None
    for data in results:
        n_absorbed = data["n_absorbed_phot"]
        n_full_edep = data["n_full_edep_1keV"]
        if n_full_edep == 0:
            ratio = float('nan')
        else:
            ratio = n_absorbed / n_full_edep
        # nan не сравнивается ни с чем: файл с нулевым n_full_edep обязан
        # попасть в нарушители, иначе K4 слепа именно к пустому прогону.
        bad = math.isnan(ratio) or not (0.90 <= ratio <= 1.10)
        if bad:
            k4_pass = False
            if worst_file is None or math.isnan(ratio) or \
               abs(ratio - 1.0) > abs(worst_ratio - 1.0):
                worst_ratio = ratio
                worst_file = data

    if k4_pass:
        print("K4 согласованность процессов: PASS")
    else:
        print("K4 согласованность процессов: FAIL")
        if worst_file is not None:
            print(f"  Худший файл {worst_file['energy_keV']:.0f} кэВ: "
                  f"отношение = {worst_ratio:.4f} (норма 0.90..1.10)")

    # None = критерий неприменим (нет данных). Он не проваливает приёмку, но и
    # не позволяет объявить её пройденной: итог тогда НЕПОЛНЫЙ, код возврата 3.
    checks = [("K1", k1_pass), ("K2", k2_pass), ("K3", k3_pass), ("K4", k4_pass)]
    failed = [t for t, v in checks if v is False]
    skipped = [t for t, v in checks if v is None]
    all_pass = not failed and not skipped

    if failed:
        print(f"ИТОГ: FAIL ({', '.join(failed)})")
        exit_code = 1
    elif skipped:
        print(f"ИТОГ: НЕПОЛНЫЙ — нет данных для {', '.join(skipped)}; "
              f"остальные критерии пройдены")
        exit_code = 3
    else:
        print("ИТОГ: PASS")
        exit_code = 0

    # JSON
    if args.json:
        report = {
            "files": [
                {
                    "energy_keV": data["energy_keV"],
                    "n_events_processed": data["n_events_processed"],
                    "n_absorbed_phot": data["n_absorbed_phot"],
                    "n_full_edep_1keV": data["n_full_edep_1keV"],
                    "n_escaped": data["n_escaped"],
                    "mean_ncompt_absorbed": data["mean_ncompt_absorbed"],
                    "sem_ncompt_absorbed": data["sem_ncompt_absorbed"]
                } for data in results
            ],
            "criteria": {
                "K1_low_energy": k1_pass,
                "K2_high_energy": k2_pass,
                "K3_monotonicity": k3_pass,
                "K4_process_energy_consistency": k4_pass
            },
            "pass": all_pass
        }
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(report, f, ensure_ascii=False, indent=1)

    sys.exit(exit_code)

if __name__ == "__main__":
    main()
