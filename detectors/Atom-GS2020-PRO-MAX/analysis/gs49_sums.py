import csv
import json
import math
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")


def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    pairs_path = os.path.join(base_dir, "../audit/gs49_pairs.csv")
    library_path = os.path.join(base_dir, "data/gs49_library_lines.csv")
    output_path = os.path.join(base_dir, "configs", "gs49_xg_sums.json")

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Load library lines
    library_lines = []
    with open(library_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            library_lines.append(line.strip())

    reader = csv.DictReader(library_lines)
    lib_data = {}  # key: (nuclide, kind), value: list of e_keV
    for row in reader:
        nuclide = row["nuclide"]
        kind = row["kind"]
        e_kev = float(row["e_keV"])
        key = (nuclide, kind)
        if key not in lib_data:
            lib_data[key] = []
        lib_data[key].append(e_kev)

    # Load pairs
    with open(pairs_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        pairs_rows = list(reader)

    total_input = len(pairs_rows)
    dropped_x_unmapped = 0
    dropped_gamma_unmapped = 0
    dropped_x_below_50 = 0
    dropped_p_threshold = 0
    dropped_e_above_2900 = 0

    # Process pairs
    merged_data = {}  # key: (nuclide, e_x_lib, e_g_lib), value: {'p_pair': float, 'p_err_sq_sum': float}

    for row in pairs_rows:
        nuclide = row["nuclide"]
        e_x_kev = float(row["e_x_keV"])
        e_gamma_kev = float(row["e_gamma_keV"])
        p_corr = float(row["p_corr"])
        p_raw_err = float(row["p_raw_err"])

        if p_corr < 1e-4:
            continue

        # Map X-ray
        x_key = (nuclide, "x")
        if x_key not in lib_data:
            dropped_x_unmapped += 1
            continue
        x_lib_energies = lib_data[x_key]
        best_x_diff = float("inf")
        best_x_energy = None
        for e in x_lib_energies:
            diff = abs(e - e_x_kev)
            if diff < best_x_diff:
                best_x_diff = diff
                best_x_energy = e

        if best_x_diff > 1.0:
            dropped_x_unmapped += 1
            continue

        # Map Gamma
        g_key = (nuclide, "g")
        if g_key not in lib_data:
            dropped_gamma_unmapped += 1
            continue
        g_lib_energies = lib_data[g_key]
        best_g_diff = float("inf")
        best_g_energy = None
        for e in g_lib_energies:
            diff = abs(e - e_gamma_kev)
            if diff < best_g_diff:
                best_g_diff = diff
                best_g_energy = e

        if best_g_diff > 0.3:
            dropped_gamma_unmapped += 1
            continue

        # Check X-ray energy threshold
        if best_x_energy < 50.0:
            dropped_x_below_50 += 1
            continue

        # Merge key
        merge_key = (nuclide, best_x_energy, best_g_energy)
        if merge_key not in merged_data:
            merged_data[merge_key] = {"p_pair": 0.0, "p_err_sq_sum": 0.0}
        merged_data[merge_key]["p_pair"] += p_corr
        merged_data[merge_key]["p_err_sq_sum"] += p_raw_err ** 2

    # Filter and format output
    output_rows = []
    for (nuclide, e_x_lib, e_g_lib), data in merged_data.items():
        e_sum = round(e_x_lib + e_g_lib, 3)
        p_pair = data["p_pair"]
        p_err = math.sqrt(data["p_err_sq_sum"])

        if e_sum > 2900.0:
            dropped_e_above_2900 += 1
            continue
        if p_pair < 2e-4:
            dropped_p_threshold += 1
            continue

        note = (
            f"#GS-49 X+gamma: {nuclide} X {e_x_lib:.3f} + g {e_g_lib:.3f}; "
            f"p={p_pair:.3e}"
        )
        output_rows.append(
            {
                "nuclide": nuclide,
                "e_x_keV": e_x_lib,
                "e_g_keV": e_g_lib,
                "e_sum_keV": e_sum,
                "p_pair": p_pair,
                "p_err": p_err,
                "note": note,
            }
        )

    # Sort by descending p_pair
    output_rows.sort(key=lambda x: x["p_pair"], reverse=True)

    # Write JSON
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(output_rows, f, ensure_ascii=False, indent=1)

    # Print statistics in Russian
    print(f"Всего входных строк: {total_input}")
    print(f"Отброшено (X не сопоставлено): {dropped_x_unmapped}")
    print(f"Отброшено (gamma не сопоставлено): {dropped_gamma_unmapped}")
    print(f"Отброшено (X < 50 кэВ): {dropped_x_below_50}")
    print(f"Отброшено (p < 2e-4): {dropped_p_threshold}")
    print(f"Отброшено (E_sum > 2900 кэВ): {dropped_e_above_2900}")
    print(f"Количество выходных строк: {len(output_rows)}")

    if output_rows:
        print("Десять крупнейших по вероятности:")
        for row in output_rows[:10]:
            print(
                f"  {row['nuclide']}, E_x={row['e_x_keV']:.3f}, "
                f"E_g={row['e_g_keV']:.3f}, E_sum={row['e_sum_keV']:.3f}, "
                f"p_pair={row['p_pair']:.3e}"
            )


if __name__ == "__main__":
    main()
