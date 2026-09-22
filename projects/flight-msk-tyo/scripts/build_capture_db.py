import sys
import math
import os
import re
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")

def product_label(symbol, A):
    return f"{A+1}-{symbol}"

def read_masses(path):
    mass_excess = {}
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if len(line) < 43 or line[0] not in " 01":
                continue
            try:
                N = int(line[4:9])
                Z = int(line[9:14])
                A = int(line[14:19])
                me = float(line[28:42].strip().replace("#", "."))
                mass_excess[(Z, A)] = me
            except:
                continue
    return mass_excess

def read_isotopes(path):
    df = pd.read_excel(path, header=None, skiprows=2)
    df[0] = df[0].ffill()
    df = df.dropna(subset=[1, 5, 8])
    df[1] = df[1].str.strip()
    df[5] = pd.to_numeric(df[5], errors="coerce")
    df[8] = pd.to_numeric(df[8], errors="coerce")
    targets = {}
    for _, row in df.iterrows():
        symbol = row[0]
        label = row[1]
        abundance = row[5]
        sigma0 = row[8]
        if not symbol or not label or math.isnan(abundance) or abundance <= 0:
            continue
        try:
            A = int(label.split("-")[0])
            Z = int(row[2])
        except:
            continue
        if symbol not in targets:
            targets[symbol] = []
        targets[symbol].append((Z, A, abundance, sigma0))
    return targets

def read_prompt(path):
    df = pd.read_excel(path, header=None, skiprows=2)
    df[1] = df[1].str.strip()
    df[5] = pd.to_numeric(df[5], errors="coerce")
    df[7] = pd.to_numeric(df[7], errors="coerce")
    df = df.dropna(subset=[1, 5, 7])
    prompt_lines = {}
    for _, row in df.iterrows():
        label = row[1]
        E = row[5]
        sigma_gamma = row[7]
        if label not in prompt_lines:
            prompt_lines[label] = []
        prompt_lines[label].append((E, sigma_gamma))
    return prompt_lines

def read_levels(path):
    levels = {}
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                idx = int(parts[0])
                flag = parts[1]
                E_level = float(parts[2])
                Jpi = parts[5]
                if not flag.replace("-", "").replace("+", "").isdigit():
                    levels[idx] = (E_level, [])
                else:
                    daughter_idx = int(parts[0])
                    Eg = float(parts[2])
                    Irel = float(parts[4])
                    alpha = float(parts[8]) if len(parts) > 8 and parts[8].replace(".", "").replace("-", "").isdigit() else 0.0
                    weight = Irel * (1 + alpha)
                    levels[idx][1].append((daughter_idx, Eg, weight, alpha))
            except:
                continue
    return levels

def merge_lines(lines):
    lines.sort(key=lambda x: x[0])
    merged = []
    i = 0
    while i < len(lines):
        E = lines[i][0]
        y_total = lines[i][1]
        j = i + 1
        while j < len(lines) and lines[j][0] - E < 0.3:
            y_total += lines[j][1]
            j += 1
        merged.append((E, y_total))
        i = j
    return merged

def classify_primary(E_lines, levels, Sn, M):
    primary = []
    secondary = []
    for i, (E, y) in enumerate(E_lines):
        Lt = Sn - E - E * E / (2 * M)
        tol = 1.5 + 3e-4 * E
        min_dist = float("inf")
        j_min = -1
        for j, (E_level, _) in levels.items():
            dist = abs(E_level - Lt)
            if dist < min_dist:
                min_dist = dist
                j_min = j
        is_primary = min_dist <= tol
        is_secondary = False
        for _, Eg, _, _ in levels[j_min][1]:
            if abs(Eg - E) <= tol:
                is_secondary = True
                break
        if (E >= 0.5 * Sn and is_primary) or (E < 0.5 * Sn and is_primary and not is_secondary):
            primary.append((i, y))
        else:
            secondary.append((i, y))
    return primary, secondary

def propagate(levels, primary_weights):
    visit_prob = {}
    for idx in levels:
        visit_prob[idx] = 0.0
    visit_prob[0] = 1.0
    for idx in sorted(levels.keys()):
        if idx == 0:
            continue
        E_level, transitions = levels[idx]
        for (daughter_idx, Eg, weight, alpha) in transitions:
            prob = weight / sum(t[2] for _, _, t in levels.items() if any(d[0] == idx for d in t))
            visit_prob[daughter_idx] += visit_prob[idx] * prob * 1.0 / (1 + alpha)
    return visit_prob

def main():
    if len(sys.argv) != 6:
        print("Usage: python build_capture_db.py <data_dir> <levels_dir> <elements> <out_db> <out_report>", file=sys.stderr)
        sys.exit(2)

    data_dir, levels_dir, elements_str, out_db, out_report = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
    elements = elements_str.split(",")

    mass_excess = read_masses(os.path.join(data_dir, "mass.mas20.txt"))
    targets = read_isotopes(os.path.join(data_dir, "isotope.xls"))
    prompt_lines = read_prompt(os.path.join(data_dir, "promptgammas.xls"))

    out_db_f = open(out_db, "w", encoding="utf-8")
    out_report_f = open(out_report, "w", encoding="utf-8")

    out_db_f.write(f"# capture_db v1 targets=<0> elements={elements_str}\n")

    written = 0
    skipped = 0
    report_lines = []
    skipped_lines = []

    for symbol in sorted(targets.keys()):
        if symbol not in elements:
            continue
        for Z, A, abundance, sigma0 in targets[symbol]:
            product_label_str = product_label(symbol, A)
            if product_label_str not in prompt_lines:
                skipped_lines.append((f"{Z}-{A}", "No PGAA lines"))
                skipped += 1
                continue

            if (Z, A) not in mass_excess or (Z, A + 1) not in mass_excess:
                skipped_lines.append((f"{Z}-{A}", "Missing mass data"))
                skipped += 1
                continue

            Sn = mass_excess[(Z, A)] + mass_excess[(0, 1)] - mass_excess[(Z, A + 1)]
            if Sn <= 0:
                skipped_lines.append((f"{Z}-{A}", "Non-positive Sn"))
                skipped += 1
                continue

            M = (A + 1) * 931494.0
            level_file = os.path.join(levels_dir, f"z{Z}.a{A+1}")
            if not os.path.exists(level_file):
                skipped_lines.append((f"{Z}-{A}", "No level file"))
                skipped += 1
                continue

            levels = read_levels(level_file)
            if not levels:
                skipped_lines.append((f"{Z}-{A}", "Empty level file"))
                skipped += 1
                continue

            E_lines = prompt_lines[product_label_str]
            merged_lines = merge_lines(E_lines)

            primary, secondary = classify_primary(merged_lines, levels, Sn, M)
            Wp = sum(y for _, y in primary)
            if Wp > 1.0:
                scale = 1.0 / Wp
                for i, _ in enumerate(primary):
                    primary[i] = (primary[i][0], primary[i][1] * scale)
                Wp = 1.0

            reachable = set()
            queue = [idx for idx, _ in primary]
            while queue:
                idx = queue.pop(0)
                if idx in reachable:
                    continue
                reachable.add(idx)
                for daughter_idx, _, _, _ in levels[idx][1]:
                    queue.append(daughter_idx)

            out_db_f.write(f"ISO {Z} {A} {Sn:.4f} {sigma0:.6g} {Wp:.6f} {len(reachable)} {len(primary)}\n")
            for idx, weight in primary:
                out_db_f.write(f"PRI {idx} {weight:.6f}\n")

            visit_prob = propagate(levels, [y for _, y in primary])
            check_lines = []
            for i, (E, y) in enumerate(merged_lines):
                if any(i == idx for idx, _ in primary):
                    continue
                if y < 0.005:
                    continue
                Lt = Sn - E - E * E / (2 * M)
                tol = 1.5 + 3e-4 * E
                min_dist = float("inf")
                j_min = -1
                for j, (E_level, _) in levels.items():
                    dist = abs(E_level - Lt)
                    if dist < min_dist:
                        min_dist = dist
                        j_min = j
                if min_dist <= tol and j_min in reachable:
                    model_yield = visit_prob[j_min] * 1.0 / (1 + levels[j_min][1][0][3]) if levels[j_min][1] else 0.0
                    if model_yield > 0 and 0.5 <= y / model_yield <= 2.0:
                        check_lines.append((E, y, model_yield))
            n_check = len(check_lines)
            n_within2 = sum(1 for _, y, model_yield in check_lines if 0.5 <= y / model_yield <= 2.0)

            for idx in sorted(reachable):
                E_level, transitions = levels[idx]
                out_db_f.write(f"LEV {idx} {E_level:.4f} {len(transitions)}\n")
                prob_sum = sum(weight for _, _, weight, _ in transitions)
                if prob_sum > 0:
                    for daughter_idx, Eg, weight, alpha in transitions:
                        prob = weight / prob_sum
                        out_db_f.write(f"TR {daughter_idx} {Eg:.4f} {prob:.6f} {alpha:.6g}\n")
                else:
                    for daughter_idx, Eg, _, alpha in transitions:
                        out_db_f.write(f"TR {daughter_idx} {Eg:.4f} 0.000000 {alpha:.6g}\n")

            out_db_f.write("END\n")

            report_lines.append((f"{Z}-{A}", f"{Sn:.3f}", f"{sigma0:.4g}", str(len(merged_lines)), str(len(primary)), f"{Wp:.4f}", f"{1-Wp:.4f}", str(n_check), str(n_within2)))

            written += 1

    out_db_f.close()

    out_report_f.write("# Capture Database Report\n\n")
    out_report_f.write("| isotope | Sn | sigma0 | PGAA lines | primaries | Wp | fallback 1-Wp | check | within2x |\n")
    out_report_f.write("|---------|----|--------|------------|-----------|----|---------------|-------|----------|\n")
    for row in report_lines:
        out_report_f.write(f"| {' | '.join(row)} |\n")

    out_report_f.write("\n## Unexplained lines (y >= 0.02)\n\n")
    for symbol in sorted(targets.keys()):
        if symbol not in elements:
            continue
        for Z, A, abundance, sigma0 in targets[symbol]:
            product_label_str = product_label(symbol, A)
            if product_label_str not in prompt_lines:
                continue
            E_lines = prompt_lines[product_label_str]
            merged_lines = merge_lines(E_lines)
            primary, _ = classify_primary(merged_lines, levels, Sn, M)
            primary_indices = {i for i, _ in primary}
            out_report_f.write(f"### {Z}-{A}\n\n")
            for i, (E, y) in enumerate(merged_lines):
                if i in primary_indices or y < 0.02:
                    continue
                out_report_f.write(f"- E = {E:.3f} keV, y = {y:.6g}\n")

    out_report_f.write("\n## Skipped isotopes\n\n")
    for reason in skipped_lines:
        out_report_f.write(f"- {reason[0]}: {reason[1]}\n")

    out_report_f.close()

    print(f"DONE isotopes=<written> skipped=<skipped>")

if __name__ == "__main__":
    main()
