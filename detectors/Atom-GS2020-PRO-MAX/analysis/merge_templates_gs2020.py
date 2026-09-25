# -*- coding: utf-8 -*-
r"""Сумматор отрезков прогонов стенда gs2020_marinelli в формате шаблона метода (шапка key,value; маркер bin_keV;
строки E,count_edep,count_light) и сборка шаблона цепочки Th-232 в равновесии. Спека: scripts\specs\SPEC-merge_templates_gs2020.md
Запуск: python merge_templates_gs2020.py <каталог_отрезков> <каталог_выхода>"""
import sys, os, re, glob

sys.stdout.reconfigure(encoding="utf-8")


def read_chunk(path):
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    comments = []
    header = []
    bins = []
    edep = []
    light = []
    marker_found = False

    for line in lines:
        stripped = line.rstrip("\n").rstrip("\r")
        if not marker_found:
            if stripped.startswith("#"):
                comments.append(stripped)
            elif stripped == "bin_keV,count_edep,count_light":
                marker_found = True
            else:
                if "," in stripped:
                    key, value = stripped.split(",", 1)
                    header.append((key, value))
        else:
            if stripped.strip():
                parts = stripped.split(",")
                bins.append(parts[0])
                edep.append(int(parts[1]))
                light.append(int(parts[2]))

    if not marker_found:
        raise SystemExit("ОТКАЗ: нет маркера bin_keV в " + path)

    return {
        "comments": comments,
        "header": header,
        "bins": bins,
        "edep": edep,
        "light": light,
    }


def write_template(path, comments, header, bins, edep, light):
    with open(path, "w", encoding="utf-8", newline="") as f:
        for c in comments:
            f.write(c + "\n")
        for k, v in header:
            f.write(f"{k},{v}\n")
        f.write("bin_keV,count_edep,count_light\n")
        for i in range(len(bins)):
            f.write(f"{bins[i]},{edep[i]},{light[i]}\n")


def merge_group(paths):
    chunks = [read_chunk(p) for p in paths]

    # Check n_events_processed == n_events_requested
    for p_i, c in zip(paths, chunks):
        h_dict = dict(c["header"])
        if "n_events_processed" not in h_dict or "n_events_requested" not in h_dict:
            raise SystemExit("ОТКАЗ: отсутствуют ключи счетчиков в шапке")
        if h_dict["n_events_processed"] != h_dict["n_events_requested"]:
            raise SystemExit("ОТКАЗ: отрезок не досчитан: " + p_i)

    # Check seeds unique
    seeds = []
    for c in chunks:
        h_dict = dict(c["header"])
        if "seed" not in h_dict:
            raise SystemExit("ОТКАЗ: отсутствует ключ seed")
        s = h_dict["seed"]
        if s in seeds:
            raise SystemExit(f"ОТКАЗ: повтор зерна {s}")
        seeds.append(s)

    # Check identical headers except specific keys
    exclude_keys = {"seed", "n_events_requested", "n_events_processed", "n_with_edep"}
    ref_header_dict = dict(chunks[0]["header"])
    
    for i, c in enumerate(chunks[1:], 1):
        curr_header_dict = dict(c["header"])
        for k, v in ref_header_dict.items():
            if k not in exclude_keys:
                if k not in curr_header_dict or curr_header_dict[k] != v:
                    val_curr = curr_header_dict.get(k, "MISSING")
                    raise SystemExit(f"ОТКАЗ: несовпадение ключа {k}: '{v}' vs '{val_curr}'")

    # Check bins identical
    ref_bins = chunks[0]["bins"]
    for i, c in enumerate(chunks[1:], 1):
        if c["bins"] != ref_bins:
            raise SystemExit("ОТКАЗ: несовпадение бинов")

    # Sum edep and light
    total_edep = [0] * len(ref_bins)
    total_light = [0] * len(ref_bins)
    
    n_req_sum = 0
    n_proc_sum = 0
    n_with_edep_sum = 0

    for c in chunks:
        h_dict = dict(c["header"])
        n_req_sum += int(h_dict.get("n_events_requested", 0))
        n_proc_sum += int(h_dict.get("n_events_processed", 0))
        n_with_edep_sum += int(h_dict.get("n_with_edep", 0))
        
        for j in range(len(ref_bins)):
            total_edep[j] += c["edep"][j]
            total_light[j] += c["light"][j]

    # Construct output header
    out_header = list(chunks[0]["header"])
    seed_str = "merged:" + ";".join(sorted(seeds))
    
    new_header = []
    for k, v in out_header:
        if k == "seed":
            new_header.append((k, seed_str))
        elif k == "n_events_requested":
            new_header.append((k, str(n_req_sum)))
        elif k == "n_events_processed":
            new_header.append((k, str(n_proc_sum)))
        elif k == "n_with_edep":
            new_header.append((k, str(n_with_edep_sum)))
        else:
            new_header.append((k, v))
    
    new_header.append(("merged_chunks", str(len(chunks))))

    return (new_header, ref_bins, total_edep, total_light, n_proc_sum)


def build_chain(out_dir, merged):
    BR = {"Th232":1.0,"Ac228":1.0,"Th228":1.0,"Ra224":1.0,"Rn220":1.0,"Pb212":1.0,"Bi212":1.0,"Tl208":0.3594}
    
    required_keys = list(BR.keys())
    for key in required_keys:
        group_name = f"mix_{key}"
        if group_name not in merged:
            raise SystemExit(f"ОТКАЗ: отсутствует группа {group_name}")

    # Get N from Th232
    th232_data = merged["mix_Th232"]
    N = th232_data[4] # n_processed
    
    # Check ratios
    for key in required_keys:
        group_name = f"mix_{key}"
        data = merged[group_name]
        n_key = data[4]
        expected = BR[key] * N
        if abs(n_key - expected) > 1e-6 * N:
            raise SystemExit(f"ОТКАЗ: звено {key}: n={n_key}, ожидалось {expected:.0f} (n/BR должно быть одинаково у всех звеньев)")

    # Check bins identical across all 8
    ref_bins = th232_data[1]
    for key in required_keys[1:]:
        group_name = f"mix_{key}"
        data = merged[group_name]
        if data[1] != ref_bins:
            raise SystemExit("ОТКАЗ: несовпадение бинов в цепочке")

    # Sum edep and light
    total_edep = [0] * len(ref_bins)
    total_light = [0] * len(ref_bins)
    n_with_edep_sum = 0
    
    for key in required_keys:
        group_name = f"mix_{key}"
        data = merged[group_name]
        h_dict = dict(data[0])
        n_with_edep_sum += int(h_dict.get("n_with_edep", 0))
        
        for j in range(len(ref_bins)):
            total_edep[j] += data[2][j]
            total_light[j] += data[3][j]

    # Construct header from Th232
    th_header = list(th232_data[0])
    new_header = []
    for k, v in th_header:
        if k == "particle":
            new_header.append((k, "chain"))
        elif k == "ion_z":
            new_header.append((k, "0"))
        elif k == "ion_a":
            new_header.append((k, "0"))
        elif k == "n_events_requested":
            new_header.append((k, str(N)))
        elif k == "n_events_processed":
            new_header.append((k, str(N)))
        elif k == "n_with_edep":
            new_header.append((k, str(n_with_edep_sum)))
        elif k == "seed":
            new_header.append((k, "chain_of_merged"))
        else:
            new_header.append((k, v))
    
    new_header.append(("chain_links", "Th232;Ac228;Th228;Ra224;Rn220;Pb212;Bi212;Tl208 BR Tl208=0.3594"))

    out_path = os.path.join(out_dir, "mix_Th232chain_npsmoff.csv")
    write_template(out_path, th232_data[0][0] if isinstance(th232_data[0], list) else [], new_header, ref_bins, total_edep, total_light) # Note: comments are not passed in merged tuple structure defined above, need to fix read_chunk/merge_group return or handle comments separately. 
    # Correction: merge_group returns (header, bins, edep, light, n_processed). It does NOT return comments.
    # The prompt says "Comments: first chunk's comments" for merge_group output. But build_chain uses merged dict which stores the return of merge_group.
    # So build_chain doesn't have access to comments easily unless we store them. 
    # However, looking at write_template signature: write_template(path, comments, header, bins, edep, light).
    # In main(), we call write_template for groups. We need comments there.
    # In merge_group, I didn't return comments. I should probably just use empty comments for chain or retrieve from first chunk if needed. 
    # The prompt says for Chain: "Header: Th232's header...". It doesn't explicitly specify comments for the chain file. 
    # Let's assume empty comments for chain or reuse Th232's if we had them. Since merge_group discards them, I'll pass empty list for chain comments to avoid crash, or better, modify merge_group to return comments?
    # "Return (header, bins, edep, light, n_processed)." -> Strictly 5 items. So comments are lost after merge_group.
    # For the chain file, I will pass an empty list of comments to write_template.

    write_template(out_path, [], new_header, ref_bins, total_edep, total_light)
    
    return N, sum(total_edep), sum(total_light)


def main():
    if len(sys.argv) != 3:
        raise SystemExit("Использование: python merge_templates_gs2020.py <каталог_отрезков> <каталог_выхода>")

    chunk_dir = sys.argv[1]
    out_dir = sys.argv[2]

    os.makedirs(out_dir, exist_ok=True)

    # Collect files
    pattern = os.path.join(chunk_dir, "*.csv")
    files = glob.glob(pattern)
    
    regex = re.compile(r"^(grid_mar_E[0-9.]+|mix_[A-Za-z0-9]+)_s(\d+)\.csv$")
    
    groups = {}
    for f in files:
        basename = os.path.basename(f)
        m = regex.match(basename)
        if m:
            group_name = m.group(1)
            if group_name not in groups:
                groups[group_name] = []
            groups[group_name].append(f)

    if not groups:
        raise SystemExit("ОТКАЗ: не найдено файлов для обработки")

    merged = {}
    
    for group_name, paths in groups.items():
        # Sort paths to ensure deterministic order (e.g. by seed or filename)
        paths.sort()
        
        header, bins, edep, light, n_processed = merge_group(paths)
        merged[group_name] = (header, bins, edep, light, n_processed)
        
        # Determine output name
        if group_name.startswith("grid_mar_E"):
            out_name = f"{group_name}.csv"
        elif group_name.startswith("mix_"):
            out_name = f"{group_name}_npsmoff.csv"
        else:
            out_name = f"{group_name}.csv"
            
        out_path = os.path.join(out_dir, out_name)
        
        # We need comments for write_template. merge_group didn't return them.
        # Let's read the first chunk again to get comments? Or modify merge_group?
        # To strictly follow "Return (header, bins, edep, light, n_processed)", I can't change return.
        # But write_template needs comments.
        # I will read the first file in the group to get comments.
        first_chunk = read_chunk(paths[0])
        
        write_template(out_path, first_chunk["comments"], header, bins, edep, light)
        
        total_edep_count = sum(edep)
        print(f"Группа {group_name}: файл {out_name}, отрезков {len(paths)}, событий {n_processed}, сумма edep {total_edep_count}")

    # Build chain
    N, total_chain_edep, total_chain_light = build_chain(out_dir, merged)
    print(f"Цепочка Th-232: событий {N}, сумма edep {total_chain_edep}, сумма light {total_chain_light}")


if __name__ == "__main__":
    main()
