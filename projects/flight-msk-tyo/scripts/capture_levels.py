import sys

def read_levels(path):
    levels = {}
    last_level_idx = None
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line: continue
            tokens = line.split()
            if len(tokens) < 2: raise ValueError(f"Недостаточно токенов в строке {i}")
            try: float(tokens[1])
            except ValueError:
                # LEVEL line
                if len(tokens) < 6: raise ValueError(f"Недостаточно токенов в строке {i}")
                try:
                    idx = int(tokens[0])
                    E = float(tokens[2])
                    halflife = float(tokens[3])
                    Jpi = tokens[4]
                    n_gamma = int(tokens[5])
                except ValueError: raise ValueError(f"Неверный формат чисел в строке {i}")
                levels[idx] = {"E": E, "trans": []}
                last_level_idx = idx
            else:
                # TRANSITION line
                if last_level_idx is None: raise ValueError(f"Переход до уровня в строке {i}")
                try:
                    daughter_idx = int(tokens[0])
                    Eg = float(tokens[1])
                    Irel = float(tokens[2])
                    multipolarity = tokens[3]
                    mixing = tokens[4]
                    alpha = float(tokens[5]) if len(tokens) > 5 else 0.0
                except ValueError: raise ValueError(f"Неверный формат чисел в строке {i}")
                weight = Irel * (1.0 + alpha)
                levels[last_level_idx]["trans"].append((daughter_idx, Eg, weight, alpha))
    return levels

def normalise(levels):
    norm_levels = {}
    for idx, data in levels.items():
        E = data["E"]
        trans = data["trans"]
        total_weight = sum(w for _, _, w, _ in trans)
        if total_weight == 0:
            norm_trans = [(d, e, 0.0, a) for d, e, _, a in trans]
        else:
            norm_trans = [(d, e, w/total_weight, a) for d, e, w, a in trans]
        norm_levels[idx] = {"E": E, "trans": norm_trans}
    return norm_levels

def main():
    if len(sys.argv) != 2:
        print("Использование: python capture_levels.py <файл>", file=sys.stderr)
        sys.exit(2)
    try:
        levels = read_levels(sys.argv[1])
        norm_levels = normalise(levels)
        n_levels = len(norm_levels)
        n_transitions = sum(len(level["trans"]) for level in norm_levels.values())
        ground_E = norm_levels[0]["E"]
        max_E = max(level["E"] for level in norm_levels.values())
        print(f"levels={n_levels} transitions={n_transitions} ground_E={ground_E:.4f} max_E={max_E:.4f}")
    except ValueError as e:
        print(str(e), file=sys.stderr)
        sys.exit(2)

if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
    main()
