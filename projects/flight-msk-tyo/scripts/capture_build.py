import sys, os
sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from capture_levels import read_levels, normalise
from capture_data import read_masses, sn_keV, read_targets, read_prompt, merge_lines, product_label
from capture_classify import classify_lines
from capture_model import reachable_levels, propagate_yields

def process_isotope(symbol, Z, A, sigma0, masses, prompt, levels_dir, tol_a, tol_b):
    Sn = sn_keV(masses, Z, A)
    if Sn is None or Sn <= 0:
        return {"skip": "no Sn", "label": f"{Z}-{A}"}
    path = os.path.join(levels_dir, f"z{Z}.a{A+1}")
    try:
        levels = normalise(read_levels(path))
    except (OSError, ValueError):
        return {"skip": "no level scheme", "label": f"{Z}-{A}"}
    label = product_label(symbol, A)
    if label not in prompt:
        return {"skip": "no PGAA lines", "label": f"{Z}-{A}"}
    if sigma0 is None or sigma0 <= 0:
        return {"skip": "no sigma0", "label": f"{Z}-{A}"}
    lines = merge_lines([(E, sg / sigma0) for E, sg in prompt[label]])
    M = (A + 1) * 931494.0
    records = classify_lines(levels, Sn, M, lines, tol_a, tol_b)
    primary = {}
    for E, y, is_p, j in records:
        if is_p:
            primary[j] = primary.get(j, 0.0) + y
    Wp = sum(primary.values())
    scaled = False
    if Wp > 1.0:
        for j in primary:
            primary[j] /= Wp
        Wp = 1.0
        scaled = True
    if not primary:
        return {"skip": "no primary lines", "label": f"{Z}-{A}"}
    reachable = reachable_levels(levels, primary)
    model = propagate_yields(levels, primary, Sn, M)
    n_check = 0
    n_within2 = 0
    unexplained = []
    for E, y, is_p, j in records:
        if not is_p and y >= 0.005:
            tol = tol_a + tol_b * E
            closest = None
            for E_m, y_m in model:
                if abs(E - E_m) <= tol:
                    if closest is None or abs(E - E_m) < abs(E - closest[0]):
                        closest = (E_m, y_m)
            n_check += 1
            if closest and 0.5 <= closest[1] / y <= 2.0:
                n_within2 += 1
            # неразобранные линии — в ТОЙ ЖЕ ветке (первая версия держала их в elif, недостижимом для y >= 0.02)
            if y >= 0.02 and (not closest or closest[1] < y / 2):
                unexplained.append((E, y, closest[1] if closest else 0.0))
    return {
        "label": label,
        "Z": Z,
        "A": A,
        "symbol": symbol,
        "Sn": Sn,
        "sigma0": sigma0,
        "Wp": Wp,
        "scaled": scaled,
        "n_lines": len(lines),
        "primary": primary,
        "reachable": reachable,
        "levels": levels,
        "n_check": n_check,
        "n_within2": n_within2,
        "unexplained": unexplained
    }

def write_db(results, elements, out_db):
    written = [r for r in results if "skip" not in r]
    with open(out_db, "w", encoding="utf-8", newline="\n") as f:
        f.write(f"# capture_db v1 targets={len(written)} elements={elements}\n")
        for r in sorted(written, key=lambda x: (x["Z"], x["A"])):
            f.write(f"ISO {r['Z']} {r['A']} {r['Sn']:.4f} {r['sigma0']:.6g} {r['Wp']:.6f} {len(r['reachable'])} {len(r['primary'])}\n")
            for j in sorted(r["primary"]):
                f.write(f"PRI {j} {r['primary'][j]:.6f}\n")
            for idx in sorted(r["reachable"]):
                level = r["levels"][idx]
                f.write(f"LEV {idx} {level['E']:.4f} {len(level['trans'])}\n")
                for d, Eg, prob, alpha in level["trans"]:
                    f.write(f"TR {d} {Eg:.4f} {prob:.6f} {alpha:.6g}\n")
            f.write("END\n")

def write_report(results, elements, out_report):
    written = [r for r in results if "skip" not in r]
    skipped = [r for r in results if "skip" in r]
    with open(out_report, "w", encoding="utf-8") as f:
        f.write(f"# Capture Report for Elements: {elements}\n\n")
        f.write("| Isotope | Sn (keV) | σ₀ (b) | PGAA lines | Primaries | Wp | Fallback 1-Wp | Check | Within2x |\n")
        f.write("|---------|----------|--------|------------|-----------|----|---------------|-------|----------|\n")
        for r in sorted(written, key=lambda x: (x["Z"], x["A"])):
            f.write(f"| {r['label']} | {r['Sn']:.3f} | {r['sigma0']:.2e} | {r['n_lines']} | {len(r['primary'])} | {r['Wp']:.4f} | {1 - r['Wp']:.4f} | {r['n_check']} | {r['n_within2']} |\n")
        if any(r["scaled"] for r in written):
            f.write("\nScaled Wp > 1.0 detected:\n\n")
            for r in written:
                if r["scaled"]:
                    f.write(f"- {r['label']}\n")
        if any(r["unexplained"] for r in written):
            f.write("\n## Unexplained lines\n\n")
            for r in written:
                if r["unexplained"]:
                    f.write(f"- **{r['label']}**:\n")
                    for E, y, y_m in r["unexplained"]:
                        f.write(f"  - {E:.2f} keV, y={y:.4f}, model={y_m:.4f}\n")
        if skipped:
            f.write("\n## Skipped isotopes\n\n")
            for r in skipped:
                f.write(f"- **{r['label']}**: {r['skip']}\n")

def main():
    if len(sys.argv) < 6:
        print("Usage: python capture_build.py <data_dir> <levels_dir> <elements> <out_db> <out_report> [tol_a=1.5] [tol_b=3.0e-4]", file=sys.stderr)
        sys.exit(2)
    data_dir = sys.argv[1]
    levels_dir = sys.argv[2]
    elements = sys.argv[3]
    out_db = sys.argv[4]
    out_report = sys.argv[5]
    tol_a = float(sys.argv[6]) if len(sys.argv) > 6 else 1.5
    tol_b = float(sys.argv[7]) if len(sys.argv) > 7 else 3e-4
    masses = read_masses(os.path.join(data_dir, "mass.mas20.txt"))
    targets = read_targets(os.path.join(data_dir, "isotope.xls"))
    prompt = read_prompt(os.path.join(data_dir, "promptgammas.xls"))
    results = []
    for symbol in elements.split(","):
        if symbol not in targets:
            print(f"Warning: element {symbol} not found in targets", file=sys.stderr)
            continue
        for Z, A, abundance, sigma0 in targets[symbol]:
            result = process_isotope(symbol, Z, A, sigma0, masses, prompt, levels_dir, tol_a, tol_b)
            results.append(result)
    written = [r for r in results if "skip" not in r]
    skipped = [r for r in results if "skip" in r]
    for r in written:
        print(f"ISO {r['Z']}-{r['A']} Sn={r['Sn']:.3f} sigma0={r['sigma0']:.4g} Wp={r['Wp']:.4f} primaries={len(r['primary'])} levels={len(r['reachable'])} check={r['n_check']} within2x={r['n_within2']}")
    if not written:
        print("ERROR no isotope processed", file=sys.stderr)
        sys.exit(3)
    write_db(results, elements, out_db)
    write_report(results, elements, out_report)
    print(f"DONE isotopes={len(written)} skipped={len(skipped)}")

if __name__ == "__main__":
    main()
