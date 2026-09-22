import sys
import math
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")

def read_masses(path):
    masses = {}
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if len(line) < 43 or line[0] not in " 01":
                continue
            try:
                N = int(line[4:9])
                Z = int(line[9:14])
                A = int(line[14:19])
                excess = float(line[28:42].strip().replace("#", "."))
                masses[(Z, A)] = excess
            except:
                continue
    return masses

def sn_keV(masses, Z, A_target):
    if (Z, A_target) not in masses or (0, 1) not in masses or (Z, A_target + 1) not in masses:
        return None
    return masses[(Z, A_target)] + masses[(0, 1)] - masses[(Z, A_target + 1)]

def read_targets(path):
    df = pd.read_excel(path, header=None, skiprows=2)
    df[0] = df[0].ffill()
    df[0] = df[0].apply(str).str.strip()
    df[2] = pd.to_numeric(df[2], errors="coerce")
    df[3] = pd.to_numeric(df[3], errors="coerce")
    df[5] = pd.to_numeric(df[5], errors="coerce")
    df[8] = pd.to_numeric(df[8], errors="coerce")
    df = df.dropna(subset=[2, 3, 5, 8])
    df = df[(df[5] > 0) & (df[8] > 0)]
    df[2] = df[2].astype(int)
    df[3] = df[3].astype(int)
    targets = {}
    for _, row in df.iterrows():
        symbol = row[0]
        Z = row[2]
        A = row[3]
        abundance = row[5]
        sigma0 = row[8]
        if symbol not in targets:
            targets[symbol] = []
        targets[symbol].append((Z, A, abundance, sigma0))
    return targets

def read_prompt(path):
    df = pd.read_excel(path, header=None, skiprows=2)
    df[1] = df[1].apply(str).str.strip()
    df[5] = pd.to_numeric(df[5], errors="coerce")
    df[7] = pd.to_numeric(df[7], errors="coerce")
    df = df.dropna(subset=[1, 5, 7])
    prompt = {}
    for _, row in df.iterrows():
        label = row[1]
        E = row[5]
        sigma_gamma = row[7]
        if label not in prompt:
            prompt[label] = []
        prompt[label].append((E, sigma_gamma))
    return prompt

def merge_lines(lines, tol=0.3):
    if not lines:
        return []
    lines.sort()
    merged = [lines[0]]
    for E, w in lines[1:]:
        last_E, last_w = merged[-1]
        if abs(E - last_E) <= tol:
            new_E = (E * w + last_E * last_w) / (w + last_w)
            merged[-1] = (new_E, w + last_w)
        else:
            merged.append((E, w))
    return merged

def product_label(symbol, A_target):
    return f"{A_target+1}-{symbol}"

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python capture_data.py <data_dir>", file=sys.stderr)
        sys.exit(2)
    data_dir = sys.argv[1]
    masses = read_masses(f"{data_dir}/mass.mas20.txt")
    print(f"masses={len(masses)}")
    print(f"Sn Al-27 = {sn_keV(masses, 13, 27):.3f}")
    print(f"Sn Fe-56 = {sn_keV(masses, 26, 56):.3f}")
    print(f"Sn H-1 = {sn_keV(masses, 1, 1):.3f}")
    targets = read_targets(f"{data_dir}/isotope.xls")
    print(f"targets_H={len(targets.get('H', []))}")
    for Z, A, abundance, sigma0 in targets.get("Al", []):
        if A == 27:
            print(f"Al27 sigma0={sigma0}")
            break
    prompt = read_prompt(f"{data_dir}/promptgammas.xls")
    print(f"prompt_labels={len(prompt)}")
    lines = prompt.get("2-H", [])
    print(f"H-2 lines={len(lines)}")
    if lines:
        E, sigma_gamma = lines[0]
        print(f"first H-2 line E={E:.3f} sigma_gamma={sigma_gamma:.4g}")
