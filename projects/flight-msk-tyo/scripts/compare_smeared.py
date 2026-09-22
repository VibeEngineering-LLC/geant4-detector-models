import sys
import math
import pandas as pd

def fmt(x):
    return "nan" if x != x else f"{x:.4g}"

def read_bench(path, Sn_keV):
    with open(path, encoding="utf-8") as f:
        line1 = f.readline().strip()
        line2 = f.readline().strip()
    
    hp_mode = None
    n_capture_gamma = None
    for token in line1[1:].split():
        k, v = token.split("=")
        if k == "hp_mode":
            hp_mode = v
        elif k == "n_capture_gamma":
            n_capture_gamma = int(v)
    
    if hp_mode is None or n_capture_gamma is None:
        sys.stderr.write("Missing hp_mode or n_capture_gamma in bench file\n")
        sys.exit(2)
    
    counts = list(map(int, line2[len("capture,"):].split(",")))
    Ncap = sum((j + 0.5) * c for j, c in enumerate(counts)) / Sn_keV
    if Ncap <= 0:
        sys.stderr.write(f"Ncap <= 0: {Ncap}\n")
        sys.exit(2)
    g = [c / Ncap for c in counts]  # на один захват (первая версия делила на Sn_keV — масштаб был бессмысленный)

    return g, Ncap, hp_mode

def read_lines(element, data_dir):
    prompt_path = f"{data_dir}/promptgammas.xls"
    isotope_path = f"{data_dir}/isotope.xls"
    
    df_prompt = pd.read_excel(prompt_path, header=None, skiprows=2)
    df_isotope = pd.read_excel(isotope_path, header=None, skiprows=2)
    
    df_prompt[1] = df_prompt[1].astype(str).str.strip()
    df_prompt[5] = pd.to_numeric(df_prompt[5], errors="coerce")
    df_prompt[7] = pd.to_numeric(df_prompt[7], errors="coerce")
    df_prompt = df_prompt.dropna(subset=[5, 7])
    
    df_isotope[0] = df_isotope[0].astype(str).ffill()
    df_isotope[1] = df_isotope[1].astype(str).str.strip()
    df_isotope[5] = pd.to_numeric(df_isotope[5], errors="coerce")
    df_isotope[8] = pd.to_numeric(df_isotope[8], errors="coerce")
    df_isotope = df_isotope.dropna(subset=[5, 8])
    
    isotope_dict = {}
    for _, row in df_isotope.iterrows():
        label = row[1]
        ab = row[5]
        sigma0 = row[8]
        isotope_dict[label] = (ab, sigma0)
    
    lines = []
    for _, row in df_prompt.iterrows():
        product_label = row[1]
        energy = row[5]
        sigma_gamma = row[7]
        
        if not product_label.endswith(f"-{element}"):
            continue
        
        Ap = int(product_label.split("-")[0])
        target_label = f"{Ap-1}-{element}"
        
        if target_label not in isotope_dict:
            continue
        
        ab, sigma0 = isotope_dict[target_label]
        num = ab / 100 * sigma_gamma
        lines.append((energy, num))
    
    # Compute denominator
    den_total = 0.0
    for _, row in df_isotope.iterrows():
        label = row[1]
        if label.endswith(f"-{element}"):
            ab = row[5]
            sigma0 = row[8]
            den_total += ab / 100 * sigma0
    
    # Normalize lines to per capture
    for i in range(len(lines)):
        energy, num = lines[i]
        y = num / den_total
        lines[i] = (energy, y)
    
    # Merge lines within 0.3 keV
    lines.sort(key=lambda x: x[0])
    merged = []
    if not lines:
        return []
    
    current_energy = lines[0][0]
    current_y = lines[0][1]
    for energy, y in lines[1:]:
        if abs(energy - current_energy) <= 0.3:
            current_y += y
        else:
            merged.append((current_energy, current_y))
            current_energy = energy
            current_y = y
    merged.append((current_energy, current_y))
    
    # Filter by energy range
    merged = [(e, y) for e, y in merged if 20 <= e <= 12000]
    
    return merged

def fold(sources, fwhm662, window_keV):
    S = [0.0] * 12000
    sigma_factor = 2.35482
    
    for E0, w in sources:
        if E0 < 1 or E0 > 12000:
            continue
        sigma = fwhm662 * math.sqrt(E0 / 661.657) / sigma_factor
        start_bin = max(0, int(E0 - 4 * sigma))
        end_bin = min(12000, int(E0 + 4 * sigma) + 1)
        
        for k in range(start_bin, end_bin):
            diff = (k + 0.5 - E0)
            exp_arg = -0.5 * (diff / sigma)**2
            weight = w * math.exp(exp_arg) / (sigma * math.sqrt(2 * math.pi))
            S[k] += weight
    
    # Window sums
    n_windows = 12000 // window_keV
    W = []
    for i in range(n_windows):
        lo = i * window_keV
        hi = (i + 1) * window_keV
        total = sum(S[lo:hi])
        W.append(total)
    
    return W

def main(bench_csv, element, Sn_keV, data_dir, fwhm662_keV=41.6, window_keV=500):
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
    
    g, Ncap, hp_mode = read_bench(bench_csv, Sn_keV)
    lines = read_lines(element, data_dir)
    
    # Model energy and multiplicity
    E_model = sum(e * y for e, y in lines)
    mult_model = sum(y for _, y in lines)
    mult_g4 = sum(g)
    E_g4 = Sn_keV
    
    # Folding
    model_sources = [(e, y) for e, y in lines]
    W_model = fold(model_sources, fwhm662_keV, window_keV)
    
    g4_sources = [(j + 0.5, g[j]) for j in range(12000) if g[j] > 0]
    W_g4 = fold(g4_sources, fwhm662_keV, window_keV)
    
    # Output
    n_lines = len(lines)
    print(f"# file={bench_csv} element={element} hp_mode={hp_mode} Ncap_est={Ncap:.6g} lines={n_lines}")
    print(f"BALANCE E_model_keV={E_model:.6g} Sn_keV={Sn_keV:.6g} mult_model={mult_model:.4g} mult_g4={mult_g4:.4g}")
    
    count = 0
    within30 = 0
    ratios = []
    
    for i in range(len(W_model)):
        lo = i * window_keV
        hi = (i + 1) * window_keV
        model = W_model[i]
        g4 = W_g4[i]
        
        if model < 1e-6:
            ratio = float('nan')
        else:
            ratio = g4 / model
        
        if model >= 0.005:
            count += 1
            if 0.7 <= ratio <= 1.3:
                within30 += 1
            ratios.append(ratio)
        
        print(f"WIN {lo}-{hi} | {fmt(model)} | {fmt(g4)} | {fmt(ratio)}")
    
    median_ratio = sorted(ratios)[len(ratios)//2] if ratios else float('nan')
    print(f"WINDOWS n={count} within30={within30} median_ratio={fmt(median_ratio)}")

if __name__ == "__main__":
    if len(sys.argv) < 5:
        sys.stderr.write("Usage: python compare_smeared.py <bench_csv> <element> <Sn_keV> <data_dir> [fwhm662_keV=41.6] [window_keV=500]\n")
        sys.exit(2)
    
    bench_csv = sys.argv[1]
    element = sys.argv[2]
    Sn_keV = float(sys.argv[3])
    data_dir = sys.argv[4]
    fwhm662_keV = float(sys.argv[5]) if len(sys.argv) > 5 else 41.6
    window_keV = int(sys.argv[6]) if len(sys.argv) > 6 else 500
    
    main(bench_csv, element, Sn_keV, data_dir, fwhm662_keV, window_keV)
