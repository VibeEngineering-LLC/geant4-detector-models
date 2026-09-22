import sys
import math
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")

def fmt(x):
    return "nan" if x != x else f"{x:.4g}"

def read_bench(path):
    with open(path, encoding="utf-8") as f:
        line1 = f.readline().strip()
        line2 = f.readline().strip()
        line3 = f.readline().strip()
    
    assert line1.startswith("# ")
    tokens = {}
    for token in line1[2:].split():
        k, v = token.split("=", 1)
        tokens[k] = v
    
    counts = list(map(int, line2[len("capture,"):].split(",")))
    return counts, tokens

def read_expected(prompt_path, isotope_path, element, Sn_keV):
    prompt = pd.read_excel(prompt_path, header=None, skiprows=2)
    prompt = prompt.iloc[:, [1, 5, 7]]
    prompt.columns = ["label", "E", "sigma_gamma"]
    prompt["label"] = prompt["label"].astype(str).str.strip()
    prompt["E"] = pd.to_numeric(prompt["E"], errors="coerce")
    prompt["sigma_gamma"] = pd.to_numeric(prompt["sigma_gamma"], errors="coerce")
    prompt = prompt.dropna(subset=["E", "sigma_gamma"])
    
    isotope = pd.read_excel(isotope_path, header=None, skiprows=2)
    isotope = isotope.iloc[:, [0, 1, 5, 8]]  # колонки сохраняют исходные номера 0, 1, 5, 8
    isotope[0] = isotope[0].ffill()
    isotope["label"] = isotope[1].astype(str).str.strip()
    isotope = isotope.dropna(subset=[0, 5, 8])
    isotope["abundance"] = pd.to_numeric(isotope[5], errors="coerce")
    isotope["sigma0"] = pd.to_numeric(isotope[8], errors="coerce")
    
    # Фильтруем по элементу
    prompt = prompt[prompt["label"].str.endswith(f"-{element}")]  # метка продукта: '<A>-<Символ>'
    
    # Строим словарь target_label -> (abundance, sigma0)
    isotope_dict = {}
    for _, row in isotope.iterrows():
        label = row["label"]
        symbol = row[0]
        ab = row["abundance"]
        sigma0 = row["sigma0"]
        # метка в isotope.xls уже ЕСТЬ метка мишени ('27-Al'): без пересчёта массового числа
        isotope_dict[label] = (ab, sigma0)
    
    # Собираем данные для каждого продукта
    expected = []
    for _, row in prompt.iterrows():
        label = row["label"]
        E = row["E"]
        sigma_gamma = row["sigma_gamma"]
        
        if "-" not in label:
            continue
            
        mass_str, sym = label.split("-", 1)
        try:
            Ap = int(mass_str.strip())
            target_label = f"{Ap-1}-{sym.strip()}"
        except:
            continue
            
        if target_label not in isotope_dict:
            continue
            
        ab, sigma0 = isotope_dict[target_label]
        
        # Вычисляем Y_el
        num = ab / 100 * sigma_gamma
        expected.append((E, num, sigma0, sigma_gamma))
    
    # Сортируем по энергии
    expected.sort(key=lambda x: x[0])
    
    # Объединяем близкие линии
    merged = []
    i = 0
    while i < len(expected):
        E0, num0, sigma00, sigma_gamma0 = expected[i]
        E_sum = E0 * num0
        num_sum = num0
        j = i + 1
        while j < len(expected) and expected[j][0] - E0 <= 0.3:
            E, num, sigma0, sigma_gamma = expected[j]
            E_sum += E * num
            num_sum += num
            j += 1
        i = j
        E_avg = E_sum / num_sum
        merged.append((E_avg, num_sum, sigma00, sigma_gamma0))
    
    # Вычисляем Y_el для каждой линии
    total_sigma0 = sum(ab / 100 * sigma0 for ab, sigma0 in isotope[isotope[0] == element][["abundance", "sigma0"]].values)
    
    lines = []
    for E, num, sigma0, sigma_gamma in merged:
        if not (50 <= E <= 11500):
            continue
        Y_exp = 100 * num / total_sigma0
        lines.append((E, Y_exp, sigma_gamma))
    
    lines.sort(key=lambda x: x[1], reverse=True)
    return lines

def line_stats(counts, E, Ncap, hw=None):
    if hw is None:
        hw = 3 + 0.0006 * E
    peak_start = int(max(0, E - hw))
    peak_end = int(min(len(counts), E + hw))
    bins = list(range(peak_start, peak_end))
    
    # Считаем sideband
    sb_start1 = max(0, peak_start - 8)
    sb_end1 = peak_start
    sb_start2 = peak_end
    sb_end2 = min(len(counts), peak_end + 8)
    
    sidebins = list(range(sb_start1, sb_end1)) + list(range(sb_start2, sb_end2))
    
    if not sidebins:
        baseline = 0
    else:
        baseline = sum(counts[i] for i in sidebins) / len(sidebins)
    
    peak_sum = sum(counts[i] for i in bins)
    area = peak_sum - baseline * len(bins)
    sigma = math.sqrt(max(peak_sum, 1))
    
    obs_per100 = 100 * area / Ncap if Ncap > 0 else 0
    
    # Центроид
    total_counts = sum(counts[i] for i in bins)
    centroid = 0
    if total_counts > 0:
        weighted_sum = sum((i + 0.5) * counts[i] for i in bins)
        centroid = weighted_sum / total_counts
    
    detected = area > 3 * sigma
    
    return {
        "area": area,
        "sigma": sigma,
        "peak_sum": peak_sum,
        "obs_per100": obs_per100,
        "centroid": centroid,
        "detected": detected
    }

def main():
    if len(sys.argv) < 5:
        print("Usage: python yield_lines.py <bench_csv> <element> <Sn_keV> <data_dir> [top_n]", file=sys.stderr)
        sys.exit(2)
    
    bench_csv = sys.argv[1]
    element = sys.argv[2]
    Sn_keV = float(sys.argv[3])
    data_dir = sys.argv[4]
    top_n = int(sys.argv[5]) if len(sys.argv) > 5 else 10
    
    counts, tokens = read_bench(bench_csv)
    
    Ncap = sum((j + 0.5) * counts[j] for j in range(len(counts))) / Sn_keV
    
    if Ncap <= 0:
        print("Error: Ncap <= 0", file=sys.stderr)
        sys.exit(2)
    
    prompt_path = f"{data_dir}/promptgammas.xls"
    isotope_path = f"{data_dir}/isotope.xls"
    
    expected_lines = read_expected(prompt_path, isotope_path, element, Sn_keV)
    
    if not expected_lines:
        print("Error: no expected lines", file=sys.stderr)
        sys.exit(2)
    
    expected_lines = expected_lines[:top_n]
    
    # Выводим заголовок
    hp_mode = tokens.get("hp_mode", "unknown")
    radius_cm = tokens.get("radius_cm", "unknown")
    n_capture_gamma = tokens.get("n_capture_gamma", "unknown")
    print(f"# file={bench_csv} element={element} hp_mode={hp_mode} n_capture_gamma={n_capture_gamma} Ncap_est={Ncap:.6g}")
    
    results = []
    for E, Y_exp, sigma_gamma in expected_lines:
        stats = line_stats(counts, E, Ncap)
        obs_over_exp = stats["obs_per100"] / Y_exp if Y_exp > 0 else float("nan")
        dE = stats["centroid"] - E if not math.isnan(stats["centroid"]) else float("nan")
        print(f"{E:.4g} | {Y_exp:.4g} | {stats['area']:.4g} | {stats['obs_per100']:.4g} | {stats['sigma']:.4g} | {fmt(obs_over_exp)} | {fmt(dE)}")
        results.append((stats["obs_per100"], Y_exp, stats["area"], stats["sigma"]))
    
    # Подсчет VERDICT
    count_ok = 0
    ratios = []
    for obs, exp, area, sigma in results:
        if not (exp > 0):
            continue
        ratio = obs / exp
        if 0.5 <= ratio <= 2.0 and area > 3 * sigma:
            count_ok += 1
        if area > 3 * sigma:
            ratios.append(ratio)
    
    median_ratio = float("nan") if not ratios else sorted(ratios)[len(ratios)//2]
    
    print(f"VERDICT lines={len(expected_lines)} ok={count_ok} median_ratio={median_ratio:.4g}")

if __name__ == "__main__":
    main()
