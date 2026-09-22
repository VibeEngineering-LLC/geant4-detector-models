import sys
import math
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")

def read_bench(path):
    with open(path, encoding="utf-8") as f:
        line1 = f.readline().strip()
        line2 = f.readline().strip()
        line3 = f.readline().strip()
    
    if not line1.startswith("#"):
        raise ValueError("Invalid header")
    
    header = {}
    for token in line1[1:].split():
        k, v = token.split("=", 1)
        header[k] = v
    
    if not line2.startswith("capture,"):
        raise ValueError("Invalid capture line")
    
    counts = list(map(int, line2[8:].split(",")))
    
    return header, counts

def read_expected(data_dir, symbols, top_n):
    prompt_path = f"{data_dir}/promptgammas.xls"
    isotope_path = f"{data_dir}/isotope.xls"
    
    df_prompt = pd.read_excel(prompt_path, header=None, skiprows=2)
    df_isotope = pd.read_excel(isotope_path, header=None, skiprows=2)
    
    df_prompt[5] = pd.to_numeric(df_prompt[5], errors="coerce")
    df_prompt[7] = pd.to_numeric(df_prompt[7], errors="coerce")
    
    df_prompt = df_prompt.dropna(subset=[5, 7])
    
    df_isotope[0] = df_isotope[0].ffill()
    df_isotope[5] = pd.to_numeric(df_isotope[5], errors="coerce")
    df_isotope = df_isotope.dropna(subset=[5])
    
    # метки в xls с пробелами ('  1-H '): без strip элементы не находились и ожидаемых линий было ноль
    isotope_abundances = dict(zip(df_isotope[1].astype(str).str.strip(), df_isotope[5]))
    
    expected_lines = {}
    
    for symbol in symbols:
        expected_lines[symbol] = []
    
    for _, row in df_prompt.iterrows():
        product_label = str(row[1]).strip()
        energy = row[5]
        sigma_gamma = row[7]
        
        if not isinstance(energy, (int, float)) or not isinstance(sigma_gamma, (int, float)):
            continue
            
        if math.isnan(energy) or math.isnan(sigma_gamma):
            continue
        
        if energy < 30 or energy > 11500:
            continue
            
        try:
            Ap, elem = product_label.split("-")
            Ap = int(Ap)
            target_label = f"{Ap-1}-{elem}"
        except:
            continue
            
        if elem not in symbols:
            continue
            
        abundance = isotope_abundances.get(target_label, None)
        if abundance is None:
            continue
            
        w = (abundance / 100) * sigma_gamma
        expected_lines[elem].append((energy, w))
    
    for symbol in symbols:
        expected_lines[symbol].sort(key=lambda x: x[1], reverse=True)
        expected_lines[symbol] = expected_lines[symbol][:top_n]
    
    return expected_lines

def merge_close_lines(lines, threshold=0.3):
    if not lines:
        return []
        
    merged = [lines[0]]
    
    for energy, weight in lines[1:]:
        last_energy, last_weight = merged[-1]
        if abs(energy - last_energy) < threshold:
            merged[-1] = ((last_energy * last_weight + energy * weight) / (last_weight + weight), last_weight + weight)
        else:
            merged.append((energy, weight))
            
    return merged

def analyse(hist, expected):
    results = []
    
    for symbol, lines in expected.items():
        if not lines:
            continue
            
        lines = merge_close_lines(lines)
        
        ref_w = max(w for _, w in lines)
        
        for energy, w in lines:
            exp_rel = w / ref_w
            
            hw = 3 + 0.0006 * energy
            lo = energy - hw
            hi = energy + hw
            
            peak_bins = range(int(math.ceil(lo - 0.5)), int(math.floor(hi - 0.5)) + 1)
            peak_bins = [i for i in peak_bins if 0 <= i < len(hist)]
            
            if not peak_bins:
                continue
                
            peak_sum = sum(hist[i] for i in peak_bins)
            
            base_lo = energy - hw - 8
            base_hi = energy + hw + 8
            
            side_lo = range(int(math.ceil(base_lo - 0.5)), int(math.floor(energy - hw - 0.5)) + 1)
            side_hi = range(int(math.ceil(energy + hw + 0.5)), int(math.floor(base_hi - 0.5)) + 1)
            
            side_bins = [i for i in side_lo if 0 <= i < len(hist)] + [i for i in side_hi if 0 <= i < len(hist)]
            
            base_per_bin = 0
            if side_bins:
                base_sum = sum(hist[i] for i in side_bins)
                base_per_bin = base_sum / len(side_bins)
                
            area = peak_sum - base_per_bin * len(peak_bins)
            
            sigma_area = math.sqrt(peak_sum + (base_per_bin ** 2 * len(peak_bins) ** 2) / max(1, len(side_bins)))
            
            centroid = 0
            total_weighted = 0
            
            for i in peak_bins:
                bin_center = i + 0.5
                count = max(0, hist[i] - base_per_bin)
                total_weighted += count * bin_center
            if total_weighted > 0:
                centroid = total_weighted / sum(max(0, hist[i] - base_per_bin) for i in peak_bins)
            else:
                centroid = float('nan')
                
            results.append({
                "element": symbol,
                "E_exp_keV": energy,
                "exp_rel": exp_rel,
                "area": area,
                "sigma_area": sigma_area,
                "obs_rel": None,
                "obs_over_exp": None,
                "centroid_keV": centroid,
                "dE_keV": centroid - energy,
                "detected": area > 3 * sigma_area
            })
    
    return results

def main():
    if len(sys.argv) < 4:
        print("Usage: python analyze_lines.py <bench_csv> <elements> <data_dir> [top_n]", file=sys.stderr)
        sys.exit(2)
        
    bench_csv = sys.argv[1]
    elements = sys.argv[2].split(",")
    data_dir = sys.argv[3]
    top_n = int(sys.argv[4]) if len(sys.argv) > 4 else 12
    
    try:
        header, hist = read_bench(bench_csv)
        expected = read_expected(data_dir, elements, top_n)
        for s in elements:
            if not expected[s]:  # тихий выход без строк — это провал, а не «чисто»
                print(f"Error: для элемента {s} не найдено ни одной ожидаемой линии", file=sys.stderr)
                sys.exit(2)
        results = analyse(hist, expected)
        
        area_ref_by_element = {}
        for r in results:
            if r["exp_rel"] == 1.0:  # сильнейшая ОЖИДАЕМАЯ линия элемента — референс по элементу
                area_ref_by_element[r["element"]] = r["area"]

        for r in results:
            ref_area = area_ref_by_element.get(r["element"], 0)
            if ref_area > 0:
                r["obs_rel"] = r["area"] / ref_area
            else:
                r["obs_rel"] = float("nan")
                
            if r["exp_rel"] > 0:
                r["obs_over_exp"] = r["obs_rel"] / r["exp_rel"]
            else:
                r["obs_over_exp"] = float("nan")
        
        print(f"# file={bench_csv} elements={','.join(elements)} hp_mode={header.get('hp_mode', 'unknown')} n_capture_gamma={header.get('n_capture_gamma', 'unknown')}")
        
        for r in results:
            print(f"{r['element']} | {r['E_exp_keV']:.4g} | {r['exp_rel']:.4g} | {r['area']:.4g} | {r['sigma_area']:.4g} | {r['obs_rel']:.4g} | {r['obs_over_exp']:.4g} | {r['centroid_keV']:.4g} | {r['dE_keV']:.4g}")
        
        for symbol in elements:
            lines = [r for r in results if r["element"] == symbol]
            if not lines:
                continue
                
            detected_lines = [r for r in lines if r["detected"]]
            count_detected = len(detected_lines)
            
            ratios = [r["obs_over_exp"] for r in detected_lines if not math.isnan(r["obs_over_exp"])]
            median_ratio = float("nan")
            if ratios:
                sorted_ratios = sorted(ratios)
                n = len(sorted_ratios)
                if n % 2 == 1:
                    median_ratio = sorted_ratios[n // 2]
                else:
                    median_ratio = (sorted_ratios[n // 2 - 1] + sorted_ratios[n // 2]) / 2
                    
            print(f"SUMMARY {symbol} lines={len(lines)} detected={count_detected} median_obs_over_exp={median_ratio:.4g}")
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(2)

if __name__ == "__main__":
    main()
