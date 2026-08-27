"""
Compare normalised spectrum SHAPES from two or more Geant4 output CSV files that
were produced with different source-disk radii, to decide whether the detector
response has saturated with respect to that radius.
"""

import csv
import math
import sys

def parse_csv(filename):
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f]
    
    # Find block headers
    metric_block_start = None
    histogram_block_start = None
    
    for i, line in enumerate(lines):
        if not line:
            continue
        parts = line.split(',', 1)
        if len(parts) < 2:
            continue
        if parts[0] == 'metric':
            metric_block_start = i
        elif parts[0] == 'bin_keV':
            histogram_block_start = i
            break
    
    if metric_block_start is None:
        raise SystemExit(f"Error: No metric block found in {filename}")
    if histogram_block_start is None:
        raise SystemExit(f"Error: No histogram block found in {filename}")
    
    # Parse metrics
    metrics = {}
    i = metric_block_start + 1
    while i < len(lines):
        line = lines[i]
        if not line or line.startswith('bin_keV'):
            break
        parts = line.split(',', 1)
        if len(parts) == 2:
            key, value = parts
            try:
                metrics[key] = float(value)
            except ValueError:
                metrics[key] = value
        i += 1
    
    # Parse histogram
    counts = {}
    per_muon = {}
    i = histogram_block_start + 1
    while i < len(lines):
        line = lines[i]
        if not line:
            break
        parts = line.split(',')
        if len(parts) >= 3:
            try:
                bin_idx = int(parts[0])
                counts[bin_idx] = int(parts[1])
                per_muon[bin_idx] = float(parts[2])
            except ValueError:
                pass
        i += 1
    
    return metrics, counts, per_muon

def compute_band_stats(counts, per_muon, start, end):
    total_counts = sum(counts.get(i, 0) for i in range(3000))
    if total_counts <= 0:
        return float('nan'), float('nan')
    
    band_counts = sum(counts.get(i, 0) for i in range(start, end + 1))
    band_per_muon = sum(per_muon.get(i, 0.0) for i in range(start, end + 1))
    
    if band_per_muon <= 0:
        frac = 0.0
    else:
        frac = band_per_muon / sum(per_muon.get(i, 0.0) for i in range(3000))
    
    if band_counts <= 0 or total_counts <= 0:
        sigma_rel = float('nan')
    else:
        f = band_counts / total_counts
        sigma_rel = math.sqrt((1 - f) / band_counts)
    
    return frac, sigma_rel

def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    
    ref_file = sys.argv[1]
    test_files = sys.argv[2:]
    
    # Parse reference file
    ref_metrics, ref_counts, ref_per_muon = parse_csv(ref_file)
    ref_r_disk = ref_metrics['r_disk_mm']
    
    # Parse test files
    test_data = []
    for test_file in test_files:
        metrics, counts, per_muon = parse_csv(test_file)
        test_data.append((metrics, counts, per_muon, test_file))
    
    # Bands
    bands = [
        ("700-1500", 700, 1500),
        ("1500-2400", 1500, 2400),
        ("2400-2999", 2400, 2999)
    ]
    
    # Print RUNS section
    print("=== RUNS ===")
    print(f"{'r_disk_mm':<10} {'n_events':<12} {'n_hits_in_crystal':<20} {'n_overflow':<12} {'sum_counts':<12} {'disk_area_cm2':<15} {'pdg_expected_per_s':<20} {'file_path'}")
    print("-" * 120)
    
    # Reference
    sum_counts_ref = sum(ref_counts.values())
    print(f"{ref_r_disk:<10} {int(ref_metrics['n_events']):<12} {int(ref_metrics['n_hits_in_crystal']):<20} {int(ref_metrics['n_overflow']):<12} {sum_counts_ref:<12} {ref_metrics['disk_area_cm2']:<15} {ref_metrics['pdg_expected_per_s']:<20} {ref_file}")
    
    # Tests
    for metrics, counts, per_muon, test_file in test_data:
        sum_counts = sum(counts.values())
        r_disk = metrics['r_disk_mm']
        n_events = int(metrics['n_events'])
        n_hits_in_crystal = int(metrics['n_hits_in_crystal'])
        n_overflow = int(metrics['n_overflow'])
        disk_area_cm2 = metrics['disk_area_cm2']
        pdg_expected_per_s = metrics['pdg_expected_per_s']
        print(f"{r_disk:<10} {n_events:<12} {n_hits_in_crystal:<20} {n_overflow:<12} {sum_counts:<12} {disk_area_cm2:<15} {pdg_expected_per_s:<20} {test_file}")
    
    # Print BAND SHAPE COMPARISON section
    print("\n=== BAND SHAPE COMPARISON ===")
    print(f"{'band_keV':<12} {'rdisk_mm':<10} {'frac_ref':<15} {'frac_test':<15} {'delta_pct':<12} {'sigma_pct':<12}")
    print("-" * 70)
    
    max_delta = 0.0
    
    for band_name, start, end in bands:
        ref_frac, ref_sigma_rel = compute_band_stats(ref_counts, ref_per_muon, start, end)
        
        for metrics, counts, per_muon, test_file in test_data:
            test_frac, test_sigma_rel = compute_band_stats(counts, per_muon, start, end)
            
            if math.isnan(ref_frac) or math.isnan(test_frac):
                delta_pct = float('nan')
            else:
                if ref_frac == 0.0:
                    delta_pct = float('nan')
                else:
                    delta_pct = (test_frac - ref_frac) / ref_frac * 100.0
            
            if math.isnan(ref_sigma_rel) or math.isnan(test_sigma_rel):
                sigma_pct = float('nan')
            else:
                sigma_pct = math.sqrt(ref_sigma_rel**2 + test_sigma_rel**2) * 100.0
            
            r_disk = metrics['r_disk_mm']
            
            print(f"{band_name:<12} {r_disk:<10} {ref_frac:<15.8f} {test_frac:<15.8f} {delta_pct:<12.3f} {sigma_pct:<12.3f}")
            
            if not math.isnan(delta_pct):
                max_delta = max(max_delta, abs(delta_pct))
    
    # Print INTERPRETATION section
    print("\n=== INTERPRETATION ===")
    
    for band_name, start, end in bands:
        ref_frac, _ = compute_band_stats(ref_counts, ref_per_muon, start, end)
        
        for metrics, counts, per_muon, test_file in test_data:
            test_frac, _ = compute_band_stats(counts, per_muon, start, end)
            
            if math.isnan(ref_frac) or math.isnan(test_frac):
                verdict = "EMPTY BAND - cannot compare"
                delta_pct = float('nan')
                sigma_pct = float('nan')
            else:
                if ref_frac == 0.0:
                    verdict = "EMPTY BAND - cannot compare"
                    delta_pct = float('nan')
                    sigma_pct = float('nan')
                else:
                    delta_pct = (test_frac - ref_frac) / ref_frac * 100.0
                    _, ref_sigma_rel = compute_band_stats(ref_counts, ref_per_muon, start, end)
                    _, test_sigma_rel = compute_band_stats(counts, per_muon, start, end)
                    
                    if math.isnan(ref_sigma_rel) or math.isnan(test_sigma_rel):
                        sigma_pct = float('nan')
                    else:
                        sigma_pct = math.sqrt(ref_sigma_rel**2 + test_sigma_rel**2) * 100.0
                    
                    if math.isnan(delta_pct):
                        verdict = "EMPTY BAND - cannot compare"
                    elif abs(delta_pct) > 5.0:
                        verdict = "ABOVE 5 PERCENT - NOT saturated"
                    elif abs(delta_pct) <= sigma_pct:
                        verdict = "within statistical error - no difference detected"
                    else:
                        verdict = "below 5 percent but above its own statistical error - small systematic shift"
            
            r_disk = metrics['r_disk_mm']
            band_counts_ref = sum(ref_counts.get(i, 0) for i in range(start, end + 1))
            band_counts_test = sum(counts.get(i, 0) for i in range(start, end + 1))
            
            print(f"{test_file} {band_name} {verdict}")
            print(f"  delta_pct={delta_pct:.3f}, sigma_pct={sigma_pct:.3f}, ref_count={band_counts_ref}, test_count={band_counts_test}")

    # Final line
    print(f"\nmaximum absolute delta_pct = {max_delta:.3f}")
    print("criterion: max |delta| > 5 percent => NOT saturated")

if __name__ == "__main__":
    main()
