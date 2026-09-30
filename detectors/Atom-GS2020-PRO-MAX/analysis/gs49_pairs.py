import sys
import os
import csv
import math

sys.stdout.reconfigure(encoding="utf-8")

def load_icc_table(filepath):
    """Load ICC correction table from CSV."""
    icc_data = []
    if not os.path.exists(filepath):
        return icc_data
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header_found = False
        for row in reader:
            if not row or row[0].startswith('#'):
                continue
            if not header_found:
                # Skip header line
                header_found = True
                continue
            try:
                dz = int(row[0])
                e_tr = float(row[1])
                ratio = float(row[2])
                icc_data.append((dz, e_tr, ratio))
            except (ValueError, IndexError):
                continue
    return icc_data

def get_daughter_z(nuclide):
    """Return list of daughter Zs for a given nuclide."""
    mapping = {
        'Ra228': [89],
        'Ac228': [90],
        'Th228': [88],
        'Ra224': [86],
        'Rn220': [84],
        'Pb212': [83],
        'Bi212': [84, 81], # Beta branch -> Po (84), Alpha branch -> Tl (81)
        'Tl208': [82]
    }
    return mapping.get(nuclide, [])

def get_k_windows(dz_list):
    """Return K X-ray energy windows for given daughter Zs."""
    windows_map = {
        81: (70.0, 90.0),
        82: (72.0, 90.0),
        83: (74.0, 93.0),
        84: (76.0, 96.0),
        86: (81.0, 101.0),
        88: (85.0, 107.0),
        89: (88.0, 110.0),
        90: (89.0, 112.0)
    }
    windows = []
    for dz in dz_list:
        if dz in windows_map:
            windows.append(windows_map[dz])
    return windows

def is_k_xray(energy, windows):
    """Check if energy falls into any K X-ray window."""
    for lo, hi in windows:
        if lo <= energy <= hi:
            return True
    return False

def get_binding_energy(z):
    """Return K binding energy for Z."""
    bk_map = {
        81: 85.530,
        82: 88.0045,
        83: 90.526,
        84: 93.105,
        86: 98.404,
        88: 103.922,
        89: 106.755,
        90: 109.651
    }
    return bk_map.get(z)

def find_icc_correction(electrons, dz_list, icc_table):
    """Find ICC correction weight for an event."""
    # If no electrons, w=1.0
    if not electrons:
        return 1.0
    
    best_match = None
    max_deviation = -1.0
    
    for e_kin in electrons:
        for dz in dz_list:
            bk = get_binding_energy(dz)
            if bk is None:
                continue
            e_tr = e_kin + bk
            
            # Find closest row in ICC table for this dz and e_tr
            candidates = []
            for row_dz, row_e_tr, row_ratio in icc_table:
                if row_dz == dz and abs(row_e_tr - e_tr) <= 0.6:
                    diff = abs(row_e_tr - e_tr)
                    candidates.append((diff, row_ratio))
            
            if candidates:
                # Sort by difference (closest first), then pick the one with largest |ratio-1| among closest?
                # "take the closest row" -> min diff. If several matches (same diff?), use largest |ratio-1|.
                # Actually, "If several matches, use the one with the largest |ratio-1|" likely refers to multiple rows within tolerance.
                # Let's find all within tolerance, pick closest. If ties in closeness, pick max deviation.
                candidates.sort(key=lambda x: (x[0], -abs(x[1]-1)))
                best_diff, best_ratio = candidates[0]
                
                dev = abs(best_ratio - 1.0)
                if dev > max_deviation:
                    max_deviation = dev
                    best_match = best_ratio
    
    return best_match if best_match is not None else 1.0

def parse_tally_file(filepath, icc_table):
    """Parse a single tally file and return accumulated pair counts."""
    raw_counts = {} # (x_energy, g_energy) -> count
    corr_counts = {} # (x_energy, g_energy) -> weighted_count
    
    n_events = 0
    data_lines = []
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Parse metadata
                if ',' in line and '=' not in line:
                    parts = line.split(',')
                    if len(parts) == 2:
                        key, val = parts[0], parts[1]
                        if key == 'n_events':
                            n_events = int(val)
                
                # Parse data lines
                if line.startswith('count;signature'):
                    continue
                
                # Data line format: <count>;P:<e1> <e2> ...;E:<e1> <e2> ...
                parts = line.split(';')
                if len(parts) == 3:
                    try:
                        count = int(parts[0])
                        p_str = parts[1]
                        e_str = parts[2]
                        
                        # Parse photons
                        photons = []
                        if p_str.startswith('P:'):
                            p_vals = p_str[2:].strip().split()
                            for pv in p_vals:
                                if pv:
                                    photons.append(float(pv))
                        
                        # Parse electrons
                        electrons = []
                        if e_str.startswith('E:'):
                            e_vals = e_str[2:].strip().split()
                            for ev in e_vals:
                                if ev:
                                    electrons.append(float(ev))
                        
                        data_lines.append((count, photons, electrons))
                    except ValueError:
                        continue
    except FileNotFoundError:
        return None, 0
    
    # Determine daughter Zs and K windows based on filename
    basename = os.path.basename(filepath)
    nuclide = basename.replace('tally_', '').replace('.txt', '')
    dz_list = get_daughter_z(nuclide)
    k_windows = get_k_windows(dz_list)
    
    if not dz_list:
        return None, 0
    
    # Process each data line
    for count, photons, electrons in data_lines:
        # Classify photons
        x_energies = []
        g_energies = []
        
        for p in photons:
            if is_k_xray(p, k_windows):
                x_energies.append(p)
            elif p >= 100.0:
                g_energies.append(p)
            # else ignored
        
        # If no X-ray or no gamma, skip
        if not x_energies or not g_energies:
            continue
        
        # Get distinct energies
        unique_x = sorted(list(set(x_energies)))
        unique_g = sorted(list(set(g_energies)))
        
        # Calculate correction weight
        w = find_icc_correction(electrons, dz_list, icc_table)
        
        # Accumulate pairs
        for x in unique_x:
            for g in unique_g:
                key = (x, g)
                raw_counts[key] = raw_counts.get(key, 0) + count
                corr_counts[key] = corr_counts.get(key, 0.0) + count * w
    
    return (raw_counts, corr_counts), n_events

def main():
    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    input_base = "C:/g4work/gs2020/gs49"
    output_dir = os.path.join(script_dir, "..", "audit")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "gs49_pairs.csv")
    
    # Load ICC table
    icc_path = os.path.join(script_dir, "data", "gs49_icc_ratio_lnhb_over_g4.csv")
    icc_table = load_icc_table(icc_path)
    
    nuclides = ['Ra228', 'Ac228', 'Th228', 'Ra224', 'Rn220', 'Pb212', 'Bi212', 'Tl208']
    
    all_rows = []
    summary_info = {} # nuclide -> {N, num_rows, top_10, xray_fraction}
    
    for nuc in nuclides:
        filepath = os.path.join(input_base, f"tally_{nuc}.txt")
        if not os.path.exists(filepath):
            continue
        
        result, n_events = parse_tally_file(filepath, icc_table)
        if result is None:
            continue
            
        raw_counts, corr_counts = result
        
        # Calculate X-ray fraction for this nuclide
        # We need to re-parse or store event info. Let's re-parse just for fraction count? 
        # Or better, modify parse_tally_file to return xray_event_count too.
        # For simplicity, let's do a quick pass again or store it.
        # Actually, I can compute it during the main loop if I expose it.
        # Let's refactor slightly: parse_tally_file returns (raw, corr, n_events_with_xray)
        
        # Re-implementing fraction calculation inside here for clarity without changing function signature too much
        # Or just trust the logic: "total fraction of events with at least one K-X-ray"
        # I need to count how many events had at least one X-ray.
        # Let's modify parse_tally_file to return this count.
        
        pass

    # Refactoring parse_tally_file to return xray_event_count
    def parse_tally_file_v2(filepath, icc_table):
        raw_counts = {}
        corr_counts = {}
        n_events = 0
        data_lines = []
        xray_event_count = 0
        
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    if ',' in line and '=' not in line:
                        parts = line.split(',')
                        if len(parts) == 2:
                            key, val = parts[0], parts[1]
                            if key == 'n_events':
                                n_events = int(val)
                    
                    if line.startswith('count;signature'):
                        continue
                    
                    parts = line.split(';')
                    if len(parts) == 3:
                        try:
                            count = int(parts[0])
                            p_str = parts[1]
                            e_str = parts[2]
                            
                            photons = []
                            if p_str.startswith('P:'):
                                p_vals = p_str[2:].strip().split()
                                for pv in p_vals:
                                    if pv:
                                        photons.append(float(pv))
                            
                            electrons = []
                            if e_str.startswith('E:'):
                                e_vals = e_str[2:].strip().split()
                                for ev in e_vals:
                                    if ev:
                                        electrons.append(float(ev))
                            
                            data_lines.append((count, photons, electrons))
                        except ValueError:
                            continue
        except FileNotFoundError:
            return None, 0, 0
        
        basename = os.path.basename(filepath)
        nuclide = basename.replace('tally_', '').replace('.txt', '')
        dz_list = get_daughter_z(nuclide)
        k_windows = get_k_windows(dz_list)
        
        if not dz_list:
            return None, 0, 0
        
        for count, photons, electrons in data_lines:
            x_energies = []
            g_energies = []
            
            has_xray = False
            for p in photons:
                if is_k_xray(p, k_windows):
                    x_energies.append(p)
                    has_xray = True
                elif p >= 100.0:
                    g_energies.append(p)
            
            if has_xray:
                xray_event_count += count
            
            if not x_energies or not g_energies:
                continue
            
            unique_x = sorted(list(set(x_energies)))
            unique_g = sorted(list(set(g_energies)))
            
            w = find_icc_correction(electrons, dz_list, icc_table)
            
            for x in unique_x:
                for g in unique_g:
                    key = (x, g)
                    raw_counts[key] = raw_counts.get(key, 0) + count
                    corr_counts[key] = corr_counts.get(key, 0.0) + count * w
        
        return (raw_counts, corr_counts), n_events, xray_event_count

    # Main loop again with v2
    all_rows = []
    
    for nuc in nuclides:
        filepath = os.path.join(input_base, f"tally_{nuc}.txt")
        if not os.path.exists(filepath):
            continue
        
        result, n_events, xray_event_count = parse_tally_file_v2(filepath, icc_table)
        if result is None:
            continue
            
        raw_counts, corr_counts = result
        
        rows_for_nuc = []
        
        for (x, g), raw_val in raw_counts.items():
            if raw_val < 10:
                continue
            
            corr_val = corr_counts.get((x, g), 0.0)
            p_raw = raw_val / n_events
            p_raw_err = math.sqrt(raw_val) / n_events
            p_corr = corr_val / n_events
            e_sum = round(x + g, 3)
            
            rows_for_nuc.append({
                'nuclide': nuc,
                'e_x_keV': x,
                'e_gamma_keV': g,
                'e_sum_keV': e_sum,
                'n_events': n_events,
                'p_raw': p_raw,
                'p_raw_err': p_raw_err,
                'p_corr': p_corr
            })
        
        # Sort by p_corr descending
        rows_for_nuc.sort(key=lambda r: r['p_corr'], reverse=True)
        
        all_rows.extend(rows_for_nuc)
        
        xray_frac = xray_event_count / n_events if n_events > 0 else 0.0
        
        summary_info[nuc] = {
            'N': n_events,
            'num_rows': len(rows_for_nuc),
            'top_10': rows_for_nuc[:10],
            'xray_frac': xray_frac
        }

    # Write CSV
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['nuclide', 'e_x_keV', 'e_gamma_keV', 'e_sum_keV', 'n_events', 'p_raw', 'p_raw_err', 'p_corr'])
        for row in all_rows:
            writer.writerow([
                row['nuclide'],
                f"{row['e_x_keV']:.3f}",
                f"{row['e_gamma_keV']:.3f}",
                f"{row['e_sum_keV']:.3f}",
                row['n_events'],
                f"{row['p_raw']:.6f}",
                f"{row['p_raw_err']:.6f}",
                f"{row['p_corr']:.6f}"
            ])

    # Print summary in Russian
    print("Результаты анализа совпадений K-рентген + гамма")
    print("-" * 50)
    
    for nuc in nuclides:
        if nuc not in summary_info:
            continue
        
        info = summary_info[nuc]
        print(f"\nНуклид: {nuc}")
        print(f"  Всего событий (N): {info['N']}")
        print(f"  Количество строк в таблице: {info['num_rows']}")
        print(f"  Доля событий с K-рентгеном: {info['xray_frac']:.4f}")
        
        print("  Топ-10 пар по p_corr:")
        for i, row in enumerate(info['top_10']):
            print(f"    {i+1}. X={row['e_x_keV']:.3f} keV, G={row['e_gamma_keV']:.3f} keV, p_corr={row['p_corr']:.6f}")

if __name__ == "__main__":
    main()
