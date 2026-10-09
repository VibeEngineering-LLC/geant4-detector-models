import sys
sys.stdout.reconfigure(encoding='utf-8')
import os
import glob
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Constants
GRID_DIRS = [
    'C:/g4work/gs2020/run_marinelli/out_v5',
    'C:/g4work/gs2020/run_marinelli/out_v5_g16_55'
]
REF_FILE = os.environ.get('GS2020_BECQMONI_EFF_CSV', 'efficiency_becqmoni_reference.csv')   # кривая BecqMoni для сравнения — локальный референс, не в репозитории
OUTDIR = 'results_eff_vs_becqmoni'
KEDGE = 33.169

def parse_grid_file(filepath):
    """
    Parses a grid CSV file.
    Returns dict with keys: E, N, counts, sample_matrix, sample_rho_g_cm3
    or None if invalid.
    """
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except Exception:
        return None

    header = {}
    data_started = False
    bins = []
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('#'):
            continue
        
        if not data_started:
            if line.startswith('bin_keV'):
                data_started = True
                continue
            # Parse key,value
            parts = line.split(',')
            if len(parts) == 2:
                key = parts[0].strip()
                val = parts[1].strip()
                header[key] = val
        else:
            # Data row: center,count_edep,count_light
            parts = line.split(',')
            if len(parts) >= 2:
                try:
                    center = float(parts[0])
                    count_edep = float(parts[1])
                    bins.append((center, count_edep))
                except ValueError:
                    continue

    if not data_started or 'energy_keV' not in header or 'n_events_processed' not in header:
        return None

    try:
        E = float(header['energy_keV'])
        N = float(header['n_events_processed'])
    except ValueError:
        return None

    if E < 5:
        return None

    # Find bin containing E
    # Bins are 1 keV wide, centered at x.5. Bin [floor, floor+1) has center floor+0.5
    # We need count_edep of the bin containing E.
    # If E is within 0.02 of integer boundary, add neighbor.
    
    # Sort bins by center to ensure order (though likely already sorted)
    bins.sort(key=lambda x: x[0])
    
    # Map center to count for quick lookup? Or just iterate.
    # Since bins are 1 keV, we can calculate expected center.
    # But safer to search.
    
    # Determine the bin index.
    # Bin center c corresponds to interval [c-0.5, c+0.5).
    # We want bin where c-0.5 <= E < c+0.5.
    # So c is in [E-0.5, E+0.5).
    # Since centers are x.5, the center is round(E - 0.5) + 0.5?
    # Let's just find the bin.
    
    counts = 0.0
    found_bin = False
    
    # Efficient search: binary search or linear. Linear is fine for small lists, but let's be safe.
    # Actually, we can just compute the expected center.
    # If E = 583.187, bin is [583, 584), center 583.5.
    # If E = 583.0, bin is [583, 584) or [582, 583)?
    # "bin [floor, floor+1)" -> floor(583.0) = 583. Bin [583, 584). Center 583.5.
    # "if E is within 0.02 keV of an integer boundary... ALSO add the count_edep of the bin on the other side"
    
    # Let's find the primary bin center.
    # floor(E) gives the lower bound. Center = floor(E) + 0.5.
    # Check if this center exists in bins.
    
    # However, floating point issues might make exact match hard.
    # Let's iterate through bins to find the one containing E.
    # Given the structure, there should be exactly one bin containing E (or two if on boundary).
    
    # Let's use a dictionary for O(1) lookup if possible, or just linear scan.
    # Linear scan is safer against missing bins.
    
    # Identify primary bin
    # The bin containing E is the one where center - 0.5 <= E < center + 0.5
    # Or simply: center = floor(E) + 0.5
    
    # Let's collect all bins that might be relevant.
    # Primary bin center:
    c_primary = math.floor(E) + 0.5
    
    # Check boundary condition
    dist_to_int = abs(E - round(E))
    is_boundary = dist_to_int < 0.02
    
    # Find count for c_primary
    count_primary = 0.0
    found_primary = False
    for c, cnt in bins:
        if abs(c - c_primary) < 1e-9:
            count_primary = cnt
            found_primary = True
            break
            
    if not found_primary:
        # Fallback: find bin containing E
        for c, cnt in bins:
            if c - 0.5 <= E < c + 0.5:
                count_primary = cnt
                found_primary = True
                break
        if not found_primary:
            return None # Bin not found
            
    counts = count_primary
    
    if is_boundary:
        # Add neighbor
        # If E is close to integer K, the other bin is on the other side of K.
        # If E < K, primary bin is [K-1, K), neighbor is [K, K+1).
        # If E > K, primary bin is [K, K+1), neighbor is [K-1, K).
        # Basically, the neighbor center is c_primary +/- 1.
        
        # Determine which side E is on relative to the integer boundary.
        # If E is slightly below integer, primary is lower bin. Neighbor is upper.
        # If E is slightly above integer, primary is upper bin. Neighbor is lower.
        
        # Let's just check both c_primary - 1 and c_primary + 1 and add the one that exists?
        # The prompt says "add the count_edep of the bin on the other side of that boundary".
        # This implies exactly one neighbor.
        
        # If E is close to integer K:
        # If E < K (e.g. 33.169 is not close to 33, but say 33.001), primary bin is [33, 34)? No.
        # floor(33.001) = 33. Bin [33, 34). Center 33.5.
        # Boundary is 33. Other side is [32, 33). Center 32.5.
        
        # If E = 33.999. floor(33.999) = 33. Bin [33, 34). Center 33.5.
        # Boundary is 34. Other side is [34, 35). Center 34.5.
        
        # So we need to know which integer boundary we are close to.
        K = round(E)
        
        # If E < K, we are in bin [K-1, K). Neighbor is [K, K+1). Center K+0.5.
        # If E > K, we are in bin [K, K+1). Neighbor is [K-1, K). Center K-0.5.
        
        if E < K:
            c_neighbor = K + 0.5
        else:
            c_neighbor = K - 0.5
            
        count_neighbor = 0.0
        for c, cnt in bins:
            if abs(c - c_neighbor) < 1e-9:
                count_neighbor = cnt
                break
                
        counts += count_neighbor

    return {
        'E': E,
        'N': N,
        'counts': counts,
        'sample_matrix': header.get('sample_matrix', 'Unknown'),
        'sample_rho_g_cm3': header.get('sample_rho_g_cm3', 'Unknown')
    }

def load_grid_data():
    data = {}
    sets = set()
    
    for d in GRID_DIRS:
        pattern = os.path.join(d, 'grid_mar_E*.csv')
        files = glob.glob(pattern)
        for f in files:
            parsed = parse_grid_file(f)
            if parsed:
                E_key = round(parsed['E'], 3)
                sets.add((parsed['sample_matrix'], parsed['sample_rho_g_cm3']))
                
                if E_key in data:
                    # Keep one with larger N
                    if parsed['N'] > data[E_key]['N']:
                        data[E_key] = parsed
                else:
                    data[E_key] = parsed
                    
    return data, sets

def load_reference():
    try:
        with open(REF_FILE, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except Exception:
        return None, None

    E_ref = []
    eps_ref = []
    err_ref = []
    
    for line in lines[1:]: # Skip header
        line = line.strip()
        if not line:
            continue
        parts = line.split(',')
        if len(parts) >= 3:
            try:
                E_ref.append(float(parts[0]))
                eps_ref.append(float(parts[1]))
                err_ref.append(float(parts[2]))
            except ValueError:
                continue
                
    if not E_ref:
        return None, None
        
    # Sort by E
    indices = np.argsort(E_ref)
    E_ref = np.array(E_ref)[indices]
    eps_ref = np.array(eps_ref)[indices]
    err_ref = np.array(err_ref)[indices]
    
    return E_ref, eps_ref, err_ref

def get_ref_value(E, E_ref, eps_ref, err_ref):
    """
    Log-log linear interpolation.
    Returns (eps, err_nearest)
    """
    if E < E_ref[0] or E > E_ref[-1]:
        return np.nan, np.nan
        
    # Find indices i such that E_ref[i] <= E <= E_ref[i+1]
    # Use searchsorted
    idx = np.searchsorted(E_ref, E, side='right') - 1
    
    if idx < 0:
        idx = 0
    if idx >= len(E_ref) - 1:
        idx = len(E_ref) - 2
        
    E0 = E_ref[idx]
    E1 = E_ref[idx+1]
    eps0 = eps_ref[idx]
    eps1 = eps_ref[idx+1]
    
    # Check discontinuity
    if E0 < KEDGE < E1:
        return np.nan, np.nan
        
    if eps0 <= 0 or eps1 <= 0:
        return np.nan, np.nan
        
    # Log-log interpolation
    # log(eps) = log(eps0) + (log(E) - log(E0)) * (log(eps1) - log(eps0)) / (log(E1) - log(E0))
    log_E = np.log(E)
    log_E0 = np.log(E0)
    log_E1 = np.log(E1)
    log_eps0 = np.log(eps0)
    log_eps1 = np.log(eps1)
    
    slope = (log_eps1 - log_eps0) / (log_E1 - log_E0)
    log_eps = log_eps0 + (log_E - log_E0) * slope
    eps = np.exp(log_eps)
    
    # Nearest reference point for error
    # "err_ref_pct = reference err_pct of the nearest reference point"
    # Nearest in E?
    if abs(E - E0) < abs(E - E1):
        err_nearest = err_ref[idx]
    else:
        err_nearest = err_ref[idx+1]
        
    return eps, err_nearest

def main():
    # Load data
    grid_data, grid_sets = load_grid_data()
    E_ref, eps_ref, err_ref = load_reference()
    
    if not grid_data:
        print("No grid data found.")
        return
        
    if E_ref is None:
        print("Reference data not found.")
        return

    # 1) Print matrix/density set
    print(f"матрица сетки: {grid_sets}")

    # Prepare results
    results = []
    for E_key in sorted(grid_data.keys()):
        d = grid_data[E_key]
        E = d['E']
        N = d['N']
        counts = d['counts']
        
        eps_our = counts / N if N > 0 else np.nan
        err_our_pct = 100.0 / np.sqrt(counts) if counts > 0 else np.nan
        
        eps_ref_val, err_ref_val = get_ref_value(E, E_ref, eps_ref, err_ref)
        
        ratio = eps_our / eps_ref_val if (eps_ref_val > 0 and not np.isnan(eps_ref_val)) else np.nan
        
        results.append({
            'E': E,
            'eps_our': eps_our,
            'err_our': err_our_pct,
            'eps_ref': eps_ref_val,
            'err_ref': err_ref_val,
            'ratio': ratio
        })

    # 2) Print table
    # Format: E_keV;eps_our;err_our_pct;eps_ref;err_ref_pct;our/ref
    # 4 significant digits for values, 3 decimals for ratio
    print("E_keV;eps_our;err_our_pct;eps_ref;err_ref_pct;our/ref")
    for r in results:
        def fmt_sig(x, sig=4):
            if np.isnan(x):
                return "nan"
            return f"{x:.{sig}g}"
            
        def fmt_dec(x, dec=3):
            if np.isnan(x):
                return "nan"
            return f"{x:.{dec}f}"
            
        print(f"{fmt_sig(r['E'])};{fmt_sig(r['eps_our'])};{fmt_sig(r['err_our'])};{fmt_sig(r['eps_ref'])};{fmt_sig(r['err_ref'])};{fmt_dec(r['ratio'])}")

    # 3) Write CSV
    os.makedirs(OUTDIR, exist_ok=True)
    csv_path = os.path.join(OUTDIR, 'fep_our_vs_becqmoni.csv')
    with open(csv_path, 'w', encoding='utf-8') as f:
        f.write("E_keV,eps_our,err_our_pct,eps_ref,err_ref_pct,ratio\n")
        for r in results:
            def fmt_csv(x):
                if np.isnan(x):
                    return ""
                return f"{x:.6g}"
            f.write(f"{fmt_csv(r['E'])},{fmt_csv(r['eps_our'])},{fmt_csv(r['err_our'])},{fmt_csv(r['eps_ref'])},{fmt_csv(r['err_ref'])},{fmt_csv(r['ratio'])}\n")

    # 4) K-edge step
    # Our grid
    E_ours = [r['E'] for r in results]
    eps_ours = [r['eps_our'] for r in results]
    
    # Largest energy below KEDGE
    below = [r for r in results if r['E'] < KEDGE and not np.isnan(r['eps_our'])]
    above = [r for r in results if r['E'] > KEDGE and not np.isnan(r['eps_our'])]
    
    if below and above:
        E_lo_our = max(r['E'] for r in below)
        E_hi_our = min(r['E'] for r in above)
        
        eps_lo_our = next(r['eps_our'] for r in below if r['E'] == E_lo_our)
        eps_hi_our = next(r['eps_our'] for r in above if r['E'] == E_hi_our)
        
        step_our = eps_hi_our / eps_lo_our
        print(f"наша ступенька: eps(E_hi)/eps(E_lo) = {step_our:.4f} ({E_lo_our:.3f}, {E_hi_our:.3f})")
    else:
        print("наша ступенька: данные отсутствуют")

    # Reference
    # Last point below KEDGE and first above
    ref_below = [(E, eps) for E, eps in zip(E_ref, eps_ref) if E < KEDGE and eps > 0]
    ref_above = [(E, eps) for E, eps in zip(E_ref, eps_ref) if E > KEDGE and eps > 0]
    
    if ref_below and ref_above:
        E_lo_ref = max(E for E, _ in ref_below)
        E_hi_ref = min(E for E, _ in ref_above)
        
        eps_lo_ref = next(eps for E, eps in ref_below if E == E_lo_ref)
        eps_hi_ref = next(eps for E, eps in ref_above if E == E_hi_ref)
        
        step_ref = eps_hi_ref / eps_lo_ref
        print(f"ступенька BecqMoni: eps(E_hi)/eps(E_lo) = {step_ref:.4f} ({E_lo_ref:.3f}, {E_hi_ref:.3f})")
    else:
        print("ступенька BecqMoni: данные отсутствуют")

    # 5) Band summary
    bands = [(16, 33.17), (33.17, 60), (60, 150), (150, 500), (500, 1500), (1500, 3000)]
    print("\nСводка по диапазонам:")
    for lo, hi in bands:
        ratios = [r['ratio'] for r in results if lo <= r['E'] < hi and not np.isnan(r['ratio'])]
        if ratios:
            med = np.median(ratios)
            mn = np.min(ratios)
            mx = np.max(ratios)
            n = len(ratios)
            print(f"Диапазон {lo}-{hi} кэВ: медиана={med:.4f}, min={mn:.4f}, max={mx:.4f}, n={n}")
        else:
            print(f"Диапазон {lo}-{hi} кэВ: нет данных")

    # 6) Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 9), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    
    # Top panel: log-log eps vs E
    # Reference
    mask_ref = (E_ref >= 10) & (E_ref <= 3000)
    E_ref_plot = E_ref[mask_ref]
    eps_ref_plot = eps_ref[mask_ref]
    
    # Split reference at KEDGE to avoid connecting across discontinuity
    # Find indices around KEDGE
    idx_kedge = np.searchsorted(E_ref_plot, KEDGE)
    
    # Plot segments
    if idx_kedge > 0:
        ax1.plot(E_ref_plot[:idx_kedge], eps_ref_plot[:idx_kedge], 'o-', markersize=3, label='BecqMoni (efficiency.csv)', color='blue')
    if idx_kedge < len(E_ref_plot):
        ax1.plot(E_ref_plot[idx_kedge:], eps_ref_plot[idx_kedge:], 'o-', markersize=3, color='blue')
        
    # Ours
    E_our_plot = [r['E'] for r in results if 10 <= r['E'] <= 3000 and not np.isnan(r['eps_our'])]
    eps_our_plot = [r['eps_our'] for r in results if 10 <= r['E'] <= 3000 and not np.isnan(r['eps_our'])]
    err_our_plot = [r['err_our'] for r in results if 10 <= r['E'] <= 3000 and not np.isnan(r['eps_our'])]
    
    # Convert err_pct to absolute error for error bars
    # err_abs = eps * (err_pct / 100)
    err_abs_plot = [e * (err / 100.0) if not np.isnan(err) else 0 for e, err in zip(eps_our_plot, err_our_plot)]
    
    ax1.errorbar(E_our_plot, eps_our_plot, yerr=err_abs_plot, fmt='o', markersize=4, capsize=2, label='Geant4, сетка монолиний', color='red')
    
    # Vertical line at KEDGE
    ax1.axvline(KEDGE, color='black', linestyle='--', linewidth=1)
    ax1.text(KEDGE, ax1.get_ylim()[1] * 0.9, 'K-край I 33,17 кэВ', rotation=90, va='top', ha='left', fontsize=8)
    
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_ylabel('ЭПП (eps)')
    ax1.set_title('ЭПП Маринелли 1 л: Geant4 против BecqMoni')
    ax1.legend(loc='best')
    ax1.grid(True, which="both", ls="-", alpha=0.5)
    
    # Bottom panel: ratio
    E_ratio_plot = [r['E'] for r in results if 10 <= r['E'] <= 3000 and not np.isnan(r['ratio'])]
    ratio_plot = [r['ratio'] for r in results if 10 <= r['E'] <= 3000 and not np.isnan(r['ratio'])]
    
    ax2.plot(E_ratio_plot, ratio_plot, 'o-', markersize=4, color='green')
    ax2.axhline(1, color='black', linestyle='-', linewidth=1)
    ax2.set_ylim(0.5, 1.5)
    ax2.set_xlabel('Энергия (кэВ)')
    ax2.set_ylabel('Отношение (Geant4 / BecqMoni)')
    ax2.grid(True, which="both", ls="-", alpha=0.5)
    
    # Set x limits
    ax1.set_xlim(10, 3000)
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, 'fep_our_vs_becqmoni.png'), dpi=130)
    plt.close()

if __name__ == '__main__':
    main()
