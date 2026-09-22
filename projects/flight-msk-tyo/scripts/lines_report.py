import os
import sys, os, math, re
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bline_table import parse_lines
NEUTRON_CAPTURE_RE = re.compile(r"\(n[,;]\s*[γg]\)")   # #REM-1 fix: «выход на 100 захватов» осмыслен ТОЛЬКО для (n,γ) —
# не для (n,n'), распада, аннигиляции и т.п.; поле "yield" в parse_lines заполнено для ВСЕХ строк таблицы, поэтому
# фильтр обязателен (было найдено на 511 кэВ/аннигиляция и 1460,82 кэВ/распад K-40 — там подпись бессмысленна)

def fwhm_keV(E):
    if os.environ.get("LINES_FWHM"):   # режим сырого (несвёрнутого) спектра: постоянная ширина окна, кэВ
        return float(os.environ["LINES_FWHM"])
    return 41.6 * math.sqrt(E / 661.657)

def net_area(E, cph, sig, E0, fwhm):
    window = np.abs(E - E0) <= 1.2 * fwhm
    if not np.any(window):
        return (0.0, float('inf'), 0.0)
    
    baseline1 = (E >= E0 - 4 * fwhm) & (E <= E0 - 2.2 * fwhm)
    baseline2 = (E >= E0 + 2.2 * fwhm) & (E <= E0 + 4 * fwhm)
    
    if np.sum(baseline1) < 2 or np.sum(baseline2) < 2:
        return (0.0, float('inf'), 0.0)
    
    # база: парабола по логарифму на обоих флангах (учёт кривизны континуума), в окне линии
    fl = baseline1 | baseline2
    ok = fl & (cph > 0)
    if np.sum(ok) < 6:
        return (0.0, float('inf'), 0.0)
    for _ in range(3):   # обрезка выбросов сверху (хвосты соседних пиков во флангах)
        coef = np.polyfit(E[ok] - E0, np.log(cph[ok]), 2)
        res = np.log(cph) - np.polyval(coef, E - E0)
        keep = ok & (res < np.median(res[ok]) + 3.0 * 1.4826 * np.median(np.abs(res[ok] - np.median(res[ok]))) + 1e-9)
        if np.sum(keep) < 6 or np.array_equal(keep, ok):
            break
        ok = keep
    basev = np.exp(np.polyval(coef, E[window] - E0))
    nwin = np.sum(window)
    net = np.sum(cph[window] - basev)
    base = float(np.mean(basev))
    sig_win = np.sum(sig[window]**2)
    sig_base = np.sum(sig[baseline1]**2) + np.sum(sig[baseline2]**2)
    nb = np.sum(baseline1) + np.sum(baseline2)
    
    sigma = math.sqrt(sig_win + nwin**2 * sig_base / nb**2) if nb > 0 else float('inf')
    
    signif = net / sigma if sigma > 0 else 0.0
    
    return (net, sigma, signif)

def attribute(cat, E0, fwhm, top=3):
    total_excess = 0
    results = []
    
    for c in cat:
        bins = c["bins"]
        w = 0
        b = 0
        
        # Signal window
        j_start = int(E0 - 1.2 * fwhm)
        j_end = int(E0 + 1.2 * fwhm)
        if j_start < 0:
            j_start = 0
        if j_end >= len(bins):
            j_end = len(bins) - 1
            
        for j in range(j_start, j_end + 1):
            w += bins[j]
            
        # Baseline windows
        baseline1_start = int(E0 - 4 * fwhm)
        baseline1_end = int(E0 - 2.2 * fwhm)
        baseline2_start = int(E0 + 2.2 * fwhm)
        baseline2_end = int(E0 + 4 * fwhm)
        
        if baseline1_start < 0:
            baseline1_start = 0
        if baseline1_end >= len(bins):
            baseline1_end = len(bins) - 1
        if baseline2_start < 0:
            baseline2_start = 0
        if baseline2_end >= len(bins):
            baseline2_end = len(bins) - 1
            
        bl1_vals = bins[baseline1_start:baseline1_end+1]
        bl2_vals = bins[baseline2_start:baseline2_end+1]
        
        if len(bl1_vals) < 2 or len(bl2_vals) < 2:
            continue
            
        b = (np.mean(bl1_vals) + np.mean(bl2_vals)) / 2 * (j_end - j_start + 1)
        
        excess = max(w - b, 0)
        total_excess += excess
        
        results.append((c["names"], excess))
    
    if total_excess == 0:
        return ""
        
    results.sort(key=lambda x: x[1], reverse=True)
    top_results = results[:top]
    
    out = []
    for name, ex in top_results:
        share = (ex / total_excess) * 100
        out.append(f"{name[0]}/{name[1]}/{name[2]} {share:.1f}%")
        
    return "; ".join(out)

def find_unlisted(E, cph, sig, listed_E, minsig=4.0):
    peaks = []
    
    for i in range(20, len(cph) - 20):
        E0 = E[i]
        fwhm = fwhm_keV(E0)
        half = max(1, int(fwhm / 2))
        lo, hi = max(0, i - half), min(len(cph), i + half + 1)
        if cph[i] < np.max(cph[lo:hi]):   # максимум в окне ±fwhm/2 (равные соседи допускаются, близкие пики сливаются ниже)
            continue
        
        fl = ((E >= E0 - 4 * fwhm) & (E <= E0 - 2.2 * fwhm)) | ((E >= E0 + 2.2 * fwhm) & (E <= E0 + 4 * fwhm))
        if np.sum(fl) < 4 or cph[i] <= 1.01 * np.mean(cph[fl]):   # пик обязан подниматься над средним флангов
            continue
        # Check distance to listed lines
        min_dist = min([abs(E0 - e) for e in listed_E])
        if min_dist <= 1.0 * fwhm:
            continue
            
        net, sigma, signif = net_area(E, cph, sig, E0, fwhm)
        
        if signif < minsig or net <= 0:
            continue
            
        peaks.append((E0, net, sigma, signif))
    
    # Merge close peaks
    merged = []
    i = 0
    while i < len(peaks):
        E0, net, sigma, signif = peaks[i]
        j = i + 1
        while j < len(peaks) and abs(peaks[j][0] - E0) <= fwhm_keV(E0):
            if peaks[j][3] > signif:
                E0, net, sigma, signif = peaks[j]
            j += 1
        merged.append((E0, net, sigma, signif))
        i = j
        
    return merged

def main():
    if len(sys.argv) == 2 and sys.argv[1] == "--selftest":
        # Test 1: Synthetic spectrum with Gaussian at 1000 keV
        E = np.arange(20.5, 3001.5)
        cph = np.full_like(E, 100.0)
        
        # Add Gaussian at 1000 keV
        sigma = fwhm_keV(1000) / 2.35482
        gaussian = 5000 * np.exp(-0.5 * ((E - 1000) / sigma)**2) / (sigma * np.sqrt(2 * np.pi))
        cph += gaussian
        
        sig = np.sqrt(cph)
        
        # Test net_area at 1000 keV
        E0 = 1000.0
        fwhm = fwhm_keV(E0)
        net, sigma, signif = net_area(E, cph, sig, E0, fwhm)
        if abs(net - 5000) / 5000 > 0.03 or signif < 10:
            print("SELFTEST FAIL: net_area test")
            return 1
            
        # Test net_area at 2000 keV (no line)
        E0 = 2000.0
        fwhm = fwhm_keV(E0)
        net, sigma, signif = net_area(E, cph, sig, E0, fwhm)
        if abs(signif) >= 3:
            print("SELFTEST FAIL: no line test")
            return 1
            
        # Test 2: attribute function
        cat = [
            {
                "names": ("gamma", "nCapture", "Skin"),
                "bins": np.zeros(3001)
            },
            {
                "names": ("neutron", "hadElastic", "Pax"),
                "bins": np.full(3001, 0.01)
            }
        ]
        
        cat[0]["bins"][998:1003] = 0.5
        
        result = attribute(cat, 1000.0, fwhm_keV(1000.0), top=3)
        if not result.startswith("gamma/nCapture/Skin") or float(result.split("%")[0].split()[-1]) < 90.0:   # доля A > 90 % (фон категории B вычитается)
            print("SELFTEST FAIL: attribute test")
            return 1
            
        # Test 3: find_unlisted
        listed_E = [500.0]
        unlisted = find_unlisted(E, cph, sig, listed_E)
        if len(unlisted) != 1 or abs(unlisted[0][0] - 1000) >= 5:
            print("SELFTEST FAIL: find_unlisted test")
            return 1
            
        listed_E = [1000.0]
        unlisted = find_unlisted(E, cph, sig, listed_E)
        if len(unlisted) != 0:
            print("SELFTEST FAIL: find_unlisted test 2")
            return 1
            
        print("SELFTEST PASS")
        return 0
        
    if len(sys.argv) != 5:
        print("Usage: lines_report.py <smeared.csv> <cat.csv> <table.md> <out_prefix>")
        return 1
        
    smeared_path = sys.argv[1]
    cat_path = sys.argv[2]
    table_path = sys.argv[3]
    out_prefix = sys.argv[4]
    
    # Load spectrum
    data = np.loadtxt(smeared_path, delimiter=",", skiprows=1)
    E = data[:, 0]
    cph = data[:, 1]
    sig = data[:, 2]
    
    # Parse table
    lines = parse_lines(table_path)
    
    # Load categories
    cat = []
    with open(cat_path, "r", encoding="utf-8") as f:
        header = f.readline().strip().split(",")
        for line in f:
            parts = line.strip().split(",")
            pg, proc, vol = parts[3], parts[4], parts[5]
            name = (pg, proc, vol)
            bins = [float(x) for x in parts[7:]]
            cat.append({"names": name, "bins": np.array(bins)})
    
    # Process lines
    found_count = 0
    found_A_B = 0
    listed_count = 0
    
    results = []
    unlisted_results = []
    
    listed_E = []
    
    all_pk = find_unlisted(E, cph, sig, [-1e9])   # пик-первым: строка таблицы засчитывается только рядом с реально найденным пиком
    for line in lines:
        E_keV = line["E_keV"]
        if not (30 <= E_keV <= 9900):
            continue
            
        listed_count += 1
        listed_E.append(E_keV)
        
        fwhm = fwhm_keV(E_keV)
        net, sigma, signif = net_area(E, cph, sig, E_keV, fwhm)
        
        pk_near = [p for p in all_pk if abs(E_keV - p[0]) <= float(os.environ.get("LINES_MATCH", 0.6 * fwhm))]
        pk_near.sort(key=lambda p: abs(E_keV - p[0]))   # ближайший пик — первым
        found = 1 if signif >= float(os.environ.get("LINES_MINSIG", "3")) and net > 0 and pk_near else 0
        origin = attribute(cat, E_keV, fwhm, top=3).replace(",", ";") if found else ""
        if found:
            found_count += 1
            text = line["text"]
            cls = text.split("|")[-1].strip()
            if cls[:1] in ("A", "B"):
                found_A_B += 1
                
        # Prepare output fields
        yld = (line.get("yield") or "").strip()   # #REM-1: вероятность захвата — «Выход» из B-gamma-lines-airframe.md (IAEA PGAA)
        is_capture = bool(NEUTRON_CAPTURE_RE.search(line["reaction"]))
        yield_lbl = "yield per 100 captures" if table_path.endswith(".en.md") else "выход на 100 захватов"
        reaction_full = line["reaction"] + (f"; {yield_lbl}: {yld}" if is_capture and yld and yld != "—" else "")
        reaction = reaction_full.replace(",", ";")
        place = line["place"].replace(",", ";")
        text = line["text"].replace(",", ";")
        
        results.append((E_keV, reaction, place + " | " + line["cls"].replace(",", ";"), net, sigma, signif, found, origin))
    
    # Find unlisted peaks
    unlisted = find_unlisted(E, cph, sig, listed_E)
    for E0, net, sigma, signif in unlisted:
        origin = attribute(cat, E0, fwhm_keV(E0), top=3)
        unlisted_results.append((E0, net, sigma, signif, origin))
    
    # Write lines output
    with open(out_prefix + "_lines.csv", "w", encoding="utf-8", newline="\n") as f:
        f.write("E_keV,reaction,place_class,net_cph,sigma_cph,signif,found,origin\n")
        for row in results:
            E0, reaction, place, net, sigma, signif, found, text = row
            f.write(f"{E0},{reaction},{place},{net},{sigma},{signif},{found},{text}\n")
            
    # Write unlisted output
    with open(out_prefix + "_unlisted.csv", "w", encoding="utf-8", newline="\n") as f:
        f.write("E_keV,net_cph,sigma_cph,signif,origin\n")
        for row in unlisted_results:
            E0, net, sigma, signif, origin = row
            f.write(f"{E0},{net},{sigma},{signif},{origin}\n")
            
    print(listed_count, found_count, found_A_B, len(unlisted_results))
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
