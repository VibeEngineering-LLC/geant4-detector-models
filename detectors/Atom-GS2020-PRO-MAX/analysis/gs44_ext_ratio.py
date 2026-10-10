# -*- coding: utf-8 -*-
# #GS-44 (29.09): r(E) = отклик внешнего фона под KCl / под водой по Geant4-прогонам внешнего поля (ext_source.hh),
# веса падающих энергий — NNLS по измеренному фону воды. Запуск: bash run_gs44_ext_ratio.sh (окружение подгонки KCl).
import os, sys, re, glob, json
import numpy as np
from scipy.optimize import nnls
sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import fit_gs2020_kcl as fk
from merge_templates_gs2020 import read_chunk
m1 = fk.m1; CAL = fk.CAL
EXT = os.environ.get("GS_EXT_DIR", r"C:\g4work\gs2020\ext")
SELFTEST = os.environ.get("GS_EXT_SELFTEST") == "1"
SAMPLE = os.environ.get("GS_EXT_SAMPLE", "kcl")   # #GS-60 (07.10): тег пробы в именах ext_<тег>_E…; kcl = KCl 740 мл (#GS-44)
TAGSFX = "" if SAMPLE == "kcl" else "_" + SAMPLE
if os.environ.get("GS_BG_R"): print("ОТКАЗ: GS_BG_R задан — r(E) считается по НЕослабленному фону воды (#GS-60)"); sys.exit(1)
FIT_LO, FIT_HI = (float(x) for x in os.environ.get("GS_EXT_FIT", "30,2800").split(","))
BANDS = [(10, 20), (20, 26), (26, 30), (30, 34), (34, 38), (38, 45), (45, 60), (60, 90), (90, 150),
         (150, 300), (300, 600), (600, 1000), (1000, 1500), (1500, 3000), (150, 3000)]

def fail(msg): print("ОТКАЗ: " + msg); sys.exit(1)

def main():
    # Step A
    s = m1.Spec(m1.bm.read(CAL.PATHS["kcl"])[0], "kcl"); b = m1.Spec(m1.bm.read(CAL.BKG_WATER_XML)[0], "bgw")
    fwhm_csv = os.path.join(EXT, "fwhm_points_gs2020.csv"); m1.write_fwhm_csv(fwhm_csv)
    bg_e = np.array([b.channel_to_energy(i) for i in range(b.n_channels)])
    r = m1.muc.unfold(s, b, [("K40", os.path.join(fk.OUT, "mix_K40_npsmoff.csv"))], fwhm_csv, lo=fk.LO, hi=fk.HI,
                    bg_energy_of_ch=bg_e, verbose=False, blur=m1.BLUR, tail=m1.TAIL)
    e = np.asarray(r["e"], float); ch_edges = r["ch_edges"]; fw = lambda E: m1.BLUR * r["fwhm"](E)
    bg = np.asarray(r["bg_scaled"], float); k = s.live_time / b.live_time

    if SAMPLE != "kcl":   # #GS-70: «нужно» — диагностика по подгонке KCl; у другой пробы её нет (nan), r(E) от неё не зависит
        je = e; jnet = jmod = jbg = np.full(len(e), np.nan)
    else:
        j = json.load(open(os.path.join(fk.OUT, "fit_kcl_bgw.json"), encoding="utf-8"))
        je = np.array(j["e"], float); jnet = np.array(j["net"], float); jmod = np.array(j["model"], float); jbg = np.array(j["bg"], float)
    if len(je) != len(e) or max(abs(je - e)) > 1e-6:
        fail("Mismatch in energy arrays from fit_kcl_bgw.json")

    raw = jnet + jbg
    print(f"k_bg={k:.4f}, каналов {len(e)}")
    mask_hi = e > FIT_HI
    mask_fit = (e >= FIT_LO) & (e <= FIT_HI)
    print(f"фон > {FIT_HI}: {bg[mask_hi].sum():.0f}")
    print(f"фон [{FIT_LO}, {FIT_HI}]: {bg[mask_fit].sum():.0f}")

    # Step B
    files = glob.glob(os.path.join(EXT, "chunks", "ext_*_E*_*.csv"))
    pat = re.compile(r"ext_([a-z0-9]+)_E([0-9.]+)_(up|down)\.csv$")
    D = {}
    energies_found = set()

    for path in files:
        basename = os.path.basename(path)
        m = pat.match(basename)
        if not m: continue
        geom, E_str, hemi = m.group(1), float(m.group(2)), m.group(3)
        energies_found.add(E_str)

        c = read_chunk(path); h = dict(c["header"])
        if h.get("src_mode") != "ext_sphere":
            fail(f"Wrong src_mode in {path}")
        n = float(h["n_events_processed"])
        hist = {float(bb): float(cc) for bb, cc in zip(c["bins"], c["edep"]) if cc > 0}

        col = m1.muc.g1s.broaden(hist, 1.0, ch_edges, fw) / n
        var = col / n

        raw1 = np.zeros(3700)
        for bb, cc in hist.items():
            idx = int(float(bb))
            if 0 <= idx < len(raw1):
                raw1[idx] += cc / n

        D[(geom, E_str, hemi)] = (col, var, raw1, n)

    energies_sorted = sorted(energies_found, key=float)
    missing = []
    for E in energies_sorted:
        for g in ("water", SAMPLE):
            for h in ("up", "down"):
                if (g, E, h) not in D:
                    missing.append(f"{g}_E{E}_{h}")

    if missing:
        fail(f"Missing files: {', '.join(missing)}")

    if SELFTEST:
        for E in energies_sorted:
            for h in ("up", "down"):
                w_data = D[("water", E, h)]
                D[(SAMPLE, E, h)] = w_data

    print(f"Files processed: {len(files)}, Energies: {len(energies_sorted)}")

    # Step C & D & E & F  (порядок кортежа D: col, var, raw1, n — шаг G чинен вручную 29.09: брал col вместо raw1)
    variants = ["4pi", "up", "down"]
    results_json = {}

    for v in variants:
        # NNLS
        F = (e >= FIT_LO) & (e <= FIT_HI)
        sig = np.sqrt(np.maximum(k * bg[F], k))
        
        X = np.zeros((F.sum(), len(energies_sorted)))
        for i, E in enumerate(energies_sorted):
            if v == "4pi":
                col_up, _, _, _ = D[("water", E, "up")]
                col_dn, _, _, _ = D[("water", E, "down")]
                R_E = 0.5 * (col_up + col_dn)
            elif v == "up":
                R_E, _, _, _ = D[("water", E, "up")]
            else:
                R_E, _, _, _ = D[("water", E, "down")]
            
            X[:, i] = R_E[F] / sig

        y = bg[F] / sig
        w, _ = nnls(X, y)

        Bw = np.zeros_like(bg)
        Vk_total = np.zeros_like(bg) # Variance of KCl model
        Vw_total = np.zeros_like(bg) # Variance of Water model
        
        for i, E in enumerate(energies_sorted):
            if w[i] == 0: continue
            
            if v == "4pi":
                col_up_w, var_up_w, _, _ = D[("water", E, "up")]
                col_dn_w, var_dn_w, _, _ = D[("water", E, "down")]
                R_E_w = 0.5 * (col_up_w + col_dn_w)
                V_E_w = 0.25 * (var_up_w + var_dn_w)

                col_up_k, var_up_k, _, _ = D[(SAMPLE, E, "up")]
                col_dn_k, var_dn_k, _, _ = D[(SAMPLE, E, "down")]
                R_E_k = 0.5 * (col_up_k + col_dn_k)
                V_E_k = 0.25 * (var_up_k + var_dn_k)
            elif v == "up":
                R_E_w, V_E_w, _, _ = D[("water", E, "up")]
                R_E_k, V_E_k, _, _ = D[(SAMPLE, E, "up")]
            else:
                R_E_w, V_E_w, _, _ = D[("water", E, "down")]
                R_E_k, V_E_k, _, _ = D[(SAMPLE, E, "down")]

            Bw += w[i] * R_E_w
            Vk_total += (w[i]**2) * V_E_k
            Vw_total += (w[i]**2) * V_E_w
            
        # KCl Model Background
        Bk = np.zeros_like(bg)
        for i, E in enumerate(energies_sorted):
            if w[i] == 0: continue
            if v == "4pi":
                col_up_k, _, _, _ = D[(SAMPLE, E, "up")]
                col_dn_k, _, _, _ = D[(SAMPLE, E, "down")]
                R_E_k = 0.5 * (col_up_k + col_dn_k)
            elif v == "up":
                R_E_k, _, _, _ = D[(SAMPLE, E, "up")]
            else:
                R_E_k, _, _, _ = D[(SAMPLE, E, "down")]
            Bk += w[i] * R_E_k

        # Chi2
        residuals = (bg[F] - Bw[F]) / sig
        chi2 = np.sum(residuals**2)
        ndof = F.sum() - np.count_nonzero(w)
        if ndof <= 0: ndof = 1 # Avoid div by zero
        chi2_ndof = chi2 / ndof

        print(f"Variant {v}: chi2/ndof = {chi2_ndof:.4f}")
        nz_weights = [(E, w[i]) for i, E in enumerate(energies_sorted) if w[i] > 0]
        for E, val in nz_weights:
            print(f"  E={E}: w={val:.6e}")

        # Band Table
        band_table = []
        for lo, hi in BANDS:
            sel = (e >= lo) & (e < hi)
            sum_Bw = Bw[sel].sum()
            sum_Bk = Bk[sel].sum()
            
            if sum_Bw == 0:
                rb = float('nan')
                drb = float('nan')
            else:
                rb = sum_Bk / sum_Bw
                term_k = Vk_total[sel].sum() / (sum_Bk**2) if sum_Bk != 0 else 0
                term_w = Vw_total[sel].sum() / (sum_Bw**2)
                drb = rb * np.sqrt(term_k + term_w)

            closure = bg[sel].sum() / sum_Bw if sum_Bw != 0 else float('nan')
            
            # Need calculation
            sum_raw = raw[sel].sum()
            sum_jmod = jmod[sel].sum()
            sum_jbg = jbg[sel].sum()
            need = (sum_raw - sum_jmod) / sum_jbg if sum_jbg != 0 else float('nan')

            print(f"{lo:4g}–{hi:4g} кэВ  нужно {need:.3f}  r {rb:.3f} ± {drb:.3f}  замыкание NNLS {closure:.3f}")
            
            band_table.append({
                "lo": lo, "hi": hi, "need": need, "r": rb, "dr": drb, "closure": closure
            })

        # Smooth Curve
        G = np.arange(10.0, 3001.0, 5.0)
        csv_rows = []
        for g in G:
            hw = max(5.0, 0.5 * fw(g))
            sel = abs(e - g) <= hw
            
            sum_Bw_g = Bw[sel].sum()
            if sum_Bw_g <= 0: continue
            
            sum_Bk_g = Bk[sel].sum()
            rb_g = sum_Bk_g / sum_Bw_g
            
            term_k_g = Vk_total[sel].sum() / (sum_Bk_g**2) if sum_Bk_g != 0 else 0
            term_w_g = Vw_total[sel].sum() / (sum_Bw_g**2)
            drb_g = rb_g * np.sqrt(term_k_g + term_w_g)
            
            csv_rows.append((g, rb_g, drb_g))

        suffix = TAGSFX + ("_selftest" if SELFTEST else "")
        csv_path = os.path.join(EXT, f"r_ext_{v}{suffix}.csv")
        with open(csv_path, "w", encoding="utf-8") as f:
            f.write("E_keV,r,dr\n")
            for row in csv_rows:
                f.write(f"{row[0]:.1f},{row[1]:.6f},{row[2]:.6f}\n")

        # Store for JSON
        weights_dict = {str(E): float(w[i]) for i, E in enumerate(energies_sorted)}
        results_json[v] = {
            "weights": weights_dict,
            "chi2_ndof": chi2_ndof,
            "bands": band_table
        }

    # Step G - Photopeak Ratio (4pi only, unbroadened)
    photopeaks = [float(x) for x in energies_sorted] if SAMPLE != "kcl" else [609.312, 1764.494, 2614.511]   # #GS-70: все узлы — для r(E) «линии+континуум» (gs70_r2.py)
    pp_results = {}
    
    for E0 in photopeaks:
        # Check if E0 is in energies_sorted (approx match)
        found_E = None
        for E_str in energies_sorted:
            if abs(float(E_str) - E0) < 1.0:
                found_E = E_str
                break
        
        if found_E is None:
            continue
            
        # Get raw1 and n for up and down for water and kcl
        _, _, r1_w_up, n_w_up = D[("water", found_E, "up")]
        _, _, r1_w_dn, n_w_dn = D[("water", found_E, "down")]
        _, _, r1_k_up, n_k_up = D[(SAMPLE, found_E, "up")]
        _, _, r1_k_dn, n_k_dn = D[(SAMPLE, found_E, "down")]

        # Peak window: int(E0) - 2 to int(E0) + 3 (exclusive of end in slice? No, inclusive usually. Python slice [start:end] excludes end. So +3 means up to index+3-1? Let's assume standard peak integration around channel.)
        # Prompt says: sum of raw1[int(E0) - 2 : int(E0) + 3]
        idx = int(E0)
        start_idx = max(0, idx - 2)
        end_idx = min(len(r1_w_up), idx + 3)
        
        # Counts per primary in peak
        cnt_w_up = r1_w_up[start_idx:end_idx].sum()
        cnt_w_dn = r1_w_dn[start_idx:end_idx].sum()
        cnt_k_up = r1_k_up[start_idx:end_idx].sum()
        cnt_k_dn = r1_k_dn[start_idx:end_idx].sum()

        # Absolute counts (per-primary * n)
        abs_w_up = cnt_w_up * n_w_up
        abs_w_dn = cnt_w_dn * n_w_dn
        abs_k_up = cnt_k_up * n_k_up
        abs_k_dn = cnt_k_dn * n_k_dn

        total_w = abs_w_up + abs_w_dn
        total_k = abs_k_up + abs_k_dn

        if total_w == 0:
            ratio_pp = float('nan')
            err_pp = float('nan')
        else:
            ratio_pp = total_k / total_w
            # Poisson error propagation for ratio K/W
            # Var(ratio) = (1/W^2)*Var(K) + (K^2/W^4)*Var(W)
            # Var(N) = N for Poisson
            var_k = total_k
            var_w = total_w
            err_pp = np.sqrt((var_k / total_w**2) + (total_k**2 * var_w / total_w**4))

        print(f"фотопик {E0:.1f}: МК KCl/вода {ratio_pp:.4f} ± {err_pp:.4f}")
        pp_results[E0] = {"ratio": ratio_pp, "error": err_pp}

    # Step H - Write JSON
    json_out = {
        "k_bg": k,
        "FIT_LO": FIT_LO,
        "FIT_HI": FIT_HI,
        "variants": results_json,
        "photopeaks": pp_results
    }
    
    suffix = TAGSFX + ("_selftest" if SELFTEST else "")
    json_path = os.path.join(EXT, f"gs44_ext_ratio{suffix}.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(json_out, f, ensure_ascii=False, indent=1)
        
    print(f"Output written to {json_path}")

if __name__ == "__main__":
    main()
