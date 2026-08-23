"""
#FIT-1 Stage-3: apparent_method1_compare.py
Proveryaet gipotezu o vliyanii LCE-tailing na 62%-deficit kontinuuma.
Chitaet boevye shablony metoda 1 i Stage-2 Y-raspredeleniya, stroit
"apparent"-versii s uchetom karty LCE(Y), primenyaet te zhe NNLS-amplitudy,
sravnivaet metriki kachestva podgonki i sohranyaet rezultaty v CSV.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8")

import io, contextlib, math, os
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import numpy as np
import fit_nuclides as fn
import fit_coverage as fc
import rcspec
import read_rcxml

LCE_MAP = {-4.5: 25.0, -3.0: 25.7, -1.5: 23.6, 0.0: 20.0,
           1.5: 18.8, 3.0: 18.2, 4.5: 15.8}
STAGE2_DIR = "D:/Claude_files/repos/geant4-detector-models/build/RadiaCode-103/_stage2_ypos"


def load_ypos(path):
    if not os.path.exists(path):
        return None, None
    y_cols = None
    dist = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(",")
            if y_cols is None:
                y_cols = [float(p[1:]) for p in parts[1:]]
                continue
            e_kev = float(parts[0])
            i = int(round(e_kev - 0.5))
            counts = np.array([float(x) for x in parts[1:]])
            dist[i] = counts
    return y_cols, dist


def global_lce_ref(all_dists):
    if not all_dists:
        return 0.0, 0
    y_cols = all_dists[0][0]
    total = np.zeros(len(y_cols))
    for _, dist in all_dists:
        for counts in dist.values():
            total += counts
    n_total = float(total.sum())
    if n_total <= 0:
        return 0.0, 0
    lce_ref = sum(total[k] * LCE_MAP[y] for k, y in enumerate(y_cols)) / n_total
    return lce_ref, int(n_total)


def build_apparent_hist(raw_hist, y_cols, dist, lce_ref):
    global_counts = None
    for counts in dist.values():
        global_counts = counts.copy() if global_counts is None else global_counts + counts
    if global_counts is not None and global_counts.sum() > 0:
        global_p = global_counts / global_counts.sum()
    else:
        global_p = None

    out = np.zeros_like(raw_hist)
    n = len(raw_hist)
    for i in range(n):
        c = raw_hist[i]
        if c <= 0:
            continue
        counts = dist.get(i)
        if counts is None or counts.sum() < 20:
            p = global_p
        else:
            p = counts / counts.sum()
        if p is None:
            out[i] += c
            continue
        e_center = i + 0.5
        for k, y in enumerate(y_cols):
            if p[k] <= 0:
                continue
            e_app = e_center * LCE_MAP[y] / lce_ref
            j = int(e_app)
            if 0 <= j < n:
                out[j] += c * p[k]
    return out


def load_column(nuc):
    wf = os.path.join(fn.BUILD, "%s_%s.csv" % (fn.WF_PREFIX, nuc))
    bg = os.path.join(fn.BG_DIR, "%s_%s.csv" % (fn.BG_PREFIX, nuc))
    if not (os.path.exists(wf) and os.path.exists(bg)):
        return None, None
    flu = fn.read_wallfield_total(wf)
    r, hz = fn.CYL["r"] / 10.0, 0.5 * (fn.CYL["z1"] - fn.CYL["z0"]) / 10.0
    area = 2 * math.pi * r * (r + 2 * hz)
    rate = flu * area / 4.0
    meta, hist = rcspec.read_spec(bg)
    t_run = float(meta["N_primaries"]) / rate
    raw_norm = hist / t_run
    return raw_norm, meta


def predict(names, A, amp_dict):
    pred = np.zeros(A.shape[0])
    for name, a in amp_dict.items():
        if name in names:
            pred += a * A[:, names.index(name)]
    return pred


def main():
    print("=== #FIT-1 Stage-3: apparent-спектр метода 1 (LCE-tailing) ===")

    stage2_have = []
    for nuc in fn.NUCS:
        y_cols, dist = load_ypos(os.path.join(STAGE2_DIR, "s2_%s_ypos.csv" % nuc))
        if dist is not None:
            stage2_have.append((nuc, y_cols, dist))

    if not stage2_have:
        print("НЕТ STAGE-2 ДАННЫХ, СТОП")
        return

    lce_ref, n_total = global_lce_ref([(yc, d) for _, yc, d in stage2_have])
    print("LCE_ref (глобальный, %d событий всех нуклидов) = %.2f %%" % (n_total, lce_ref))

    names_raw, cols_baseline, cols_apparent = [], [], []
    for nuc in fn.NUCS:
        raw, meta = load_column(nuc)
        if raw is None:
            print("[--] %s: нет шаблона" % nuc)
            continue
        found = next((t for t in stage2_have if t[0] == nuc), None)
        if found is not None:
            _, y_cols, dist = found
            apparent_raw = build_apparent_hist(raw, y_cols, dist, lce_ref)
        else:
            apparent_raw = raw.copy()
            print("[--] %s: нет Stage-2 Y-данных, apparent = baseline" % nuc)
        cps_baseline = rcspec.fold(raw, "103")
        cps_apparent = rcspec.fold(apparent_raw, "103")
        names_raw.append(nuc)
        cols_baseline.append(cps_baseline)
        cols_apparent.append(cps_apparent)

    names_b, cols_b = fn.merge_by_chain(names_raw, cols_baseline)
    names_a, cols_a = fn.merge_by_chain(names_raw, cols_apparent)
    if names_b != names_a:
        raise SystemExit("names_b != names_a - rasoshlis spiski nuklidov")

    mu, pdg = fn.load_muons()
    if mu is not None:
        names_b.append("mu")
        cols_b.append(mu)
        names_a.append("mu")
        cols_a.append(mu)

    smp = read_rcxml.read(fn.MEASURED)[0]
    cnt = smp.counts[:-1].astype(float)
    ch = np.arange(len(cnt))
    e_meas = np.asarray(sum(c * ch ** i for i, c in enumerate(fn.CAL_ROOM)))[:len(cnt)]
    cps_meas = cnt / smp.live

    A_b = np.zeros((len(e_meas), len(cols_b)))
    for k, c in enumerate(cols_b):
        A_b[:, k] = fn.fl.rebin_model_to_meas(np.arange(len(c)) + 0.5, c, e_meas)
    A_a = np.zeros((len(e_meas), len(cols_a)))
    for k, c in enumerate(cols_a):
        A_a[:, k] = fn.fl.rebin_model_to_meas(np.arange(len(c)) + 0.5, c, e_meas)

    mu_idx = names_b.index("mu") if "mu" in names_b else None
    with contextlib.redirect_stdout(io.StringIO()):
        amp_dict, a_mu = fn.fit_by_lines(names_b, A_b, e_meas, cps_meas, smp.live, mu_idx, pdg)

    pred_baseline = predict(names_b, A_b, amp_dict)
    pred_apparent = predict(names_a, A_a, amp_dict)

    print("")
    print("=== SRAVNENIE (te zhe NNLS-amplitudy, raznye shablony baseline/apparent) ===")
    print("metrika                    baseline      apparent")
    print("form_residual_pct %%      %10.2f    %10.2f" % (
        fc.form_residual_pct(pred_baseline, cps_meas),
        fc.form_residual_pct(pred_apparent, cps_meas)))
    print("fraction_covered          %10.4f    %10.4f" % (
        fc.fraction_covered(pred_baseline, cps_meas),
        fc.fraction_covered(pred_apparent, cps_meas)))
    chi2_b, ndf_b = fc.chi2_of(pred_baseline, cps_meas, smp.live, e_meas)
    chi2_a, ndf_a = fc.chi2_of(pred_apparent, cps_meas, smp.live, e_meas)
    print("chi2/ndf                  %10.2f    %10.2f" % (chi2_b / ndf_b, chi2_a / ndf_a))

    print("")
    print("=== PO POLOSAM ===")
    print("%-12s %10s %10s %10s %10s" % ("polosa,keV", "izmereno", "baseline", "apparent", "apparent/izm"))
    for lo, hi in ((20, 100), (100, 300), (300, 700), (700, 1500),
                   (1500, 2000), (2000, 2400), (2400, 2830)):
        m = (e_meas >= lo) & (e_meas < hi)
        ym = cps_meas[m].sum()
        bm = pred_baseline[m].sum()
        am = pred_apparent[m].sum()
        print("%5d-%-6d %10.5f %10.5f %10.5f %10.3f" % (lo, hi, ym, bm, am, am / ym if ym else float("nan")))

    print("")
    print("=== TREBUET TOLKOVANIYA (par.31.A #SA-4) ===")
    print("1. LCE_ref = globalnoe vzveshennoe srednee po VSEM sobytiyam vseh")
    print("   nuklidov - ne edinstvenno vozmozhnyj vybor (alternativa: moda,")
    print("   maksimum LCE kak 'istinnaya' kalibrovochnaya tochka). Smena referensa")
    print("   sdvigaet ABSOLYUTNUYU shkalu apparent-spektra, forma razmazyvaniya")
    print("   mezhdu Y-tochkami ne menyaetsya.")
    print("2. Odni i te zhe NNLS-amplitudy primeneny k raznym shablonam (baseline")
    print("   i apparent) - eto NE polnyj re-fit s apparent-shablonami. Esli")
    print("   raznica metrik sushchestvenna, sleduyushchij shag - pereschitat NNLS na")
    print("   apparent-shablonah napryamuyu (mozhet dat drugie amplitudy).")
    print("3. rcspec.fold() primenyaet PASPORTNOE FWHM(E), otkalibrovannoe na")
    print("   realnom pribore - ono MOZHET uzhe chastichno vklyuchat vklad LCE-")
    print("   dispersii v shirinu pikov. Dvojnogo ucheta tut staraemsya izbezhat")
    print("   tem, chto LCE-kernel primenyaetsya k RAW (do fold) energovydeleniyu,")
    print("   a ne poverh uzhe svyornutogo spektra - no eto ne strogoe")
    print("   dokazatelstvo otsutstviya peresecheniya effektov.")
    print("4. Stage-2 statistika (N=1e6 na nuklid) daet shumnoe P(y|E) na redkih")
    print("   kanalah - fallback na globalnoe P(y) nuklida (porog 20 sobytij)")
    print("   sglazhivaet eto, no ostaetsya istochnikom neopredelennosti.")

    out_dir = os.path.join(_HERE, "results")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "apparent_vs_baseline_20260823.csv")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("E_keV,measured,baseline,apparent\n")
        for i in range(len(e_meas)):
            f.write("%.2f,%.6e,%.6e,%.6e\n" % (e_meas[i], cps_meas[i], pred_baseline[i], pred_apparent[i]))
    print("")
    print("Сохранено: %s" % out_path)


if __name__ == "__main__":
    main()
