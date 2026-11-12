import json
import math
import sys
import pathlib
import argparse

sys.stdout.reconfigure(encoding="utf-8")

PREFIXES = ["k24g1t", "dtF", "dtP", "dtH", "dtE1", "dtE2", "dtE3", "dtE4"]
VARIANTS = ["dtP", "dtH", "dtE1", "dtE2", "dtE3", "dtE4"]
CONTROL = "k24g1t"
CHAINS = ["Th228chain", "Ra228chain"]
VM = "marinelli_OISN16_2016"
PEAKS = {"Th228chain": [238.6, 583.2, 727.3, 2614.5], "Ra228chain": [338.3, 911.2, 969.0]}
BANDS = [[50, 150], [150, 260], [300, 500], [700, 1500], [1600, 2600]]

_EDGES = list(range(20, 301, 5)) + list(range(310, 1001, 10)) + list(range(1020, 2801, 20))

_out_dir = None


def load(prefix, chain):
    p = pathlib.Path(_out_dir) / f"{prefix}_{chain}_{VM}.csv"
    if not p.is_file():
        raise SystemExit(f"Missing template CSV: {p}")
    n = None
    spectrum = {}
    in_spec = False
    with open(p, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            first = parts[0].strip()
            if first == "n_events_processed":
                if len(parts) >= 2:
                    n = float(parts[1])
            elif first == "bin_keV":
                in_spec = True
            elif in_spec and len(parts) >= 2:
                spectrum[float(parts[0])] = float(parts[1])
    if n is None:
        raise SystemExit(f"n_events_processed not found in {p}")
    return spectrum, n


def area(h, lo, hi):
    total = 0.0
    for e, c in h.items():
        if lo <= e < hi:
            total += c
    return total


def clean(obj):
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
        return obj
    if isinstance(obj, dict):
        return {k: clean(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [clean(v) for v in obj]
    return obj


def main():
    ap = argparse.ArgumentParser(description="Export detail data as data.js")
    ap.add_argument("--out-dir", required=True, help="Directory with template CSVs")
    ap.add_argument("--fit-dir", required=True, help="Directory with fit_<prefix>.json")
    ap.add_argument("--dest", required=True, help="Path to output data.js")
    ap.add_argument("--text", required=True, help="page_text.json (page prose, embedded as window.G1SD_TEXT)")
    args = ap.parse_args()

    global _out_dir
    _out_dir = args.out_dir
    fit_dir = args.fit_dir
    dest = pathlib.Path(args.dest)

    all_prefixes = ["dtF"] + VARIANTS + [CONTROL]

    # Pre-load all spectra
    loaded = {}
    for chain in CHAINS:
        loaded[chain] = {}
        for P in all_prefixes:
            loaded[chain][P] = load(P, chain)

    # --- Peaks ---
    peaks = {}
    for chain in CHAINS:
        peaks[chain] = {}
        F, nF = loaded[chain]["dtF"]
        for X in VARIANTS + [CONTROL]:
            hX, nX = loaded[chain][X]
            pk_list = []
            for E in PEAKS[chain]:
                cf = area(F, E - 2, E + 2)
                cx = area(hX, E - 2, E + 2)
                r = (cx / nX) / (cf / nF)
                s = r * math.sqrt(1.0 / max(cf, 1) + 1.0 / max(cx, 1))
                pk_list.append({"E": E, "r": round(r, 4), "s": round(s, 4)})
            bd_list = []
            for lo_b, hi_b in BANDS:
                rb = (area(hX, lo_b, hi_b) / nX) / (area(F, lo_b, hi_b) / nF)
                bd_list.append({"lo": lo_b, "hi": hi_b, "r": round(rb, 4)})
            peaks[chain][X] = {"peaks": pk_list, "bands": bd_list}

    # --- Spectra ---
    n_bins = len(_EDGES) - 1
    centres = [(_EDGES[i] + _EDGES[i + 1]) / 2.0 for i in range(n_bins)]
    widths = [_EDGES[i + 1] - _EDGES[i] for i in range(n_bins)]

    spectra = {}
    for chain in CHAINS:
        F, nF = loaded[chain]["dtF"]
        y_dict = {}
        for P in all_prefixes:
            hP, nP = loaded[chain][P]
            y_list = []
            for i in range(n_bins):
                a = area(hP, _EDGES[i], _EDGES[i + 1])
                y_i = a / nP / widths[i]
                y_list.append(float(f"{y_i:.6g}"))
            y_dict[P] = y_list

        ratio_dict = {}
        for P in VARIANTS + [CONTROL]:
            hP, nP = loaded[chain][P]
            r_list = []
            s_list = []
            for i in range(n_bins):
                aF = area(F, _EDGES[i], _EDGES[i + 1])
                aP = area(hP, _EDGES[i], _EDGES[i + 1])
                if aF <= 0 or aP <= 0:
                    r_list.append(None)
                    s_list.append(None)
                else:
                    rr = (aP / nP) / (aF / nF)
                    rs = rr * math.sqrt(1.0 / max(aF, 1) + 1.0 / max(aP, 1))
                    r_list.append(round(rr, 4))
                    s_list.append(round(rs, 4))
            ratio_dict[P] = {"r": r_list, "s": s_list}

        spectra[chain] = {
            "e": [round(c, 1) for c in centres],
            "y": y_dict,
            "ratio": ratio_dict,
        }

    # --- Fits ---
    fits = {}
    fit_meta = None
    for P in PREFIXES:
        fp = pathlib.Path(fit_dir) / f"fit_{P}.json"
        if not fp.is_file():
            raise SystemExit(f"Missing fit file: {fp}")
        with open(fp, encoding="utf-8") as fh:
            d = json.load(fh)
        res = d["results"][0]
        g = res["groups"]
        fits[P] = {
            "ra": round(g["Ra228chain"]["A_over_passport"], 4),
            "th": round(g["Th228chain"]["A_over_passport"], 4),
            "dra": round(g["Ra228chain"]["dA_Bq"] / g["Ra228chain"]["A_Bq"], 4),
            "dth": round(g["Th228chain"]["dA_Bq"] / g["Th228chain"]["A_Bq"], 4),
            "chi2": round(res["chi2_ref_nu"], 4),
            "passport_bq": round(res["passport_bq"], 1),
            "unc_pct": round(res["unc_pct"], 4),
            "date_meas": res["date_meas"],
            "live_s": round(res["live_s"], 1),
        }
        if P == "dtF":
            fit_meta = {
                "lo_keV": d["lo_keV"],
                "hi_keV": d["hi_keV"],
                "fwhm_source": d["fwhm_source"],
            }

    # --- Assemble and write ---
    data = {
        "generated": "2026-11-12",
        "ref": "dtF",
        "variants": VARIANTS,
        "control": CONTROL,
        "chains": CHAINS,
        "fit_meta": fit_meta,
        "fits": fits,
        "peaks": peaks,
        "spectra": spectra,
    }
    data = clean(data)

    dest.parent.mkdir(parents=True, exist_ok=True)
    text = "window.G1SD = " + json.dumps(data, ensure_ascii=False, separators=(",", ":")) + ";\n"
    page_text = json.loads(pathlib.Path(args.text).read_text(encoding="utf-8"))
    text += "window.G1SD_TEXT = " + json.dumps(page_text, ensure_ascii=False, separators=(",", ":")) + ";\n"
    with open(dest, "w", encoding="utf-8") as fh:
        fh.write(text)

    # --- Summary ---
    size = dest.stat().st_size
    for chain in CHAINS:
        dtP_ratios = " ".join(f"{p['r']:.4f}" for p in peaks[chain]["dtP"]["peaks"])
        print(f"{chain}, {n_bins} bins, dtP peak ratios [{dtP_ratios}], {size} bytes")


if __name__ == "__main__":
    main()
