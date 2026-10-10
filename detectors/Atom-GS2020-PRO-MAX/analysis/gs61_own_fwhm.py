# -*- coding: utf-8 -*-
import os, sys, re, csv
sys.stdout.reconfigure(encoding="utf-8")
if "SPECTRAVIBE_ROOT" not in os.environ: raise SystemExit("ОТКАЗ: задайте SPECTRAVIBE_ROOT")
sys.path.insert(0, os.path.join(os.environ["SPECTRAVIBE_ROOT"], "scripts"))
from gamma.calibration.fwhm_measure import measure_fwhm, FwhmMeasurementError
import numpy as np
import fit_gs2020_th232_m1 as m1

def fmt(v, d=2): return "—" if v is None else f"{v:.{d}f}"
def get_bm(path):
    try: txt = open(path, encoding="utf-8", errors="replace").read()
    except: return None, None
    m = re.search(r"<PowerFwhmCalibration>.*?<Coefficient>([\d.eE+-]+)</Coefficient>.*?<Coefficient>([\d.eE+-]+)</Coefficient>.*?</PowerFwhmCalibration>", txt, re.S)
    return (float(m.group(1)), float(m.group(2))) if m else (None, None)

tags = sys.argv[1:] if len(sys.argv) > 1 else ["kcl", "sample"]
lines = {"kcl": [238.632, 583.187, 609.312, 911.204, 1120.287, 1460.822, 1764.494, 2614.511],
         "sample": [238.632, 338.32, 583.187, 727.33, 911.204, 1460.822, 2614.511]}
out_dir = os.path.join("..", "results", "gs61_own_fwhm")
os.makedirs(out_dir, exist_ok=True)

for tag in tags:
    path = m1.CAL.PATHS[tag]
    s = m1.Spec(m1.bm.read(path)[0], tag)
    counts = np.array(s.counts, float)
    e = np.asarray(s._e, float)
    ch = np.arange(len(e), dtype=float)
    a, b = get_bm(path)
    cal = list(np.polyfit(ch, e, 4)[::-1])
    print(f"== {tag}: {os.path.basename(path)}")
    rows = []
    for E in lines.get(tag, lines["sample"]):
        c0 = float(np.interp(E, e, ch))
        seed = a * c0**b if a is not None else 0.05 * c0
        status, reason, c, w, fw, cen, unc, spr, bm_keV = "отказ", "", None, None, None, None, None, None, None
        try:
            m = measure_fwhm(counts, peak_channel=c0, energy_cal=cal, seed_fwhm_channels=seed)
            if m.passed and m.fwhm_channels is not None:
                c, w = m.centroid_channel, m.fwhm_channels
                fw = np.interp(c + w/2, ch, e) - np.interp(c - w/2, ch, e)
                cen = np.interp(c, ch, e)
                unc = (m.fwhm_uncertainty_keV * fw / m.fwhm_keV) if (m.fwhm_uncertainty_keV is not None and m.fwhm_keV is not None) else None
                spr = m.method_spread_pct
                status = "успех"
            else:
                reason = m.reason
        except FwhmMeasurementError as ex:
            reason = str(ex)
        cc = c if c is not None else c0
        bm_keV = (np.interp(cc + a*cc**b/2, ch, e) - np.interp(cc - a*cc**b/2, ch, e)) if a is not None else None
        sl = m1.FWHM_SL.get(E)
        print(f"  E={E:.1f}  центр {fmt(cen)}  ПШПВ {fmt(w)} кан = {fmt(fw)} ± {fmt(unc)} кэВ (разброс методов {fmt(spr)} %) | BecqMoni {fmt(bm_keV)} | СпектраЛайн {fmt(sl)} | {status} {reason}")
        rows.append([E, cen, w, fw, unc, spr, bm_keV, sl, status, reason])
    with open(os.path.join(out_dir, f"{tag}.csv"), "w", newline="", encoding="utf-8") as f:
        wtr = csv.writer(f)
        wtr.writerow(["E_nominal","centroid_keV","fwhm_ch","fwhm_keV","unc_keV","spread_pct","becqmoni_keV","sl_keV","status","reason"])
        for r in rows: wtr.writerow(["" if x is None else x for x in r])
