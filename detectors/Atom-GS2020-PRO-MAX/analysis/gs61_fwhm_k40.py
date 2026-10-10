# -*- coding: utf-8 -*-
import sys
sys.stdout.reconfigure(encoding="utf-8")
import json, numpy as np
from scipy.optimize import curve_fit

path = sys.argv[1] if len(sys.argv) > 1 else r"web\gs2020-th232-page\dist\data-k40.js"
with open(path, encoding="utf-8") as f:
    txt = f.read()
d = json.loads(txt[txt.index("{"):txt.rindex("}") + 1])

E = np.array(d["spectrum"]["e_of_ch"])
arrays = {
    "counts": np.array(d["spectrum"]["counts"]),
    "model_counts": np.array(d["spectrum"]["model_counts"]),
    "model2_counts": np.array(d["spectrum"]["model2_counts"]),
    "model2_full_counts": np.array(d["spectrum"]["model2_full_counts"]),
}
k, p = d["fwhm_cal"]["k"], d["fwhm_cal"]["p"]
pts = [(round(pt["E_nominal"], 1), round(pt["fwhm_keV"], 2)) for pt in d["fwhm_cal"]["points"]]

def gauss(x, a, mu, s, b0, b1):
    return a * np.exp(-0.5 * ((x - mu) / s) ** 2) + b0 + b1 * (x - 1460)

def fwhm_b(x, y):
    n = y - np.interp(x, [x[0], x[-1]], [np.mean(y[:5]), np.mean(y[-5:])])
    k_idx = np.argmax(n)
    h = n[k_idx] / 2
    L = np.interp(h, n[:k_idx + 1], x[:k_idx + 1])
    R = np.interp(-h, -n[k_idx:], x[k_idx:])
    return x[k_idx], R - L

for lo, hi in [(1370, 1560), (1340, 1590), (1390, 1535)]:
    mask = (E > lo) & (E < hi)
    x, y_all = E[mask], {name: arr[mask] for name, arr in arrays.items()}
    print(f"окно {lo}–{hi} кэВ:")
    for name, y in y_all.items():
        try:
            popt, _ = curve_fit(gauss, x, y, p0=[np.max(y), 1460, 32, np.min(y), 0], sigma=np.sqrt(np.maximum(y, 1)))
            c_a, f_a = popt[1], 2.3548 * abs(popt[2])
        except Exception:
            c_a, f_a = np.nan, np.nan
        c_b, f_b = fwhm_b(x, y)
        print(f"  {name}: центр {c_a:.1f}, ПШПВ гаусс {f_a:.1f} / полувысота {f_b:.1f}")

print(f"закон страницы: k={k:.4f} p={p:.4f} → ПШПВ(1460,8)={k * 1460.822**p:.2f}; точки подгонки: {pts}")
