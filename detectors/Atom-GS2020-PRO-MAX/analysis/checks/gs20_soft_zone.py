# -*- coding: utf-8 -*-
"""GS-20: мягкая зона GS2020 Th-232 — нетто против модели М1 по полосам и пикам. Запуск: python gs20_soft_zone.py <out_dir> [json ...]"""
import sys, os, json; import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
BANDS = [(20, 40), (40, 60), (60, 90), (90, 120), (120, 150), (150, 200), (200, 300), (300, 600), (600, 1500)]
PEAKS = [("K-рентген Pb/Bi 72–90", (66, 96), (56, 64), (98, 106)), ("238,6 Pb-212", (218, 260), (200, 212), (262, 268)),
         ("583,2 Tl-208", (550, 616), (530, 545), (620, 635)), ("911,2 Ac-228", (870, 950), (845, 862), (990, 1005))]
def band(e, v, a, b):
    mask = (e >= a) & (e < b); return float(v[mask].sum())
def peak_area(e, v, win, lb, rb):
    yl = float(v[(e >= lb[0]) & (e < lb[1])].mean()); xl = (lb[0] + lb[1]) / 2
    yr = float(v[(e >= rb[0]) & (e < rb[1])].mean()); xr = (rb[0] + rb[1]) / 2
    mask = (e >= win[0]) & (e < win[1]); base = yl + (yr - yl) * (e[mask] - xl) / (xr - xl)
    return float((v[mask] - base).sum())
def main():
    out = sys.argv[1]; files = sys.argv[2:] if len(sys.argv) > 2 else sorted(f for f in os.listdir(out) if f.startswith("fit_m1") and f.endswith(".json"))
    for f in files:
        d = json.load(open(os.path.join(out, f), encoding="utf-8")); e = np.asarray(d["e"], float); net = np.asarray(d["net"], float); mod = np.asarray(d["chain"]["model"], float)
        print(f"== {f}  (k_bg {d['k_bg']:.4f}, A {d['chain']['A_Bq']:.1f} Бк)")
        print("полоса, кэВ | нетто | модель | модель/нетто | избыток модели в долях нетто")
        for a, b in BANDS:
            n = band(e, net, a, b); mo = band(e, mod, a, b)
            if n == 0: print(f"{a:>4}–{b:<5} | {n:12.0f} | {mo:12.0f} | нет отсчётов")
            else: print(f"{a:>4}–{b:<5} | {n:12.0f} | {mo:12.0f} | {mo / n:6.3f} | {(mo - n) / n:+.3f}")
        print("пик | площадь нетто | площадь модели | модель/нетто")
        for name, win, lb, rb in PEAKS:
            pn = peak_area(e, net, win, lb, rb); pm = peak_area(e, mod, win, lb, rb)
            if pn == 0: print(f"{name} | {pn:12.0f} | {pm:12.0f} | нет отсчётов")
            else: print(f"{name} | {pn:12.0f} | {pm:12.0f} | {pm / pn:6.3f}")
        print()
if __name__ == "__main__": main()
