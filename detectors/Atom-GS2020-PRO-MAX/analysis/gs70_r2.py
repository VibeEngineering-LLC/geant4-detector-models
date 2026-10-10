# -*- coding: utf-8 -*-
# #GS-70 (оператор 10.10 «делай»): r(E) фона «линии + континуум». r по NNLS-сетке у лёгкой пробы «звенит» на узлах сетки (фотопик ≈1,06,
# континуум ≈0,93), у KCl ≈ воды — гладко. Здесь: фон воды = пики своих линий (гаусс ПШПВ карточки + линейная подложка, МНК) + континуум;
# r_eff = r_c + (r_p − r_c)·пик/фон; r_c — полосы 4π (gs44_ext_ratio JSON), r_p(E) — фотопики МК по всем узлам. argv: ratio.json fit.json fwhm.csv out.csv
import sys, json, numpy as np; sys.stdout.reconfigure(encoding="utf-8")
R, F = (json.load(open(p, encoding="utf-8")) for p in sys.argv[1:3]); fw = np.genfromtxt(sys.argv[3], delimiter=",", invalid_raise=False)
B = [b for b in R["variants"]["4pi"]["bands"] if b["hi"] - b["lo"] < 1000]; bc = np.array([np.sqrt(b["lo"] * b["hi"]) for b in B])
rc = lambda E: np.interp(np.log(E), np.log(bc), [b["r"] for b in B])
P = sorted((float(k), v["ratio"]) for k, v in R["photopeaks"].items()); rp = lambda E: np.interp(E, [p[0] for p in P], [p[1] for p in P])
e, bg = np.array(F["e"]), np.array(F["bg"]); pk = np.zeros_like(bg)
for E0 in (238.632, 295.224, 351.932, 583.187, 609.312, 911.204, 1120.287, 1460.822, 1764.494, 2614.511):
    s = np.interp(E0, fw[:, 0], fw[:, 1]) / 2.3548; w = np.abs(e - E0) < 3 * s
    g = np.exp(-0.5 * ((e[w] - E0) / s) ** 2); M = np.vstack([g, np.ones(w.sum()), e[w] - E0]).T
    a = np.linalg.lstsq(M, bg[w], rcond=None)[0]; pk[w] += max(a[0], 0) * g
    print(f"линия {E0:8.2f}: пик фона {max(a[0], 0) * g.sum():9.0f} отсч. ({100 * max(a[0], 0) * g.sum() / bg[w].sum():4.1f} % окна), r_p {rp(E0):.4f}, r_c {rc(E0):.4f}")
fr = np.clip(pk / np.maximum(bg, 1e-9), 0, 1); Eg = np.arange(10.0, 3001.0, 5.0)
reff = rc(Eg) + (rp(Eg) - rc(Eg)) * np.interp(Eg, e, fr)
with open(sys.argv[4], "w", encoding="utf-8") as f:
    f.write("E_keV,r,dr\n"); [f.write(f"{x:.1f},{y:.6f},0.003\n") for x, y in zip(Eg, reff)]
print("r_eff в точках: " + ", ".join(f"{x:g}→{np.interp(x, Eg, reff):.4f}" for x in (100, 200, 238.6, 400, 609.3, 1000, 1460.8, 2000, 2614.5)))
