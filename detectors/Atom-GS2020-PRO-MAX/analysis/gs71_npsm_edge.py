# -*- coding: utf-8 -*-
# #GS-71 этап 2: сдвиг края Комптона «свет − депозит», который даёт NPSM на монолинии (столбцы count_edep/count_light одного прогона).
# Свет растянут по вершине фотопика (как gs68_npsm_compare.py); оба свёрнуты гауссом ПШПВ карточки; ступенька — как gs71_edge.py. argv: csv E0 [csv E0 ...]
import sys, json, numpy as np; from scipy.optimize import curve_fit; from scipy.special import erfc; sys.stdout.reconfigure(encoding="utf-8")
C = json.load(open("<REPOS>/geant4-detector-models/detectors/Atom-GS2020-PRO-MAX/dataset/detector_card.json", encoding="utf-8"))["fwhm"]
P = np.array(C["points_sl"]); fw = lambda e: np.exp(np.interp(np.log(np.maximum(e, 1)), np.log(P[:, 0]), np.log(P[:, 1])))
def step(E, h, c, s, a, b): return h * erfc((E - c) / (np.sqrt(2) * abs(s))) / 2 + a + b * (E - c)
def fold(x, y):
    sg = fw(x) / 2.3548; G = np.exp(-0.5 * ((x[None, :] - x[:, None]) / sg[:, None]) ** 2); G /= G.sum(axis=1, keepdims=True); return G.T @ y
for path, E0 in zip(sys.argv[1::2], map(float, sys.argv[2::2])):
    t = open(path, encoding="utf-8").read().splitlines(); i = t.index("bin_keV,count_edep,count_light"); h = dict(l.split(",", 1) for l in t[:i] if "," in l)
    a = np.array([[float(v) for v in r.split(",")] for r in t[i + 1:]]); x, ed, li = a[:, 0], a[:, 1], a[:, 2]; m = x < 1.2 * E0; x, ed, li = x[m], ed[m], li[m]
    r = (x > 0.3 * E0) & (x < 1.05 * E0); p0 = x[r][np.argmax(li[r])]; w = np.abs(x - p0) < 0.03 * p0; k = E0 / ((x[w] * li[w]).sum() / li[w].sum())
    lis = np.diff(np.interp(np.concatenate([x - 0.5, [x[-1] + 0.5]]), np.concatenate([x - 0.5, [x[-1] + 0.5]]) * k, np.concatenate([[0], np.cumsum(li)])))
    Ec = E0 * (1 - 1 / (1 + 2 * E0 / 510.999)); res = []
    for y in (fold(x, ed), fold(x, lis)):
        s = (x > 0.75 * Ec) & (x < min(1.25 * Ec, 0.91 * E0)); q, cv = curve_fit(step, x[s], y[s], p0=[y[s][0] - y[s][-1], Ec, 0.03 * Ec, y[s][-1], 0], maxfev=20000); res.append((q[1], np.sqrt(cv[1, 1]), abs(q[2])))   # W-174: окно ниже фотопика
    print(f"{E0:7.1f} кэВ η={h.get('npsm_eta', '?')}: край депозит {res[0][0]:.1f}±{res[0][1]:.1f}, свет {res[1][0]:.1f}±{res[1][1]:.1f} → сдвиг {res[1][0] - res[0][0]:+.1f} кэВ "
          f"({100 * (res[1][0] - res[0][0]) / Ec:+.2f} %), ширина ×{res[1][2] / res[0][2]:.2f}; k {k:.4f}")
