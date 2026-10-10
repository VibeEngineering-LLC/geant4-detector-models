# -*- coding: utf-8 -*-
# #GS-71 (план 10.10, этап 1): край Комптона изм против модели одной метрикой на трёх спектрах GS2020 — в окне Ec·(1±0,25) подгонка
# ступеньки h·erfc((E−c)/(√2·s))/2 + a + b·(E−Ec) отдельно к net и к model; печать c, s, h и отношений. argv: метка json Ec [метка json Ec ...]
import sys, json, numpy as np; from scipy.optimize import curve_fit; from scipy.special import erfc
sys.stdout.reconfigure(encoding="utf-8")
def step(E, h, c, s, a, b): return h * erfc((E - c) / (np.sqrt(2) * abs(s))) / 2 + a + b * (E - c)
print("спектр            | край c, кэВ изм / мод (Δ) | ширина s, кэВ изм / мод (×) | высота h изм / мод (×)")
for nm, p, Ec in zip(sys.argv[1::3], sys.argv[2::3], map(float, sys.argv[3::3])):
    j = json.load(open(p, encoding="utf-8")); e = np.array(j["e"]); E0 = (Ec + np.sqrt(Ec * Ec + 2 * 510.999 * Ec)) / 2; w = (e > 0.75 * Ec) & (e < min(1.25 * Ec, 0.91 * E0)); x = e[w]   # W-174: окно не захватывает фотопик
    res = []
    for y in (np.array(j["net"])[w], np.array(j["model"])[w]):
        p0 = [y[x < Ec].mean() - y[x > Ec].mean(), Ec, 0.03 * Ec, y[x > Ec].mean(), 0.0]
        q, cv = curve_fit(step, x, y, p0=p0, sigma=np.sqrt(np.maximum(np.abs(y), 1)), maxfev=20000); res.append((q, np.sqrt(np.diag(cv))))
    (qd, ed), (qm, em) = res
    print(f"{nm:<17} | {qd[1]:7.1f}±{ed[1]:.1f} / {qm[1]:7.1f} ({qd[1]-qm[1]:+.1f}) | {abs(qd[2]):5.1f}±{ed[2]:.1f} / {abs(qm[2]):5.1f} (×{abs(qd[2]/qm[2]):.2f})"
          f" | {qd[0]:9.0f} / {qm[0]:9.0f} (×{qd[0]/qm[0]:.3f})")
