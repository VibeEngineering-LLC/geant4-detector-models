# -*- coding: utf-8 -*-
# #GS-44 (29.09): измеренное ослабление внешнего фона пробой — площади линий фона в сыром спектре KCl (raw = net + bg)
# против фона воды x k (bg), обе в шкале KCl и в одном живом времени. Независимо от МК. Запуск: python gs44_bg_lines.py [папка]
import sys, json, numpy as np
sys.stdout.reconfigure(encoding="utf-8")
D = sys.argv[1] if len(sys.argv) > 1 else r"C:\g4work\gs2020\kcl_v4w85_83"
j = json.load(open(D + r"\fit_kcl_bgw.json", encoding="utf-8"))
e, net, bg, k = np.asarray(j["e"], float), np.asarray(j["net"], float), np.asarray(j["bg"], float), float(j["k_bg"])
raw = net + bg
def area(y, var, E0, w):   # площадь над линейной подложкой (боковые окна 2w..3w), дисперсия
    c, l, r = np.abs(e - E0) < w, (e >= E0 - 3 * w) & (e < E0 - 2 * w), (e > E0 + 2 * w) & (e <= E0 + 3 * w)
    nc, nl, nr = c.sum(), l.sum(), r.sum()
    a = y[c].sum() - nc * (y[l].mean() + y[r].mean()) / 2
    return a, var[c].sum() + (nc / 2) ** 2 * (var[l].sum() / nl ** 2 + var[r].sum() / nr ** 2)
print("k_bg=%.4f, живое KCl %.0f с" % (k, j["live_s"]))
for name, E0, w in (("Pb-212 239", 238.632, 14), ("Pb-214 352", 351.93, 18), ("Bi-214 609", 609.312, 28), ("Bi-214 1120", 1120.29, 40), ("Bi-214 1764", 1764.49, 50), ("Tl-208 2614", 2614.51, 65)):
    ar, vr = area(raw, raw, E0, w); ab, vb = area(bg, k * bg, E0, w)
    rt = ar / ab; drt = abs(rt) * np.sqrt(vr / ar ** 2 + vb / ab ** 2) if ar > 0 and ab > 0 else np.nan
    print("%-12s KCl %9.0f ± %6.0f   вода×k %8.0f ± %5.0f   KCl/вода %.3f ± %.3f" % (name, ar, np.sqrt(vr), ab, np.sqrt(vb), rt, drt))
