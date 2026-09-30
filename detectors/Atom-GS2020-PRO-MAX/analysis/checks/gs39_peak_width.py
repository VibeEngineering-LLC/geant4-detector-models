# -*- coding: utf-8 -*-
# #GS-39 (оператор 29.09 «у тория пшпв 2614 теперь меньше чем нужно»): ПШПВ пиков в (изм.−фон) и в модели М1 на оси страницы.
# Гаусс + линия, центр и ширина свободны (curve_fit), окно ±1·ПШПВ закона. Печать: ширина изм., модели, отношение, центроиды.
import sys, json, numpy as np
from scipy.optimize import curve_fit
sys.stdout.reconfigure(encoding="utf-8")
P = r"D:\cloud-folder\work-folder\GEANT4\web\gs2020-th232-page\%s"
FW = lambda E: 0.7880 * E ** 0.6335
g = lambda x, a, mu, s, c, d: a * np.exp(-0.5 * ((x - mu) / s) ** 2) + c + d * (x - mu)
def fit(e, v, E0):
    s = np.abs(e - E0) < FW(E0); x, y = e[s], v[s]
    p, cv = curve_fit(g, x, y, p0=[y.max() - y.min(), E0, FW(E0) / 2.355, y.min(), 0.0], sigma=np.sqrt(np.maximum(np.abs(y), 1)))
    return p[1], 2.3548 * abs(p[2]), 2.3548 * np.sqrt(cv[2, 2])
for fn, tag in (("gs2020_th232_data.json", "торий"), ("gs2020_k40_data.json", "KCl")):
    d = json.load(open(P % fn, encoding="utf-8"))
    for key, sp in (("spectrum", d["spectrum"]), ("cs", d.get("cs", {}).get("spectrum"))):
        if not sp or "counts" not in sp: continue   # cs.spectrum несёт только модели второго фона
        e = np.asarray(sp.get("e_of_ch", d["spectrum"]["e_of_ch"])); net = np.asarray(sp["counts"], float) - np.asarray(sp["bg_counts"], float); m = np.asarray(sp["model_counts"], float)
        for E0 in ((238.63, 583.19, 911.2, 2614.51) if tag == "торий" else (1460.82,)):
            mn, wn, dw = fit(e, net, E0); mm, wm, _ = fit(e, m, E0)
            print("%-6s %-8s %7.1f: ПШПВ изм. %6.1f ± %4.1f, модель %6.1f, модель/изм. %.3f; центроид изм. %7.1f, модель %7.1f" % (tag, key, E0, wn, dw, wm, wm / wn, mn, mm))
