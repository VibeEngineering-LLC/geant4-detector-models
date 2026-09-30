# -*- coding: utf-8 -*-
# #GS-38 (оператор 29.09 «что опять с калибровкой?», «фон ты к времени набора образца привел?»): линии ФОНА по отдельности в KCl и
# в фоне·k на ОДНОЙ оси страницы (фон переведён на ось KCl по энергии). Разные центроиды = разные шкалы; разные площади = разный фон.
# Подгонка каждого спектра: a·G(E−μ) + c + d·(E−E0), ширина по закону ПШПВ страницы, μ свободен (сетка ±8 кэВ, шаг 0,05).
import sys, json, numpy as np
sys.stdout.reconfigure(encoding="utf-8")
d = json.load(open(r"D:\cloud-folder\work-folder\GEANT4\web\gs2020-th232-page\gs2020_k40_data.json", encoding="utf-8"))
sp = d["spectrum"]; e = np.asarray(sp["e_of_ch"]); y = np.asarray(sp["counts"], float); bg = np.asarray(sp["bg_counts"], float)
FW = lambda E: 0.7880 * E ** 0.6335
def fit(v, E0):
    s = np.abs(e - E0) < 1.5 * FW(E0); x = e[s]; sg = FW(E0) / 2.3548; w = 1 / np.maximum(v[s], 1.0); best = None
    for mu in np.arange(E0 - 8, E0 + 8.001, 0.05):
        A = np.vstack([np.exp(-0.5 * ((x - mu) / sg) ** 2), np.ones_like(x), x - E0]).T
        p = np.linalg.solve(A.T @ (A * w[:, None]), A.T @ (v[s] * w)); c2 = np.sum(w * (v[s] - A @ p) ** 2)
        if best is None or c2 < best[0]: best = (c2, mu, p[0] * sg * 2.5066)
    return best[1], best[2]
print("линия      центроид KCl  центроид фона  разница   площадь KCl  площадь фона·k  отношение")
for E0 in (238.63, 295.22, 351.93, 583.19, 609.31, 911.2, 1120.29, 1764.49, 2614.51):
    mk, ak = fit(y, E0); mb, ab = fit(bg, E0)
    print("%7.1f   %10.2f   %12.2f   %+6.2f   %11.0f   %13.0f   %8.3f" % (E0, mk, mb, mk - mb, ak, ab, ak / ab if ab else float("nan")))
