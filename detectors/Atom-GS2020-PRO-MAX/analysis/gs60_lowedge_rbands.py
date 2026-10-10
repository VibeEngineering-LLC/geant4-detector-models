# #GS-60: (м1 + r(E)·фон)/изм по полосам от 15 кэВ — фон воды × r(E) внешнего поля (gs44_ext_ratio.py, GS_EXT_SAMPLE=<тег>).
# Подгонка не переделывается (м1 со страницы): выше 150 кэВ r ≈ 1, поэтому оценка мягкой зоны прямая.
import json, sys, numpy as np
sys.stdout.reconfigure(encoding="utf-8")
D = r"<WORKDIR>\GEANT4\web\gs2020-th232-page\dist"
X = r"C:\g4work\gs2020\ext"
BANDS = [(15, 25), (25, 35), (35, 50), (50, 70), (70, 100), (100, 150), (150, 200), (200, 300), (300, 500), (15, 150)]
for fn, pre, tags in (("data.js", "window.G1S=", ["th"]), ("data-k40.js", "window.GS_K40=", ["kcl1l", "kcl"])):
    s = json.loads(open(D + "\\" + fn, encoding="utf-8").read().strip()[len(pre):].rstrip(";"))["spectrum"]
    E = np.array(s["e_of_ch"], float); c, bg, m1 = (np.array(s[k], float) for k in ("counts", "bg_counts", "model_counts"))
    R = {}
    for t in tags:
        a = np.loadtxt(X + "\\r_ext_4pi" + ("" if t == "kcl" else "_" + t) + ".csv", delimiter=",", skiprows=1)
        R[t] = np.interp(E, a[:, 0], a[:, 1])
    print(f"== {fn}: (м1 + фон·r)/изм; без r — столбец «r=1»")
    print("  полоса кэВ     r=1   " + "".join(f"{t:>9}" for t in tags) + "   нетто изм/м1 при r(" + tags[0] + ")")
    for a, b in BANDS:
        k = (E >= a) & (E < b); C = c[k].sum(); M = m1[k].sum()
        row = "".join(f"{(M + (bg[k] * R[t][k]).sum()) / C:9.3f}" for t in tags)
        net = (C - (bg[k] * R[tags[0]][k]).sum()) / M
        print(f"  {a:5}-{b:<5} {(M + bg[k].sum()) / C:7.3f}  {row}   {net:8.3f}")
