# Дрейф KCl: полосы модель/измерение (как gs60_lowedge_bands.py) по выгрузкам K-40 нескольких частей рядом. Аргументы: метка=путь_к_json
import json, sys, numpy as np
sys.stdout.reconfigure(encoding="utf-8")
BANDS = [(15, 25), (25, 35), (35, 50), (50, 70), (70, 100), (100, 150), (150, 200), (200, 300), (300, 500), (500, 1000), (1600, 2500)]   # = gs60_lowedge_bands.py:6
cols = []
for a in sys.argv[1:]:
    tag, path = a.split("=", 1)
    s = json.load(open(path, encoding="utf-8"))["spectrum"]
    E, c, bg, m1, m2 = (np.array(s[k], float) for k in ("e_of_ch", "counts", "bg_counts", "model_counts", "model2_counts"))
    cols.append((tag, [((m1[k].sum() + bg[k].sum()) / c[k].sum(), (m2[k].sum() + bg[k].sum()) / c[k].sum()) for k in ((E >= lo) & (E < hi) for lo, hi in BANDS)]))
print(" полоса кэВ  " + "".join("%-18s" % ("%s м1/м2" % t) for t, _ in cols) + "   (модель+фон)/изм")
for i, (lo, hi) in enumerate(BANDS):
    print(" %5d-%-5d " % (lo, hi) + "".join("%7.3f %7.3f    " % v[i] for _, v in cols))
