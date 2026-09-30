# -*- coding: utf-8 -*-
# #GS-44 (30.09): чувствительность r к спектру падающего поля — отношение депозитов KCl/вода по полосам для КАЖДОЙ монолинии E
# (4π = up+down, без уширения, прямо из chunks/ext_*.csv). Запуск: python gs44_ext_perE.py [папка ext]
import sys, os, glob, re, numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from merge_templates_gs2020 import read_chunk
X = sys.argv[1] if len(sys.argv) > 1 else r"C:\g4work\gs2020\ext"
BANDS = [(20, 45), (45, 60), (60, 90), (90, 150), (150, 300), (300, 1000), (1000, 3000)]
H = {}
for p in glob.glob(os.path.join(X, "chunks", "ext_*_E*_*.csv")):
    g, E, h = re.match(r"ext_(water|kcl)_E([0-9.]+)_(up|down)\.csv$", os.path.basename(p)).groups()
    c = read_chunk(p); a = np.zeros(3700); a[[int(float(b)) for b in c["bins"]]] = c["edep"]
    H[(g, float(E))] = H.get((g, float(E)), 0) + a   # up + down = 4π, одинаковое число историй
print("E_пад, кэВ | " + " | ".join("%d–%d" % b for b in BANDS) + "   (KCl/вода по депозитам; «—» если в воде < 200 отсчётов)")
for E in sorted({e for _, e in H}):
    w, k = H[("water", E)], H[("kcl", E)]
    cells = ["%.3f±%.3f" % (k[l:h].sum() / w[l:h].sum(), k[l:h].sum() / w[l:h].sum() * np.sqrt(1 / k[l:h].sum() + 1 / w[l:h].sum())) if w[l:h].sum() >= 200 and k[l:h].sum() > 0 else "—" for l, h in BANDS]
    print("%8.1f | %s" % (E, " | ".join(cells)))
