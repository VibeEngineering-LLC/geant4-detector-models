# #GS-60: модель/измерение по полосам от 15 кэВ на данных страницы GS2020 (Th-232 и KCl).
# T = (изм − модель)/фон — пропускание фона, при котором модель сходится с измерением.
import json, sys, numpy as np
sys.stdout.reconfigure(encoding="utf-8")
D = r"<WORKDIR>\GEANT4\web\gs2020-th232-page\dist"
BANDS = [(15, 25), (25, 35), (35, 50), (50, 70), (70, 100), (100, 150), (150, 200), (200, 300), (300, 500), (500, 1000), (1600, 2500)]
for fn, pre in (("data.js", "window.G1S="), ("data-k40.js", "window.GS_K40=")):
    d = json.loads(open(D + "\\" + fn, encoding="utf-8").read().strip()[len(pre):].rstrip(";"))
    s = d["spectrum"]; E = np.array(s["e_of_ch"], float)
    c, bg, m1, m2 = (np.array(s[k], float) for k in ("counts", "bg_counts", "model_counts", "model2_counts"))
    print(f"== {fn}: {d['meta'].get('template_source')}, матрица {d['meta'].get('matrix_name')} ρ {d['meta'].get('matrix_density_g_cm3')}")
    print(" полоса кэВ    изм        фон     изм-фон   (м1+фон)/изм  м1/(изм-фон)  м2/(изм-фон)  T(м1)  T(м2)")
    for a, b in BANDS:
        k = (E >= a) & (E < b); C, B, M, M2 = c[k].sum(), bg[k].sum(), m1[k].sum(), m2[k].sum()
        print(f" {a:5}-{b:<5} {C:10.0f} {B:10.0f} {C-B:10.0f}   {(M+B)/C:8.3f}    {M/(C-B):9.2f}    {M2/(C-B):9.2f}   {(C-M)/B:6.3f} {(C-M2)/B:6.3f}")
