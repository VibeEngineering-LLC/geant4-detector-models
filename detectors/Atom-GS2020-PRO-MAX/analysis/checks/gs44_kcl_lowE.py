# -*- coding: utf-8 -*-
# #GS-44 (оператор 29.09: «у калия после вычета 32 вылез»): нетто KCl в мягкой зоне по массивам подгонки (e, net, model, bg),
# поиск пика 20–45 кэВ и проверка линии Cs-137 662 кэВ (если 32 кэВ — рентген Ba). Запуск: python gs44_kcl_lowE.py [папка]
import sys, json, numpy as np
sys.stdout.reconfigure(encoding="utf-8")
D = sys.argv[1] if len(sys.argv) > 1 else r"C:\g4work\gs2020\kcl_v3w85"
J = sys.argv[2] if len(sys.argv) > 2 else "fit_kcl_bgw.json"   # 30.09: имя JSON вторым аргументом (fit_kcl_bgw_r.json — фон × r(E), #GS-44)
j = json.load(open(D + "\\" + J, encoding="utf-8")); print("JSON:", D + "\\" + J)
e, net, mod, bg = (np.asarray(j[k], float) for k in ("e", "net", "model", "bg"))
raw = net + bg
print("каналов %d; e[0..3]=%s; окно подгонки %g–%g кэВ" % (len(e), np.round(e[:4], 2), j["lo"], j["hi"]))
for lo, hi in ((10, 20), (20, 26), (26, 30), (30, 34), (34, 38), (38, 45), (45, 60), (60, 90), (90, 150)):
    s = (e >= lo) & (e < hi)
    print("%3g–%3g кэВ: сырой %9.0f  фон×k %9.0f  нетто %9.0f  модель %9.0f  нетто/модель %6.3f" % (lo, hi, raw[s].sum(), bg[s].sum(), net[s].sum(), mod[s].sum(), net[s].sum() / mod[s].sum() if mod[s].sum() else np.nan))
s = (e >= 20) & (e < 45); i = np.flatnonzero(s)[np.argmax(net[s])]
print("максимум нетто 20–45 кэВ: канал %d, E=%.2f кэВ, нетто %.0f, модель %.0f" % (i, e[i], net[i], mod[i]))
def peak(E0, w):   # площадь над линейной подложкой по боковым окнам, ± погрешность
    c, l, r = (np.abs(e - E0) < w), (e >= E0 - 3 * w) & (e < E0 - 2 * w), (e > E0 + 2 * w) & (e <= E0 + 3 * w)
    base = (net[l].mean() + net[r].mean()) / 2 * c.sum(); a = net[c].sum() - base
    return a, np.sqrt(raw[c].sum() + bg[c].sum() + raw[l | r].sum())
for name, E0, w in (("Cs-137 662", 661.657, 35.0), ("Ba Kα 32", 32.0, 3.5), ("I Kα 28,5", 28.5, 3.0), ("Pb Kα 75", 75.0, 4.0)):
    a, da = peak(E0, w); print("%-12s площадь над подложкой в нетто: %9.0f ± %6.0f (%.1fσ)" % (name, a, da, a / da))
