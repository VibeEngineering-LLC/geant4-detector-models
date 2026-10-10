# #GS-68: остатки изм/мод вокруг пика по JSON подгонки (e, net, model, var). argv: json E0 [полуокно=300] [шаг=20]
import sys, json, numpy as np
sys.stdout.reconfigure(encoding="utf-8")
d = json.load(open(sys.argv[1], encoding="utf-8")); E0 = float(sys.argv[2])
hw = float(sys.argv[3]) if len(sys.argv) > 3 else 300.0; st = float(sys.argv[4]) if len(sys.argv) > 4 else 20.0
e, n, m = (np.array(d[k], float) for k in ("e", "net", "model"))
v = np.array(d["var"], float) if "var" in d else np.abs(n) + 1.0  # нет var — σ по нетто (занижена: без фона)
print(f"{sys.argv[1]}  E0={E0}  A/ожид={d.get('ratio', float('nan'))}")
print("полоса, кэВ от E0 | изм | мод | изм/мод | (изм-мод)/σ")
for lo in np.arange(E0 - hw, E0 + hw, st):
    k = (e >= lo) & (e < lo + st); N, M, V = n[k].sum(), m[k].sum(), v[k].sum()
    print(f"{lo - E0:+5.0f}..{lo + st - E0:+5.0f} | {N:.0f} | {M:.0f} | {N / M if M > 0 else float('nan'):.3f} | {(N - M) / V ** 0.5 if V > 0 else 0:.1f}")
