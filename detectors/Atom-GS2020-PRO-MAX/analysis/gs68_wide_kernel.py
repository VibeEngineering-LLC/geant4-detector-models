# #GS-68: досвёртка модели ядром (1-w)·δ + w·Gauss(σ_w) и масштаб a в окне; подгонка w, σ_w, a к net. argv: json lo hi
import sys, json, numpy as np
from scipy.optimize import least_squares
sys.stdout.reconfigure(encoding="utf-8")
d = json.load(open(sys.argv[1], encoding="utf-8")); lo, hi = float(sys.argv[2]), float(sys.argv[3])
e, n, m = (np.array(d[k], float) for k in ("e", "net", "model"))
v = np.array(d["var"], float) if "var" in d else np.abs(n) + 1.0
big = (e > lo - 400) & (e < hi + 400); eb, mb = e[big], m[big]; win = (eb >= lo) & (eb <= hi)
nw, vw = n[big][win], v[big][win]
def conv(w, s):
    G = np.exp(-0.5 * ((eb[:, None] - eb[None, :]) / s) ** 2); G /= G.sum(axis=1, keepdims=True)
    return (1 - w) * mb + w * (G.T @ mb)  # перенос отсчётов из бина j по гауссу σ_w (отсчёты сохраняются)
def res(p): return ((p[0] * conv(p[1], p[2]))[win] - nw) / vw ** 0.5
r0 = res([1.0, 0.0, 30.0]); a0 = least_squares(lambda p: res([p[0], 0.0, 30.0]), [1.0]).x[0]
c0 = (res([a0, 0.0, 30.0]) ** 2).sum(); fit = least_squares(res, [a0, 0.1, 40.0], bounds=([0.5, 0, 3], [2, 1, 300]))
c1 = (fit.fun ** 2).sum(); k = win.sum()
print(f"{lo:.0f}-{hi:.0f}: без досвёртки a {a0:.4f} χ²/ν {c0/(k-1):.2f} | с досвёрткой a {fit.x[0]:.4f} w {fit.x[1]:.3f} σ_w {fit.x[2]:.1f} кэВ χ²/ν {c1/(k-3):.2f}")
for b in np.arange(lo, hi, 40):
    s = (eb[win] >= b) & (eb[win] < b + 40); M1 = (fit.x[0] * conv(fit.x[1], fit.x[2]))[win][s].sum()
    print(f"  {b:.0f}-{b+40:.0f} изм/мод_до {nw[s].sum()/(a0*mb[win][s].sum()):.3f}  изм/мод_после {nw[s].sum()/M1:.3f}")
