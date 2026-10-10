# #GS-68: досвёртка модели разбросом усиления (1-w)·δ + w·N(0; s·E) СОВМЕСТНО по окнам (свой масштаб a у окна).
# argv: json w0 s0 lo:hi ... [fix] — с «fix» w,s фиксированы (проверка), подгоняются только a. Печать χ²/ν без/с по окнам.
import sys, json, numpy as np
from scipy.optimize import least_squares
sys.stdout.reconfigure(encoding="utf-8")
d = json.load(open(sys.argv[1], encoding="utf-8")); fix = "fix" in sys.argv; W = [tuple(map(float, a.split(":"))) for a in sys.argv[4:] if ":" in a]
e, n, m = (np.array(d[k], float) for k in ("e", "net", "model")); v = np.array(d["var"], float) if "var" in d else np.abs(n) + 1.0
S = []
for lo, hi in W:
    big = (e > lo - 400) & (e < hi + 400); eb, mb = e[big], m[big]; win = (eb >= lo) & (eb <= hi)
    S.append((eb, mb, win, n[big][win], v[big][win], (eb[None, :] - eb[:, None]) ** 2))
def conv(s, w, sr):  # строка j — перенос из бина j по гауссу σ = sr·E_j (отсчёты сохраняются)
    eb, mb, win, _, _, D2 = s; G = np.exp(-0.5 * D2 / (sr * eb[:, None]) ** 2); G /= G.sum(axis=1, keepdims=True)
    return ((1 - w) * mb + w * (G.T @ mb))[win]
def res(p): return np.concatenate([(p[k] * conv(s, p[-2], p[-1]) - s[3]) / s[4] ** 0.5 for k, s in enumerate(S)])
nW = len(S); ws = [float(sys.argv[2]), float(sys.argv[3])]; lb, ub = np.r_[np.full(nW, .5), 0, .002], np.r_[np.full(nW, 2), 1, .2]
f0 = least_squares(lambda a: res(np.r_[a, 0.0, 0.03]), np.ones(nW)).fun
if fix: f1 = least_squares(lambda a: res(np.r_[a, ws]), np.ones(nW)).fun
else: fit = least_squares(res, np.r_[np.ones(nW), ws], bounds=(lb, ub)); f1 = fit.fun; print(f"СОВМЕСТНО: w {fit.x[-2]:.4f}  s {100*fit.x[-1]:.3f} % от E")
i = 0
for (lo, hi), s in zip(W, S):
    k = s[2].sum(); print(f"{lo:.0f}-{hi:.0f}: χ²/ν без {(f0[i:i+k]**2).sum()/k:.2f} → с {(f1[i:i+k]**2).sum()/k:.2f}"); i += k
