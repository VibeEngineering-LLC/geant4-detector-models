# -*- coding: utf-8 -*-
"""Мягкая зона GS2020: какой сдвиг δ (кэВ) шкалы фона и масштаб a лучше описывают измерение: изм(E) ≈ a·фон(E−δ)·t_изм/t_фон
+ b + c·E в окне [lo, hi]. a≈1 при δ≠0 — сдвиг шкал; δ≈0 при a<1 — ослабление фона пробой. Запуск: python gs20_lowE_shift.py kcl bgw 35 110"""
import os, sys; import numpy as np
sys.stdout.reconfigure(encoding="utf-8"); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cal"))
import gs2020_calib as C
def spec(tag):   # плотность счёта на кэВ по своей шкале спектра и живое время
    d = C.read(tag); c = np.asarray(d["counts"], float); e = C.energy_axis(tag, len(c)); return e, c / np.gradient(e), d["live_time"]
(em, ym, tm), (eb, yb, tb) = spec(sys.argv[1]), spec(sys.argv[2]); lo, hi = float(sys.argv[3]), float(sys.argv[4])
w = (em > lo) & (em < hi); x, y = em[w], ym[w]; sig = np.sqrt(np.maximum(C.read(sys.argv[1])["counts"][w], 1)) / np.gradient(em)[w]
res = []
for d in np.arange(-15, 15.001, 0.25):
    A = np.vstack([np.interp(x - d, eb, yb) * tm / tb, np.ones_like(x), x - x.mean()]).T / sig[:, None]
    k, *_ = np.linalg.lstsq(A, y / sig, rcond=None); res.append((float(np.sum((A @ k - y / sig) ** 2) / (w.sum() - 3)), float(d), float(k[0])))
at0 = [r for r in res if r[1] == 0.0][0]; best = min(res)
print("%s против %s, окно %g–%g кэВ: δ=0 → χ²/ν %.2f, a %.3f | лучший δ %+.2f кэВ → χ²/ν %.2f, a %.3f"
      % (sys.argv[1], sys.argv[2], lo, hi, at0[0], at0[2], best[1], best[0], best[2]))
