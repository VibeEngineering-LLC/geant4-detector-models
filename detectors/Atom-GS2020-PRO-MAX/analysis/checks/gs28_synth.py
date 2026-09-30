# -*- coding: utf-8 -*-
"""#GS-28 #DBG-1: синтетика — табличные линии групп + сдвиг SHIFT + континуум, пуассон; скан δ обязан вернуть SHIFT."""
import sys, os, numpy as np
sys.stdout.reconfigure(encoding="utf-8"); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "cal"))
import gs2020_calib_mplet as M
SHIFT = float(sys.argv[1]) if len(sys.argv) > 1 else 3.0
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1); x_all = M.C.energy_axis("bg")
for g in M.GROUPS_DEF:
    gl = [l for l in (M.get_lines_in_window(*c, g["lo"], g["hi"]) for c in g["comp"]) if l]
    m = (x_all >= g["lo"]) & (x_all <= g["hi"]); x = x_all[m]; dE = np.gradient(x)
    mu = 2000.0 * (1 - 0.3 * (x - x.mean()) / (x[-1] - x[0]))
    for lines in gl:
        for E, I in lines: mu += 3000.0 * I * dE * np.exp(-0.5 * ((x - E - SHIFT) / M.sigma(E)) ** 2) / (M.sigma(E) * M.CONST_2SQRTPI)
    y = rng.poisson(mu).astype(float); ds = np.arange(-8, 8.05, 0.05)
    c2 = [(M.solve_for_delta(x, y, gl, g["cont"], d)[1] or np.inf) for d in ds]
    print(f"{g['name']}: δ̂ = {ds[int(np.argmin(c2))]:+.2f} (задано {SHIFT:+.2f}), χ²/ν(min) {min(c2):.3f}, компонент {len(gl)}")
