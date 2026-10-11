# -*- coding: utf-8 -*-
# #GS-87: σ амплитуд метода 1 черники с дисперсией приведённого фона (сэндвич (XᵀWX)⁻¹XᵀW·V·WX(XᵀWX)⁻¹, V = var + bg·k_bg), окно 150–3000 кэВ.
# Запуск: python gs87_k40_bgnoise.py <fit_berry_bgw_*.json> <живое время фона, с>
import sys, json, numpy as np
j = json.load(open(sys.argv[1], encoding="utf-8")); Tb = float(sys.argv[2]); kb = j["live_s"] / Tb
e = np.array(j["e"]); sel = (e >= 150) & (e <= 3000); v = np.array(j["var"])[sel]; b = np.array(j["bg"])[sel]
names = list(j["nuclides"]); A = np.array([j["nuclides"][n]["A_Bq"] for n in names]); cols = np.array(j["cols"])[:, sel]
X = (cols / np.where(A == 0, 1, A)[:, None]).T; W = 1 / v
M = np.linalg.inv(X.T @ (X * W[:, None])); V2 = v + b * kb
C1 = M; C2 = M @ (X.T @ (X * (W * W * V2)[:, None])) @ M
for i, n in enumerate(names): print(f"{n}: A {A[i]:.2f}  σ(вес) {np.sqrt(C1[i,i]):.2f}  σ(+фон) {np.sqrt(C2[i,i]):.2f}  стат(JSON) {j['nuclides'][n]['dA_stat_Bq']:.2f}")
print("k_bg", round(kb, 4), "var включает фон?", np.allclose(v, np.maximum(np.array(j['net'])[sel] + b, 1)) )
