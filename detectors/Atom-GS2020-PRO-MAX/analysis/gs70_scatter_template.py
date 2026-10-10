# -*- coding: utf-8 -*-
# #GS-70 (оператор 10.10: «пузо пока выдели отдельным шаблоном»): мешающий шаблон однократно рассеянных γ E0 вне модели —
# спектр E'(θ) по Клейну–Нишине (изотропно по Ω) × отклик прибора по сетке монолиний grid_mar_E*.csv (узлы E0/(1+2ε)…E0).
# Запуск: python gs70_scatter_template.py <E0 кэВ> <папка сетки> <выходной csv>. Амплитуда — мешающая, не активность.
import sys, os, glob, numpy as np; sys.stdout.reconfigure(encoding="utf-8"); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import merge_templates_gs2020 as mt
E0, GD, OUT = float(sys.argv[1]), sys.argv[2], sys.argv[3]; eps = E0 / 510.999; Emin = E0 / (1 + 2 * eps)
P = {float(os.path.basename(p)[10:-4]): p for p in glob.glob(os.path.join(GD, "grid_mar_E*.csv"))}; nodes = sorted(P)
nodes = sorted({x for x in nodes if Emin <= x <= E0 - 1.5}); Emin = min(Emin, nodes[0]) - 1e-6
c = np.linspace(-1, 1, 200001); Ep = E0 / (1 + eps * (1 - c)); q = Ep / E0; w = q * q * (q + 1 / q - (1 - c * c))   # dσ/dΩ (Клейн–Нишина)
mid = np.concatenate([[Emin], 0.5 * (np.array(nodes[1:]) + np.array(nodes[:-1])), [E0]]); W = np.histogram(Ep, mid, weights=w)[0]; W /= W.sum()
acc = None; NORM = 10 ** 9
for x, wi in zip(nodes, W):
    ch = mt.read_chunk(P[x]); n = float(dict(ch["header"])["n_events_processed"])
    v = np.array(ch["edep"], float) * wi / n; acc = v if acc is None else acc + v; bins, hdr, com = ch["bins"], ch["header"], ch["comments"]
hdr = [(k, (str(NORM) if k == "n_events_processed" else v)) for k, v in hdr] + [("gs70_scatter", f"KN single scatter E0={E0} nodes={len(nodes)}")]
cnt = np.rint(acc * NORM).astype(int).tolist(); mt.write_template(OUT, com, hdr, bins, cnt, cnt)
print(f"узлов {len(nodes)} ({nodes[0]:.1f}–{nodes[-1]:.1f} кэВ), доля веса ниже 300 кэВ {W[np.array(nodes) < 300].sum():.3f}; отсчётов на 1e9 рассеяний {sum(cnt)}")
