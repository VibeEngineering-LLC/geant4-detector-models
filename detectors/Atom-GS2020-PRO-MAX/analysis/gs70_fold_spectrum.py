# -*- coding: utf-8 -*-
# #GS-70 (оператор 10.10: «тормозное отдельно»): шаблон М1 из фотонного спектра (k_keV, dN/dk на распад; напр. IB KUB ib_spectrum_kub.py)
# × отклик прибора по сетке монолиний grid_mar_E*.csv: вес узла = ∫dN/dk по его интервалу (середины между узлами); n = 1e9 распадов.
# Запуск: python gs70_fold_spectrum.py <спектр.csv> <папка сетки> <выходной csv> [kmin кэВ=20]
import sys, os, glob, numpy as np; sys.stdout.reconfigure(encoding="utf-8"); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import merge_templates_gs2020 as mt
SP, GD, OUT = sys.argv[1:4]; KMIN = float(sys.argv[4]) if len(sys.argv) > 4 else 20.0
k, y = np.array([l.split(",")[:2] for l in open(SP, encoding="utf-8") if l[:1].isdigit()], float).T
P = {float(os.path.basename(p)[10:-4]): p for p in glob.glob(os.path.join(GD, "grid_mar_E*.csv"))}
nodes = np.array(sorted({x for x in P if KMIN <= x <= k.max()})); mid = np.concatenate([[KMIN], 0.5 * (nodes[1:] + nodes[:-1]), [k.max() + 0.5]])
W = np.array([y[(k >= a) & (k < b)].sum() for a, b in zip(mid[:-1], mid[1:])])   # фотонов на распад в интервале узла (бины 1 кэВ)
acc = None; NORM = 10 ** 9
for x, wi in zip(nodes, W):
    ch = mt.read_chunk(P[x]); n = float(dict(ch["header"])["n_events_processed"])
    v = np.array(ch["edep"], float) * wi / n; acc = v if acc is None else acc + v; bins, hdr, com = ch["bins"], ch["header"], ch["comments"]
hdr = [(kk, (str(NORM) if kk == "n_events_processed" else vv)) for kk, vv in hdr] + [("gs70_folded", f"{os.path.basename(SP)} kmin={KMIN} nodes={len(nodes)}")]
cnt = np.rint(acc * NORM).astype(int).tolist(); mt.write_template(OUT, com, hdr, bins, cnt, cnt)
print(f"узлов {len(nodes)} ({nodes[0]:.1f}–{nodes[-1]:.1f} кэВ); фотонов на распад ≥{KMIN:g} кэВ {W.sum():.4e} (в спектре {y[k >= KMIN].sum():.4e}); отсчётов на 1e9 распадов {sum(cnt)}")
