# -*- coding: utf-8 -*-
# #GS-70 (#CFG-1 «Да, по таблице» 10.10): шаблон М1 в шкале СВЕТА (NPSM Payne, GS2020_NPSM=1) — гистограмма count_light растянута
# k = E0/вершина света фотопика (как scripts\gs68_npsm_compare.py), перебинирована в 1 кэВ и записана в count_edep. argv: вход выход E0
import sys, os, numpy as np; sys.stdout.reconfigure(encoding="utf-8"); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import merge_templates_gs2020 as mt
IN, OUT, E0 = sys.argv[1], sys.argv[2], float(sys.argv[3]); ch = mt.read_chunk(IN)
x = np.array(ch["bins"], float); li = np.array(ch["light"], float); ed = np.array(ch["edep"], float)
r = (x > 0.3 * E0) & (x < 1.05 * E0); p0 = x[r][np.argmax(li[r])]; w = np.abs(x - p0) < 0.03 * p0
pl = (x[w] * li[w]).sum() / li[w].sum(); k = E0 / pl; si = np.sqrt((li[w] * (x[w] - pl) ** 2).sum() / li[w].sum()) * k
cum = np.concatenate([[0.0], np.cumsum(li)]); edges = np.concatenate([x - 0.5, [x[-1] + 0.5]])   # бины света [x-0,5; x+0,5)
new = np.diff(np.interp(edges, edges * k, cum)); cnt = np.rint(new).astype(int).tolist()
hdr = ch["header"] + [("gs70_light_scale_k", f"{k:.6f}"), ("gs70_light_sigma_intr_keV", f"{si:.3f}")]
mt.write_template(OUT, ch["comments"], hdr, ch["bins"], cnt, cnt)   # обе колонки = свет в шкале энергии: read_template при npsm_enabled=1 берёт count_light
print(f"вершина света {pl:.2f} → {E0} (k {k:.5f}); σ света в пике {si:.2f} кэВ; сумма свет {li.sum():.0f} → {sum(cnt)}, депозит {ed.sum():.0f}")
for lo, hi in ((150, 250), (250, 450), (450, 500), (500, 550), (600, 700)):
    m = (x >= lo) & (x < hi); print(f"  {lo}-{hi}: свет/депозит {np.array(cnt)[m].sum() / max(ed[m].sum(), 1):.3f}")
