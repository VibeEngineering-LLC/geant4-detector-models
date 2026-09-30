# -*- coding: utf-8 -*-
r"""#GS-42: влияние приближения донора add_ib_template (дисперсия IB-части шаблона М1 = col/N_ion) на A метода 1 KCl —
та же подгонка критерием A2V с ТОЧНОЙ дисперсией col_mix/N_ion + col_IB/n_eff; контроль — A2V с приближённой дисперсией
обязан дать A2 побитово. Только чтение, JSON не пишет. Окружение — как у run_gs2020_kcl.sh (GS_OUT KCl, GS_BG_WATER=1)."""
import os, sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import fit_gs2020_kcl as fk
import gs2020_extra_components as gx
import crit_bench_m2 as cbm2   # путь analysis/ донора уже в sys.path (fit_gs2020_th232_m1)
m1 = fk.m1
s, b = m1.Spec(m1.bm.read(m1.CAL.PATHS["kcl"])[0], "kcl"), m1.Spec(m1.bm.read(m1.CAL.BKG_WATER_XML)[0], "bgw")
tpl, fw = os.path.join(fk.OUT, "mix_K40_npsmoff.csv"), os.path.join(fk.OUT, "fwhm_points_gs2020.csv")
comp = gx.load("K40", fk.OUT, beta=False)
bg_e = np.array([b.channel_to_energy(i) for i in range(b.n_channels)])
with gx.m1_ib(m1.muc.g1s, {tpl: [(comp, 1.0)]}):
    r = m1.muc.unfold(s, b, [("K40", tpl)], fw, lo=fk.LO, hi=fk.HI, bg_energy_of_ch=bg_e, verbose=False, blur=m1.BLUR, tail=m1.TAIL)
ib = gx.fold(comp, m1.muc.g1s.broaden, r["ch_edges"], lambda E: m1.BLUR * r["fwhm"](E))["ib"][0]
n, k, T = r["n_events"][0], s.live_time / b.live_time, s.live_time
A = {}
for lab, V in (("A2V приближённая (контроль)", r["cols"][0] / n), ("A2V точная", (r["cols"][0] - ib) / n + ib / comp["ib"]["n_eff"])):
    A[lab] = cbm2.fit_A2V(r["cols"], r["counts"], r["bg_scaled"], k, V[None, :], r["sel"])[0][0] / T
    print("%s: A %.4f Бк; относительно A2 подгонки %+.3e" % (lab, A[lab], A[lab] / r["activities"][0] - 1))
