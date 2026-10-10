# -*- coding: utf-8 -*-
# #GS-70 черника (оператор 10.10: «Определяй активности. Это первый боевой разбор»): метод 1 на 3 шаблона — Cs-137 (до Ba-137m), K-40,
# Sr-90+Y-90 (1:1, одна амплитуда, gs70_sry_template.py) + фон воды × r(E) черники (GS_BG_R). Ширина — карточка (#DET-1, strict).
# Долив β K-40 (#GS-45, +3,6 % β) не применён: β-шаблона в геометрии черники нет. Запуск: GS_OUT=<шаблоны> GS_CAL_SHAPE=1 GS_BG_R=<csv> python fit_gs2020_berry.py
import os, sys, json; import numpy as np; sys.stdout.reconfigure(encoding="utf-8"); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fit_gs2020_th232_m1 as m1, fit_gs2020_kcl as fk; CAL = m1.CAL; fk.LO = float(os.environ.get("GS_LO", fk.LO)); TAG = os.environ.get("GS_BERRY_TAG", "")
MASS_KG = float(os.environ.get("GS_BERRY_MASS", "482")) / 1000; CS_SCALE = 85.131 / 85.10   # Iγ 662 в шаблоне Geant4 11.2 (RadioactiveDecay5.6 z55.a137 94,699 %; PhotonEvaporation5.7 α 0,1124) / LNHB 85,10
NUC = [x for x in [("Cs137", CS_SCALE), ("Cs137L", CS_SCALE), ("K40", fk.K40_SCALE), ("SrY90", 1.0), ("Scat662", 1.0), ("IBSrY90", 1.0)] if x[0] in os.environ.get("GS_BERRY_NUC", "Cs137,K40,SrY90").split(",")]   # SrY90: A — распадов Sr-90 в с
s = m1.Spec(m1.bm.read(CAL.PATHS["berry"])[0], "berry"); b = m1.Spec(m1.bm.read(CAL.BKG_WATER_XML)[0], "bgw")
fw = os.path.join(m1.OUT, "fwhm_points_gs2020.csv"); m1.write_fwhm_csv(fw, strict=True)
r = m1.muc.unfold(s, b, [(n, os.path.join(m1.OUT, f"mix_{n}_npsmoff.csv")) for n, _ in NUC], fw, lo=fk.LO, hi=fk.HI, verbose=False, blur=m1.BLUR,
                  tail=m1.TAIL, tvar=m1.TVAR, bg_energy_of_ch=np.array([b.channel_to_energy(i) for i in range(b.n_channels)]))
model = (r["cols"] * r["coef"][:, None]).sum(axis=0); sel, net, var, e, L = r["sel"], r["net"], r["var"], r["e"], s.live_time
chi2 = float(np.sum((model[sel] - net[sel]) ** 2 / var[sel])); nd = int(sel.sum()) - len(NUC); fd = np.genfromtxt(fw, delimiter=",", invalid_raise=False)
FW = lambda E: float(np.interp(E, fd[:, 0], fd[:, 1])); m1.SHAPE_PEAKS = [661.657, 1460.822]
shp, shp3 = m1.shape_residual(e, net, var, model, FW)[0], m1.shape_residual(e, net, var, model, FW, 3.0)[0]
print("ФОРМА ПИКОВ (#SHAPE-1, χ²/ν ±1,5 / ±3 ПШПВ): " + "; ".join(f"{p:.1f} → {shp[p]:.2f} / {shp3[p]:.2f}" for p in m1.SHAPE_PEAKS))
res = {n: {"A_Bq": float(r["activities"][i]) * k, "dA_stat_Bq": float(r["sd"][i]) / L * k} for i, (n, k) in enumerate(NUC)}
for n, v in res.items(): print(f"{n} (метод 1): A {v['A_Bq']:.2f} ± {v['dA_stat_Bq']:.2f} (стат) Бк; {v['A_Bq'] / MASS_KG:.1f} ± {v['dA_stat_Bq'] / MASS_KG:.1f} Бк/кг")
print(f"χ²/ν {chi2 / nd:.3f} (ν {nd}), окно {fk.LO:g}–{fk.HI:g} кэВ, живое {L:.1f} с, k_bg {L / b.live_time:.4f}, GS_BG_R {os.environ.get('GS_BG_R')}")
json.dump({"nuclides": res, "mass_kg": MASS_KG, "chi2": chi2, "ndof": nd, "shape": {str(p): [shp[p], shp3[p]] for p in shp}, "live_s": L, "e": np.asarray(e).tolist(),
           "net": np.asarray(net).tolist(), "model": model.tolist(), "bg": np.asarray(r["bg_scaled"]).tolist(), "var": np.asarray(var).tolist(),
           "cols": (r["cols"] * r["coef"][:, None]).tolist()}, open(os.path.join(m1.OUT, f"fit_berry_bgw{TAG}.json"), "w", encoding="utf-8"), ensure_ascii=False)
