# -*- coding: utf-8 -*-
"""GS-24: KCl 740 мл 829 г — метод 1 (шаблон полного распада K-40) + фон Маринелли+вода. Запуск: python fit_gs2020_kcl.py (GS_OUT — папка с mix_K40_npsmoff.csv)"""
import os, sys, json, math; import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fit_gs2020_th232_m1 as m1; CAL = m1.CAL
import gs2020_extra_components as gx   # #GS-42: IB K-40 в шаблон М1 (β уже в распаде иона); GS_EXTRA=0 — прежнее поведение
MASS_G = 829.0; VOL_ML = 740.0
K_FRAC = 39.0983 / 74.551
# #GS-45 (оператор 30.09.2026: база выходов LNHB-DDEP 2025, LNHB K-40_tables.pdf с.1-2). GS_K40_DB=ensdf — прежние константы.
K40_DB = os.environ.get("GS_K40_DB", "lnhb")
IG_TEMPLATE = 10.66   # выход γ 1460 внутри шаблона Geant4 (ENSDF: 8,1377+2,2291+0,29223 %, z19.a40)
IG_K40, ABUND_K40, T12_K40_Y = (10.34, 1.1668e-4, 1.2522e9) if K40_DB == "lnhb" else (10.66, 1.17e-4, 1.248e9)
K40_SCALE = IG_TEMPLATE / IG_K40   # A_истинная = A_подгонки · 10,66/Iγ: шаблон М1 несёт 10,66 % внутри
K40_BQ_PER_G_K = 6.02214076e23 / 39.0983 * ABUND_K40 * math.log(2) / (T12_K40_Y * 365.25 * 86400)
A_EXP = MASS_G * K_FRAC * K40_BQ_PER_G_K; LO, HI = 150.0, 3000.0; OUT = m1.OUT
# #GS-44 (29.09): GS_BG_R=<csv E_keV,r,dr> — фон воды × r(E), отношение откликов внешнего поля KCl/вода (gs44_ext_ratio.py);
# JSON с суффиксом _r, основной не перезаписывается. Дисперсия фона остаётся k·bg (при r<1 — завышена в 1/r, консервативно).
BG_R = os.environ.get("GS_BG_R"); SUF = "_r" if BG_R else ""
def apply_bg_r(b):
    if not BG_R: return
    t = np.genfromtxt(BG_R, delimiter=",", names=True); rr = lambda E: float(np.interp(E, t["E_keV"], t["r"]))
    b.counts = [c * rr(b.channel_to_energy(i)) for i, c in enumerate(b.counts)]
    print("ФОН × r(E) (#GS-44): %s; r(30)=%.3f r(60)=%.3f r(200)=%.3f r(1460.8)=%.3f r(2614.5)=%.3f" % (BG_R, *(rr(x) for x in (30, 60, 200, 1460.8, 2614.5))))
def main():
    s = m1.Spec(m1.bm.read(CAL.PATHS["kcl"])[0], "kcl")
    b = m1.Spec(m1.bm.read(CAL.BKG_WATER_XML)[0], "bgw"); apply_bg_r(b)
    fwhm_csv = os.path.join(OUT, "fwhm_points_gs2020.csv"); m1.write_fwhm_csv(fwhm_csv)
    tpl = os.path.join(OUT, "mix_K40_npsmoff.csv")
    if not os.path.exists(tpl): raise SystemExit("ОТКАЗ: нет шаблона " + tpl)
    bg_e = np.array([b.channel_to_energy(i) for i in range(b.n_channels)])
    comp = gx.load("K40", OUT, beta=False) if gx.ENABLED else None
    if not gx.ENABLED: print("ФИЗИКА (#GS-42): GS_EXTRA=0 — IB выключен (прежнее поведение)")
    with gx.m1_ib(m1.muc.g1s, {tpl: [(comp, 1.0)]} if comp else {}):
        r = m1.muc.unfold(s, b, [("K40", tpl)], fwhm_csv, lo=LO, hi=HI, bg_energy_of_ch=bg_e, verbose=False, blur=m1.BLUR, tail=m1.TAIL)
    A = float(r["activities"][0]) * K40_SCALE; dA = float(r["sd"][0]) / s.live_time * K40_SCALE
    model = (r["cols"] * r["coef"][:, None]).sum(axis=0); sel = r["sel"]; net = r["net"]; var = r["var"]; e = r["e"]
    extra = None
    if comp:   # вклад IB в модель (тот же столбец, что вошёл в шаблон через m1_ib) и оценка приближения дисперсии донора
        ib = gx.fold(comp, m1.muc.g1s.broaden, r["ch_edges"], lambda E: m1.BLUR * r["fwhm"](E))["ib"][0]
        exc = gx.m1_var_excess([(1.0, ib, comp["ib"]["n_eff"])], r["n_events"][0], r["coef"][0], var, sel)
        # ib_col: поканальный столбец IB K-40 на сетке e/model (та же величина, что фолдится в band_counts ниже) —
        # для слоя IB на странице (GS-42, п.2); БЕЗ округления, чтобы Σ(ib_col[sel]) == band_counts["window"].
        ib_col = ib * r["coef"][0]
        extra = {"K40": {"br": 1.0, "ib": gx.band_counts(ib_col, e, sel), "ib_col": [float(x) for x in ib_col],
                         "ib_file": comp["ib"]["file"],
                         "Y": comp["ib"]["Y"], "drawn": comp["ib"]["drawn"]}, "amplitude": "K40, метод 1", "var_excess_max": exc}
        print("IB K-40 (#GS-42, метод 1): в окне %.0f отсчётов (%.2f %% модели в окне); приближение дисперсии донора ≤ %.1e дисперсии канала"
              % (extra["K40"]["ib"]["window"], 100 * extra["K40"]["ib"]["window"] / model[sel].sum(), exc))
    chi2 = float(np.sum((model[sel] - net[sel]) ** 2 / var[sel])); ndof = int(sel.sum()) - 1
    fwhm_data = np.genfromtxt(fwhm_csv, delimiter=',', invalid_raise=False)
    if fwhm_data.ndim == 1: fwhm_data = np.array([[fwhm_data[0], fwhm_data[1]]])
    def FWHM(E): return float(np.interp(E, fwhm_data[:, 0], fwhm_data[:, 1]))
    m1.SHAPE_PEAKS = [1460.822]; shp, _ = m1.shape_residual(e, net, var, model, FWHM)
    print(f"ФОРМА ПИКОВ (#SHAPE-1, χ²/ν в окне ±1,5 ПШПВ): 1460.8 → {shp[1460.822]:.2f}")
    print(f"K-40 (метод 1): A {A:.1f} ± {dA:.1f} (стат) Бк; ожидается по массе {A_EXP:.1f} Бк; отношение {A / A_EXP:.4f}; χ²/ν {chi2 / ndof:.3f}")
    # #GS-31 (оператор 28.09 «отклонился от МЕТОДА. никаких площадей, только метод 1 и 2»): пересчёт по площади пика снят
    with open(os.path.join(OUT, "fit_kcl_bgw%s.json" % SUF), "w", encoding="utf-8") as f:
        json.dump(dict({"extra_components": extra} if extra else {}, **{"A_Bq": A, "dA_stat_Bq": dA, "A_expected_Bq": A_EXP, "ratio": A / A_EXP, "chi2": chi2, "ndof": ndof, "shape_1460": shp[1460.822], "live_s": s.live_time, "k_bg": s.live_time / b.live_time, "lo": LO, "hi": HI, "e": e.tolist(), "net": net.tolist(), "model": model.tolist(), "bg": np.asarray(r["bg_scaled"]).tolist(), "var": np.asarray(var).tolist()}), f, ensure_ascii=False)   # var — дисперсия критерия A2 (с дисперсией шаблона), для пересчёта #SA-10
if __name__ == "__main__": main()
