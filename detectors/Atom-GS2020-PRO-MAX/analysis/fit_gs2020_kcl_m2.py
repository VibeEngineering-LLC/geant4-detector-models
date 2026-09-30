# -*- coding: utf-8 -*-
r"""#GS-24: метод 2 для K-40 в KCl (линия библиотеки × прямой моно-γ прогон grid_mar_E1460.822.csv) + фон Маринелли+вода.
Механика донора Gamma-1S/web-th232 (grid_response, run_method2, fit_A2V). Запуск: GS_OUT=C:\g4work\gs2020\kcl python fit_gs2020_kcl_m2.py"""
import sys, os, json, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
DONOR = r"D:\repos-folder\repos\geant4-detector-models\detectors\Gamma-1S\web-th232"
os.environ.setdefault("G4MODELS_SOURCE_CONFIG", os.path.join(HERE, "configs", "th232_gs2020_full_xray.yaml"))  # нужен донору при импорте
os.environ.setdefault("SPECTRAVIBE_ROOT", r"D:\cloud-folder\Дозиметрия\ИИ\1 Скилы\0_Work\gamma-spectrum-analysis")
sys.path.insert(0, HERE)
import fit_gs2020_th232_m1 as m1
import fit_gs2020_kcl as fk
sys.path.insert(0, DONOR)
import export_amticseu_data as eam
import export_data as ed
import export_ra226_data as erd
import mix_unfold_g1s as g1s
import crit_bench_m2 as cbm2
if m1.TAIL is not None:
    g1s.TAIL_T = m1.TAIL
eam.BUILD_OUT = m1.OUT
ed.E_FIT_HI = fk.HI

import gs2020_extra_components as gx   # #GS-42: β + IB K-40 в модель нуклида; GS_EXTRA=0 — прежнее поведение
PARTS = os.environ.get("GS_EXTRA_PARTS", "beta,ib,xr").split(",")   # #GS-45: какие добавки включены (при GS_EXTRA≠0)
IB_SCALE = float(os.environ.get("GS_IB_SCALE", "1"))                # #GS-45: множитель IB (вилка неопределённости)

KEYS = ["K40"]
LIB = [(1460.82, fk.IG_K40, "K40", "LNHB-DDEP 2025 (K-40_tables.pdf с.2): P(γ+ce) 10,34(7) %" if fk.K40_DB == "lnhb"
        else "ENSDF, IAEA LiveChart API (J. Chen 2015): 10,66(13) %")]   # #GS-45: база выходов — GS_K40_DB
# #XR-1 (#GS-42): K/L-рентген Ar — дочернее ядро K-40 при электронном захвате; на 100 распадов, ENSDF через IAEA LiveChart
# API (J. Chen 2015). Узлы сетки grid_mar_E0.265/2.956/2.958/3.19.csv в папке KCl (прямые прогоны).
LIB_AR_XR = [(0.265, 0.02363, "K40", "ENSDF, IAEA LiveChart API (J. Chen 2015): L-рентген Ar (ЭЗ K-40), #XR-1"),
             (2.956, 0.29564, "K40", "ENSDF, IAEA LiveChart API (J. Chen 2015): KA2 Ar (ЭЗ K-40), #XR-1"),
             (2.958, 0.58555, "K40", "ENSDF, IAEA LiveChart API (J. Chen 2015): KA1 Ar (ЭЗ K-40), #XR-1"),
             (3.19, 0.09508, "K40", "ENSDF, IAEA LiveChart API (J. Chen 2015): KB Ar (ЭЗ K-40), #XR-1")]
IG_LNHB = 10.34

def main():
    s = m1.Spec(m1.bm.read(m1.CAL.PATHS["kcl"])[0], "kcl")
    b = m1.Spec(m1.bm.read(m1.CAL.BKG_WATER_XML)[0], "bgw"); fk.apply_bg_r(b)   # #GS-44: GS_BG_R — фон × r(E)

    fwhm_csv = os.path.join(m1.OUT, "fwhm_points_gs2020.csv")
    m1.write_fwhm_csv(fwhm_csv)
    fwhm = g1s.make_fwhm(fwhm_csv)

    tpl = os.path.join(m1.OUT, "mix_K40_npsmoff.csv")
    if not os.path.exists(tpl):
        raise SystemExit("ОТКАЗ: нет шаблона " + tpl)

    bg_e = np.array([b.channel_to_energy(i) for i in range(b.n_channels)])
    r = m1.muc.unfold(s, b, [("K40", tpl)], fwhm_csv, lo=fk.LO, hi=fk.HI, bg_energy_of_ch=bg_e, verbose=False, blur=m1.BLUR, tail=m1.TAIL)

    live = s.live_time
    k_bg = s.live_time / b.live_time

    resp, n_nodes = eam.grid_response(r["ch_edges"], fwhm, e=r["e"], file_e=lambda c: float(s.channel_to_energy(int(c))))
    print(f"Число узлов сетки: {n_nodes}")

    def n_of(E):
        keys = list(resp.n_map.keys())
        closest = min(keys, key=lambda k: abs(k - E))
        if abs(closest - E) > 0.01:
            raise SystemExit(f"ОТКАЗ: нет узла сетки {E:.3f} кэВ")
        return resp.n_map[closest]

    lib = LIB + (LIB_AR_XR if gx.ENABLED and "xr" in PARTS else [])
    print(f"Строк-рентгена в библиотеке (#XR-1): {sum('#XR-1' in ln[3] for ln in lib)}")
    var2 = {}
    _, by_nuc_w, lines_m2, n_sum = erd.run_method2(lib, [], resp, r["e"], r["ch_edges"], KEYS, var_acc=var2, n_of=n_of)

    W2 = np.array([by_nuc_w["K40"]])
    V2 = np.array([var2["K40"]])
    cnt = np.asarray(s.counts, dtype=float)
    sel = r["sel"]
    # #GS-42: β + IB K-40 (на распад) — в модель нуклида с той же амплитудой, что его линии, после той же свёртки (М1)
    xc = gx.fold(gx.load("K40", m1.OUT), g1s.broaden, r["ch_edges"], lambda E: m1.BLUR * fwhm(E)) if gx.ENABLED else None
    if xc:
        wb, wi = float("beta" in PARTS), IB_SCALE * float("ib" in PARTS)   # #GS-45: разложение по компонентам
        W2[0] += wb * xc["beta"][0] + wi * xc["ib"][0]
        V2[0] += wb * xc["beta"][1] + wi * wi * xc["ib"][1]
    else:
        print("ФИЗИКА (#GS-42): GS_EXTRA=0 — β, IB и рентген Ar выключены (прежнее поведение)")

    a2, sd2, _ = cbm2.fit_A2V(W2, cnt, r["bg_scaled"], k_bg, V2, sel)

    # #GS-26: lines[] в JSON нужен вкладке K-40 (таблица линий метода 2, как у тория th232_m2:195) — здесь нет
    # цепочки распада (K-40 не звено, а конечная точка), поэтому BR тривиальный и predicted_net = weight_per_branch*a2
    lines_out = [dict(ln, predicted_net=float(ln.get("weight_per_branch", 0.0) * a2[0])) for ln in lines_m2]

    model = W2.T @ a2
    net = cnt - r["bg_scaled"]
    var_full = np.maximum(cnt + k_bg * r["bg_scaled"], 1.0) + V2.T @ a2 ** 2
    chi2 = float(np.sum((model[sel] - net[sel]) ** 2 / var_full[sel]))
    ndof = int(sel.sum()) - 1

    A = float(a2[0] / live)
    dA = float(sd2[0] / live)
    A_lnhb = A * LIB[0][1] / IG_LNHB

    FWHM = lambda E: float(fwhm(E))
    m1.SHAPE_PEAKS = [1460.822]
    shp, _ = m1.shape_residual(r["e"], net, var_full, model, FWHM)
    shape_val = shp[1460.822]

    print(f"ФОРМА ПИКОВ (#SHAPE-1, χ²/ν в окне ±1,5 ПШПВ): 1460.8 → {shape_val:.2f}")
    print(f"K-40 (метод 2, Iγ {LIB[0][1]:.2f} % {fk.K40_DB}): A {A:.1f} ± {dA:.1f} (стат) Бк; ожидается по массе {fk.A_EXP:.1f} Бк; отношение {A/fk.A_EXP:.4f}; χ²/ν {chi2/ndof:.3f}")
    print(f"Чувствительность к ядерным данным: с Iγ 10,34 % (LNHB-DDEP 2025) A {A_lnhb:.1f} Бк; отношение {A_lnhb/fk.A_EXP:.4f}")

    extra = None
    if xc:   # #GS-42: предсказанные отсчёты β и IB (амплитуда K-40 метода 2) в окне и полосах <150, 150–400, >400 кэВ
        # *_col: поканальные столбцы β/IB (та же величина, что фолдится в band_counts ниже) — для слоёв BETA/IB
        # на странице (GS-42, п.2); БЕЗ округления, чтобы Σ(col[sel]) == band_counts["window"].
        beta_col = xc["beta"][0] * a2[0]
        ib_col = xc["ib"][0] * a2[0]
        extra = {"amplitude": "K40, метод 2", "bands_keV": "lt150: E<150; 150_400: 150≤E<400; gt400: E≥400 (вся шкала)",
                 "K40": {"br": 1.0, "beta": gx.band_counts(beta_col, r["e"], sel),
                         "ib": gx.band_counts(ib_col, r["e"], sel),
                         "beta_col": [float(x) for x in beta_col], "ib_col": [float(x) for x in ib_col]}}
        m_win = float(model[sel].sum())
        print("β/IB K-40 (#GS-42, метод 2): в окне β %.0f, IB %.0f отсчётов (%.2f %% и %.2f %% модели в окне)"
              % (extra["K40"]["beta"]["window"], extra["K40"]["ib"]["window"],
                 100 * extra["K40"]["beta"]["window"] / m_win, 100 * extra["K40"]["ib"]["window"] / m_win))
    out_json = os.path.join(m1.OUT, "fit_kcl_m2_bgw%s.json" % fk.SUF)
    data = dict({"extra_components": extra} if extra else {}, **{
        "A_Bq": A,
        "dA_stat_Bq": dA,
        "A_expected_Bq": fk.A_EXP,
        "ratio": A / fk.A_EXP,
        "chi2": chi2,
        "ndof": ndof,
        "shape_1460": shape_val,
        "Igamma_pct": LIB[0][1],
        "Igamma_src": LIB[0][3],
        "A_lnhb_Bq": A_lnhb,
        "Igamma_lnhb_pct": IG_LNHB,
        "n_nodes": n_nodes,
        "n_channels_fit": int(sel.sum()),
        "lines": lines_out,
        "live_s": live,
        "k_bg": k_bg,
        "lo": fk.LO,
        "hi": fk.HI,
        "e": r["e"].tolist(),
        "net": net.tolist(),
        "model": model.tolist(),
        "bg": np.asarray(r["bg_scaled"]).tolist(),
        "var": var_full.tolist()   # дисперсия критерия A2 (с дисперсией узла сетки), для пересчёта #SA-10
    })
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)

if __name__ == "__main__":
    main()
