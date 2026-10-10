# -*- coding: utf-8 -*-
r"""#GS-74 (оператор 10.10 «а метод 2?», «только метод 1 и 2»): метод 2 для черники — линии библиотеки × прямые моно-γ прогоны
grid_mar_E*.csv в геометрии черники (run_gs74_berry_m2_runs.sh, CFG1-2026-10-10-berry-m2.md) + β/e⁻ Cs-137 и K-40 (beta_merged, GS2020_BETA_ONLY)
с амплитудой своего нуклида + β Sr-90+Y-90 (шаблон метода 1 mix_SrY90_npsmoff.csv, мешающий параметр) + фон воды × r(E).
Запуск: GS_OUT=C:\g4work\gs2020\berry GS_CAL_SHAPE=1 GS_BG_WATER=1 GS_BG_R=<csv> python fit_gs2020_berry_m2.py"""
import sys, os, json, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
DONOR = r"<REPOS>\geant4-detector-models\detectors\Gamma-1S\web-th232"
os.environ.setdefault("G4MODELS_SOURCE_CONFIG", os.path.join(HERE, "configs", "th232_gs2020_full_xray.yaml"))  # нужен донору при импорте
os.environ.setdefault("SPECTRAVIBE_ROOT", r"<DOSIM>\ИИ\1 Скилы\0_Work\gamma-spectrum-analysis")
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

MASS_KG = float(os.environ.get("GS_BERRY_MASS", "482")) / 1000
KEYS = ["Cs137", "K40"]
CS_SRC = "LNHB-DDEP 2025: Iγ 661,657 кэВ 85,10 %"
XR_SRC = "IAEA decay_rads rad_types=x (библиотека донора amticseu.yaml:256-258); Kβ 1,350 % разделена по относительным выходам ЛСРМ (xray_lines_lsrm.csv: Kb3 4,47, Kb1 8,63, Kb2 2,73), #XR-1"
LIB = [(661.657, 85.10, "Cs137", CS_SRC),
       (31.817, 1.990, "Cs137", XR_SRC), (32.194, 3.667, "Cs137", XR_SRC),
       (36.304, 1.350 * 4.47 / 15.83, "Cs137", XR_SRC), (36.378, 1.350 * 8.63 / 15.83, "Cs137", XR_SRC), (37.255, 1.350 * 2.73 / 15.83, "Cs137", XR_SRC),
       (1460.822, fk.IG_K40, "K40", "LNHB-DDEP 2025: P(γ+ce) 10,34 %")]
LIB_AR_XR = [(0.265, 0.02363, "K40", "ENSDF, IAEA LiveChart API (J. Chen 2015): L-рентген Ar (ЭЗ K-40), #XR-1"),
             (2.956, 0.29564, "K40", "ENSDF, IAEA LiveChart API (J. Chen 2015): KA2 Ar (ЭЗ K-40), #XR-1"),
             (2.958, 0.58555, "K40", "ENSDF, IAEA LiveChart API (J. Chen 2015): KA1 Ar (ЭЗ K-40), #XR-1"),
             (3.19, 0.09508, "K40", "ENSDF, IAEA LiveChart API (J. Chen 2015): KB Ar (ЭЗ K-40), #XR-1")]

def main():
    s = m1.Spec(m1.bm.read(m1.CAL.PATHS["berry"])[0], "berry")
    b = m1.Spec(m1.bm.read(m1.CAL.BKG_WATER_XML)[0], "bgw"); fk.apply_bg_r(b)

    fwhm_csv = os.path.join(m1.OUT, "fwhm_points_gs2020.csv")
    m1.write_fwhm_csv(fwhm_csv, strict=True)
    fwhm = g1s.make_fwhm(fwhm_csv)

    tpl_sr = os.path.join(m1.OUT, "mix_SrY90_npsmoff.csv")
    if not os.path.exists(tpl_sr):
        raise SystemExit("ОТКАЗ: нет шаблона " + tpl_sr)

    bg_e = np.array([b.channel_to_energy(i) for i in range(b.n_channels)])
    r = m1.muc.unfold(s, b, [("SrY90", tpl_sr)], fwhm_csv, lo=fk.LO, hi=fk.HI, bg_energy_of_ch=bg_e, verbose=False, blur=m1.BLUR, tail=m1.TAIL, tvar=m1.TVAR)

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

    lib = LIB + LIB_AR_XR
    print(f"Строк-рентгена в библиотеке (#XR-1): {sum('#XR-1' in ln[3] for ln in lib)}")

    var2 = {}
    _, by_nuc_w, lines_m2, n_sum = erd.run_method2(lib, [], resp, r["e"], r["ch_edges"], KEYS, var_acc=var2, n_of=n_of, var_of=gx.node_var_of(eam.BUILD_OUT, r["ch_edges"], fwhm))

    bfold = {}
    for nuc in KEYS:
        comp = {"beta": gx._load_beta(nuc, m1.OUT, gx._ref(m1.OUT, nuc)), "ib": None}
        bfold[nuc] = gx.fold(comp, g1s.broaden, r["ch_edges"], lambda E: m1.BLUR * fwhm(E))["beta"]

    W2 = np.array([by_nuc_w["Cs137"] + bfold["Cs137"][0], by_nuc_w["K40"] + bfold["K40"][0], np.asarray(r["cols"][0], dtype=float)])
    V2 = np.array([var2["Cs137"] + bfold["Cs137"][1], var2["K40"] + bfold["K40"][1], np.zeros(len(r["e"]))])
    print("β SrY90: дисперсия шаблона не учитывается (10^8 распадов)")

    cnt = np.asarray(s.counts, dtype=float)
    sel = r["sel"]
    a2, sd2, _ = cbm2.fit_A2V(W2, cnt, r["bg_scaled"], k_bg, V2, sel)

    model = W2.T @ a2
    net = cnt - r["bg_scaled"]
    var_full = np.maximum(cnt + k_bg * r["bg_scaled"], 1.0) + V2.T @ a2 ** 2
    chi2 = float(np.sum((model[sel] - net[sel]) ** 2 / var_full[sel]))
    ndof = int(sel.sum()) - 3

    A_Cs137 = float(a2[0] / live)
    dA_Cs137 = float(sd2[0] / live)
    A_K40 = float(a2[1] / live)
    dA_K40 = float(sd2[1] / live)
    A_SrY90 = float(a2[2] / live)
    dA_SrY90 = float(sd2[2] / live)
    print("SrY90: мешающий параметр, не публикуется")

    FWHM = lambda E: float(fwhm(E))
    m1.SHAPE_PEAKS = [661.657, 1460.822]
    shp = m1.shape_residual(r["e"], net, var_full, model, FWHM)[0]
    shp3 = m1.shape_residual(r["e"], net, var_full, model, FWHM, 3.0)[0]
    print("ФОРМА ПИКОВ (#SHAPE-1, χ²/ν ±1,5 / ±3 ПШПВ): " + "; ".join(f"{p:.1f} → {shp[p]:.2f} / {shp3[p]:.2f}" for p in m1.SHAPE_PEAKS))

    print(f"Cs137 (метод 2): A {A_Cs137:.2f} ± {dA_Cs137:.2f} (стат) Бк; {A_Cs137/MASS_KG:.1f} ± {dA_Cs137/MASS_KG:.1f} Бк/кг")
    print(f"K40 (метод 2): A {A_K40:.2f} ± {dA_K40:.2f} (стат) Бк; {A_K40/MASS_KG:.1f} ± {dA_K40/MASS_KG:.1f} Бк/кг")
    print(f"χ²/ν {chi2/ndof:.3f} (ν {ndof}), окно {fk.LO:g}–{fk.HI:g} кэВ, узлов {n_nodes}")

    lines_out = []
    for ln in lines_m2:
        d = dict(ln)
        nuc = d.get("nuclide", d.get("nuc"))
        if nuc is not None and nuc in KEYS:
            d["predicted_net"] = float(d.get("weight_per_branch", 0.0) * a2[KEYS.index(nuc)])
        else:
            d["predicted_net"] = 0.0
        lines_out.append(d)

    cols = [list(map(float, W2[i] * a2[i])) for i in range(3)]
    beta_cols = {n: list(map(float, bfold[n][0] * a2[KEYS.index(n)])) for n in KEYS}

    out_json = os.path.join(m1.OUT, "fit_berry_m2_bgw.json")
    data = {
        "nuclides": {
            "Cs137": {"A_Bq": A_Cs137, "dA_stat_Bq": dA_Cs137},
            "K40": {"A_Bq": A_K40, "dA_stat_Bq": dA_K40},
            "SrY90": {"A_Bq": A_SrY90, "dA_stat_Bq": dA_SrY90}
        },
        "mass_kg": MASS_KG,
        "chi2": chi2,
        "ndof": ndof,
        "shape": {str(p): [shp[p], shp3[p]] for p in m1.SHAPE_PEAKS},
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
        "var": var_full.tolist(),
        "cols": cols,
        "col_names": ["Cs137", "K40", "SrY90"],
        "beta_cols": beta_cols,
        "library": [list(x) for x in lib]
    }
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)

if __name__ == "__main__":
    main()
