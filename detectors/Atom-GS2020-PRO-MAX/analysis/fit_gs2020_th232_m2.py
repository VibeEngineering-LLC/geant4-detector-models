# -*- coding: utf-8 -*-
r"""Метод 2 (линии библиотеки × прямые моно-γ прогоны сетки grid_mar_E*.csv, сумм-пики с F_B) для КИ Th-232 в
Маринелли на GS2020. Импорт донора Gamma-1S/web-th232 (grid_response, run_method2, fit_A2V, fit_E1).
#XR-1 (оператор 25.09, заменяет прежнее #M2-1): рентген K/L учитывается ВСЕГДА — библиотека `th232_gs2020_lib05.yaml`
несёт K/L-линии всех звеньев (`make_th232_lib05.add_xray_lines`), узлы сетки под них — jobs_xray.txt (42 энергии).
Спека: scripts\specs\SPEC-fit_gs2020_th232_m2.md"""

import sys, os, json, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))   # detectors/Atom-GS2020-PRO-MAX/analysis → корень репо
DONOR = os.path.join(REPO_ROOT, "detectors", "Gamma-1S", "web-th232")
# GS_M2_CONFIG — своя библиотека (th232_gs2020_lib05.yaml, порог 0,5 %); по умолчанию конфиг донора (порог 2 %)
FULL_XRAY_CFG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "configs", "th232_gs2020_full_xray.yaml")
os.environ["G4MODELS_SOURCE_CONFIG"] = os.environ.get("GS_M2_CONFIG", FULL_XRAY_CFG)   # #XR-1: 2%-библиотека + рентген, не голый донор
if not os.environ.get("SPECTRAVIBE_ROOT"):
    raise RuntimeError("Переменная окружения SPECTRAVIBE_ROOT не установлена")
sys.path.insert(0, HERE)
import fit_gs2020_th232_m1 as m1          # Spec, true_energy, write_fwhm_csv, XML_SAMPLE, XML_BG, OUT, LO, HI, PASSPORT_BQ, PASSPORT_UNC, CHAIN, bm, muc
sys.path.insert(0, DONOR)
import export_amticseu_data as eam        # grid_response; reads module global BUILD_OUT at call time
import export_data as ed                  # ed.GAMMA_LIBRARY: list of (E_keV, I_pct, nuclide_key, note); ed.SUM_PEAKS
import export_ra226_data as erd           # erd.run_method2
import mix_unfold_g1s as g1s              # g1s.make_fwhm(csv_path) -> fwhm(E) function
import crit_bench as crb                  # crb.fit_E1
import crit_bench_m2 as cbm2              # cbm2.fit_A2V
eam.BUILD_OUT = m1.OUT                    # our grid runs grid_mar_E*.csv live here
ed.E_FIT_HI = m1.HI                       # #SUM-1: run_method2 отбрасывает суммы выше E_FIT_HI донора (его конфиг), у нас окно до m1.HI

KEYS = ["Ac228", "Ra224", "Pb212", "Bi212", "Tl208"]
BR = {"Ac228":1.0,"Ra224":1.0,"Pb212":1.0,"Bi212":1.0,"Tl208":0.3594}

def main():
    s = m1.Spec(m1.bm.read(m1.XML_SAMPLE)[0])
    b = m1.Spec(m1.bm.read(m1.XML_BG)[0], m1.BG_TAG)   # #CAL-2: своя шкала фона по его реперам
    if m1.BG_T:  # тот же режим ослабления фона сосудом, что в методе 1 (GS_BG_T="f,rhot")
        b.counts = [c * m1.bg_transmission(b.channel_to_energy(i)) for i, c in enumerate(b.counts)]
        print("ФОН ОСЛАБЛЕН сосудом: f=%g, ρt=%g г/см²" % tuple(m1.BG_T))

    fwhm_csv = os.path.join(m1.OUT, "fwhm_points_gs2020.csv")
    m1.write_fwhm_csv(fwhm_csv)
    fwhm = g1s.make_fwhm(fwhm_csv)

    r = m1.muc.unfold(s, b, [("Tl208", os.path.join(m1.OUT, "mix_Tl208_npsmoff.csv"))], fwhm_csv, lo=m1.LO, hi=m1.HI, bg_energy_of_ch=np.array([b.channel_to_energy(i) for i in range(b.n_channels)]), verbose=False)
    live = s.live_time
    k_bg = s.live_time / b.live_time

    file_e = lambda c: float(s.channel_to_energy(int(c)))
    resp, n_nodes = eam.grid_response(r["ch_edges"], fwhm, e=r["e"], file_e=file_e)
    print(f"Узлов сетки: {n_nodes}")

    lib_m2 = [(E, I, nk, note) for (E, I, nk, note) in ed.GAMMA_LIBRARY if nk in KEYS]
    if len(lib_m2) == 0:
        raise SystemExit("ОТКАЗ: библиотека пуста")
    
    counts_per_key = {k: 0 for k in KEYS}
    for _, _, nk, _ in lib_m2:
        counts_per_key[nk] += 1
    print("Строк библиотеки по ключам:")
    for k in KEYS:
        print(f"  {k}: {counts_per_key[k]}")

    def n_of(E):
        Ek = min(resp.n_map, key=lambda x: abs(x - E))
        if abs(Ek - E) > 0.01:
            raise SystemExit("ОТКАЗ: нет узла сетки %.3f кэВ" % E)
        return resp.n_map[Ek]

    var2 = {}
    _, by_nuc_w, lines_m2, n_sum = erd.run_method2(lib_m2, ed.SUM_PEAKS, resp, r["e"], r["ch_edges"], KEYS, var_acc=var2, n_of=n_of)
    print(f"Суммарных пиков: {n_sum}")

    # #XR-1 (25.09, «отдельный видимый слой»): второй вызов run_method2 ТОЛЬКО по строкам библиотеки
    # с меткой рентгена — аддитивная разность даёт честный вклад K/L-рентгена на графике, не меняя W2/подгонку.
    lib_m2_xray = [(E, I, nk, note) for (E, I, nk, note) in lib_m2 if "#XR-1" in note]
    by_nuc_w_x = {k: np.zeros_like(r["e"]) for k in KEYS}
    if lib_m2_xray:
        _, by_nuc_w_x, _, _ = erd.run_method2(lib_m2_xray, [], resp, r["e"], r["ch_edges"], KEYS, var_acc={}, n_of=n_of)
    print(f"Строк-рентгена в библиотеке (#XR-1): {len(lib_m2_xray)}")

    W2 = np.array([by_nuc_w[k] for k in KEYS])
    W2x = np.array([by_nuc_w_x[k] for k in KEYS])
    V2 = np.array([var2[k] for k in KEYS])
    cnt = np.asarray(s.counts, dtype=float)
    sel = r["sel"]
    # Диагностика (оператор 25.09 «да»): GS_PEAKWIN=k — подгонка только в окнах ±k·ПШПВ вокруг линий библиотеки,
    # континуум между пиками исключён. Проверяет гипотезу «М2 завышает, добирая континуум, которого нет в модели».
    PW = float(os.environ.get("GS_PEAKWIN", "0"))
    if PW > 0:
        w = np.zeros(len(r["e"]), dtype=bool)
        for E, _, _, _ in lib_m2:
            w |= np.abs(r["e"] - E) <= PW * fwhm(E)
        sel = sel & w
        print("ОКНА ПИКОВ: ±%.2f·ПШПВ, каналов в подгонке %d" % (PW, int(sel.sum())))

    a2, sd2, _ = cbm2.fit_A2V(W2, cnt, r["bg_scaled"], k_bg, V2, sel)
    a2e, sd2e, ex2e = crb.fit_E1(W2, cnt, r["bg_scaled"], k_bg, np.ones(len(KEYS)), sel)

    model = W2.T @ a2
    net = cnt - r["bg_scaled"]
    var_crit = np.maximum(cnt[sel] + k_bg * r["bg_scaled"][sel], 1.0) + V2[:, sel].T @ a2 ** 2
    chi2 = float(np.sum((model[sel] - net[sel]) ** 2 / var_crit))
    ndof = int(sel.sum()) - len(KEYS)
    birge = math.sqrt(max(chi2 / ndof, 1.0))

    A = a2 / live
    dA_stat = sd2 / live
    dA = dA_stat * birge
    E1_Bq = a2e / live
    E1_sd_Bq = sd2e / live
    
    chain_eq = {}
    for i, k in enumerate(KEYS):
        chain_eq[k] = A[i] / BR[k]

    print("\nРезультаты Метода 2:")
    print(f"{'Ключ':<8} {'A2 Bq':>10} {'±dA_stat':>10} {'±dA (Birge)':>12} {'E1 Bq':>10} {'Chain-eq Bq':>12} {'Ratio/Pass':>10}")
    print("-" * 85)
    
    valid_chain = []
    for i, k in enumerate(KEYS):
        ratio = chain_eq[k] / m1.PASSPORT_BQ if m1.PASSPORT_BQ > 0 else float('nan')
        print(f"{k:<8} {A[i]:>10.2f} {dA_stat[i]:>10.2f} {dA[i]:>12.2f} {E1_Bq[i]:>10.2f} {chain_eq[k]:>12.2f} {ratio:>10.4f}")
        
        if not math.isinf(dA[i]) and dA[i] > 0:
            rel_err = dA[i] / A[i] if A[i] != 0 else float('inf')
            if rel_err < 0.2:
                valid_chain.append((chain_eq[k], dA[i] / BR[k]))

    if valid_chain:
        weights = [1.0 / (da/BR[k])**2 for k, (_, da) in zip(KEYS, valid_chain)] # Note: BR is constant per key, but we need to map back correctly. 
        # Actually, the prompt says "weights 1/(dA/BR)^2". Since dA is uncertainty of A, and chain_eq = A/BR, uncertainty of chain_eq is dA/BR.
        # So weight is 1 / (unc_chain)**2.
        
        vals = [c for c, _ in valid_chain]
        unc = [u for _, u in valid_chain]
        w = [1.0/u**2 for u in unc]
        w_sum = sum(w)
        weighted_mean = sum(v * wi for v, wi in zip(vals, w)) / w_sum
        weighted_var = 1.0 / w_sum
        weighted_sd = math.sqrt(weighted_var)
        
        print(f"\nВзвешенное среднее эквивалента цепи: {weighted_mean:.2f} ± {weighted_sd:.2f} Bq")
        ratio_pass = weighted_mean / m1.PASSPORT_BQ if m1.PASSPORT_BQ > 0 else float('nan')
        print(f"Отношение к паспорту: {ratio_pass:.4f}")
    else:
        weighted_mean = float('nan')
        weighted_sd = float('nan')

    print(f"\nChi2: {chi2:.2f}, ndof: {ndof}, Chi2/ndof: {chi2/ndof:.2f}, Birge: {birge:.2f}")

    print("\nТРЕБУЕТ ТОЛКОВАНИЯ:")
    issues = []
    for i, k in enumerate(KEYS):
        if A[i] == 0:
            issues.append(f"  - Активность {k} равна нулю.")
    
    if chi2 / ndof > 2:
        issues.append(f"  - Chi2/ndof ({chi2/ndof:.2f}) превышает 2.0.")
        
    if not math.isnan(weighted_mean) and m1.PASSPORT_BQ > 0:
        combined_sigma = math.sqrt(weighted_sd**2 + (0.06 * m1.PASSPORT_BQ)**2)
        diff = abs(weighted_mean - m1.PASSPORT_BQ)
        if diff > 2 * combined_sigma:
            issues.append(f"  - Взвешенное среднее отличается от паспорта более чем на 2 сигмы.")
            
    issues.append(f"  - метод 2: рентген K/L учтён в библиотеке (#XR-1), окно от {m1.LO:.0f} кэВ")
    
    for issue in issues:
        print(issue)

    model_by_key = {}
    for i, k in enumerate(KEYS):
        model_by_key[k] = [round(float(x), 3) for x in W2[i] * a2[i]]

    # ГЛАВНЫЙ результат (оператор 25.09: «цепочка в равновесии»): активность звена = BR·A цепочки,
    # столбец цепочки W = Σ BR·W_k, дисперсия узлов V = Σ BR²·V_k. Поузловая подгонка выше — диагностика.
    br = np.array([BR[k] for k in KEYS])
    Wc, Vc = (br[:, None] * W2).sum(axis=0)[None, :], (br[:, None] ** 2 * V2).sum(axis=0)[None, :]
    ac, sdc, _ = cbm2.fit_A2V(Wc, cnt, r["bg_scaled"], k_bg, Vc, sel)
    ace, sdce, _ = crb.fit_E1(Wc, cnt, r["bg_scaled"], k_bg, np.ones(1), sel)
    mc = (Wc.T @ ac)
    vc = np.maximum(cnt[sel] + k_bg * r["bg_scaled"][sel], 1.0) + Vc[:, sel].T @ ac ** 2
    chi2c = float(np.sum((mc[sel] - net[sel]) ** 2 / vc)); ndofc = int(sel.sum()) - 1
    birgec = math.sqrt(max(chi2c / ndofc, 1.0))
    Ac, dAc = float(ac[0] / live), float(sdc[0] / live)
    print("\nЦЕПОЧКА В РАВНОВЕСИИ (метод 2): A2 %.1f ± %.1f (стат) ± %.1f (Бирге %.2f) Бк; E1 %.1f Бк; отношение к паспорту %.4f; χ²/ν %.3f"
          % (Ac, dAc, dAc * birgec, birgec, float(ace[0] / live), Ac / m1.PASSPORT_BQ, chi2c / ndofc))
    xray_layer = (br[:, None] * W2x).sum(axis=0) * ac[0]   # честный слой К-рентгена (#XR-1), из stack[k] ниже вычтен
    chain_fit = {"A_Bq": Ac, "dA_stat_Bq": dAc, "dA_Bq": dAc * birgec, "birge": birgec, "E1_Bq": float(ace[0] / live),
                 "chi2": chi2c, "ndof": ndofc, "model": [round(float(x), 3) for x in mc],
                 "stack": {k: [round(float(x), 3) for x in (W2[i] - W2x[i]) * BR[k] * ac[0]] for i, k in enumerate(KEYS)},
                 "xray_stack": [round(float(x), 3) for x in xray_layer],
                 "lines": [dict(ln, predicted_net=ln.get("weight_per_branch", 0.0) * BR.get(ln.get("nuclide"), 1.0) * ac[0])
                           for ln in lines_m2]}

    out_data = {
        "keys": KEYS,
        "A_Bq": [float(x) for x in A],
        "dA_stat_Bq": [float(x) for x in dA_stat],
        "dA_Bq": [float(x) for x in dA],
        "E1_Bq": [float(x) for x in E1_Bq],
        "E1_sd_Bq": [float(x) for x in E1_sd_Bq],
        "chain_eq": {k: float(v) for k, v in chain_eq.items()},
        "weighted_chain_Bq": float(weighted_mean) if not math.isnan(weighted_mean) else None,
        "weighted_chain_sd": float(weighted_sd) if not math.isnan(weighted_sd) else None,
        "passport_Bq": m1.PASSPORT_BQ,
        "chi2": chi2,
        "ndof": ndof,
        "birge": birge,
        "n_nodes": n_nodes,
        "n_sum": n_sum,
        "e": [round(float(x), 3) for x in r["e"]],
        "net": [float(x) for x in net],
        "model": [float(x) for x in model],
        "model_by_key": model_by_key,
        "sel": [bool(x) for x in sel],
        "chain": chain_fit
    }

    with open(os.path.join(m1.OUT, "fit_m2%s%s%s.json" % ("_lib05" if "lib05" in os.environ.get("GS_M2_CONFIG", "") else "", "_pw%g" % PW if PW > 0 else "", "_bgT%g_%g" % tuple(m1.BG_T) if m1.BG_T else ("_bgw" if m1.BG_WATER else ""))), "w", encoding="utf-8") as f:
        json.dump(out_data, f, ensure_ascii=False, indent=1)

if __name__ == "__main__":
    main()
