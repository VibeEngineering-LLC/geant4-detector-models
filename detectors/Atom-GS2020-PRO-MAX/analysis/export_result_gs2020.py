# -*- coding: utf-8 -*-
r"""Выгрузка итогов разложения КИ Th-232 (GS2020, Маринелли) для страницы web/gs2020-th232: result.js.
Источник — JSON подгонок out_v5 (fit_gs2020_th232_m1.py / _m2.py). Числа на странице — только отсюда (#PUB-1).
Спека: scripts\specs\SPEC-export_result_gs2020.md. Запуск: python export_result_gs2020.py"""

import sys
import os
import json
import datetime
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

OUT = r"C:\g4work\gs2020\run_marinelli\out_v5"
DST = r"D:\cloud-folder\work-folder\GEANT4\web\gs2020-th232\result.js"
FITS = [
    ("m1", "fit_m1.json", "Метод 1", "шаблоны полного распада звеньев (RDM)"),
    ("m2", "fit_m2_lib05.json", "Метод 2", "моно-γ сетка × библиотека ENSDF ≥ 0,5 % (42 линии)"),
    ("m2_lib2", "fit_m2.json", "Метод 2, библиотека 2 %", "моно-γ сетка × библиотека донора ≥ 2 % (20 линий)"),
    ("m1_bgT06", "fit_m1_bgT0.6_6.02.json", "Метод 1, фон ослаблен сосудом (f 0,6)", "фон × T(E), T(1461) = 0,828 по K-40; 60 % фона через сосуд"),
    ("m1_bgT1", "fit_m1_bgT1_3.365.json", "Метод 1, фон ослаблен сосудом (f 1)", "фон × T(E), T(1461) = 0,828 по K-40; весь фон через сосуд")
]
LO, HI = 150.0, 2900.0
RBIN = 25.0

def load(name):
    with open(os.path.join(OUT, name), encoding="utf-8") as f:
        return json.load(f)

def main():
    fits = {tag: load(fname) for tag, fname, _, _ in FITS}
    
    passport = fits["m1"]["passport_Bq"]
    for tag, fname, _, _ in FITS[1:]:
        if abs(fits[tag]["passport_Bq"] - passport) > 1e-6:
            raise SystemExit("ОТКАЗ: паспорт в JSON разный")

    e_m1 = np.array(fits["m1"]["e"])
    net_m1 = np.array(fits["m1"]["net"])
    
    for tag, fname, _, _ in FITS[1:]:
        e_cur = np.array(fits[tag]["e"])
        net_cur = np.array(fits[tag]["net"])
        if not np.allclose(e_m1, e_cur, atol=1e-6):
            raise SystemExit("ОТКАЗ: разные e/net в JSON подгонок")
        if not tag.startswith("m1_bgT") and not np.allclose(net_m1, net_cur, rtol=0, atol=1e-3):  # М1 пишет net с округлением до 3 знаков; у bgT фон ослаблен — net другой по построению
            raise SystemExit("ОТКАЗ: разные e/net в JSON подгонок")

    rows = []
    for tag, fname, title, what in FITS:
        c = fits[tag]["chain"]
        row = {
            "tag": tag,
            "title": title,
            "what": what,
            "A": c["A_Bq"],
            "dA": c["dA_Bq"],
            "dA_stat": c["dA_stat_Bq"],
            "birge": c["birge"],
            "E1": c["E1_Bq"],
            "ratio": c["A_Bq"] / passport,
            "chi2_ndof": c["chi2"] / c["ndof"],
            "ndof": c["ndof"]
        }
        rows.append(row)

    links = []
    # M1 links
    m1_data = fits["m1"]
    names = m1_data["names"]
    activities_A2 = m1_data["activities_A2"]
    sd_A2_Bq = m1_data["sd_A2_Bq"]
    chain_eq_m1 = m1_data["chain_eq"]
    
    for i, n in enumerate(names):
        a = activities_A2[i]
        if a > 0:
            sd_val = sd_A2_Bq[i]
            if np.isnan(sd_val):
                sd_out = None
            else:
                sd_out = sd_val
            links.append({
                "method": "m1",
                "link": n,
                "A": a,
                "dA": sd_out,
                "chain_eq": chain_eq_m1[i]
            })

    # M2 (lib05) links
    m2_data = fits["m2"]
    keys = m2_data["keys"]
    A_Bq_m2 = m2_data["A_Bq"]
    dA_Bq_m2 = m2_data["dA_Bq"]
    chain_eq_m2 = m2_data["chain_eq"]
    
    for i, k in enumerate(keys):
        links.append({
            "method": "m2",
            "link": k,
            "A": A_Bq_m2[i],
            "dA": dA_Bq_m2[i],
            "chain_eq": chain_eq_m2[k]
        })

    # Spectrum arrays
    e = np.array(fits["m1"]["e"])
    net = np.array(fits["m1"]["net"])
    de = np.gradient(e)
    
    mask = (e >= LO) & (e <= HI)
    
    e_filtered = np.round(e[mask], 2).tolist()
    net_filtered = np.round(net[mask] / de[mask], 2).tolist()
    
    model_m1 = np.array(fits["m1"]["chain"]["model"])
    m1_filtered = np.round(model_m1[mask] / de[mask], 2).tolist()
    
    model_m2 = np.array(fits["m2"]["chain"]["model"])
    m2_filtered = np.round(model_m2[mask] / de[mask], 2).tolist()

    # Ratio net/model in bins
    bin_centers = []
    ratio_m1 = []
    sigma_m1 = []
    ratio_m2 = []
    sigma_m2 = []
    
    b = LO
    while b < HI:
        b_end = b + RBIN
        bin_mask = (e >= b) & (e < b_end)
        
        s_net = np.sum(net[bin_mask])
        s_mod_m1 = np.sum(model_m1[bin_mask])
        s_mod_m2 = np.sum(model_m2[bin_mask])
        
        center = b + RBIN / 2.0
        
        if s_mod_m1 > 0:
            r1 = s_net / s_mod_m1
            sig1 = r1 / np.sqrt(s_mod_m1)
        else:
            r1 = None
            sig1 = None
            
        if s_mod_m2 > 0:
            r2 = s_net / s_mod_m2
            sig2 = r2 / np.sqrt(s_mod_m2)
        else:
            r2 = None
            sig2 = None
            
        bin_centers.append(round(center, 4))
        ratio_m1.append(round(r1, 4) if r1 is not None else None)
        sigma_m1.append(round(sig1, 4) if sig1 is not None else None)
        ratio_m2.append(round(r2, 4) if r2 is not None else None)
        sigma_m2.append(round(sig2, 4) if sig2 is not None else None)
        
        b = b_end

    data = {
        "generated": datetime.datetime.now().isoformat(timespec='seconds'),
        "passport_Bq": passport,
        "passport_unc_pct": 6.0,
        "lo": LO,
        "hi": HI,
        "rows": rows,
        "links": links,
        "fit": {
            "E_keV": e_filtered,
            "net": net_filtered,
            "m1": m1_filtered,
            "m2": m2_filtered
        },
        "ratio": {
            "E_keV": bin_centers,
            "m1": ratio_m1,
            "m1_sd": sigma_m1,
            "m2": ratio_m2,
            "m2_sd": sigma_m2
        }
    }

    with open(DST, "w", encoding="utf-8") as f:
        f.write("window.GS2020_RESULT = " + json.dumps(data, ensure_ascii=False) + ";\n")

    print(f"Файл записан: {DST}")
    
    for row in rows:
        print(f"{row['title']}: A2 {row['A']:.1f} ± {row['dA']:.1f} Бк, E1 {row['E1']:.1f}, отношение {row['ratio']:.4f}, χ²/ν {row['chi2_ndof']:.3f}")
        
    print(f"Точек спектра: {len(e_filtered)}")
    print(f"Бинов отношения: {len(bin_centers)}")

    print("\nТРЕБУЕТ ТОЛКОВАНИЯ:")
    
    issues_found = False
    
    # Rows with |ratio - 1| > 0.06
    for row in rows:
        if abs(row['ratio'] - 1) > 0.06:
            print(f"  {row['title']}: отношение {row['ratio']:.4f} (вне паспорта ±6%)")
            issues_found = True
            
    # Rows with chi2_ndof > 2
    for row in rows:
        if row['chi2_ndof'] > 2:
            print(f"  {row['title']}: χ²/ν = {row['chi2_ndof']:.3f} (> 2)")
            issues_found = True
            
    # Ratio bins where m1 ratio deviates from 1 by more than 0.10
    count_lines = 0
    for i, center in enumerate(bin_centers):
        r = ratio_m1[i]
        if r is not None and abs(r - 1) > 0.10:
            print(f"  Бин {center:.1f} кэВ: отношение m1 = {r:.4f}")
            count_lines += 1
            issues_found = True
            if count_lines >= 15:
                break
                
    if not issues_found:
        print("  нет")

if __name__ == "__main__":
    main()
