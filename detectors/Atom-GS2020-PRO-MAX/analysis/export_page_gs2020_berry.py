# -*- coding: utf-8 -*-
r"""#GS-72: выгрузка данных вкладки «Черника» страницы GS2020 в КОНТРАКТЕ донорской страницы (как export_page_gs2020_k40.py);
на странице объект называется window.GS_BERRY (build_page_gs2020_berry.py).
Источник чисел: подгонка метода 1 fit_gs2020_berry.py — fit_berry_bgw_R2_g014_beta.json (принятая, 3 шаблона: Cs-137, K-40,
Sr-90+Y-90) и её вариант _beta_ib.json (+ слой внутреннего тормозного KUB; амплитуда 0, модель побитово та же — проверяется).
Метод 2 (#GS-74): fit_berry_m2_bgw.json, блок строит export_berry_m2.build(); слои — линии Cs-137/K-40, тормозное β/e⁻ Cs-137/K-40, β Sr-90+Y-90, фон.
Строки сравнения (площадь пика × ε, стерильные пересчёты) — из export_sources_berry.berry_source() (§33, не копия).
Спека: scripts\specs\SPEC-export_page_gs2020_berry.md
Запуск: GS_OUT=C:\g4work\gs2020\berry GS_CAL_SHAPE=1 PYTHONIOENCODING=utf-8 python export_page_gs2020_berry.py"""
import os, sys, json, math, datetime
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
os.environ.setdefault("GS_OUT", r"C:\g4work\gs2020\berry")
os.environ.setdefault("GS_CAL_SHAPE", "1")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import export_page_gs2020_k40 as ek      # tpl_header, gdml_material, xml_times, r4, rnum, load, SL_CENTROID, N_EFF_MIN, PAGE (§33: берём готовое)
import export_sources_berry as esb       # berry_source(): строки сравнения черники
import fit_gs2020_th232_m1 as m1
import fit_gs2020_kcl as fk
import export_berry_m2 as eb2        # #GS-74: блок метода 2 (сгенерирован по SPEC-export_berry_m2.md)

BERRY = r"C:\g4work\gs2020\berry"
PAGE = ek.PAGE
DST = os.path.join(PAGE, "gs2020_berry_data.json")
F_M1 = os.path.join(BERRY, "fit_berry_bgw_R2_g014_beta.json")
F_IB = os.path.join(BERRY, "fit_berry_bgw_R2_g014_beta_ib.json")
F_M2 = os.path.join(BERRY, "fit_berry_m2_bgw.json")   # #GS-74: метод 2 (fit_gs2020_berry_m2.py)
F_FWHM = os.path.join(BERRY, "fwhm_points_gs2020.csv")
N_EFF_MIN = ek.N_EFF_MIN
MASS_KG = 0.482
CS_SCALE = 85.131 / 85.10      # fit_gs2020_berry.py:6 — множитель Iγ 662, тот же, что в подгонке
BG_COLOR = "#b8b2a2"
SL_CS_PER_KG_2024, SL_CS_UNC_PCT, SL_K_PER_KG, SL_K_UNC_PCT = 4190.0, 3.0, 264.0, 9.0   # results\ref_beta1s\README.md
SL_YEARS = 2.018       # от 01.10.2024 до 09.10.2026 (README.md, шапка таблицы)
T_HALF_CS = 30.08      # лет, LNHB-DDEP
FACTS = {"edge_662_pct": 2.5, "edge_1461_pct": 1.75, "resid_pct": 14.0, "resid_lo_keV": 450, "resid_hi_keV": 550,
         "cs_sys_pct": 4.0, "k_noise_Bq": 14.0, "k_sys_Bq": 30.0, "variants_cs_pct": 2.0, "gain_w": 0.14, "gain_s_pct": 2.8,
         "level_mm": -16.5, "k_frac_model_pct": 0.6, "k_frac_found_pct": 0.9}   # results\RESULT-2026-10-10-berry-gs70.md, «Оговорки»

NUC = [
    ("CS137", "Cs137", "mix_Cs137_npsmoff.csv", CS_SCALE, "Cs-137", "#b8860b", "МК-шаблон полного распада Cs-137 (до Ba-137m) в объёме пробы; K-рентгеновское излучение Ba — внутри шаблона"),
    ("K40", "K40", "mix_K40_npsmoff.csv", fk.K40_SCALE, "K-40", "#2e7d32", "МК-шаблон полного распада K-40 в объёме пробы; предварительно: малая разность пика пробы и фонового калия, зависит от r(E)"),
    ("SRY90", "SrY90", "mix_SrY90_npsmoff.csv", 1.0, "внешнее тормозное излучение β Sr-90+Y-90 (Geant4)", "#7b1fa2", "β-электроны Sr-90 и Y-90 (1:1) тормозятся в пробе, стенке сосуда и корпусе; мешающий параметр, активность Sr-90 не публикуется"),
]   # #GS-76 (оператор 10.10 «а зачем оно? нигде не учитывали»): слой внутреннего тормозного KUB снят (амплитуда 0, у K-40 выключен #GS-45)

def refuse(msg):
    raise SystemExit("ОТКАЗ: " + msg)

def crop_to(n_show):
    def crop(obj, n_full):
        if isinstance(obj, list):
            return obj[:n_show]
        if isinstance(obj, dict):
            return {k: crop(v, n_full) for k, v in obj.items()}
        return obj
    return crop

def fwhm_cal():
    conv = {}
    with open(F_FWHM, encoding="utf-8") as f:
        lines = f.readlines()
    for line in lines[1:]:
        parts = line.strip().split(",")
        if len(parts) != 2:
            continue
        try:
            e = float(parts[0])
            fw = float(parts[1])
            conv[e] = fw
        except ValueError:
            continue
    if len(conv) < 3:
        refuse("fwhm_cal: меньше 3 точек в " + F_FWHM)
    points = []
    used_E = []
    used_F = []
    scales = []
    for E in sorted(m1.FWHM_SL):
        sl = m1.FWHM_SL[E]
        key = None
        for k in conv:
            if abs(k - E) < 1e-3:
                key = k
                break
        used = key is not None
        fw = conv[key] if used else sl
        sc = fw / sl if used else None
        if used:
            used_E.append(E)
            used_F.append(fw)
            scales.append(sc)
        Ec = ek.SL_CENTROID.get(E, E)
        pt = {
            "E_nominal": E,
            "E_centroid": Ec,
            "fwhm_keV": fw,
            "d_fwhm_keV": 0.0,
            "res_pct": 100 * fw / Ec,
            "shift_keV": Ec - E,
            "fwhm_sl_keV": sl,
            "own_keV": None,
            "own_unc_keV": None,
            "becqmoni_keV": None,
            "own_used": False,
            "own_reason": "",
            "scale": sc,
            "used": used,
            "n_lines_window": 1,
            "fwhm_model_keV": None,
            "dev_pct": None,
            "reject": ""
        }
        points.append(pt)
    c = np.polyfit(np.log(used_E), np.log(used_F), 1)
    p = float(c[0])
    k = math.exp(c[1])
    rms = 100 * math.sqrt(np.mean([(k * E**p / f - 1)**2 for E, f in zip(used_E, used_F)]))
    for pt in points:
        E = pt["E_nominal"]
        pt["fwhm_model_keV"] = k * E**p
        if pt["used"]:
            pt["dev_pct"] = 100 * (pt["fwhm_keV"] / pt["fwhm_model_keV"] - 1)
        else:
            pt["reject"] = "точки нет в наборе данных детектора (референсный спектр тория), в свёртку не входит: СпектраЛайн " + ek.rnum(pt["fwhm_sl_keV"], 2) + " кэВ"
    return {
        "source": "набор данных детектора #DET-1: " + os.path.basename(F_FWHM) + " (референсный спектр тория) подгонки черники",
        "k": k,
        "p": p,
        "rms_dev_pct": rms,
        "n_used": len(used_E),
        "n_anchors": len(points),
        "fwhm662_law": k * 661.657**p,
        "fwhm662_cs": k * 661.657**p,
        "res662_pct": 100 * k * 661.657**p / 661.657,
        "scale": np.mean(scales),
        "scale_uniform": (max(scales) - min(scales)) <= 1e-6,
        "k40_sl_keV": m1.FWHM_SL.get(1460.822, 0.0),
        "k40_conv_keV": float(m1.muc.g1s.make_fwhm(F_FWHM)(1460.822)),   # интерполяция в логарифмах по точкам свёртки (как в export_page_gs2020_k40)
        "k40_law_keV": k * 1460.822**p,
        "points": points
    }

def main():
    jm = ek.load(F_M1)
    jb = ek.load(F_IB)
    if list(jm["nuclides"]) != ["Cs137", "K40", "SrY90"]:
        refuse("порядок/состав nuclides в " + F_M1 + " изменился: " + str(list(jm["nuclides"])))
    if list(jb["nuclides"]) != ["Cs137", "K40", "SrY90", "IBSrY90"]:
        refuse("порядок/состав nuclides в " + F_IB + " изменился: " + str(list(jb["nuclides"])))
    for name in jm["nuclides"]:
        a1 = jm["nuclides"][name]["A_Bq"]
        a2 = jb["nuclides"][name]["A_Bq"]
        if abs(a1 - a2) / max(abs(a1), 1e-30) > 1e-9:
            refuse("A_Bq " + name + " различается между " + F_M1 + " и " + F_IB)
    if jb["nuclides"]["IBSrY90"]["A_Bq"] != 0.0:
        refuse("амплитуда внутреннего тормозного не нулевая — слой нельзя выдавать как «≈ 0»")
    if np.max(np.abs(np.asarray(jm["model"]) - np.asarray(jb["model"]))) > 1e-6:
        refuse("модель различается между " + F_M1 + " и " + F_IB)
    if abs(jm["chi2"] - jb["chi2"]) > 1e-6 or jm["live_s"] != jb["live_s"]:
        refuse("chi2/live_s различаются между " + F_M1 + " и " + F_IB)

    sp = m1.bm.read(m1.CAL.PATHS["berry"])[0]
    s = m1.Spec(sp, "berry")
    bsp = m1.bm.read(m1.CAL.BKG_WATER_XML)[0]
    b = m1.Spec(bsp, "bgw")
    e_of_ch = np.array([s.channel_to_energy(i) for i in range(s.n_channels)])
    if not np.allclose(e_of_ch, jm["e"], atol=1e-6):
        refuse("энергии каналов не совпадают с fit JSON")
    live = s.live_time
    if abs(live - jm["live_s"]) > 1e-6:
        refuse("live_time не совпадает с fit JSON")
    counts = np.asarray(s.counts, float)
    bg = np.asarray(jm["bg"], float)
    if not np.allclose(counts, np.asarray(jm["net"]) + bg, atol=1e-6):
        refuse("counts != net + bg")
    cols = np.asarray(jb["cols"], float)
    if cols.shape[0] != 4:
        refuse("cols.shape[0] != 4")
    if np.max(np.abs(cols.sum(axis=0) - np.asarray(jb["model"]))) > 1e-6:
        refuse("sum(cols) != model")
    lo, hi = fk.LO, fk.HI
    sel = [(lo <= x <= hi) for x in jm["e"]]

    stack = {}
    trusted = {}
    noise = {}
    n_events = {}
    for i, (key, jn, tpl, k, lab, col, note) in enumerate(NUC):
        stack[key] = ek.r4(cols[i])
        if tpl is None:
            trusted[key] = [True] * len(cols[i])
            noise[key] = 0.0
            continue
        hdr = ek.tpl_header(os.path.join(BERRY, tpl))
        n_events[key] = int(hdr["n_events_processed"])
        act = jb["nuclides"][jn]["A_Bq"] / k
        scale = n_events[key] / (act * live)
        trusted[key] = [bool(v * scale >= N_EFF_MIN) for v in cols[i]]
        total = sum(cols[i])
        if total > 0:
            noise[key] = sum(cols[i][j] for j in range(len(cols[i])) if not trusted[key][j]) / total
        else:
            noise[key] = 0.0
    stack["BG"] = ek.r4(bg)
    trusted["BG"] = [True] * len(bg)
    noise["BG"] = 0.0
    sel_idx = [i for i, v in enumerate(sel) if v]
    sum_stack = sum(sum(stack[key][i] for i in sel_idx) for key in stack)
    sum_model_bg = sum(np.asarray(jm["model"])[i] for i in sel_idx) + sum(bg[i] for i in sel_idx)
    if abs(sum_stack - sum_model_bg) / max(abs(sum_model_bg), 1e-30) > 1e-3:
        refuse("сумма слоёв в окне подгонки не сходится с model+bg")
    model = np.asarray(jm["model"], float)
    sb = {
        "model_counts": ek.r4(model),
        "model2_counts": ek.r4(model),
        "model2_full_counts": ek.r4(model),
        "stack": stack,
        "trusted": trusted,
        "n_eff_min": N_EFF_MIN,
        "noise_frac": noise,
        "stack2": stack,
        "stack2_full": stack,
        "stack2_chan": {},
        "stack2_chan_full": {}
    }

    chi2_ndof = jm["chi2"] / jm["ndof"]
    birge = math.sqrt(max(chi2_ndof, 1.0))
    dec = 2.0 ** (-SL_YEARS / T_HALF_CS)
    sl_cs_kg = SL_CS_PER_KG_2024 * dec
    if abs(sl_cs_kg - 3999.6) > 0.15:
        refuse("sl_cs_kg не совпадает с README")
    per = {}
    for key, jn, tpl, k, lab, col, note in NUC:
        v = jb["nuclides"][jn]
        A = v["A_Bq"]
        dst = v["dA_stat_Bq"]
        ent = {
            "A_Bq": A,
            "dA_Bq": dst * birge if dst == dst else None,
            "share": None,
            "nuisance": (jn == "SrY90"),
            "ref_ratio": None,
            "per_kg": A / MASS_KG
        }
        if jn == "SrY90":
            ent["A_Bq"] = None
            ent["dA_Bq"] = None
            ent["per_kg"] = None
        if jn == "Cs137":
            ent["ref_ratio"] = A / MASS_KG / sl_cs_kg
        if jn == "K40":
            ent["ref_ratio"] = A / MASS_KG / SL_K_PER_KG
        per[key] = ent
    shp = jm["shape"]
    M1 = {
        "A_Bq": per["CS137"]["A_Bq"],
        "dA_Bq": per["CS137"]["dA_Bq"],
        "dA_stat_Bq": jm["nuclides"]["Cs137"]["dA_stat_Bq"],
        "birge": birge,
        "E1_Bq": None,
        "bg_amplitude": 1.0,
        "d_bg_amplitude": 0.0,
        "chi2": jm["chi2"],
        "ndof": jm["ndof"],
        "chi2_ndof": chi2_ndof,
        "E_fit_lo": lo,
        "E_fit_hi": hi,
        "sys_floor": 0.0,
        "xray_total_per_branch_pct": 0.0,
        "ratio_to_passport": per["CS137"]["ref_ratio"],
        "d_ratio": 0.0,
        "shape_662": shp["661.657"][0],
        "shape_1460": shp["1460.822"][0],
        "per_nuclide": per
    }
    jm2 = ek.load(F_M2)     # #GS-74: метод 2 — тот же фон, живое время, окно
    if not np.allclose(np.asarray(jm2["bg"], float), bg, atol=1e-6) or abs(jm2["live_s"] - live) > 1e-6 or (jm2["lo"], jm2["hi"]) != (lo, hi):
        refuse("метод 2: фон/живое время/окно не совпадают с методом 1")
    b2 = eb2.build(jm2, bg, len(e_of_ch), MASS_KG, sl_cs_kg, SL_K_PER_KG)
    M2 = b2["method2"]
    sb["stack2"] = b2["stack2"]
    sb["stack2_full"] = b2["stack2"]
    sb["model2_counts"] = ek.r4(jm2["model"])
    sb["model2_full_counts"] = ek.r4(jm2["model"])

    hdr = ek.tpl_header(os.path.join(BERRY, "mix_Cs137_npsmoff.csv"))
    gdml = hdr["_title"].split("GDML ")[1].strip()
    mat = hdr["sample_matrix"]
    rho = float(hdr["sample_rho_g_cm3"])
    dens, comp = ek.gdml_material(gdml, mat)
    if abs(dens - rho) > 1e-3:
        refuse("плотность GDML не совпадает с sample_rho")
    t0, t1 = ek.xml_times(m1.CAL.PATHS["berry"])
    cal = m1.CAL_OWN
    if "berry" not in cal or "bgw" not in cal:
        refuse("cal не содержит berry/bgw")
    k_bg = s.live_time / b.live_time
    if abs(k_bg / 0.6075 - 1) > 1e-3:
        refuse("k_bg не совпадает с 0.6075")

    rows = esb.berry_source()["results"]
    pick = {r["lab"]: r for r in rows}
    cs_a = sl_cs_kg * MASS_KG
    cs_da = cs_a * SL_CS_UNC_PCT / 100
    k_a = SL_K_PER_KG * MASS_KG
    k_da = k_a * SL_K_UNC_PCT / 100
    cmp = {
        "mass_kg": MASS_KG,
        "cs": [
            # оператор 10.10 («удельная другая была» → «сделай пометку что пересчитано на дату»): в подписи — дата пересчёта и исходное значение
            {"lab": "Бета-1С (СпектраЛайн): " + ek.rnum(SL_CS_PER_KG_2024, 0) + " Бк/кг на дату его измерения 01.10.2024; пересчитано на 09.10.2026", "A": cs_a, "dA": cs_da, "kind": "ref", "note": "здесь — с распадом Cs-137 за " + ek.rnum(SL_YEARS, 3) + " года, на дату нашего измерения, и на массу этой пробы; другая навеска и другой прибор; параметры модели не подбирались"},
            {"lab": "метод 1", "A": per["CS137"]["A_Bq"], "dA": per["CS137"]["dA_Bq"], "kind": "m1", "note": "подгонка полного спектра шаблонами Geant4"},
            b2["cmp"]["cs"],
        ],   # #GS-74 (оператор 10.10 «это мы не используем. только метод 1 и 2»): площадь×ε и стерильный пересчёт — только в отчёте
        "k": [
            {"lab": "Бета-1С (СпектраЛайн): " + ek.rnum(SL_K_PER_KG, 0) + " Бк/кг на дату его измерения 01.10.2024", "A": k_a, "dA": k_da, "kind": "ref", "note": "пересчитано на массу этой пробы; у K-40 распад за 2 года пренебрежимо мал; другая навеска и другой прибор"},
            {"lab": "метод 1", "A": per["K40"]["A_Bq"], "dA": per["K40"]["dA_Bq"], "kind": "m1", "note": "предварительно: целиком зависит от r(E) фона"},
            b2["cmp"]["k"]
        ]
    }
    beta1s = {
        "cs_2024_per_kg": SL_CS_PER_KG_2024,
        "cs_unc_pct": SL_CS_UNC_PCT,
        "cs_now_per_kg": sl_cs_kg,
        "k_per_kg": SL_K_PER_KG,
        "k_unc_pct": SL_K_UNC_PCT,
        "years": SL_YEARS,
        "t_half_cs": T_HALF_CS,
        "date_ref": "01.10.2024",
        "date_now": "09.10.2026",
        "sr90_per_kg_2024": 100.0,
        "img": "img/spectraline_beta1s_berry_2024-10-01.jpg"
    }
    passport = {
        "A_Bq": cs_a,
        "dA_Bq": cs_da,
        "Bq_per_kg": sl_cs_kg,
        "unc_pct": SL_CS_UNC_PCT,
        "mass_g": MASS_KG * 1000,
        "vol_ml": 1000.0,
        "date_certified": "опорное: Бета-1С, 01.10.2024; пересчёт на 09.10.2026",
        "date_measured": t0[:10] + " — " + t1[:10],
        "decay_factor": dec,
        "unc_components_pct": {},
        "center_note": "опорное значение другого прибора; аттестации у пробы нет"
    }

    meta = {
        "live_s": live,
        "real_s": float(sp.real),
        "bg_live_s": b.live_time,
        "bg_real_s": float(bsp.real),
        "bg_scale_time": k_bg,
        "cal_sample": {
            "coefs": cal["berry"]["coeffs"],
            "order": len(cal["berry"]["coeffs"]) - 1,
            "n_channels": s.n_channels,
            "repr_max_dev_keV": cal["berry"].get("repr_max_dev_keV")
        },
        "cal_bg": {
            "coefs": cal["bgw"]["coeffs"],
            "order": len(cal["bgw"]["coeffs"]) - 1,
            "n_channels": b.n_channels,
            "repr_max_dev_keV": cal["bgw"].get("repr_max_dev_keV")
        },
        "cal_file": os.path.basename(m1.CAL.OUT_JSON),
        "sys_floor_pct": 0.0,
        "xray_span_lo_keV": 0.0,
        "xray_span_hi_keV": 0.0,
        "template_decays": [
            {"nuclide": "Cs-137", "n": n_events["CS137"]},
            {"nuclide": "K-40", "n": n_events["K40"]},
            {"nuclide": "Sr-90+Y-90", "n": n_events["SRY90"]}
        ],
        "nuclide_list_ru": "Cs-137, K-40, Sr-90+Y-90",
        "matrix_name": "сушёная черника",
        "gdml_material": mat,
        "matrix_density_g_cm3": rho,
        "matrix_composition": comp,
        "template_source": "Geant4, стенд gs2020_marinelli, " + os.path.basename(BERRY),
        "date_start": t0,
        "date_end": t1
    }

    nuclides = []
    for key, jn, tpl, k, lab, col, note in NUC:
        nuclides.append({
            "key": key,
            "label_ru": lab,
            "label_en": lab,
            "color": col,
            "note": note,
            "branching": 1.0
        })
    nuclides.append({
        "key": "BG",
        "label_ru": "фон (приведён)",
        "label_en": "background",
        "color": BG_COLOR,
        "branching": 1.0,
        "note": "фон — сосуд Маринелли 1 л с дистиллированной водой, своя шкала энергии, × отношение живых времён × r(E) черники (отдельно для линий фона и для континуума, Geant4, #GS-60); не подгоняется"
    })
    contrib = {key: sum(stack[key]) for key in stack}
    nuclides.sort(key=lambda x: (x["key"] in ("SRY90", "IBSRY90"), -contrib[x["key"]]))
    nuclides.extend(b2["nuclides_extra"])     # слои метода 2 — после отсортированных слоёв метода 1
    per_total = sum(contrib.values())
    for key in per:
        per[key]["share"] = contrib[key] / per_total

    n_full = len(e_of_ch)
    n_show = int(np.searchsorted(e_of_ch, 3000.0, side="right"))
    crop = crop_to(n_show)
    data = {
        "meta": meta,
        "fwhm_cal": fwhm_cal(),
        "passport": passport,
        "nuclides": nuclides,
        "channels": [],
        "spectrum": crop(dict({"e_of_ch": ek.r4(e_of_ch), "counts": [int(round(c)) for c in counts], "bg_counts": ek.r4(bg)}, **sb), n_full),
        "cs": {"method1": M1, "method2": M2, "method2_full": M2, "spectrum": crop(sb, n_full)},
        "method1": M1,
        "method2": M2,
        "method2_full": M2,
        "library": {"i_threshold_pct": 0.0, "fixed_n": 0, "full_n": 0, "full_threshold_pct": 0.0},
        "reference_lines": [[661.657, "Cs-137", "Cs-137 661,7"], [1460.822, "K-40", "K-40 1460,8"], [2614.511, "Tl-208 (фон)", "Tl-208 2614,5"]],
        "detector_lines": [],
        "cmp": cmp,
        "beta1s": beta1s,
        "facts": FACTS,
        "provenance": {"fit": os.path.basename(F_M1), "fit_ib": os.path.basename(F_IB), "fit_m2": os.path.basename(F_M2), "fit_log": "fit_R2_g014_beta.log", "built": datetime.date.today().isoformat()}
    }
    with open(DST, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, separators=(",", ":"), allow_nan=False)

    size_kb = os.path.getsize(DST) / 1024
    print("Сохранено: " + DST + " (" + ek.rnum(size_kb, 1) + " КБ)")
    print("Шкала: " + m1.CAL.OUT_JSON + ", каналов " + str(n_show) + " из " + str(n_full) + ", окно подгонки " + ek.rnum(lo, 1) + "–" + ek.rnum(hi, 1) + " кэВ")
    for key in ("CS137", "K40"):
        e = per[key]
        print(key + ": A=" + ek.rnum(e["A_Bq"], 1) + " Бк, dA=" + ek.rnum(e["dA_Bq"], 1) + " Бк, " + ek.rnum(e["per_kg"], 1) + " Бк/кг, ref_ratio=" + ek.rnum(e["ref_ratio"], 3) + ", share=" + ek.rnum(e["share"] * 100, 1) + " %")
    print("χ²/ν=" + ek.rnum(chi2_ndof, 3) + ", ν=" + str(jm["ndof"]) + ", Birge=" + ek.rnum(birge, 3) + ", shape_662=" + ek.rnum(shp["661.657"][0], 4) + ", shape_1460=" + ek.rnum(shp["1460.822"][0], 4))
    print("Доля шума: " + ", ".join(key + "=" + ek.rnum(noise[key], 4) for key in noise))
    print("Бета-1С Cs-137 на 09.10.2026: " + ek.rnum(sl_cs_kg, 1) + " Бк/кг (распад ×" + ek.rnum(dec, 4) + ")")
    print("GDML материал: " + mat + " (" + ", ".join("%s %.1f %%" % (c["element"], 100 * c["mass_fraction"]) for c in comp) + ")")
    print("ТРЕБУЕТ ТОЛКОВАНИЯ:")
    print("  амплитуда внутреннего тормозного KUB = 0.0 (модель побитово совпадает с принятой), срез нулевой и помечен на странице")
    print("  метод 2: Cs-137 " + ek.rnum(M2["A_Bq"], 1) + " Бк, K-40 " + ek.rnum(M2["per_nuclide"]["K40"]["A_Bq"], 1) + " Бк, χ²/ν " + ek.rnum(M2["chi2_ndof"], 3) + "; расхождение с методом 1 по Cs-137 " + ek.rnum(100 * (M2["A_Bq"] / per["CS137"]["A_Bq"] - 1), 2) + " %")
    print("  K-40 предварительно: стерильный пересчёт 158 Бк против метода 1 (" + ek.rnum(per["K40"]["A_Bq"], 0) + " Бк)")

if __name__ == "__main__":
    main()
