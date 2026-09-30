# -*- coding: utf-8 -*-
r"""#GS-26: выгрузка данных вкладки K-40 (KCl) страницы GS2020 в КОНТРАКТЕ донорской страницы Гамма-1С Th-232
(window.G1S, эталон — export_page_gs2020.py); на странице объект называется window.GS_K40 (build_page_gs2020_k40.py).
У K-40 нет цепочки (один нуклид), один вариант фона (Маринелли+вода) и одна линия — поэтому cs.* и method2_full
ДУБЛИРУЮТ method1/method2. Источник чисел — JSON подгонок KCl (fit_gs2020_kcl.py, fit_gs2020_kcl_m2.py) в GS_OUT.
Спека: scripts\specs\SPEC-export_page_gs2020_k40.md.
Запуск: GS_OUT=C:\g4work\gs2020\kcl GS_CAL_SHAPE=1 PYTHONIOENCODING=utf-8 python export_page_gs2020_k40.py"""

import sys
import os
import re
import json
import math
import datetime
import xml.etree.ElementTree as ET
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fit_gs2020_th232_m1 as m1      # bm, Spec, CAL, CAL_OWN, OUT, FWHM_SL, muc (muc.g1s.make_fwhm — ядро свёртки)
import fit_gs2020_kcl as fk           # MASS_G, VOL_ML, K_FRAC, K40_BQ_PER_G_K, A_EXP

OUT = m1.OUT
PAGE = r"D:\cloud-folder\work-folder\GEANT4\web\gs2020-th232-page"
DST = os.path.join(PAGE, "gs2020_k40_data.json")

F_M1 = os.path.join(OUT, "fit_kcl_bgw.json")
F_M2 = os.path.join(OUT, "fit_kcl_m2_bgw.json")
F_TPL = os.path.join(OUT, "mix_K40_npsmoff.csv")
F_GRID = os.path.join(OUT, "grid_mar_E1460.822.csv")
F_FWHM = os.path.join(OUT, "fwhm_points_gs2020.csv")

E_K40 = 1460.822
N_EFF_MIN = 4.0
COLOR_K40 = "#1b8a8f"
COLOR_BG = "#b8b2a2"

# audit/GS-26-k40-expected-unc.md
U_THETA = 0.000008 / 0.011668   # #GS-45 (30.09): база LNHB-DDEP 2025 (K-40_tables.pdf с.1): 0,011668 (8) %
U_THALF = 0.0027 / 1.2522       # T½ 1,2522 (27)·10⁹ лет
PUR_LO = 0.998 * (1.0 - 0.008)
PUR_HI = 1.0
U_PUR = (PUR_HI - PUR_LO) / math.sqrt(12.0) / ((PUR_HI + PUR_LO) / 2.0)
U_EXP = math.sqrt(U_THETA ** 2 + U_THALF ** 2 + U_PUR ** 2)

SL_CENTROID = {238.632: 238.735, 583.187: 583.747, 911.204: 910.859, 964.766: 964.420, 968.971: 968.625, 1460.822: 1460.866, 2614.511: 2614.506}
# центроиды таблицы «Параметры пиков» СпектраЛайн (спектр Th-232), README-референсы.md, строки 134-141

REF_NUC = {238.632: "Pb-212 (фон)", 609.312: "Bi-214 (фон)", 1460.822: "K-40", 2614.511: "Tl-208 (фон)"}
# реперы шкалы KCl — gs2020_calib.py REFS["kcl"]

def refuse(msg):
    raise SystemExit("ОТКАЗ: " + msg)

def load(p):
    if not os.path.exists(p):
        refuse("нет файла " + p)
    return json.load(open(p, encoding="utf-8"))

def r4(a):
    return [round(float(x), 4) for x in a]

def window_sum(col, sel):
    return float(sum(v for v, s in zip(col, sel) if s))

def check_layer_sum(name, model, stack, sel, tol=1e-3):
    """GS-42 п.2: Σ ВСЕХ слоёв стека (K40+BETA+IB+BG) в окне sel = model + слой BG (аналог проверки
    в export_page_gs2020.py, здесь её раньше не было)."""
    m = window_sum(model, sel) + window_sum(stack["BG"], sel)
    tot = sum(window_sum(col, sel) for col in stack.values())
    if m > 0 and abs(m - tot) / m > tol:
        refuse("Σ слоёв страницы (BETA/IB/BG) ≠ модель+фон в %s (%.3f против %.3f)" % (name, m, tot))

def rnum(x, d):
    return ("%.*f" % (d, x)).replace(".", ",")

def tpl_header(path):
    hdr = {}
    with open(path, encoding="utf-8") as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("bin_keV"):
            break
        if "," in line and not line.startswith("#"):
            key, val = line.split(",", 1)
            hdr[key.strip()] = val.strip()
        elif line.startswith("#"):
            hdr["_title"] = line.strip()
    return hdr

def gdml_material(gdml_path, name):
    with open(gdml_path, encoding="utf-8", errors="replace") as f:
        text = f.read()
    pat = r'<material name="' + re.escape(name) + r'"[^>]*>(.*?)</material>'
    m = re.search(pat, text, re.S)
    if not m:
        refuse("в GDML нет материала " + name)
    body = m.group(1)
    dm = re.search(r'<D value="([0-9.eE+-]+)"', body)
    if not dm:
        refuse("нет плотности в материале " + name)
    density = float(dm.group(1))
    fracs = []
    for fm in re.finditer(r'<fraction n="([0-9.eE+-]+)" ref="([^"]+)"', body):
        n = float(fm.group(1))
        ref = fm.group(2)
        el = ref.replace("G4_", "") if ref.startswith("G4_") else ref
        fracs.append({"element": el, "mass_fraction": n})
    return density, fracs

def xml_times(path):
    root = ET.parse(path).getroot()
    st = root.find(".//StartTime")
    en = root.find(".//EndTime")
    return (st.text if st is not None else None, en.text if en is not None else None)

def detector_lines():
    out = []
    E = E_K40
    nk = "K-40"
    grid_path = F_GRID
    if not os.path.exists(grid_path):
        refuse("нет файла " + grid_path)
    e, c = [], []
    for l in open(grid_path, encoding="utf-8"):
        p = l.split(",")
        try:
            e.append(float(p[0]))
            c.append(float(p[1]))
        except ValueError:
            pass
    e, c = np.array(e), np.array(c)
    pk = c[np.floor(e) == math.floor(E)].sum()
    if pk <= 0:
        refuse("pk <= 0 в grid_mar_E%s.csv" % E)
    bg = c[((e >= E - 38) & (e < E - 35)) | ((e >= E - 25) & (e < E - 22))].mean()
    exc = lambda a, b: float((c[(e >= E - b) & (e < E - a)] - bg).sum()) / pk
    out.append({"E_keV": E, "nuclide": nk, "E_esc_ka": round(E - 28.5, 1), "frac_ka": round(exc(27.4, 29.6), 5),
                "E_esc_kb": round(E - 32.3, 1), "frac_kb": round(exc(31.2, 33.4), 5)})
    return out

def fwhm_cal():
    if not os.path.exists(F_FWHM):
        refuse("нет файла " + F_FWHM)
    conv = {}
    with open(F_FWHM, encoding="utf-8") as f:
        lines = f.readlines()
    for l in lines:
        l = l.strip()
        if not l or l.startswith("E_keV"):
            continue
        parts = l.split(",")
        try:
            E_val = float(parts[0])
            fwhm_val = float(parts[1])
            conv[E_val] = fwhm_val
        except (ValueError, IndexError):
            pass

    law_conv = m1.muc.g1s.make_fwhm(F_FWHM)

    points = []
    used_points_E = []
    used_points_fwhm = []
    scales = []

    for E in sorted(m1.FWHM_SL.keys()):
        sl = m1.FWHM_SL[E]
        found_key = None
        for k in conv:
            if abs(k - E) < 1e-3:
                found_key = k
                break
        
        if found_key is not None:
            used = True
            fwhm = conv[found_key]
            scale = fwhm / sl
            scales.append(scale)
            used_points_E.append(E)
            used_points_fwhm.append(fwhm)
        else:
            used = False
            fwhm = sl
            scale = None
        
        Ec = SL_CENTROID.get(E, E)
        
        # Placeholder for model values, will be filled after polyfit if needed, 
        # but spec says calculate law first. We need to store points then fill model/dev later?
        # Spec: "points: for each E ... a dict ... fwhm_model_keV: k*E**p"
        # So we must compute k, p first.
        
        points.append({
            "E_nominal": E, 
            "E_centroid": Ec, 
            "fwhm_keV": fwhm, 
            "d_fwhm_keV": 0.0, 
            "res_pct": 100 * fwhm / Ec, 
            "shift_keV": Ec - E,
            "fwhm_sl_keV": sl, 
            "scale": scale, 
            "used": used, 
            "n_lines_window": 1,
            "fwhm_model_keV": None, # To be filled
            "dev_pct": None, # To be filled
            "reject": "" if used else None # To be filled for unused
        })

    if len(scales) < 3:
        refuse("меньше 3 использованных точек для ПШПВ")

    scale_uniform = (max(scales) - min(scales)) <= 1e-6
    fw_scale = sum(scales) / len(scales)

    # Power law fit
    log_E = [math.log(E) for E in used_points_E]
    log_fwhm = [math.log(f) for f in used_points_fwhm]
    coeffs = np.polyfit(log_E, log_fwhm, 1)
    p = coeffs[0]
    k = math.exp(coeffs[1])

    # RMS dev
    rms_dev_pct = 0.0
    if len(used_points_E) > 0:
        err_sq = []
        for i, E in enumerate(used_points_E):
            fwhm_model = k * (E ** p)
            err_sq.append((fwhm_model / used_points_fwhm[i] - 1) ** 2)
        rms_dev_pct = 100 * math.sqrt(sum(err_sq) / len(err_sq))

    # Fill points with model values
    k40_sl = m1.FWHM_SL.get(E_K40, None)
    k40_conv = float(law_conv(E_K40)) if law_conv else 0.0
    k40_law = k * (E_K40 ** p)

    for pt in points:
        E = pt["E_nominal"]
        fwhm_model = k * (E ** p)
        pt["fwhm_model_keV"] = fwhm_model
        if pt["used"]:
            pt["dev_pct"] = 100 * (pt["fwhm_keV"] / fwhm_model - 1)
            pt["reject"] = ""
        else:
            # Unused point logic
            if k40_sl is not None:
                reject_text = "контрольная точка K-40, в свёртку не входит: СпектраЛайн " + rnum(k40_sl, 2) + " кэВ (× " + rnum(fw_scale, 2) + " = " + rnum(k40_sl * fw_scale, 2) + "), свёртка модели " + rnum(k40_conv, 2) + ", закон k·E^p " + rnum(k40_law, 2) + " кэВ"
            else:
                reject_text = "контрольная точка K-40, в свёртку не входит"
            pt["reject"] = reject_text

    fwhm662_law = k * (661.657 ** p)

    return {
        "source": "таблица «Параметры пиков» СпектраЛайн (спектр Th-232 того же прибора, README-референсы.md) × множитель свёртки; точки — " + os.path.basename(F_FWHM) + " подгонки KCl",
        "k": k, 
        "p": p, 
        "rms_dev_pct": rms_dev_pct, 
        "n_used": len(used_points_E), 
        "n_anchors": len(points),
        "fwhm662_law": fwhm662_law, 
        "fwhm662_cs": fwhm662_law, 
        "res662_pct": 100 * fwhm662_law / 661.657,
        "scale": fw_scale, 
        "scale_uniform": scale_uniform, 
        "k40_sl_keV": k40_sl if k40_sl else 0.0, 
        "k40_conv_keV": k40_conv, 
        "k40_law_keV": k40_law,
        "points": points
    }

def method1_block(j, A_exp):
    chi2_ndof = j["chi2"] / j["ndof"]
    birge = math.sqrt(max(chi2_ndof, 1.0))
    dA = j["dA_stat_Bq"] * birge
    A = j["A_Bq"]
    return {
        "A_Bq": A, 
        "dA_Bq": dA, 
        "dA_stat_Bq": j["dA_stat_Bq"], 
        "birge": birge, 
        "E1_Bq": None, 
        "bg_amplitude": 1.0,
        "d_bg_amplitude": 0.0, 
        "chi2": j["chi2"], 
        "ndof": j["ndof"], 
        "chi2_ndof": chi2_ndof, 
        "E_fit_lo": j["lo"], 
        "E_fit_hi": j["hi"],
        "sys_floor": 0.0, 
        "xray_total_per_branch_pct": 0.0, 
        "ratio_to_passport": A / A_exp, 
        "d_ratio": dA / A_exp,
        "shape_1460": j["shape_1460"], 
        "per_nuclide": {"K40": {"A_Bq": A, "dA_Bq": dA, "share": 1.0}}
    }

def method2_block(j, A_exp):
    if "lines" not in j or not j["lines"]:
        refuse("method2: нет линий")
    
    chi2_ndof = j["chi2"] / j["ndof"]
    birge = math.sqrt(max(chi2_ndof, 1.0))
    dA = j["dA_stat_Bq"] * birge
    A = j["A_Bq"]

    lines = []
    for ln in j["lines"]:
        lines.append({
            "E_keV": ln["E_keV"], 
            "nuclide": ln["nuclide"], 
            "I_gamma_pct": ln.get("I_pct"), 
            "branch": 1.0, 
            "eps_peak": ln["eps_peak"],
            "weight_per_branch": ln["weight_per_branch"], 
            "predicted_net": ln["predicted_net"], 
            "kind": ln["kind"], 
            "note": ln.get("note", ""),
            "E1_keV": None, 
            "E2_keV": None, 
            "I1_pct": None, 
            "I2_pct": None
        })

    return {
        "A_Bq": A, 
        "dA_Bq": dA, 
        "dA_stat_Bq": j["dA_stat_Bq"], 
        "birge": birge, 
        "E1_Bq": None, 
        "bg_amplitude": 1.0,
        "chi2": j["chi2"], 
        "ndof": j["ndof"], 
        "chi2_ndof": chi2_ndof, 
        "n_lines": len(lines), 
        "n_channels_fit": j["n_channels_fit"],
        "n_sum_peaks": 0, 
        "n_sum_peaks_total": 0, 
        "n_xray_energies": 0, 
        "n_nodes": j["n_nodes"], 
        "ratio_to_passport": A / A_exp,
        "d_ratio": dA / A_exp, 
        "shape_1460": j["shape_1460"], 
        "Igamma_src": j["Igamma_src"], 
        "Igamma_lnhb_pct": j["Igamma_lnhb_pct"],
        "A_lnhb_Bq": j["A_lnhb_Bq"], 
        "lines": lines
    }

def extra_layers(j, n):
    """GS-42 п.2: слои BETA/IB подгонки K-40 j (JSON плоский, extra_components — верхнего уровня) — поканальные
    столбцы *_col нуклида K40 на сетке model (GS_EXTRA=0 — extra_components отсутствует, слои нулевые)."""
    d = (j.get("extra_components") or {}).get("K40") or {}
    beta = np.asarray(d["beta_col"], dtype=float) if "beta_col" in d else np.zeros(n)
    ib = np.asarray(d["ib_col"], dtype=float) if "ib_col" in d else np.zeros(n)
    return beta, ib

def spectrum_block(jm1, jm2, bg, live, n_events):
    A = jm1["A_Bq"]
    scale = n_events / (A * live)
    mod1 = np.asarray(jm1["model"], float)
    mod2 = np.asarray(jm2["model"], float)
    
    trusted_k = [bool(v * scale >= N_EFF_MIN) for v in mod1]
    
    sum_mod1 = mod1.sum()
    if sum_mod1 > 0:
        noise_sum = sum(v for v, t in zip(mod1, trusted_k) if not t)
        noise_k = float(noise_sum / sum_mod1)
    else:
        noise_k = 0.0

    # GS-42 п.2: слои BETA/IB — М1 несёт только IB (β уже внутри шаблона распада иона, не отделим), М2 — оба;
    # вклад вычитается из слоя K40, чтобы не считать его дважды (тот же приём, что у страницы Th-232).
    beta1, ib1 = extra_layers(jm1, len(mod1))
    beta2, ib2 = extra_layers(jm2, len(mod2))
    k40_1 = mod1 - ib1
    k40_2 = mod2 - beta2 - ib2

    return {
        "model_counts": r4(mod1),
        "model2_counts": r4(mod2),
        "model2_full_counts": r4(mod2),
        "stack": {"K40": r4(k40_1), "BETA": [0.0] * len(bg), "IB": r4(ib1), "BG": r4(bg)},
        "trusted": {"K40": trusted_k, "BETA": [False] * len(bg), "IB": [False] * len(bg), "BG": [True] * len(bg)},
        "n_eff_min": N_EFF_MIN,
        "noise_frac": {"K40": noise_k, "BETA": 0.0, "IB": 0.0, "BG": 0.0},
        "stack2": {"K40": r4(k40_2), "BETA": r4(beta2), "IB": r4(ib2), "BG": r4(bg)},
        "stack2_full": {"K40": r4(k40_2), "BETA": r4(beta2), "IB": r4(ib2), "BG": r4(bg)},
        "stack2_chan": {},
        "stack2_chan_full": {}
    }

def main():
    jm1 = load(F_M1)
    jm2 = load(F_M2)

    if abs(jm1["A_expected_Bq"] - jm2["A_expected_Bq"]) > 1e-6:
        refuse("разные ожидаемые активности в JSON подгонок")
    
    if not np.allclose(jm1["e"], jm2["e"], atol=1e-6):
        refuse("разные оси энергии в JSON подгонок")

    sp = m1.bm.read(m1.CAL.PATHS["kcl"])[0]
    s = m1.Spec(sp, "kcl")
    bsp = m1.bm.read(m1.CAL.BKG_WATER_XML)[0]
    b = m1.Spec(bsp, "bgw")

    if not hasattr(sp, "real"):
        refuse("sp не имеет атрибута real")
    if not hasattr(bsp, "real"):
        refuse("bsp не имеет атрибута real")

    e_of_ch = np.array([s.channel_to_energy(i) for i in range(s.n_channels)])
    
    if len(e_of_ch) != len(jm1["e"]):
        refuse("шкала KCl не совпадает с подгонкой (нужен тот же GS_CAL_SHAPE, что у fit_gs2020_kcl.py)")
    if not np.allclose(e_of_ch, jm1["e"], atol=1e-6):
        refuse("шкала KCl не совпадает с подгонкой (нужен тот же GS_CAL_SHAPE, что у fit_gs2020_kcl.py)")

    print("ШКАЛА KCl: " + m1.CAL.OUT_JSON + ", калибровка по форме: " + str(m1.CAL.SHAPE))

    live = s.live_time
    if abs(live - jm1["live_s"]) > 1e-6:
        refuse("разное живое время")
    
    k_bg = s.live_time / b.live_time
    if abs(k_bg / jm1["k_bg"] - 1) > 1e-9:
        refuse("разный коэффициент фона")

    counts = np.asarray(s.counts, float)
    bg = np.asarray(jm1["bg"], float)

    if not np.allclose(counts, np.asarray(jm1["net"]) + bg, atol=1e-6):
        refuse("net+bg подгонки М1 ≠ счёт файла KCl")
    
    if not np.allclose(jm2["bg"], bg, atol=1e-6):
        refuse("фон М2 не совпадает с фоном М1")

    hdr = tpl_header(F_TPL)
    n_events = int(hdr["n_events_processed"])
    sb = spectrum_block(jm1, jm2, bg, live, n_events)   # GS-42 п.2: считаем один раз, переиспользуем ниже
    sel_win = [(jm1["lo"] <= x <= jm1["hi"]) for x in jm1["e"]]
    check_layer_sum("method1", sb["model_counts"], sb["stack"], sel_win)
    check_layer_sum("method2", sb["model2_counts"], sb["stack2"], sel_win)
    rho = float(hdr["sample_rho_g_cm3"])
    mat = hdr["sample_matrix"]
    
    title_line = hdr.get("_title", "")
    gdml_path = ""
    if "GDML " in title_line:
        gdml_path = title_line.split("GDML ")[1].strip()
    
    if not gdml_path or not os.path.exists(gdml_path):
        refuse("не найден GDML файл: " + gdml_path)

    dens, comp = gdml_material(gdml_path, mat)
    if abs(dens - rho) > 1e-3:
        refuse("плотность в GDML не совпадает с заголовком шаблона")

    A_exp = jm1["A_expected_Bq"]
    if abs(A_exp / fk.A_EXP - 1) > 1e-9:
        refuse("ожидаемая активность не совпадает с константой fit_gs2020_kcl")

    t0, t1 = xml_times(m1.CAL.PATHS["kcl"])

    passport = {
        "A_Bq": A_exp, 
        "dA_Bq": A_exp * U_EXP, 
        "Bq_per_kg": A_exp / (fk.MASS_G / 1000.0), 
        "unc_pct": 100 * U_EXP,
        "mass_g": fk.MASS_G, 
        "vol_ml": fk.VOL_ML, 
        "date_certified": "расчёт по массе и ядерным данным, без аттестации",
        "date_measured": (t0[:10] + " — " + t1[:10]) if (t0 and t1) else "не найдено в файле", 
        "decay_factor": 1.0,
        "unc_components_pct": {"abundance": 100 * U_THETA, "half_life": 100 * U_THALF, "purity": 100 * U_PUR},
        "center_note": "чистота KCl в центральном значении принята полной (audit/GS-26-k40-expected-unc.md: вопрос о поправке 0,995 открыт)"
    }

    cal = m1.CAL_OWN
    if cal is None or "kcl" not in cal or "bgw" not in cal:
        refuse("нет калибровки KCl или BGW в CAL_OWN")

    meta = {
        "live_s": live, 
        "real_s": float(sp.real), 
        "bg_live_s": b.live_time, 
        "bg_real_s": float(bsp.real), 
        "bg_scale_time": k_bg,
        "cal_sample": {"coefs": cal["kcl"]["coeffs"], "order": len(cal["kcl"]["coeffs"]) - 1, "n_channels": s.n_channels,
                       "repr_max_dev_keV": cal["kcl"].get("repr_max_dev_keV")},
        "cal_bg": {"coefs": cal["bgw"]["coeffs"], "order": len(cal["bgw"]["coeffs"]) - 1, "n_channels": b.n_channels,
                   "repr_max_dev_keV": cal["bgw"].get("repr_max_dev_keV")},
        "cal_file": os.path.basename(m1.CAL.OUT_JSON), 
        "sys_floor_pct": 0.0, 
        "xray_span_lo_keV": 0.0, 
        "xray_span_hi_keV": 0.0,
        "template_decays": [{"nuclide": "K-40", "n": n_events}], 
        "nuclide_list_ru": "K-40", 
        "matrix_name": "KCl",
        "gdml_material": mat, 
        "matrix_density_g_cm3": rho, 
        "matrix_composition": comp,
        "template_source": "Geant4, стенд gs2020_marinelli, " + os.path.basename(OUT.rstrip("\\/")),
        "date_start": t0, 
        "date_end": t1
    }

    nuclides = [
        {"key": "K40", "label_ru": "K-40", "label_en": "K-40", "color": COLOR_K40,
         "note": "МК-шаблон полного распада K-40 (Geant4, ионный источник Z = 19, A = 40) в объёме пробы KCl", "branching": 1.0},
        {"key": "BG", "label_ru": "фон (приведён)", "label_en": "background", "color": COLOR_BG,
         "note": "фон — сосуд Маринелли 1 л с дистиллированной водой, своя шкала энергии, × отношение живых времён; не подгоняется", "branching": 1.0},
        # GS-42 п.2: тормозное β и внутреннее тормозное (IB) — отдельные слои (см. export_page_gs2020.py)
        {"key": "BETA", "label_ru": "тормозное β", "label_en": "β bremsstrahlung", "color": "#2b6cb0",
         "note": "в методе 1 отдельно не выделяется (уже внутри шаблона распада иона); в методе 2 — тормозное излучение электронов β-распада (Geant4)", "branching": 1.0},
        {"key": "IB", "label_ru": "внутреннее тормозное (IB)", "label_en": "internal bremsstrahlung (IB)", "color": "#c0392b",
         "note": "фотоны внутреннего тормозного при β-распаде по таблице KUB; Geant4 их не рождает, добавлены отдельно (#GS-42)", "branching": 1.0}
    ]

    contrib = {"K40": sum(sb["stack"]["K40"]), "BG": sum(bg), "BETA": sum(sb["stack2"]["BETA"]), "IB": sum(sb["stack"]["IB"])}
    if contrib["IB"] == 0:   # 01.10 (оператор): IB K-40 убран из расчёта до числового подтверждения — пустой слой в легенде не показываем
        nuclides = [n for n in nuclides if n["key"] != "IB"]
    nuclides.sort(key=lambda n: (n["key"] in ("BETA", "IB"), -contrib[n["key"]]))

    reference_lines = []
    for E, _hw in m1.CAL.REFS["kcl"]:
        if E not in REF_NUC:
            refuse("репер %s не найден в REF_NUC" % E)
        nuc_name = REF_NUC[E]
        short_name = nuc_name.split(" ")[0] + " " + rnum(E, 1)
        reference_lines.append([E, nuc_name, short_name])

    n_full = len(e_of_ch)
    n_show = int(np.searchsorted(e_of_ch, 3000.0, side="right"))

    def crop(obj):
        if isinstance(obj, list) and len(obj) == n_full:
            return obj[:n_show]
        elif isinstance(obj, dict):
            return {k: crop(v) for k, v in obj.items()}
        else:
            return obj

    M1 = method1_block(jm1, A_exp)
    M2 = method2_block(jm2, A_exp)
    fw = fwhm_cal()

    data = {
        "meta": meta, 
        "fwhm_cal": fw, 
        "passport": passport, 
        "nuclides": nuclides, 
        "channels": [],
        "spectrum": crop(dict({"e_of_ch": r4(e_of_ch), "counts": [int(round(c)) for c in counts], "bg_counts": r4(bg)}, **sb)),
        "cs": {"method1": M1, "method2": M2, "method2_full": M2, "spectrum": crop(sb)},
        "method1": M1, 
        "method2": M2, 
        "method2_full": M2,
        "library": {"i_threshold_pct": 0.0, "fixed_n": M2["n_lines"], "full_n": M2["n_lines"], "full_threshold_pct": 0.0},
        "reference_lines": reference_lines, 
        "detector_lines": detector_lines()
    }

    os.makedirs(PAGE, exist_ok=True)
    with open(DST, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    
    print("Файл сохранен: %s (%.1f КБ)" % (DST, os.path.getsize(DST) / 1024))

    print("ожидаемая по массе: %.1f ± %.1f Бк (%.2f %%: распространённость %.2f, T½ %.2f, чистота %.2f)" % (
        A_exp, A_exp * U_EXP, 100 * U_EXP, 100 * U_THETA, 100 * U_THALF, 100 * U_PUR))

    for name, blk in (("method1", M1), ("method2", M2)):
        print("%s: A = %.1f ± %.1f Бк (стат %.1f, Бирге %.3f), к ожидаемой %.4f, χ²/ν %.3f, форма 1460 %.2f" % (
            name, blk["A_Bq"], blk["dA_Bq"], blk["dA_stat_Bq"], blk["birge"], 
            blk["ratio_to_passport"], blk["chi2_ndof"], blk["shape_1460"]))

    print("ПШПВ: точки свёртки %d из %d, множитель %.4f (одинаковый: %s); закон k=%.4f p=%.4f, СКО %.2f %%; K-40: СпектраЛайн %.2f, свёртка %.2f, закон %.2f кэВ" % (
        fw["n_used"], fw["n_anchors"], fw["scale"], fw["scale_uniform"], 
        fw["k"], fw["p"], fw["rms_dev_pct"], 
        fw["k40_sl_keV"], fw["k40_conv_keV"], fw["k40_law_keV"]))

    det_lines = detector_lines()
    if det_lines:
        dl = det_lines[0]
        print("пики вылета 1460,8: Kα %.3f %%, Kβ %.3f %%" % (100 * dl["frac_ka"], 100 * dl["frac_kb"]))

    noise_k = sb["noise_frac"]["K40"]
    print("шаблон K-40: %d распадов, n_eff_min %.0f, доля шума слоя %.4f; каналов на странице %d из %d" % (
        n_events, N_EFF_MIN, noise_k, n_show, n_full))

    issues = []
    for name, blk in (("method1", M1), ("method2", M2)):
        sig = math.sqrt(blk["dA_Bq"]**2 + passport["dA_Bq"]**2)
        if abs(blk["A_Bq"] - A_exp) > 2 * sig:
            issues.append("%s: расхождение с ожидаемой %.1f Бк = %.1f σ (σ = ошибка метода ⊕ неопределённость ожидаемой)" % (
                name, blk["A_Bq"] - A_exp, abs(blk["A_Bq"] - A_exp) / sig))

    if not fw["scale_uniform"]:
        issues.append("множитель ПШПВ свёртки не одинаков по точкам")

    if noise_k > 0.5:
        issues.append("доля шума шаблона K-40 > 50 %")

    issues.append("точки ПШПВ — файл " + os.path.basename(F_FWHM) + ", перезаписывается каждой подгонкой; соответствие именно подгонке М1 — по порядку run_gs2020_kcl.sh (tail → main → m2), не по содержимому JSON")

    print("ТРЕБУЕТ ТОЛКОВАНИЯ:")
    if issues:
        for iss in issues:
            print("  " + iss)
    else:
        print("  нет")

if __name__ == "__main__":
    main()
