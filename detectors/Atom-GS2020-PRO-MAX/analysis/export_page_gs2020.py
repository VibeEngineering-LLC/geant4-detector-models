# -*- coding: utf-8 -*-
r"""Выгрузка данных страницы GS2020 Th-232 в КОНТРАКТЕ донорской страницы Гамма-1С Th-232 (web-th232: g1s-th232.js,
build_page.py) — оформление и интерактив берутся у донора без правок (оператор 25.09: «всё в едином согласованном стиле»).
Источник чисел — JSON подгонок out_v5. Переключатель донора «закон ПШПВ» (lines | cs) здесь означает фон:
lines — фон S31_18 как снят (без сосуда), cs — фон реального измерения «Маринелли 1 л + дист. вода» (оператор 26.09,
GS_BG_WATER=1, промежуточный замер 11,7 ч) — заменил прежнее модельное ослабление сосудом (GS_BG_T), которое
хуже описывало форму (χ²/ν хуже во всех трёх подгонках).
Спека: scripts\specs\SPEC-export_page_gs2020.md. Запуск: python export_page_gs2020.py"""

import sys
import os
import json
import math
import numpy as np
import yaml

sys.stdout.reconfigure(encoding="utf-8")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fit_gs2020_th232_m1 as m1      # bm, Spec, true_energy, XML_SAMPLE, XML_BG, PEAK_TABLE, CHAIN, PASSPORT_BQ, PASSPORT_UNC, LO, HI

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
OUT = m1.OUT
PAGE = os.environ.get("GS2020_PAGE_WORK", os.path.join(REPO_ROOT, ".work", "gs2020-th232-page"))
DST = os.path.join(PAGE, "gs2020_th232_data.json")
DONOR_CFG = os.path.join(REPO_ROOT, "detectors", "Gamma-1S", "web-th232", "configs", "th232.yaml")
LIB2_CFG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "configs", "th232_gs2020_full_xray.yaml")   # #XR-1: донор + рентген, не голый DONOR_CFG
LIB05_CFG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "configs", "th232_gs2020_lib05.yaml")
FILES = {"m1": "fit_m1.json", "m2": "fit_m2.json", "m2f": "fit_m2_lib05.json"}
N_PER_BR = 5.5e7
N_EFF_MIN = 4.0
PASSPORT = {"A_Bq": m1.PASSPORT_BQ, "dA_Bq": m1.PASSPORT_BQ * m1.PASSPORT_UNC, "Bq_per_kg": 910.0, "unc_pct": 6.0, "mass_g": 1052.0, "date_certified": "образец с известной активностью (дата не указана)", "date_measured": "2026-09-25", "decay_factor": 1.0}

import fit_gs2020_th232_m2 as m2    # ed.SUM_PEAKS: (E1, E2, нуклид, I1 %, I2 %, примечание, fb %) — выходы пар сумм-пиков
SUM_I = {(round(t[0], 3), round(t[1], 3)): (t[3], t[4]) for t in m2.ed.SUM_PEAKS}


def load(name):
    return json.load(open(os.path.join(OUT, name), encoding="utf-8"))

def bgw(name):
    return name.replace(".json", "_bgw.json")

def r4(a):
    return [round(float(x), 4) for x in a]

def method1_block(j):
    c = j["chain"]
    A = c["A_Bq"]
    T = live
    links = [k for k, _ in m1.CHAIN]
    BR = dict(m1.CHAIN)
    per_nuclide = {k: {"A_Bq": A * BR[k], "dA_Bq": c["dA_Bq"] * BR[k], "share": float(sum(c["stack"][k])) / float(sum(c["model"]))} for k in links}
    return {"A_Bq": A, "dA_Bq": c["dA_Bq"], "dA_stat_Bq": c["dA_stat_Bq"], "birge": c["birge"], "E1_Bq": c["E1_Bq"],
            "bg_amplitude": 1.0, "d_bg_amplitude": 0.0, "chi2": c["chi2"], "ndof": c["ndof"], "chi2_ndof": c["chi2"] / c["ndof"],
            "E_fit_lo": m1.LO, "E_fit_hi": m1.HI, "sys_floor": 0.0, "xray_total_per_branch_pct": 0.0,
            "ratio_to_passport": A / PASSPORT["A_Bq"], "d_ratio": c["dA_Bq"] / PASSPORT["A_Bq"], "per_nuclide": per_nuclide}

def method2_block(j, n_lines):
    c = j["chain"]
    A = c["A_Bq"]
    BR = dict(m1.CHAIN)
    lines = [{"E_keV": ln["E_keV"], "nuclide": ln["nuclide"], "I_gamma_pct": ln.get("I_pct"),
              "branch": BR.get(ln["nuclide"], 1.0), "eps_peak": ln["eps_peak"], "weight_per_branch": ln["weight_per_branch"],
              "predicted_net": ln["predicted_net"], "kind": ln["kind"], "note": ln.get("note", ""),
              "E1_keV": ln.get("E1_keV"), "E2_keV": ln.get("E2_keV"),
              "I1_pct": SUM_I.get((round(ln.get("E1_keV") or 0, 3), round(ln.get("E2_keV") or 0, 3)), (None, None))[0],
              "I2_pct": SUM_I.get((round(ln.get("E1_keV") or 0, 3), round(ln.get("E2_keV") or 0, 3)), (None, None))[1]} for ln in c["lines"]]
    return {"A_Bq": A, "dA_Bq": c["dA_Bq"], "dA_stat_Bq": c["dA_stat_Bq"], "birge": c["birge"], "E1_Bq": c["E1_Bq"],
            "bg_amplitude": 1.0, "chi2": c["chi2"], "ndof": c["ndof"], "chi2_ndof": c["chi2"] / c["ndof"],
            "n_lines": n_lines, "n_channels_fit": int(sum(j["sel"])), "n_sum_peaks": j["n_sum"], "n_sum_peaks_total": j["n_sum"], "n_xray_energies": 0,
            "n_nodes": j["n_nodes"], "ratio_to_passport": A / PASSPORT["A_Bq"], "d_ratio": c["dA_Bq"] / PASSPORT["A_Bq"], "lines": lines}

def spectrum_block(jm1, jm2, jm2f, bg):
    A = jm1["chain"]["A_Bq"]
    stack = jm1["chain"]["stack"]
    scale = N_PER_BR / (A * live)
    trusted = {k: [bool(float(v) * scale >= N_EFF_MIN) for v in stack[k]] for k in stack}
    noise_frac = {k: (sum(v for v, t in zip(stack[k], trusted[k]) if not t) / sum(stack[k]) if sum(stack[k]) > 0 else 0.0) for k in stack}
    zero = [0.0] * len(jm1["e"])
    # #XR-1 (25.09, «отдельный видимый слой»): в М1 рентген размазан внутри общего Geant4-шаблона
    # каждого звена (fluo=1) и постфактум не отделим без нового прогона — слой остаётся нулевой
    # заглушкой (см. note нуклида XRAY). В М2 рентген — отдельные строки библиотеки, слой честный
    # (fit_gs2020_th232_m2.py: xray_stack = разница полного и без-рентгеновского вызова run_method2).
    xray2 = r4(jm2["chain"].get("xray_stack", zero))
    xray2f = r4(jm2f["chain"].get("xray_stack", zero))
    return {"model_counts": r4(jm1["chain"]["model"]), "model2_counts": r4(jm2["chain"]["model"]), "model2_full_counts": r4(jm2f["chain"]["model"]),
            # BG — приведённый фон отдельным слоем (оператор 25.09 «а почему фон не вычтен?»): фон ~половина спектра,
            # без слоя разрыв «измерение − модель» читается как невычтенный фон. Сумма слоёв = модель + фон.
            "stack": dict({k: r4(v) for k, v in stack.items()}, XRAY=zero, BG=r4(bg)),
            "trusted": dict(trusted, XRAY=[False] * len(zero), BG=[True] * len(zero)), "n_eff_min": N_EFF_MIN,
            "noise_frac": dict(noise_frac, XRAY=0.0, BG=0.0),
            "stack2": dict({k: r4(v) for k, v in jm2["chain"]["stack"].items()}, XRAY=xray2, BG=r4(bg)),
            "stack2_full": dict({k: r4(v) for k, v in jm2f["chain"]["stack"].items()}, XRAY=xray2f, BG=r4(bg)),
            "stack2_chan": {}, "stack2_chan_full": {}}

def fwhm_cal():
    pts = []
    used_pts = []
    for row in m1.PEAK_TABLE:
        e_file, e_lib, fwhm = row
        if e_lib == 1460.822:
            pts.append({"E_nominal": e_lib, "E_centroid": e_file, "fwhm_keV": fwhm, "d_fwhm_keV": 0.0, "res_pct": 100*fwhm/e_file,
                        "shift_keV": e_file - e_lib, "used": False, "reject": "фоновая линия K-40"})
        else:
            pts.append({"E_nominal": e_lib, "E_centroid": e_file, "fwhm_keV": fwhm, "d_fwhm_keV": 0.0, "res_pct": 100*fwhm/e_file,
                        "shift_keV": e_file - e_lib, "used": True, "reject": "", "n_lines_window": 1})
            used_pts.append((e_lib, fwhm))

    log_e = [math.log(e) for e, _ in used_pts]
    log_f = [math.log(f) for _, f in used_pts]
    coeffs = np.polyfit(log_e, log_f, 1)
    p = coeffs[0]
    log_k = coeffs[1]
    k = math.exp(log_k)

    rms_dev_pct = 0.0
    if len(used_pts) > 0:
        residuals = []
        for e, f in used_pts:
            pred = k * (e ** p)
            residuals.append((pred / f) - 1.0)
        rms_dev_pct = 100.0 * math.sqrt(sum(r**2 for r in residuals) / len(residuals))

    for q in pts:  # закон в точке и отклонение — поля таблицы ширин донора
        q["fwhm_model_keV"] = k * q["E_nominal"] ** p
        q["dev_pct"] = 100.0 * (q["fwhm_keV"] / q["fwhm_model_keV"] - 1.0)
    fwhm662_law = k * (661.657 ** p)
    res662_pct = 100.0 * fwhm662_law / 661.657

    return {"source": "таблица пиков прибора (AtomSpectra), линии Th-232 этого спектра", "k": k, "p": p, "rms_dev_pct": rms_dev_pct, "n_used": len(used_pts), "n_anchors": len(pts),
            "fwhm662_law": fwhm662_law, "fwhm662_cs": fwhm662_law, "res662_pct": res662_pct, "points": pts}

def lib_lines(cfg_path):
    cfg = yaml.safe_load(open(cfg_path, encoding="utf-8"))
    return cfg["library"]["lines"]

live = None

def main():
    global live

    jm1 = load(FILES["m1"])
    jm2 = load(FILES["m2"])
    jm2f = load(FILES["m2f"])
    cm1 = load(bgw(FILES["m1"]))
    cm2 = load(bgw(FILES["m2"]))
    cm2f = load(bgw(FILES["m2f"]))

    all_jsons = [jm1, jm2, jm2f, cm1, cm2, cm2f]
    ref_e = np.array(jm1["e"])
    for j in all_jsons[1:]:
        if not np.allclose(ref_e, np.array(j["e"]), atol=1e-6):
            raise SystemExit("ОТКАЗ: разные оси энергии в JSON подгонок")

    for name, j in [("jm1", jm1), ("jm2", jm2), ("jm2f", jm2f), ("cm1", cm1), ("cm2", cm2), ("cm2f", cm2f)]:
        c = j["chain"]
        if "stack" not in c:
            raise SystemExit(f"ОТКАЗ: stack отсутствует в {name}")
        sel = np.array(j["sel"], dtype=bool)
        model_sum = sum(c["model"][i] for i, s in enumerate(sel) if s)
        stack_sum = 0.0
        for k, v in c["stack"].items():
            stack_sum += sum(v[i] for i, s in enumerate(sel) if s)
        if model_sum > 0:
            rel_diff = abs(model_sum - stack_sum) / model_sum
            if rel_diff > 1e-3:
                raise SystemExit(f"ОТКАЗ: Σ вкладов звеньев ≠ модели цепочки в {name}")

    for name, j in [("jm2", jm2), ("jm2f", jm2f), ("cm2", cm2), ("cm2f", cm2f)]:
        if "lines" not in j["chain"]:
            raise SystemExit(f"ОТКАЗ: chain['lines'] отсутствует в {name}")

    s = m1.Spec(m1.bm.read(m1.XML_SAMPLE)[0])
    bsp = m1.bm.read(m1.XML_BG)[0]
    b = m1.Spec(bsp, "bg")   # #CAL-2: своя шкала фона по его реперам
    sp = m1.bm.read(m1.XML_SAMPLE)[0]

    if not hasattr(s, 'live_time'): raise SystemExit("Отсутствует атрибут live_time у sample")
    if not hasattr(sp, 'real'): raise SystemExit("Отсутствует атрибут real у sample spectrum object")
    if not hasattr(sp, 'cal'): raise SystemExit("Отсутствует атрибут cal у sample spectrum object")
    if not hasattr(b, 'live_time'): raise SystemExit("Отсутствует атрибут live_time у background")
    if not hasattr(bsp, 'real'): raise SystemExit("Отсутствует атрибут real у background spectrum object")
    if not hasattr(bsp, 'cal'): raise SystemExit("Отсутствует атрибут cal у background spectrum object")

    live = s.live_time
    k_bg = s.live_time / b.live_time

    e_of_ch = jm1["e"]
    counts = [int(c) for c in s.counts]
    bg_counts = r4(np.asarray(b.counts, float) * k_bg)

    assert len(s.counts) == len(b.counts) == len(e_of_ch), "Не совпадают длины массивов данных"

    n2 = len(lib_lines(LIB2_CFG))
    n05 = len(lib_lines(LIB05_CFG))

    donor_cfg_data = yaml.safe_load(open(DONOR_CFG, encoding="utf-8"))
    nuclides = []
    for n in donor_cfg_data["nuclides"]:
        nuclides.append({"key": n["key"], "label_ru": n["label_ru"], "label_en": n["label_en"], "color": n["color"], "note": n.get("note_ru", ""), "branching": n["br"]})
    nuclides.append({"key": "BG", "label_ru": "фон (приведён)", "label_en": "background", "color": "#b8b2a2", "note": "фон лаборатории × отношение живых времён; в альтернативном режиме — прямое измерение фона в сосуде Маринелли с дист. водой (промежуточный замер 11,7 ч), × отношение живых времён", "branching": 1.0})
    nuclides.append({"key": "XRAY", "label_ru": "K-рентген", "label_en": "K X-rays", "color": "#6b5f4a",
                      "note": "в методе 1 отдельно не выделяется (рождается внутри общего шаблона звена, "
                              "не отделим без нового прогона); в методе 2 — сумма строк библиотеки #XR-1", "branching": 1.0})

    reference_lines = [[float(l["e_kev"]), l["nuclide"], "%s %.1f" % (l["nuclide"], float(l["e_kev"]))] for l in lib_lines(LIB2_CFG)]

    meta = {"live_s": s.live_time, "real_s": float(sp.real), "bg_live_s": b.live_time, "bg_real_s": float(bsp.real), "bg_scale_time": k_bg,
            "cal_sample": {"coefs": m1.CAL_OWN["sample"]["coeffs"], "order": len(m1.CAL_OWN["sample"]["coeffs"]) - 1,   # #CAL-2: своя шкала, не полином файла
                           "n_channels": len(s.counts)},
            "cal_bg": {"coefs": m1.CAL_OWN["bg"]["coeffs"], "order": len(m1.CAL_OWN["bg"]["coeffs"]) - 1, "n_channels": len(b.counts)},
            "sys_floor_pct": 0.0, "xray_span_lo_keV": 0.0, "xray_span_hi_keV": 0.0,
            "template_decays": [{"nuclide": k, "n": int(round(N_PER_BR * br))} for k, br in m1.CHAIN],
            "nuclide_list_ru": ", ".join(k for k, _ in m1.CHAIN), "matrix_name": "Epoxy_crumb", "matrix_density_g_cm3": 0.7987,
            "matrix_composition": [], "template_source": "Geant4, стенд gs2020_marinelli, out_v5"}

    # #CAL-1: у фона своя шкала — плотность отсчётов фона на кэВ переносится на сетку образца по энергии
    eb, es = np.array([b.channel_to_energy(i) for i in range(len(b.counts))]), np.asarray(e_of_ch, float)
    bg_arr = np.interp(es, eb, np.asarray(b.counts, float) / np.gradient(eb)) * np.gradient(es) * k_bg
    bg_counts = r4(bg_arr)
    # #PAGE-4 (оператор 25.09 «шаблоны отсортируй по вкладу»): легенда — по убыванию вклада слоя в модель М1 (фон — тоже слой),
    # K-рентген (в этой модели нулевой) — последним
    contrib = {k: float(np.sum(v)) for k, v in jm1["chain"]["stack"].items()}
    contrib["BG"] = float(np.sum(bg_arr))
    nuclides.sort(key=lambda n: (n["key"] == "XRAY", -contrib.get(n["key"], 0.0)))
    # Панель «cs»: фон — реальное измерение «Маринелли 1 л + дист. вода» (оператор 26.09), не модельное ослабление.
    bw = m1.Spec(m1.bm.read(m1.CAL.BKG_WATER_XML)[0], "bgw")
    k_bg_w = s.live_time / bw.live_time
    ebw = np.array([bw.channel_to_energy(i) for i in range(len(bw.counts))])
    bg_water_arr = np.interp(es, ebw, np.asarray(bw.counts, float) / np.gradient(ebw)) * np.gradient(es) * k_bg_w
    fw = fwhm_cal()

    # Оператор 26.09 «спектр на всех картинках на 3000 обрежь»: массивы по каналам (длина n_full, реально до
    # ~4774 кэВ) обрезаются для ОТОБРАЖЕНИЯ; окно подгонки LO..HI=150..3600 (#SUM-1) не меняется — обрезка
    # чисто визуальная, χ² считается по полному окну.
    n_full = len(es)
    n_show = int(np.searchsorted(es, 3000.0, side="right"))

    def crop_spec(obj):
        if isinstance(obj, list):
            return obj[:n_show] if len(obj) == n_full else obj
        if isinstance(obj, dict):
            return {k: crop_spec(v) for k, v in obj.items()}
        return obj

    data = {
        "meta": meta,
        "fwhm_cal": fw,
        "passport": PASSPORT,
        "nuclides": nuclides,
        "channels": [],
        "spectrum": crop_spec(dict({"e_of_ch": r4(e_of_ch), "counts": counts, "bg_counts": bg_counts}, **spectrum_block(jm1, jm2, jm2f, bg_arr))),
        "cs": {
            "method1": method1_block(cm1),
            "method2": method2_block(cm2, n2),
            "method2_full": method2_block(cm2f, n05),
            "spectrum": crop_spec(spectrum_block(cm1, cm2, cm2f, bg_water_arr))
        },
        "method1": method1_block(jm1),
        "method2": method2_block(jm2, n2),
        "method2_full": method2_block(jm2f, n05),
        "library": {"i_threshold_pct": 2.0, "fixed_n": n2, "full_n": n05, "full_threshold_pct": 0.5},
        "reference_lines": reference_lines
    }

    os.makedirs(PAGE, exist_ok=True)
    with open(DST, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, separators=(",", ":"), allow_nan=False)

    size_kb = os.path.getsize(DST) / 1024.0
    print(f"Файл сохранен: {DST} ({size_kb:.1f} КБ)")

    blocks_to_check = [
        ("method1 (измеренный)", data["method1"]),
        ("method2 (измеренный)", data["method2"]),
        ("method2_full (измеренный)", data["method2_full"]),
        ("method1 (ослабленный фон)", data["cs"]["method1"]),
        ("method2 (ослабленный фон)", data["cs"]["method2"]),
        ("method2_full (ослабленный фон)", data["cs"]["method2_full"])
    ]

    for name, block in blocks_to_check:
        ratio = block["ratio_to_passport"]
        chi2_ndof = block["chi2_ndof"]
        print(f"{name}: A = {block['A_Bq']:.1f} ± {block['dA_Bq']:.1f} Бк, к паспорту {ratio:.3f}, χ²/ν {chi2_ndof:.2f}")

    print(f"Закон ПШПВ: k={fw['k']:.4f}, p={fw['p']:.4f}, rms={fw['rms_dev_pct']:.2f}%")

    issues = []
    for name, block in blocks_to_check:
        if abs(block["ratio_to_passport"] - 1.0) > 0.06:
            issues.append(f"  {name}: отклонение от паспорта > 6%")

    spec_blocks = [data["spectrum"], data["cs"]["spectrum"]]
    for i, sb in enumerate(spec_blocks):
        label = "измеренный фон" if i == 0 else "ослабленный фон"
        nf = sb.get("noise_frac", {})
        for k, val in nf.items():
            if val > 0.5:
                issues.append(f"  {label}, нуклид {k}: доля шума > 50%")

    print("ТРЕБУЕТ ТОЛКОВАНИЯ:")
    if not issues:
        print("  нет")
    else:
        for issue in issues:
            print(issue)

if __name__ == "__main__":
    main()
