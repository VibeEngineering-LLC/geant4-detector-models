# -*- coding: utf-8 -*-
r"""#GS-74: данные подвкладки «метод 2» черники (страница GS2020): блок method2, стек слоёв метода 2, строки сравнения.
Источник — fit_berry_m2_bgw.json (scripts\fit_gs2020_berry_m2.py). Чистая функция build(): файлов не читает и не пишет;
вызывается из export_page_gs2020_berry.py. Слои метода 2 в сумме с фоном дают модель+фон побитово — проверяется.
Спека: scripts\specs\SPEC-export_berry_m2.md"""
import math
import numpy as np


def refuse(msg):
    raise SystemExit("ОТКАЗ: " + msg)


def r4(a):
    return [round(float(x), 4) for x in a]


KEYS = {"Cs137": "CS137", "K40": "K40", "SrY90": "SRY90"}
ELEM = {"CS137": "Ba", "K40": "Ar"}
EXTRA_NUCLIDES = [
    {"key": "M2_CS137", "label_ru": "Cs-137: γ-линии и K-рентгеновское излучение Ba", "label_en": "Cs-137: gamma and Ba X-ray lines", "color": "#b8860b", "branching": 1.0,
     "note": "отклик Geant4 на линии библиотеки Cs-137: γ 661,7 кэВ и K-рентгеновское излучение Ba (#XR-1), с собственными линиями детектора"},
    {"key": "M2_K40", "label_ru": "K-40: линия 1460,8 кэВ и рентгеновское излучение Ar", "label_en": "K-40: 1460.8 keV line and Ar X-ray", "color": "#2e7d32", "branching": 1.0,
     "note": "отклик Geant4 на линию 1460,8 кэВ и рентгеновское излучение Ar (#XR-1); рентгеновское излучение Ar до кристалла не доходит"},
    {"key": "M2_BCS137", "label_ru": "внешнее тормозное излучение (β/e⁻) Cs-137", "label_en": "external bremsstrahlung (beta/e-) of Cs-137", "color": "#e3c76f", "branching": 1.0,
     "note": "тормозное излучение β-электронов и электронов конверсии Cs-137 в пробе, сосуде и корпусе (Geant4)"},
    {"key": "M2_BK40", "label_ru": "внешнее тормозное излучение (β/e⁻) K-40", "label_en": "external bremsstrahlung (beta/e-) of K-40", "color": "#8fc28f", "branching": 1.0,
     "note": "тормозное излучение β-электронов K-40 в пробе, сосуде и корпусе (Geant4)"},
]
NOTE_M2 = {"CS137": "γ-линия и K-рентгеновское излучение Ba из библиотеки, отклик Geant4, плюс внешнее тормозное излучение (β/e⁻); амплитуда одна на нуклид",
           "K40": "предварительно: линия 1460,8 кэВ из библиотеки, отклик Geant4, плюс внешнее тормозное излучение β; целиком зависит от r(E) фона",
           "SRY90": "мешающий параметр: внешнее тормозное излучение β Sr-90+Y-90 (шаблон Geant4); активность Sr-90 не публикуется"}


def build(jm2, bg, n_show, mass_kg, ref_cs_per_kg, ref_k_per_kg):
    if jm2["col_names"] != ["Cs137", "K40", "SrY90"]:
        refuse("col_names не совпадает с ожидаемым порядком")
    if list(jm2["nuclides"]) != ["Cs137", "K40", "SrY90"]:
        refuse("порядок nuclides не совпадает")
    if abs(jm2["mass_kg"] - mass_kg) >= 1e-9:
        refuse("масса не совпадает")
    if len(bg) != len(jm2["model"]):
        refuse("длина фона не совпадает с моделью")

    cols = np.asarray(jm2["cols"], float)
    bc = {"Cs137": np.asarray(jm2["beta_cols"]["Cs137"], float), "K40": np.asarray(jm2["beta_cols"]["K40"], float)}
    bgv = np.asarray(bg, float)
    model = np.asarray(jm2["model"], float)

    col_sum = np.sum(cols, axis=0)
    if np.max(np.abs(col_sum - model)) > 1e-6:
        refuse("сумма столбцов не равна модели")

    lines_cs = cols[0] - bc["Cs137"]
    lines_k = cols[1] - bc["K40"]
    if lines_cs.min() < -1e-6 or lines_k.min() < -1e-6:
        refuse("слой линий отрицателен")

    full = {
        "M2_CS137": lines_cs,
        "M2_K40": lines_k,
        "M2_BCS137": bc["Cs137"],
        "M2_BK40": bc["K40"],
        "SRY90": cols[2],
        "BG": bgv
    }

    sum_layers = np.sum(np.sum(np.array(list(full.values())), axis=0))
    model_plus_bg = float(model.sum() + bgv.sum())
    if abs(sum_layers - model_plus_bg) / model_plus_bg > 1e-9:
        refuse("сумма слоёв не равна модели+фону")

    non_bg_sum = np.sum(np.array([lines_cs, lines_k, bc["Cs137"], bc["K40"], cols[2]]), axis=0)
    if np.max(np.abs(non_bg_sum - model)) > 1e-6:
        refuse("сумма слоёв без фона не равна модели")

    stack2 = {key: r4(arr[:n_show]) for key, arr in full.items()}

    chi2_ndof = jm2["chi2"] / jm2["ndof"]
    birge = math.sqrt(max(chi2_ndof, 1.0))
    total_all = float(model.sum() + bgv.sum())

    per = {}
    fit_keys = ["Cs137", "K40", "SrY90"]
    page_keys = ["CS137", "K40", "SRY90"]
    for i, fk in enumerate(fit_keys):
        pk = page_keys[i]
        A = jm2["nuclides"][fk]["A_Bq"]
        dst = jm2["nuclides"][fk]["dA_stat_Bq"]
        nuclide_total = float(cols[i].sum())
        share = nuclide_total / total_all

        if fk == "Cs137":
            beta_frac = float(bc["Cs137"].sum()) / nuclide_total
            ref_ratio = (A / mass_kg) / ref_cs_per_kg
            nuisance = False
        elif fk == "K40":
            beta_frac = float(bc["K40"].sum()) / nuclide_total
            ref_ratio = (A / mass_kg) / ref_k_per_kg
            nuisance = False
        else:
            beta_frac = None
            ref_ratio = None
            nuisance = True

        entry = {
            "A_Bq": A,
            "dA_Bq": dst * birge,
            "dA_stat_Bq": dst,
            "per_kg": A / mass_kg,
            "ref_ratio": ref_ratio,
            "share": share,
            "beta_frac": beta_frac,
            "nuisance": nuisance,
            "note_ru": NOTE_M2[pk]
        }

        if nuisance:
            entry["A_Bq"] = None
            entry["dA_Bq"] = None
            entry["dA_stat_Bq"] = None
            entry["per_kg"] = None

        per[pk] = entry

    lines = []
    n_xray = 0
    n_below = 0
    for ln in jm2["lines"]:
        page_key = KEYS[ln["nuclide"]]
        if ln["E_keV"] < 40.0:
            n_xray += 1   # рентген в библиотеке (все, включая Ar ниже 25 кэВ)
        if ln["E_keV"] < 25.0:   # #GS-81 (оператор 10.10 «ниже 25 кэв не указывай линии»): в подгонке остаются, в таблице нет
            n_below += 1
            continue
        if ln["E_keV"] < 40.0:
            text = "рентгеновское излучение " + ELEM[page_key] + ": " + ln["note"]
        else:
            text = ln["note"]
        if ln["eps_peak"] == 0:
            text = text + "; до кристалла не доходит (ε = 0)"

        lines.append({
            "E_keV": ln["E_keV"],
            "nuclide": page_key,
            "I_gamma_pct": ln["I_pct"],
            "branch": 1.0,
            "eps_peak": ln["eps_peak"],
            "weight_per_branch": ln["weight_per_branch"],
            "predicted_net": ln["predicted_net"],
            "kind": "line",
            "note": text,
            "E1_keV": None,
            "E2_keV": None,
            "I1_pct": None,
            "I2_pct": None
        })

    method2 = {
        "A_Bq": per["CS137"]["A_Bq"],
        "dA_Bq": per["CS137"]["dA_Bq"],
        "dA_stat_Bq": per["CS137"]["dA_stat_Bq"],
        "birge": birge,
        "E1_Bq": None,
        "bg_amplitude": 1.0,
        "chi2": jm2["chi2"],
        "ndof": jm2["ndof"],
        "chi2_ndof": chi2_ndof,
        "n_lines": len(lines) + n_below,
        "n_lines_shown": len(lines),
        "n_below_25": n_below,
        "n_xray_shown": n_xray - n_below,
        "n_channels_fit": jm2["n_channels_fit"],
        "n_sum_peaks": 0,
        "n_sum_peaks_total": 0,
        "n_xray_energies": n_xray,
        "n_nodes": jm2["n_nodes"],
        "ratio_to_passport": per["CS137"]["ref_ratio"],
        "d_ratio": 0.0,
        "shape_662": jm2["shape"]["661.657"][0],
        "shape_1460": jm2["shape"]["1460.822"][0],
        "E_fit_lo": jm2["lo"],
        "E_fit_hi": jm2["hi"],
        "per_nuclide": per,
        "lines": lines,
        "computed": True
    }

    cmp = {
        "cs": {
            "lab": "метод 2",
            "A": per["CS137"]["A_Bq"],
            "dA": per["CS137"]["dA_Bq"],
            "kind": "m2",
            "note": "функция полного поглощения × библиотека линий, отклик Geant4"
        },
        "k": {
            "lab": "метод 2",
            "A": per["K40"]["A_Bq"],
            "dA": per["K40"]["dA_Bq"],
            "kind": "m2",
            "note": "предварительно: целиком зависит от r(E) фона"
        }
    }

    nuclides_extra = [dict(x) for x in EXTRA_NUCLIDES]
    share_check = {"sum_layers": float(sum_layers), "model_plus_bg": model_plus_bg}

    return {
        "method2": method2,
        "stack2": stack2,
        "nuclides_extra": nuclides_extra,
        "cmp": cmp,
        "share_check": share_check
    }
