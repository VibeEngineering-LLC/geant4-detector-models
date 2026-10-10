# -*- coding: utf-8 -*-
# #GS-70 (оператор 10.10 «вкладку на вебке оформляй»): источник «черника» для gs2020_sources.json — из подгонки метода 1
# fit_berry_bgw_R2_g014_beta.json (fit_gs2020_berry.py, фон × r2, β Sr+Y мешающий). Вторые пути — results\RESULT-2026-10-10-berry-gs70.md.
import json, math
FIT = r"C:\g4work\gs2020\berry\fit_berry_bgw_R2_g014_beta.json"
def berry_source():
    j = json.load(open(FIT, encoding="utf-8")); N = j["nuclides"]; kb = 72547.5 / 119423.3
    idx = [i for i, x in enumerate(j["e"]) if 20 <= x <= 3000]; g = lambda a, d=1: [round(a[i], d) for i in idx]
    res = [{"lab": "метод 1: Cs-137", "A": N["Cs137"]["A_Bq"], "dA": N["Cs137"]["dA_stat_Bq"], "ratio": None,
            "note": ("%.0f Бк/кг сушёной пробы; систематическая погрешность ±4 %%; χ²/ν %.3f, форма пика 661,7 кэВ χ²/ν %.2f" % (N["Cs137"]["A_Bq"] / j["mass_kg"], j["chi2"] / j["ndof"], j["shape"]["661.657"][0])).replace(".", ",")},
           {"lab": "Cs-137: площадь пика × ε набора", "A": 1952.5, "dA": 6.6, "ratio": None, "note": "площадь СпектраЛайна 1 632 436 × ε 1,1965·10⁻² × f(ρ) 1,132 (экстраполяция ниже 0,70 г/см³)"},
           {"lab": "Cs-137: независимый пересчёт", "A": 1907.0, "dA": None, "ratio": None, "note": "собственная подгонка пика, неопределённость ±4–5 %"},
           {"lab": "метод 1: K-40", "A": N["K40"]["A_Bq"], "dA": N["K40"]["dA_stat_Bq"], "ratio": None,
            "note": "%.0f Бк/кг (%.0f Бк); погрешность вычитания фона ±14 Бк, систематическая погрешность ±30 Бк (в абсолютных значениях активности пробы; зависит от поправки r(E) на ослабление фона); предварительно" % (N["K40"]["A_Bq"] / j["mass_kg"], N["K40"]["A_Bq"])},
           {"lab": "K-40: независимый пересчёт", "A": 158.0, "dA": None, "ratio": None, "note": "155–162 Бк: нетто пика 1461 минус фоновый пик × r 1,058–1,063"}]
    n, m, b = j["net"], j["model"], j["bg"]
    spec = {"ch0": idx[0], "e": g(j["e"], 2), "meas": [round(n[i] + b[i], 1) for i in idx], "bg": g(b), "model": g(m), "lo": 150.0, "hi": 3000.0, "view": [25, 3000],
            "resid": [round((n[i] - m[i]) / math.sqrt(max(n[i] + b[i] * (1 + kb), 1.0)), 2) for i in idx],
            "comps": [{"label": l, "color": c, "data": g(j["cols"][k], 2)} for k, (l, c) in enumerate((("Cs-137", "#b8860b"), ("K-40", "#2e7d32"), ("β Sr-90+Y-90", "#7b1fa2")))],
            "anchors": [[661.657, "Cs-137 661,7"], [1460.822, "K-40 1460,8"]]}
    return {"label": "Черника: Cs-137 + K-40", "status": "ready", "results": res, "spectrum": spec}
