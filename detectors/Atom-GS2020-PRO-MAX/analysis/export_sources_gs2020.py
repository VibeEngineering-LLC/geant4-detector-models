# -*- coding: utf-8 -*-
#GS-24: данные вкладок источников страницы GS2020 (K-40 — из fit_kcl_bgw.json) → gs2020_sources.json
import sys
import os
import json
import math

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import fit_gs2020_kcl as fk

KCL = r"C:\g4work\gs2020\kcl_v4w85_83\fit_kcl_bgw.json"   # Маринелли 2 v4, конус колодца (29.09)
OUT = r"D:\cloud-folder\work-folder\GEANT4\web\gs2020-th232-page\gs2020_sources.json"

def f(x, d):
    """Форматирование числа: группировка целой части через U+202F, запятая вместо точки."""
    if x is None:
        return "None"
    
    # Формируем строку с плавающей точкой
    raw = "%.*f" % (d, x)
    
    parts = raw.split('.')
    integer_part = parts[0]
    decimal_part = parts[1] if len(parts) > 1 else ""
    
    # Обработка знака
    sign = ""
    if integer_part.startswith('-'):
        sign = "-"
        integer_part = integer_part[1:]
    elif integer_part.startswith('+'):
        sign = "+"
        integer_part = integer_part[1:]
        
    # Группировка целой части по 3 цифры с конца
    # Используем узкий неразрывный пробел U+202F
    formatted_int = ""
    count = 0
    for char in reversed(integer_part):
        if count > 0 and count % 3 == 0:
            formatted_int = "\u202f" + formatted_int
        formatted_int = char + formatted_int
        count += 1
        
    result = sign + formatted_int
    if decimal_part:
        result += "," + decimal_part
        
    return result

# Загрузка данных
try:
    with open(KCL, 'r', encoding='utf-8') as f_in:
        data = json.load(f_in)
except Exception as e:
    raise SystemExit(f"Ошибка чтения файла {KCL}: {e}")

# Извлечение полей
A_Bq = data['A_Bq']
dA_stat_Bq = data['dA_stat_Bq']
A_expected_Bq = data['A_expected_Bq']
chi2 = data['chi2']
ndof = data['ndof']
shape_1460 = data['shape_1460']
live_s = data['live_s']
k_bg = data['k_bg']
lo = data['lo']
hi = data['hi']

e_arr = data['e']
net_arr = data['net']
model_arr = data['model']
bg_arr = data['bg']

# Проверка длин массивов
len_e = len(e_arr)
len_net = len(net_arr)
len_model = len(model_arr)
len_bg = len(bg_arr)

if not (len_e == len_net == len_model == len_bg):
    raise SystemExit("Длины массивов e, net, model, bg различаются")

# Фильтрация индексов: 20 <= e <= 3000
idx = [i for i, val in enumerate(e_arr) if 20 <= val <= 3000]

if not idx:
    raise SystemExit("Массив idx пуст (нет каналов в диапазоне 20-3000 кэВ)")

# Вычисление ratio
ratio = A_Bq / A_expected_Bq

# Формирование спектра
spectrum_e = [round(e_arr[i], 2) for i in idx]
spectrum_meas = [net_arr[i] + bg_arr[i] for i in idx]
spectrum_bg = [bg_arr[i] for i in idx]
spectrum_model = [model_arr[i] for i in idx]

# resid = (net - model) / sqrt(max(net + bg * (1 + k_bg), 1.0))
spectrum_resid = []
for i in idx:
    net_val = net_arr[i]
    model_val = model_arr[i]
    bg_val = bg_arr[i]
    
    denominator_sq = max(net_val + bg_val * (1 + k_bg), 1.0)
    denom = math.sqrt(denominator_sq)
    resid_val = (net_val - model_val) / denom
    spectrum_resid.append(round(resid_val, 2))

spectrum = {
    "ch0": idx[0],
    "e": spectrum_e,
    "meas": [round(v, 1) for v in spectrum_meas],
    "bg": [round(v, 1) for v in spectrum_bg],
    "model": [round(v, 1) for v in spectrum_model],
    "resid": spectrum_resid,
    "lo": lo,
    "hi": hi,
    "view": [25, 3000],
    "anchors": [[1460.822, "K-40 1460,8"]]
}

# Формирование fill
k40_mass = f(fk.MASS_G, 0) + " г"
k40_vol = f(fk.VOL_ML, 0) + " мл"
k40_rho = f(fk.MASS_G / fk.VOL_ML, 2) + " г/см³"
k40_kfrac = f(fk.K_FRAC, 5)
k40_spec = f(fk.K40_BQ_PER_G_K, 2) + " Бк/г K"
k40_aexp = f(A_expected_Bq, 0) + " Бк"
k40_live = f(live_s / 3600, 1) + " ч"
k40_a = f(A_Bq, 0) + " ± " + f(dA_stat_Bq, 0) + " Бк"
k40_ratio = f(ratio, 4)
k40_chi2 = f(chi2 / ndof, 3)
k40_shape = f(shape_1460, 2)

fill = {
    "k40_mass": k40_mass,
    "k40_vol": k40_vol,
    "k40_rho": k40_rho,
    "k40_kfrac": k40_kfrac,
    "k40_spec": k40_spec,
    "k40_aexp": k40_aexp,
    "k40_live": k40_live,
    "k40_a": k40_a,
    "k40_ratio": k40_ratio,
    "k40_chi2": k40_chi2,
    "k40_shape": k40_shape
}

# Формирование results
note1 = f"{k40_mass} × {k40_kfrac} (доля K в KCl) × {k40_spec}; чистота «ч» принята полной"
note2 = f"шаблон полного распада K-40 (Geant4), окно {f(lo, 0)}–{f(hi, 0)} кэВ, фон Маринелли+вода; χ²/ν = {k40_chi2}; форма пика 1460,8 кэВ: χ²/ν = {k40_shape}"
note2 = "ядро как у тория (без хвоста, ПШПВ ×1,05); " + note2
# W-155: вариант с донорским хвостом ядра (run_gs2020_kcl.sh → fit_kcl_bgw_tail.json) — строка чувствительности
tl = json.load(open(KCL.replace(".json", "_tail.json"), encoding="utf-8"))
row_tail = {"lab": "метод 1: ядро с хвостом (донор Гамма-1С)", "A": tl["A_Bq"], "dA": tl["dA_stat_Bq"], "ratio": tl["A_Bq"] / tl["A_expected_Bq"],
            "note": "χ²/ν = %s; форма пика 1460,8 кэВ: χ²/ν = %s" % (f(tl["chi2"] / tl["ndof"], 3), f(tl["shape_1460"], 2))}
fill["k40_tail_ratio"] = f(row_tail["ratio"], 4)
# #GS-24 (оператор 28.09 «только метод 1 и 2»): метод 2 — fit_gs2020_kcl_m2.py → fit_kcl_m2_bgw.json
m2 = json.load(open(KCL.replace("fit_kcl_bgw.json", "fit_kcl_m2_bgw.json"), encoding="utf-8"))
row_m2 = {"lab": "метод 2: линия × отклик Geant4", "A": m2["A_Bq"], "dA": m2["dA_stat_Bq"], "ratio": m2["ratio"],
          "note": "γ 1460,82 кэВ, выход %s %% (ENSDF, API IAEA); χ²/ν = %s; форма пика: χ²/ν = %s; с выходом %s %% (LNHB-DDEP 2025) — %s от ожидаемой"
                  % (f(m2["Igamma_pct"], 2), f(m2["chi2"] / m2["ndof"], 3), f(m2["shape_1460"], 2), f(m2["Igamma_lnhb_pct"], 2), f(m2["A_lnhb_Bq"] / m2["A_expected_Bq"], 4))}
fill.update({"k40_m2_a": f(m2["A_Bq"], 0) + " ± " + f(m2["dA_stat_Bq"], 0) + " Бк", "k40_m2_ratio": f(m2["ratio"], 4),
             "k40_m2_chi2": f(m2["chi2"] / m2["ndof"], 2), "k40_ig": f(m2["Igamma_pct"], 2), "k40_ig_lnhb": f(m2["Igamma_lnhb_pct"], 2)})

results = [
    {
        "lab": "ожидаемая по массе",
        "A": A_expected_Bq,
        "dA": None,
        "ratio": 1.0,
        "note": note1
    },
    {
        "lab": "метод 1: полный спектр",
        "A": A_Bq,
        "dA": dA_stat_Bq,
        "ratio": ratio,
        "note": note2
    },
    row_m2,
    row_tail
]

# Итоговый объект
output_obj = {
    "fill": fill,
    "order": ["th", "k40", "cs137", "ra226"],
    "sources": {
        "th": {"label": "Th-232", "status": "ready"},   # #GS-29 (оператор 28.09): слово «КИ» по GS2020 не писать нигде
        "k40": {
            "label": "K-40 (KCl)",
            "status": "ready",
            "results": results,
            "spectrum": spectrum
        },
        "cs137": {"label": "Cs-137", "status": "soon"},
        "ra226": {"label": "Ra-226", "status": "soon"}
    }
}

# Запись в файл
with open(OUT, 'w', encoding='utf-8') as f_out:
    json.dump(output_obj, f_out, ensure_ascii=False, separators=(",", ":"))

# Вычисление размера файла в КБ
file_size_bytes = os.path.getsize(OUT)
file_size_kb = file_size_bytes / 1024.0

# Печать итоговой строки
print(f"gs2020_sources.json: каналов {len(idx)} (с {idx[0]}), A {k40_a}, к ожидаемой {k40_ratio}, {f(file_size_kb, 0)} КБ")
