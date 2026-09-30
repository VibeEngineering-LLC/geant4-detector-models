# -*- coding: utf-8 -*-
r"""Выгрузка данных страницы хода работ GS2020 Th-232 (web/gs2020-th232/data.js). План:
D:\cloud-folder\work-folder\GEANT4\PLAN-2026-09-25-gs2020-th232.md. Запуск: python export_progress_gs2020.py"""

import sys, os, json, csv, math, datetime
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"D:\repos-folder\repos\geant4-detector-models\common\py")
import becqmoni as bm

REF = r"D:\cloud-folder\Дозиметрия\Спектры\Atom GS2020 PRO MAX\Референсы"
SPE_TH = os.path.join(REF, "Калибровка Th-232 (без вычета фона).xml")
SPE_BG = os.path.join(REF, "Фон лаба S31_18.xml")
OUT_DIR = r"C:\g4work\gs2020\run_marinelli\out_v5"  # ред. 4 геометрии, сетка 2e7 (оператор 25.09)
DST = r"D:\cloud-folder\work-folder\GEANT4\web\gs2020-th232\data.js"

import glob, yaml
TH232_CFG = r"D:\repos-folder\repos\geant4-detector-models\detectors\Gamma-1S\web-th232\configs\th232.yaml"
with open(TH232_CFG, encoding="utf-8") as _f:
    LINE_NUC = {float(l["e_kev"]): l["nuclide"] for l in yaml.safe_load(_f)["library"]["lines"]}

PASSPORT = {"bq_per_kg": 910.0, "unc_pct": 6.0, "mass_g": 1052.0}
E_MAX = 3000
STEP = 2

STAGES = [
    {"n": 1, "name": "Геометрия сосуда и источника в GDML", "status": "done", "note": "ред. 4 одобрена оператором: PP-вкладыш колодца 2 мм, корпус 109 мм, ρ 0,7987; дно колодца на верх прибора z=41,5 (W-133…W-135)"},
    {"n": 2, "name": "Стенд gs2020_marinelli, приёмка #SA-6 (K-40 1e5)", "status": "done", "note": "ред. 4: 1e5 за 6 с, 100000/100000 событий, пик 1461 — 0,532 %, GeomNav1002 ×4"},
    {"n": 3, "name": "Моно-прогоны 6 линий по 2e6 (ред. 1)", "status": "rejected", "note": "недействительны: сосуд стоял на торце (W-133), не было PP колодца"},
    {"n": 4, "name": "Предварительная NNLS, 4 линии, фон свободный", "status": "rejected", "note": "χ²/ν 180, фон 1,38 — не метод"},
    {"n": 5, "name": "Формат вывода стенда = формат шаблона метода + ионный режим (метод 1)", "status": "done", "note": "v2: шапка read_template, ion:Z:A + nucleusLimits; тест Tl-208 2e4 OK"},
    {"n": 6, "name": "Метод 2: прямые прогоны на 20 линиях библиотеки", "status": "done", "note": "28 энергий (20 линий + 8 для сумм-пиков) × 2e7, отрезками по 5e6"},
    {"n": 7, "name": "Метод 1: полный распад 8 звеньев цепочки", "status": "done", "note": "8 звеньев, n/BR = 5,5e7 (Tl208 1,977e7), отрезками по 5e6; шаблон цепочки = сумма звеньев"},
    {"n": 8, "name": "Подгонка A2 (фон фиксирован k) + E1, сверка с паспортом", "status": "done", "note": "Цепочка в равновесии, паспорт 957,3 ± 57,4 Бк. М1: 985,0 ± 1,6 Бк (×1,029), E1 978,7, χ²/ν 4,35. М2 (библиотека 2 %, 20 линий): 1192,4 ± 2,0 Бк (×1,246), χ²/ν 5,04; М2 (библиотека 0,5 % ENSDF, 42 линии, сетка 47 × 2e7): 1096,6 ± 1,6 Бк (×1,146), E1 1095,0, χ²/ν 4,12 — расхождение с М1 21 % → 11 %. Причина: модель М2 на 1 Бк на 18 % ниже М1 (1000–1500 кэВ — 0,58: линии Ac228 <2 % не в библиотеке; 300–1000 — 0,75…0,80: их комптон и тормозное β; >2750 — 0,09: суммирование). С библиотекой 0,5 %: модель М2/М1 на 1 Бк 0,894 (1000–1500: 0,77; >2750: 0,09 — суммирование, тормозное β не моделируется). В окнах пиков М2 те же ×1,25 — подложка. Шаблоны n/BR 5,5e7, сетка 28 × 2e7"},
    {"n": 9, "name": "Итоговый раздел страницы + независимый пересчёт (#SA-10)", "status": "work", "note": "раздел «Итог» опубликован как рабочий. #SA-10: метод площадей смещён (суррогат 0,51–0,92); с ослабленным по K-40 фоном пики и континуум сходятся, М1 1243–1261 Бк; площади прибора ≤ 900–1016 Бк. Итог не установлен — ждёт замеров фона с водой и KCl"}
]

def count_geomnav(path):
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8", errors="replace") as f:
        content = f.read()
    return content.count("GeomNav1002")

def density(spec):
    """Отсчёты на кэВ в родных каналах файла (шкала своя у каждого файла, #CAL-0).
    Прежняя перебиновка на 2 кэВ давала пилу: канал ~2,4 кэВ шире шага."""
    ch = np.arange(len(spec.n))
    lo, hi = spec.energy(ch - 0.5), spec.energy(ch + 0.5)
    return (lo + hi) / 2, np.asarray(spec.n, dtype=float) / (hi - lo)

os.environ.setdefault("SPECTRAVIBE_ROOT", r"D:\cloud-folder\Дозиметрия\ИИ\1 Скилы\0_Work\gamma-spectrum-analysis")
sys.path.insert(0, r"D:\repos-folder\repos\geant4-detector-models\detectors\Gamma-1S\analysis")
from mix_unfold_g1s import read_template  # формат шаблона метода, импорт (§33)

eff_list = []
for csv_path in sorted(glob.glob(os.path.join(OUT_DIR, "grid_mar_E*.csv")),
                       key=lambda p: float(os.path.basename(p)[10:-4])):
    E0 = float(os.path.basename(csv_path)[10:-4])
    hist, n, _ = read_template(csv_path)
    key = LINE_NUC.get(E0, "сумм.")  # энергии компонент и crossover сумм-пиков
    log_name = os.path.basename(csv_path)[:-4] + ".log"
    # пик — строго полное поглощение: бин 1 кэВ с E0 и соседний снизу (как grid_response донора)
    b0 = math.floor(E0) + 0.5
    peak = hist.get(b0, 0.0) + hist.get(b0 - 1.0, 0.0)
    total = float(sum(v for k, v in hist.items() if k >= 1.0))
    eff_peak = peak / n
    eff_total = total / n
    d_eff_peak = math.sqrt(peak) / n
    geomnav = count_geomnav(os.path.join(OUT_DIR, log_name))
    eff_list.append({
        "nuclide": key,
        "E_keV": E0,
        "n_events": n,
        "eff_peak": eff_peak,
        "eff_total": eff_total,
        "d_eff_peak": d_eff_peak,
        "geomnav": geomnav
    })

sp, _ = bm.read(SPE_TH)
bg, _ = bm.read(SPE_BG)

k_bg = sp.live / bg.live
E_keV_centers, counts = density(sp)
E_bg, d_bg = density(bg)
bg_scaled = np.interp(E_keV_centers, E_bg, d_bg) * k_bg
sel = (E_keV_centers > 0) & (E_keV_centers <= E_MAX)
E_keV_centers, counts, bg_scaled = E_keV_centers[sel], counts[sel], bg_scaled[sel]
net = counts - bg_scaled

data = {
    "generated": datetime.datetime.now().isoformat(timespec="seconds"),
    "stages": STAGES,
    "eff": eff_list,
    "passport": PASSPORT,
    "activity_passport_bq": float(PASSPORT["bq_per_kg"] * PASSPORT["mass_g"] / 1000),
    "live_s": float(sp.live),
    "live_bg_s": float(bg.live),
    "k_bg": k_bg,
    "spectrum": {
        "E_keV": [round(float(x), 1) for x in E_keV_centers],
        "counts": [round(float(x), 1) for x in counts],
        "bg": [round(float(x), 2) for x in bg_scaled],
        "net": [round(float(x), 2) for x in net]
    }
}

os.makedirs(os.path.dirname(DST), exist_ok=True)
with open(DST, "w", encoding="utf-8") as f:
    f.write("window.GS2020_PROGRESS = " + json.dumps(data, ensure_ascii=False) + ";\n")

print(DST)
print(len(eff_list))
print(f"{k_bg:.5f}")

mask = (E_keV_centers >= 150) & (E_keV_centers <= 3000)
sum_counts = int(np.sum(counts[mask]))
sum_bg = int(np.sum(bg_scaled[mask]))
print(sum_counts, sum_bg)

for item in eff_list:
    gn = str(item["geomnav"]) if item["geomnav"] is not None else "None"
    print(f"{item['nuclide']:<6s} {item['E_keV']:8.3f} eff_peak={item['eff_peak']*100:.5f}% ± {item['d_eff_peak']*100:.5f}% eff_total={item['eff_total']*100:.4f}% geomnav={gn}")
