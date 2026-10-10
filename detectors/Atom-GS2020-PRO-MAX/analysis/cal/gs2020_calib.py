import sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sa10_th232 as s

REFS = {
    "sample": [(238.632, 30), (338.320, 30), (583.187, 45), (727.330, 45), (911.204, 55), (1460.822, 80), (2614.511, 150)],   # K-40 (фон во время замера) — единственный репер между 911 и 2614; Ac-228 1459/1496 слабее на порядок
    "bg":     [(238.632, 30), (351.932, 25), (609.312, 40), (1460.822, 80), (2614.511, 150)],   # Bi-214 1120/1764 слиты с соседями Bi-214 — не реперы
}
# bgw (Маринелли+вода, промежуточный замер 11,7 ч) — СВОЯ калибровка по пикам невозможна (W-147, оператор 26.09
# «911, 1461 не на месте»): за 11,7 ч в спектре нет разрешимых пиков вообще, фит по REFS подгонял бы шум/континуум
# (проверено — добавление репера 911 дало невязку +15,7 кэВ и испортило соседний 609). Оператор (26.09, п.1) выбрал:
# наследовать калибровку от "bg" (7 сут, тот же прибор без сосуда) — заводские file_coeffs у bg и bgw СОВПАДАЮТ
# (проверка в calibrate_bgw_from_bg), т.е. шкала прибора между замерами не менялась.
def calibrate_bgw_from_bg():
    d_bg = s.read_atomspectra_xml(PATHS["bg"]); d_bgw = s.read_atomspectra_xml(PATHS["bgw"])
    if list(d_bg["coeffs"]) != list(d_bgw["coeffs"]):
        raise SystemExit("ОТКАЗ: file_coeffs bg %s != bgw %s — наследование калибровки запрещено, нужна своя (оператор 26.09 п.1 предполагал совпадение)"
                          % (d_bg["coeffs"], d_bgw["coeffs"]))
    r = dict(calibrate("bg")); r["tag"] = "bgw"; r["inherited_from"] = "bg"
    # #GS-23 (оператор 28.09 «надо откалибровать сдвиг»): замер воды 33,2 ч уплыл относительно bg в мягкой зоне
    # (238,6: −4,1±1,0 канала; 1460: −0,5±0,6) — сдвиг bgw→bg по общим реперам, за крайними — постоянный
    e_file = s.channel_to_energy(np.arange(len(d_bgw["counts"]), dtype=float), d_bgw["coeffs"])
    sh = []
    for E0, hw in REFS["bg"]:
        fb = s.fit_peak(e_file, d_bg["counts"], E0, hw); fw = s.fit_peak(e_file, d_bgw["counts"], E0, hw)
        if not fb or not fw or "mu" not in fb or "mu" not in fw:
            raise SystemExit("ОТКАЗ: bgw/bg, репер %.3f кэВ — фит не сошёлся" % E0)
        sh.append([float(fw["mu"]), float(fb["mu"] - fw["mu"]), float(math.hypot(fb["mu_err"], fw["mu_err"]))])
    r["shift_nodes"] = sorted(sh)
    return r
# Фон с водой: путь из GS2020_BG_WATER_XML, по умолчанию — файл в референсах (оператор 26.09: «образец с торием считать с фоном с водой»)
BKG_WATER_XML = os.environ.get("GS2020_BG_WATER_XML", os.path.join(os.path.dirname(s.BKG_XML), "Фон Маринелли 1 л вода дист (33,2 ч, 25-26.09).xml"))   # оператор 27.09 «обнови фон» — то же измерение, длиннее (33,2ч, конец 26.09 23:34); прежние — 22ч, 18,2 ч (26.09) и 11,7 ч (24.09)
PATHS = {"sample": s.SAMPLE_XML, "bg": s.BKG_XML, "bgw": BKG_WATER_XML,
         # #GS-47 (оператор 07.10: «KCl 1 литр готов. 1085 г. Прими и замени спектр для калия»); прежний — GS2020_KCL_XML=<путь>
         # «KCl ч 740 мл 829 г (9,3 ч, 27-28.09).xml». В новом файле шкала ЕСТЬ (как у S31_18) — #CAL-0 сверить с реперами.
         # 09.10 (оператор «да, пересчитать»): рабочий — накопление того же замера 88,3 ч; 41,5 ч — GS2020_KCL_XML=<путь>
         # #GS-62 (оператор 09.10 «Сумма с выравниванием»): сумма 88,3 ч без выравнивания смещена (усиление части 2 +0,2 %) —
         # рабочий спектр = часть 1 + часть 2, перебинованная на шкалу части 1 (scripts\gs2020_kcl_align.py, audit\GS-62-kcl-drift.md)
         "kcl": os.environ.get("GS2020_KCL_XML", r"<WORKDIR>\GEANT4\results\kcl_parts\KCl ч 1 л 1085 г сумма выровненная (88,3 ч).xml")}
# #GS-24: в файле KCl калибровки НЕТ (ни в xml, ни в spe) — начальная шкала от файла bgw (тот же прибор), дальше — свои
# реперы спектра KCl (583/911/1120/1764 слиты или в шуме за 9,3 ч — не реперы, проверено 28.09)
REFS["kcl"] = [(238.632, 30), (609.312, 40), (1460.822, 80), (2614.511, 150)]
# #GS-23, ред. 2 (28.09): замер воды 33,2 ч калибруется ПО СВОИМ реперам (все 5 в пороге, rms 1,56 кэВ) — проще
# наследования от bg + сдвига (calibrate_bgw_from_bg) и равнозначно ему: центры пиков по обеим шкалам совпали до 0,01 кэВ
# (W-154: «двойной сдвиг» был ложной диагностикой — сравнивались центроиды разных методов). W-147 (11,7 ч) к 33,2 ч не относится.
REFS["bgw"] = REFS["bg"]
# #GS-70 (оператор 10.10: «Черника с цезием и стронцием. Маринелли, 1 литр 482 г»): шкала в файле есть (#CAL-0 — scripts\gs70_cal_check.py)
PATHS["berry"] = os.environ.get("GS2020_BERRY_XML", os.path.join(os.path.dirname(s.BKG_XML), "Черника Cs+Sr Маринелли 1 л 482 г (20,2 ч, 09-10.10).xml"))
REFS["berry"] = [(661.657, 45), (1460.822, 80), (2614.511, 150)]
def read(tag):
    if tag != "kcl": return s.read_atomspectra_xml(PATHS[tag])
    import xml.etree.ElementTree as ET
    es = ET.parse(PATHS["kcl"]).getroot().find(".//EnergySpectrum")
    return {"live_time": float(es.find("LiveTime").text), "coeffs": s.read_atomspectra_xml(PATHS["bgw"])["coeffs"],
            "counts": np.array([float(x.text) for x in es.find("Spectrum").findall("DataPoint")]), "path": PATHS["kcl"]}
REPR_ORDER = 4
OUT_JSON = r"C:\g4work\gs2020\run_marinelli\out_v5\cal_own.json"
# #GS-28 совместная калибровка (ЛСРМ): GS_CAL_MPLET=1 — узлы из подгонки ФОРМЫ групп на шкале файла
# (gs2020_calib_mplet.py при GS28_SCALE=file → cal_mplet_file.json) вместо одиночных гауссов реперов (609 фона — смесь
# с Tl-208 583, смещение ~8 кэВ; 911 не был узлом). Результат — в cal_own_mplet.json; подгонки берут его при
# GS_CAL_OWN_JSON=<путь>. Узел годен: не СЛАБО, не КРАЙ, χ²/ν ≤ 3, σ_δ ≤ 5 % ПШПВ.
MPLET = os.environ.get("GS_CAL_MPLET") == "1"
if MPLET: OUT_JSON = OUT_JSON.replace("cal_own.json", "cal_own_mplet.json")
# GS_CAL_HYBRID=1 (оператор 28.09 «да»): у bg/bgw репер 609 (смесь с Tl-208 583) заменён узлом формы группы 583+609,
# добавлен узел формы 911 (у фонов его не было); остальные реперы — одиночные гауссы, как прежде → cal_own_hybrid.json
HYBRID = os.environ.get("GS_CAL_HYBRID") == "1"
if HYBRID: OUT_JSON = OUT_JSON.replace("cal_own.json", "cal_own_hybrid.json")
OUT_JSON = os.environ.get("GS_CAL_OWN_JSON", OUT_JSON)
# #GS-37: GS_CAL_SHAPE=1 — шкала из калибровки по ФОРМЕ (gs2020_calib_shape.py → cal_shape.json): d(E_файла) кусочно-
# линейная по узлам значимых групп спектра (либо полином), отклики Geant4, совместно по всем группам
SHAPE_JSON = r"C:\g4work\gs2020\run_marinelli\out_v5\cal_shape.json"
SHAPE = os.environ.get("GS_CAL_SHAPE") == "1"
if SHAPE: OUT_JSON = SHAPE_JSON
def mplet_nodes(tag, groups=None):
    with open(os.path.join(os.path.dirname(OUT_JSON), "cal_mplet_file.json"), encoding="utf-8") as f: g = json.load(f)[tag]
    return [(v["E_node"], v["delta_keV"], v["sigma_keV"]) for k, v in g.items() if groups is None or k in groups
            if not v["weak"] and not v["edge"] and v["chi2_nu"] <= 3 and "E_node" in v
            and v["sigma_keV"] <= 0.05 * 0.6061 * v["E_node"] ** 0.669]   # σ_δ ≤ 5 % ПШПВ (порог невязки #CAL-0 — 25 %)

def calibrate(tag):
    d = read(tag)
    ch = np.arange(len(d["counts"]), dtype=float)
    e_file = s.channel_to_energy(ch, d["coeffs"])
    mu_files, dEs, sigs, fwhms, E0s = [], [], [], [], []
    hyb = mplet_nodes(tag, ("609", "911")) if HYBRID and tag in ("bg", "bgw") else []
    if hyb: print("   %s: гибрид — узлы формы %s" % (tag, ", ".join("%.1f (δ %+.2f ± %.2f)" % n for n in hyb)))
    for E0, dl, sd in (mplet_nodes(tag) if MPLET else hyb):
        mu_files.append(E0 + dl); sigs.append(sd); fwhms.append(0.6061 * E0 ** 0.669); E0s.append(E0)
    for E0, hw in ([] if MPLET else [r for r in REFS[tag] if not (hyb and abs(r[0] - 609.312) < 1)]):
        f = s.fit_peak(e_file, d["counts"], E0, hw)
        if f is None or "error" in f or "mu" not in f:
            raise SystemExit("ОТКАЗ: %s, репер %.3f кэВ — фит не сошёлся" % (tag, E0))
        mu_files.append(f["mu"]); sigs.append(max(f["mu_err"], 0.05 * f["fwhm"])); fwhms.append(f["fwhm"]); E0s.append(E0)   # пол 5 % ПШПВ: смеси в NaI сдвигают центроид сильнее статистики
    # Поправка dE(E_файла) — интерполяция через СВОИ реперы спектра, за крайними — постоянная. Полиномы 2-й и 3-й степени
    # отвергнуты 25.09 (K-40 образца −8 кэВ; 3-я расходится на 74 кэВ за 2614), см. W-140.
    xs = np.array(mu_files); ys = np.array(E0s) - xs; o = np.argsort(xs); xs, ys = xs[o], ys[o]
    p = [[float(a), float(b)] for a, b in zip(xs, ys)]
    e_own = e_file + np.interp(e_file, xs, ys)
    c_repr = np.polynomial.polynomial.polyfit(ch, e_own, REPR_ORDER)
    repr_dev = float(np.max(np.abs(np.polynomial.polynomial.polyval(ch, c_repr) - e_own)))
    refs_out = []
    for i, E0 in enumerate(E0s):
        mu_file = mu_files[i]; fwhm = fwhms[i]
        m = np.abs(xs - mu_file) > 1e-9      # проверка с исключением репера: шкала без него предсказывает его энергию
        E_own = mu_file + np.interp(mu_file, xs[m], ys[m]); resid = E_own - E0; limit = 0.25 * fwhm; ok = abs(resid) <= limit
        refs_out.append({"E_lib": float(E0), "mu_file": float(mu_file), "mu_err": float(sigs[i]), "fwhm": float(fwhm),
                         "E_own": float(E_own), "resid_keV": float(resid), "limit_keV": float(limit), "ok": bool(ok)})
    rms = float(math.sqrt(np.mean([r["resid_keV"]**2 for r in refs_out])))
    return {"tag": tag, "file_coeffs": list(d["coeffs"]), "corr_nodes": p, "coeffs": list(c_repr),
            "repr_max_dev_keV": repr_dev, "refs": refs_out, "rms_keV": rms}

def energy_axis(tag, n_channels=None):
    if not os.path.exists(OUT_JSON): raise SystemExit("ОТКАЗ: нет %s — сначала python gs2020_calib.py" % OUT_JSON)
    with open(OUT_JSON, encoding="utf-8") as f: data = json.load(f)
    r = data[tag]
    n = n_channels if n_channels is not None else len(read(tag)["counts"])
    e_file = s.channel_to_energy(np.arange(n), r["file_coeffs"])
    if SHAPE:   # тот же базис, что gs2020_calib_shape.basis: узлы — interp с постоянной за краями, иначе полином по u
        th = np.asarray(r["theta"], float)
        if r.get("knots"): return e_file + np.interp(e_file, r["knots"], th)
        u = (e_file - 1000) / 1000.0
        return e_file + sum(t * u ** j for j, t in enumerate(th))
    nodes = np.array(r["corr_nodes"])
    if r.get("shift_nodes"):   # #GS-23: bgw приводится к шкале файла bg до поправки по реперам bg
        sn = np.array(r["shift_nodes"]); e_file = e_file + np.interp(e_file, sn[:, 0], sn[:, 1])
    return e_file + np.interp(e_file, nodes[:, 0], nodes[:, 1])

if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    if SHAPE: raise SystemExit("ОТКАЗ: GS_CAL_SHAPE=1 — cal_shape.json пишет gs2020_calib_shape.py, не этот скрипт")
    res = {tag: calibrate(tag) for tag in ("sample", "bg")}
    res["bgw"] = calibrate("bgw")   # W-154: своя шкала, не calibrate_bgw_from_bg()
    res["kcl"] = calibrate("kcl")
    any_fail = False
    for tag, r in res.items():
        print("== %s: коэффициенты файла %s" % (tag, ", ".join("%.6g" % c for c in r["file_coeffs"])))
        print("   свои E(ch), степень 4: %s; представление ≤ %.3f кэВ" % (", ".join("%.6g" % c for c in r["coeffs"]), r["repr_max_dev_keV"]))
        for ref in r["refs"]:
            status = "ok" if ref["ok"] else "ПРЕВЫШЕН"
            print("%9.3f  μ_файл %9.3f ± %.2f  ПШПВ %6.1f  свой %9.3f  невязка без него %+6.2f кэВ  (порог ±%.1f) %s" %
                  (ref["E_lib"], ref["mu_file"], ref["mu_err"], ref["fwhm"], ref["E_own"], ref["resid_keV"], ref["limit_keV"], status))
            if not ref["ok"]: any_fail = True
        print("   rms %.2f кэВ" % r["rms_keV"])
    with open(OUT_JSON, "w", encoding="utf-8") as f: json.dump(res, f, ensure_ascii=False, indent=1)
    e_s = energy_axis("sample"); ch_s = np.arange(len(e_s)); e_b = energy_axis("bg")
    for E in [100, 238.632, 1460.822, 2614.511, 3500]:
        c = np.interp(E, e_s, ch_s); Es = np.interp(c, ch_s, e_s); Eb = np.interp(c, np.arange(len(e_b)), e_b)
        print("канал %.1f: образец %.2f кэВ, фон %.2f кэВ, разница %+.2f" % (c, Es, Eb, Es - Eb))
    if any_fail: print("ТРЕБУЕТ ТОЛКОВАНИЯ: невязки выше 0,25·ПШПВ"); sys.exit(1)
