import sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sa10_th232 as s

REFS = {
    "sample": [(238.632, 30), (338.320, 30), (583.187, 45), (727.330, 45), (911.204, 55), (1460.822, 80), (2614.511, 150)],   # K-40 (фон во время замера) — единственный репер между 911 и 2614; Ac-228 1459/1496 слабее на порядок
    "bg":     [(238.632, 30), (351.932, 25), (609.312, 40), (1460.822, 80), (2614.511, 150)],   # Bi-214 1120/1764 слиты с соседями Bi-214 — не реперы
}
# bgw (Маринелли+вода, промежуточный замер 11,7 ч) — СВОЯ калибровка по пикам невозможна (W-147, оператор 26.09
# «911, 1461 не на месте»): за 11,7 ч в спектре нет разрешимых пиков вообще, фит по REFS подгонял бы шум/континуум.
# Оператор (26.09, п.1) выбрал: наследовать калибровку от "bg" (7 сут, тот же прибор без сосуда) — заводские
# file_coeffs у bg и bgw СОВПАДАЮТ (проверка в calibrate_bgw_from_bg), т.е. шкала прибора между замерами не менялась.
def calibrate_bgw_from_bg():
    d_bg = s.read_atomspectra_xml(PATHS["bg"]); d_bgw = s.read_atomspectra_xml(PATHS["bgw"])
    if list(d_bg["coeffs"]) != list(d_bgw["coeffs"]):
        raise SystemExit("ОТКАЗ: file_coeffs bg %s != bgw %s — наследование калибровки запрещено, нужна своя"
                          % (d_bg["coeffs"], d_bgw["coeffs"]))
    r = dict(calibrate("bg")); r["tag"] = "bgw"; r["inherited_from"] = "bg"
    return r
BKG_WATER_XML = os.environ.get("GS2020_BG_WATER_XML", os.path.join(os.path.dirname(s.BKG_XML), "Фон Маринелли 1 л вода дист (промеж 11,7 ч, 24.09).xml"))
PATHS = {"sample": s.SAMPLE_XML, "bg": s.BKG_XML, "bgw": BKG_WATER_XML}
REPR_ORDER = 4
OUT_JSON = os.path.join(os.environ.get("GS2020_OUT", r"C:\g4work\gs2020\run_marinelli\out_v5"), "cal_own.json")

def calibrate(tag):
    d = s.read_atomspectra_xml(PATHS[tag])
    ch = np.arange(len(d["counts"]), dtype=float)
    e_file = s.channel_to_energy(ch, d["coeffs"])
    mu_files, dEs, sigs, fwhms, E0s = [], [], [], [], []
    for E0, hw in REFS[tag]:
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
    r = data[tag]; nodes = np.array(r["corr_nodes"])
    n = n_channels if n_channels is not None else len(s.read_atomspectra_xml(PATHS[tag])["counts"])
    e_file = s.channel_to_energy(np.arange(n), r["file_coeffs"])
    return e_file + np.interp(e_file, nodes[:, 0], nodes[:, 1])

if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    res = {tag: calibrate(tag) for tag in ("sample", "bg")}
    res["bgw"] = calibrate_bgw_from_bg()
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
