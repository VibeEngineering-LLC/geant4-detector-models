# -*- coding: utf-8 -*-
"""SA-10 stub"""
import sys
sys.stdout.reconfigure(encoding="utf-8")
import os
import re
import math
import csv
import xml.etree.ElementTree as ET
import numpy as np
from scipy.optimize import curve_fit

if not os.environ.get("GS2020_REF") or not os.environ.get("GS2020_OUT"):
    raise RuntimeError("Переменные окружения GS2020_REF и GS2020_OUT не установлены (см. README.md)")
_REF = os.environ["GS2020_REF"]
SAMPLE_XML = os.path.join(_REF, "Калибровка Th-232 (без вычета фона).xml")
BKG_XML = os.path.join(_REF, "Фон лаба S31_18.xml")
GRID_DIR = os.environ["GS2020_OUT"]
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
ENSDF_CSV = os.path.join(_REPO_ROOT, "detectors", "Gamma-1S", "web-th232", "data", "ensdf_th232_chain_lines.csv")
OUT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_TXT = os.path.join(OUT_DIR, "sa10_th232_result.txt")
BR_BI212_TO_TL208 = 0.3594
MATCH_TOL_KEV = 0.05
LINES = [
    ("Ac228", 911.20, "Ac-228"),
    ("Ac228", 968.97, "Ac-228"),
    ("Pb212", 238.63, "Pb-212"),
    ("Bi212", 727.33, "Bi-212"),
    ("Tl208", 583.19, "Tl-208"),
    ("Tl208", 2614.51, "Tl-208"),
]

def read_atomspectra_xml(path):
    tree = ET.parse(path)
    root = tree.getroot()
    es = root.find(".//EnergySpectrum")
    live_time = float(es.find("LiveTime").text)
    coeffs_node = es.find("EnergyCalibration/Coefficients")
    coeffs = [float(c.text) for c in coeffs_node.findall("Coefficient")]
    counts = [float(dp.text) for dp in es.find("Spectrum").findall("DataPoint")]
    counts = np.array(counts, dtype=float)
    return {"live_time": live_time, "coeffs": coeffs, "counts": counts, "path": path}
def channel_to_energy(ch, coeffs):
    ch = np.asarray(ch, dtype=float)
    e = np.zeros_like(ch)
    for i, c in enumerate(coeffs):
        e += c * ch ** i
    return e

def energy_to_channel_grid(coeffs, n_channels):
    ch = np.arange(n_channels)
    return channel_to_energy(ch, coeffs)

def gauss_lin(x, a, mu, sigma, b0, b1):
    return a * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + b0 + b1 * (x - mu)
def fit_peak(energies, counts, e_nominal, half_window_keV, mu_slack=0.4, sigma_max_frac=0.6):
    mask = (energies >= e_nominal - half_window_keV) & (energies <= e_nominal + half_window_keV)
    x = energies[mask]
    y = counts[mask]
    if len(x) < 8:
        return None
    bkg0 = np.median(np.concatenate([y[: max(2, len(y) // 6)], y[-max(2, len(y) // 6):]]))
    a0 = max(y.max() - bkg0, 1.0)
    mu0 = e_nominal
    sigma0 = half_window_keV / 3.0
    p0 = [a0, mu0, sigma0, bkg0, 0.0]
    mu_lo, mu_hi = e_nominal - mu_slack * half_window_keV, e_nominal + mu_slack * half_window_keV
    sigma_hi = max(sigma_max_frac * half_window_keV, 1.0)
    bounds = ([0.0, mu_lo, 0.3, 0.0, -np.inf], [np.inf, mu_hi, sigma_hi, np.inf, np.inf])
    try:
        popt, pcov = curve_fit(gauss_lin, x, y, p0=p0, bounds=bounds,
            sigma=np.sqrt(np.maximum(y, 1.0)), absolute_sigma=True, maxfev=40000)
    except Exception as exc:
        return {"error": str(exc)}
    a, mu, sigma, b0, b1 = popt
    perr = np.sqrt(np.abs(np.diag(pcov)))
    fwhm = abs(sigma) * 2.3548200450309493
    fwhm_err = perr[2] * 2.3548200450309493
    return {"mu": mu, "mu_err": perr[1], "fwhm": fwhm, "fwhm_err": fwhm_err,
            "sigma": abs(sigma), "amp": a, "b0": b0, "b1": b1, "x": x, "y": y}
def net_peak_area(energies, counts, mu, fwhm, side_frac=1.0):
    half = 1.5 * fwhm
    lo, hi = mu - half, mu + half
    peak_mask = (energies >= lo) & (energies <= hi)
    side_w = side_frac * fwhm
    left_mask = (energies >= lo - side_w) & (energies < lo)
    right_mask = (energies > hi) & (energies <= hi + side_w)
    if peak_mask.sum() < 3 or left_mask.sum() < 2 or right_mask.sum() < 2:
        return None
    gross = counts[peak_mask].sum()
    n_ch = peak_mask.sum()
    left_mean = counts[left_mask].mean()
    right_mean = counts[right_mask].mean()
    left_e_mean = energies[left_mask].mean()
    right_e_mean = energies[right_mask].mean()
    slope = (right_mean - left_mean) / (right_e_mean - left_e_mean)
    bkg_per_channel = left_mean + slope * (energies[peak_mask] - left_e_mean)
    bkg_sum = bkg_per_channel.sum()
    net = gross - bkg_sum
    var_gross = gross
    var_left_mean = counts[left_mask].sum() / (left_mask.sum() ** 2)
    var_right_mean = counts[right_mask].sum() / (right_mask.sum() ** 2)
    var_bkg_per_ch = 0.25 * var_left_mean + 0.25 * var_right_mean
    var_bkg_sum = var_bkg_per_ch * (n_ch ** 2)
    sigma_net = math.sqrt(max(var_gross, 0.0) + max(var_bkg_sum, 0.0))
    return {"net": net, "sigma": sigma_net, "gross": gross, "bkg": bkg_sum,
            "n_ch": int(n_ch), "lo": lo, "hi": hi}
def read_grid_efficiency(path, e_line):
    header_meta = {}
    rows = []
    with open(path, encoding="utf-8") as f:
        reader = csv.reader(f)
        header_found = False
        for row in reader:
            if not row:
                continue
            if not header_found:
                if row[0] == "bin_keV":
                    header_found = True
                    continue
                if len(row) >= 2:
                    header_meta[row[0]] = row[1]
                continue
            rows.append(row)
    bins = np.array([float(r[0]) for r in rows])
    counts = np.array([float(r[1]) for r in rows])
    n_events = float(header_meta.get("n_events_processed", "nan"))
    grid_energy = float(header_meta.get("energy_keV", "nan"))
    bin_center = math.floor(e_line) + 0.5
    idx = np.where(np.isclose(bins, bin_center))[0]
    if len(idx) == 0:
        return None
    i = idx[0]
    c_line = counts[i]
    c_below = counts[i - 1] if i - 1 >= 0 else 0.0
    total = c_line + c_below
    eff = total / n_events
    eff_err = math.sqrt(max(total, 1.0)) / n_events
    return {"eff": eff, "eff_err": eff_err, "n_events": n_events,
        "grid_energy": grid_energy, "c_line": c_line, "c_below": c_below,
        "bin_line": bins[i], "bin_below": bins[i - 1] if i - 1 >= 0 else None, "path": path}

def find_grid_file(grid_dir, e_line, tol=MATCH_TOL_KEV):
    best = None
    best_d = None
    for fn in os.listdir(grid_dir):
        m = re.match(r"grid_mar_E([\d.]+)\.csv$", fn)
        if not m:
            continue
        e = float(m.group(1))
        d = abs(e - e_line)
        if d <= tol and (best_d is None or d < best_d):
            best = os.path.join(grid_dir, fn)
            best_d = d
    return best, best_d
def read_ensdf_yield(path, nuclide, e_nominal, tol_kev=0.05):
    candidates = []
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["nuclide"] != nuclide:
                continue
            if row.get("line_type", "gamma") != "gamma":
                continue
            e = float(row["E_keV"])
            if abs(e - e_nominal) <= tol_kev:
                candidates.append(row)
    if not candidates:
        with open(path, encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row["nuclide"] != nuclide:
                    continue
                if row.get("line_type", "gamma") != "gamma":
                    continue
                e = float(row["E_keV"])
                if abs(e - e_nominal) <= 0.5:
                    candidates.append(row)
    if not candidates:
        return None
    best = min(candidates, key=lambda r: abs(float(r["E_keV"]) - e_nominal))
    i_percent = float(best["I_percent"])
    unc_percent = float(best.get("unc_I_percent", "0") or 0.0)
    return {"E_keV": float(best["E_keV"]), "I_percent": i_percent,
            "unc_I_percent": unc_percent, "row": best}
def main():
    log_lines = []
    def log(*args):
        s = " ".join(str(a) for a in args)
        print(s)
        log_lines.append(s)

    log("=== SA-10 независимый пересчёт: активность цепочки Th-232 (Маринелли) ===")
    log("Скрипт:", os.path.abspath(__file__))
    log()
    log("--- Чтение спектра образца ---")
    sample = read_atomspectra_xml(SAMPLE_XML)
    log("Файл:", SAMPLE_XML)
    log("Живое время образца, с:", sample["live_time"])
    log("Коэффициенты калибровки образца:", sample["coeffs"])
    n_ch_sample = len(sample["counts"])
    e_sample = energy_to_channel_grid(sample["coeffs"], n_ch_sample)
    log("Число каналов:", n_ch_sample, " диапазон энергий:", e_sample[0], "-", e_sample[-1], "кэВ")
    log()
    log("--- Чтение спектра фона ---")
    bkg = read_atomspectra_xml(BKG_XML)
    log("Файл:", BKG_XML)
    log("Живое время фона, с:", bkg["live_time"])
    log("Коэффициенты калибровки фона:", bkg["coeffs"])
    n_ch_bkg = len(bkg["counts"])
    e_bkg = energy_to_channel_grid(bkg["coeffs"], n_ch_bkg)
    log("Число каналов:", n_ch_bkg, " диапазон энергий:", e_bkg[0], "-", e_bkg[-1], "кэВ")
    log()
    calib_equal = (sample["coeffs"] == bkg["coeffs"])
    log("Калибровки образца и фона совпадают побитово:", calib_equal,
        "(проверено по факту в файлах, не предполагалось заранее — #CAL-0/LESSONS брифа)")
    log()
    scale_bkg = sample["live_time"] / bkg["live_time"]
    log("Масштаб фона по живому времени (t_sample/t_bkg):", scale_bkg)
    log()
    def bkg_counts_on_sample_grid():
        if n_ch_sample == n_ch_bkg and calib_equal:
            return bkg["counts"]
        return np.interp(e_sample, e_bkg, bkg["counts"], left=0.0, right=0.0)
    bkg_on_sample = bkg_counts_on_sample_grid()

    results = []
    assumptions = []
    not_done = []

    log("--- Предкалибровка разрешения детектора по двум изолированным сильным линиям ---")
    fit238 = fit_peak(e_sample, sample["counts"], 238.63, 45.0, mu_slack=0.3, sigma_max_frac=0.5)
    fit2614 = fit_peak(e_sample, sample["counts"], 2614.51, 90.0, mu_slack=0.3, sigma_max_frac=0.5)
    ok238 = fit238 and "error" not in fit238
    ok2614 = fit2614 and "error" not in fit2614
    if ok238 and ok2614:
        e1, f1 = fit238["mu"], fit238["fwhm"]
        e2, f2 = fit2614["mu"], fit2614["fwhm"]
        b_res = (f2 ** 2 - f1 ** 2) / (e2 - e1)
        a_res = f1 ** 2 - b_res * e1
        log("  238.63 кэВ (окно +/-45): центроид=%.2f ПШПВ=%.2f кэВ" % (e1, f1))
        log("  2614.51 кэВ (окно +/-90): центроид=%.2f ПШПВ=%.2f кэВ" % (e2, f2))
        log("  Принятый закон разрешения: ПШПВ(E)^2 = %.2f + %.5f*E, E в кэВ" % (a_res, b_res))

        def fwhm_law(E):
            v = a_res + b_res * E
            return math.sqrt(v) if v > 1.0 else 1.0
    else:
        assumptions.append("Предкалибровка разрешения по 238.63/2614.51 кэВ не удалась — использован дефолт 7.5%*sqrt(662*E)")
        log("  ПРЕДКАЛИБРОВКА НЕ УДАЛАСЬ, используется дефолт закон")

        def fwhm_law(E):
            return 0.075 * math.sqrt(661.7 * E)
    log()

    log("--- Пики: фит, площадь, эффективность, выход, активность ---")
    for nuclide, e_nom, label in LINES:
        log()
        log("* Линия", label, e_nom, "кэВ")
        grid_path, d_e = find_grid_file(GRID_DIR, e_nom, tol=MATCH_TOL_KEV)
        if grid_path is None:
            log("  НЕТ файла сетки эффективности в допуске +/-", MATCH_TOL_KEV, "кэВ — линия ПРОПУЩЕНА.")
            not_done.append(label + " " + str(e_nom) + " кэВ: нет файла grid_mar_E*.csv в допуске")
            continue
        log("  Файл сетки:", grid_path, "(|dE|=%.4f кэВ)" % d_e)
        eff_info = read_grid_efficiency(grid_path, e_nom)
        if eff_info is None:
            log("  Не найден бин сетки, покрывающий энергию", e_nom, "кэВ — линия ПРОПУЩЕНА.")
            not_done.append(label + " " + str(e_nom) + " кэВ: не найден бин в сетке")
            continue
        log("  n_events_processed=%.0f, бин линии=%s (счёт %.0f), бин ниже=%s (счёт %.0f)" % (
            eff_info["n_events"], eff_info["bin_line"], eff_info["c_line"],
            eff_info["bin_below"], eff_info["c_below"]))
        log("  Эффективность пика = %.6e +/- %.2e" % (eff_info["eff"], eff_info["eff_err"]))
        yld = read_ensdf_yield(ENSDF_CSV, nuclide, e_nom)
        if yld is None:
            log("  Не найдена линия", nuclide, e_nom, "кэВ в ENSDF CSV — линия ПРОПУЩЕНА.")
            not_done.append(label + " " + str(e_nom) + " кэВ: нет строки в ENSDF CSV")
            continue
        i_percent = yld["I_percent"]
        i_unc_percent = yld["unc_I_percent"]
        yield_per_decay = i_percent / 100.0
        yield_unc = i_unc_percent / 100.0
        branch_note = ""
        if nuclide == "Tl208":
            branch_note = " x BR(Bi212->Tl208)=%.4f" % BR_BI212_TO_TL208
            yield_per_decay = yield_per_decay * BR_BI212_TO_TL208
            yield_unc = yield_unc * BR_BI212_TO_TL208
        log("  ENSDF: E=%.3f кэВ, I=%.4f%% +/- %.4f%%%s" % (
            yld["E_keV"], i_percent, i_unc_percent, branch_note))
        log("  Выход на распад цепочки = %.6e +/- %.2e" % (yield_per_decay, yield_unc))
        fwhm_guess = fwhm_law(e_nom)
        other_e = [e2 for _, e2, _ in LINES if abs(e2 - e_nom) > 0.01]
        min_dist = min(abs(e2 - e_nom) for e2 in other_e) if other_e else 1e9
        half_window = min(max(1.6 * fwhm_guess, 15.0), 0.45 * min_dist)
        if half_window < 15.0:
            half_window = 15.0
            assumptions.append(f"{label} {e_nom} кэВ: соседняя запрошенная линия ближе, чем разрешение позволяет "
                                f"выбрать окно >=15 кэВ без риска перекрытия (ближайшая на {min_dist:.1f} кэВ)")
        fit = fit_peak(e_sample, sample["counts"], e_nom, half_window, mu_slack=0.35, sigma_max_frac=0.55)
        if fit is None or "error" in fit:
            log("  Фит пика не сошёлся — линия ПРОПУЩЕНА.", fit.get("error") if fit else "мало точек")
            not_done.append(label + " " + str(e_nom) + " кэВ: фит гауссианы не сошёлся")
            continue
        mu = fit["mu"]
        fwhm = fit["fwhm"]
        log("  Фит (окно +/-%.1f кэВ вокруг %.2f): центроид=%.3f +/- %.3f кэВ, ПШПВ=%.3f +/- %.3f кэВ" % (
            half_window, e_nom, mu, fit["mu_err"], fwhm, fit["fwhm_err"]))
        if abs(mu - e_nom) > 3.0:
            assumptions.append(label + " " + str(e_nom) + " кэВ: фитированный центроид %.2f кэВ отличается от номинала более чем на 3 кэВ" % mu)
        for r_prev in results:
            if abs(r_prev["mu"] - mu) < 0.6 * min(fwhm, r_prev["fwhm"]):
                msg = ("линии %s %.2f и %s %.2f кэВ: фиты сошлись к практически одному и тому же центроиду "
                       "(%.1f и %.1f кэВ) — детектор НЕ разрешает их по отдельности, площади дублируют один блок") % (
                    r_prev["label"], r_prev["E_keV"], label, e_nom, r_prev["mu"], mu)
                log("  ВНИМАНИЕ (слияние линий):", msg)
                assumptions.append(msg)
        area_sample = net_peak_area(e_sample, sample["counts"], mu, fwhm)
        if area_sample is None:
            log("  Недостаточно каналов для площади/подложки в образце — линия ПРОПУЩЕНА.")
            not_done.append(label + " " + str(e_nom) + " кэВ: мало каналов для площади в образце")
            continue
        area_bkg = net_peak_area(e_sample, bkg_on_sample, mu, fwhm)
        if area_bkg is None:
            assumptions.append(label + " " + str(e_nom) + " кэВ: не удалось выделить окно в фоне тем же методом — фон НЕ вычтен")
            net_after_bkg = area_sample["net"]
            sigma_after_bkg = area_sample["sigma"]
            bkg_net_scaled = 0.0
        else:
            bkg_net_scaled = area_bkg["net"] * scale_bkg
            net_after_bkg = area_sample["net"] - bkg_net_scaled
            sigma_after_bkg = math.sqrt(area_sample["sigma"] ** 2 + (area_bkg["sigma"] * scale_bkg) ** 2)
        log("  Площадь брутто-подложка (образец): %.1f +/- %.1f (окно %d каналов, %.2f-%.2f кэВ)" % (
            area_sample["net"], area_sample["sigma"], area_sample["n_ch"], area_sample["lo"], area_sample["hi"]))
        if area_bkg is not None:
            log("  Площадь брутто-подложка (фон, то же окно): %.1f +/- %.1f, после масштаба %.6f: %.1f" % (
                area_bkg["net"], area_bkg["sigma"], scale_bkg, bkg_net_scaled))
        log("  Чистая площадь (образец - фон*масштаб): %.1f +/- %.1f" % (net_after_bkg, sigma_after_bkg))
        interference = None
        if abs(e_nom - 583.19) < 0.01:
            interference = "рядом линия ~580 кэВ (Tl-208/Bi-212) — не разделена"
        if abs(e_nom - 911.20) < 0.01:
            interference = "рядом линия ~904 кэВ (Ac-228) — не разделена"
        if interference:
            log("  ВНИМАНИЕ (наложение):", interference)
            assumptions.append(label + " " + str(e_nom) + " кэВ: " + interference)
        live_time = sample["live_time"]
        eff = eff_info["eff"]
        eff_err = eff_info["eff_err"]
        activity = net_after_bkg / (live_time * eff * yield_per_decay)
        rel_area = sigma_after_bkg / net_after_bkg if net_after_bkg != 0 else float("inf")
        rel_eff = eff_err / eff if eff != 0 else 0.0
        rel_yield = yield_unc / yield_per_decay if yield_per_decay != 0 else 0.0
        rel_tot = math.sqrt(rel_area ** 2 + rel_eff ** 2 + rel_yield ** 2)
        activity_sigma = abs(activity) * rel_tot
        log("  Активность = чистая_площадь / (t_live * eff * yield) = %.3f +/- %.3f Бк" % (activity, activity_sigma))
        log("    (rel: площадь %.1f%%, эфф-ть %.1f%%, выход %.1f%%)" % (rel_area*100, rel_eff*100, rel_yield*100))
        results.append({"label": label, "E_keV": e_nom, "mu": mu, "fwhm": fwhm,
            "net_area": net_after_bkg, "net_area_sigma": sigma_after_bkg,
            "eff": eff, "eff_err": eff_err, "yield": yield_per_decay, "yield_err": yield_unc,
            "activity": activity, "activity_sigma": activity_sigma})

    log()
    log("=== Итоговая таблица ===")
    header = "%-12s%13s%10s%18s%12s%10s%16s" % (
        "линия", "центроид,кэВ", "ПШПВ,кэВ", "площадь+-s", "эфф-ть", "выход", "A+-s,Бк")
    log(header)
    for r in results:
        area_str = "%.0f+-%.0f" % (r["net_area"], r["net_area_sigma"])
        act_str = "%.2f+-%.2f" % (r["activity"], r["activity_sigma"])
        log("%-12s%13.2f%10.2f%18s%12.4e%10.4e%16s" % (
            r["label"], r["mu"], r["fwhm"], area_str, r["eff"], r["yield"], act_str))
    log()
    if results:
        weights = np.array([1.0 / (r["activity_sigma"] ** 2) for r in results if r["activity_sigma"] > 0])
        vals = np.array([r["activity"] for r in results if r["activity_sigma"] > 0])
        if len(weights) > 0 and weights.sum() > 0:
            a_weighted = float(np.sum(weights * vals) / np.sum(weights))
            a_weighted_sigma = float(1.0 / math.sqrt(np.sum(weights)))
            log("Взвешенное среднее по всем линиям: %.3f +/- %.3f Бк (N=%d линий, веса 1/s^2)" % (
                a_weighted, a_weighted_sigma, len(vals)))
            if len(vals) > 1:
                scatter = float(np.std(vals, ddof=1))
                log("Стандартное отклонение активностей между линиями (внешняя оценка): %.3f Бк" % scatter)
        else:
            log("Взвешенное среднее не посчитано — нет линий с определённой σ.")
    else:
        log("Взвешенное среднее не посчитано — нет ни одной успешной линии.")
    log()
    log("=== Допущения, которые могли исказить число ===")
    default_assumptions = [
        "ПШПВ на каждой линии получена подгонкой гауссианы+линейного фона к самому спектру образца "
        "(начальное приближение окна из эмпирического закона ~8.5% при 662 кэВ, масштаб sqrt(E); "
        "итоговый ПШПВ — из фита, т.к. в файле нет таблицы пиков прибора).",
        "Подложка под пиком — линейная по средним значениям боковых окон шириной 1 ПШПВ; "
        "реальная форма континуума (комптоновские ступени соседних линий) может отличаться от прямой.",
        "Фон вычитается с фиксированным масштабом по отношению живых времён (без подгонки формы/интенсивности); "
        "статистика фона на слабых линиях и на 2614.51 кэВ может быть низкой.",
    ]
    default_assumptions.append(
        "Эффективность взята как сумма count_edep в бине линии и бине ниже, делённая на n_events_processed, "
        "БЕЗ размытия аппаратным разрешением (сетка Geant4 — выделенная энергия шагом 1 кэВ без свёртки с ПШПВ). "
        "Это может НЕ учитывать перераспределение отсчётов размытием между окном интегрирования измеренного пика "
        "(+-1.5 ПШПВ, десятки кэВ) и двумя узкими бинами эффективности (2 кэВ) — систематика неизвестного знака "
        "и величины, не измерена в этом прогоне.")
    default_assumptions.append(
        "Выход линии на распад цепочки взят из ENSDF CSV проекта без поправки на суммирование совпадений "
        "(true coincidence summing) в геометрии Маринелли для каскадных переходов Ac-228/Tl-208.")
    default_assumptions.append(
        "Предполагается вековое равновесие всей цепочки Th-232 (одна активность на все линии); "
        "расхождение между линиями по нуклидам, если системное, может означать нарушение равновесия "
        "(например, эманирование Rn-220), а не только случайную ошибку.")
    default_assumptions.append(
        "КРИТИЧНО для Ac-228 911.20 и 968.97 кэВ: закон разрешения (по двум опорным точкам) предсказывает "
        "ПШПВ~90-110 кэВ в этой области, но окно фита принудительно сужено до +/-26 кэВ, чтобы линии 911.20/968.97 "
        "не слились в одну — это ЗАНИЖАЕТ их площадь (интегрируется только часть истинного пика) и делает эти две "
        "активности МЕНЕЕ надёжными, чем Pb-212(238.63)/Tl-208(583.19, 2614.51), которые дали окна, согласованные "
        "с законом разрешения.")
    default_assumptions.append(
        "Закон разрешения ПШПВ(E)^2=a+b*E определён всего по ДВУМ опорным точкам (238.63 и 2614.51 кэВ) — "
        "экстраполяция/интерполяция в область 580-970 кэВ не проверена независимо и может быть неточной.")
    for a in default_assumptions:
        log("-", a)
    for a in assumptions:
        log("-", a)

    log()
    log("=== Что НЕ удалось сделать ===")
    if not_done:
        for n in not_done:
            log("-", n)
    else:
        log("- Все запрошенные линии обработаны (файлы сеток эффективности найдены в допуске +/-%s кэВ)." % MATCH_TOL_KEV)

    log()
    log("Файл вывода:", OUT_TXT)
    with open(OUT_TXT, "w", encoding="utf-8") as f:
        f.write("\n".join(log_lines) + "\n")

if __name__ == "__main__":
    main()

