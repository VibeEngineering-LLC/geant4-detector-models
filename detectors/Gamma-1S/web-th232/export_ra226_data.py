# -*- coding: utf-8 -*-
"""Экспорт данных для лёгкой страницы Ra-226 -- ТОЛЬКО метод 2 (по прямому
указанию оператора 09.08.2026), в двух вариантах библиотеки (отобранная
I>=2% / полная цепочка), сумм-пики те же в обоих (физика каскада не
зависит от того, какие ОДИНОЧНЫЕ линии показаны отдельно).

Метод 1 здесь нет: iso_*.csv по звеньям Ra-226 не прогнаны (задача #182,
следующий шаг). Схема JSON -- СВОЯ, короче g1s_th232_data.json (нет
метода 1, нет варианта "cs" по отдельной калибровке цезия, нет масок
достоверности МК-статистики) -- под лёгкий фронтенд ra226.js, не
g1s-th232.js.

Запуск:
    python export_ra226_data.py
Переменные окружения: G4MODELS_BUILD_GAMMA_1S (сетка отклика, та же, что
и у Th-232 -- тот же сосуд/матрица/детектор), G4MODELS_RA226_BG_SPE
(приватный фон, в репозитории его нет).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault("G4MODELS_SOURCE_CONFIG",
                      os.path.join(HERE, "configs", "ra226.yaml"))
os.environ.setdefault("PYTHONIOENCODING", "utf-8")

sys.path.insert(0, HERE)
import export_data as ed  # noqa: E402

SPECTRAVIBE_ROOT = (r"C:\Users\Дмитрий\Мой диск\Дозиметрия\ИИ\1 Скилы"
                    r"\0_Work\gamma-spectrum-analysis")
sys.path.insert(0, os.path.join(SPECTRAVIBE_ROOT, "scripts"))
from gamma.io.lsrm_spe import read_lsrm_spe  # noqa: E402


def read_pair():
    cfg = ed._CFG
    sample_rel = cfg["source"]["measured_sample_spe_rel"]
    sample_path = os.path.join(str(ed.KIT), *sample_rel.split("/"))
    bg_env = cfg["source"]["measured_background_spe_env"]
    bg_path = os.environ.get(bg_env)
    if not bg_path:
        raise SystemExit("Задайте %s -- путь к фоновому .spe." % bg_env)
    s = read_lsrm_spe(sample_path)
    b = read_lsrm_spe(bg_path)
    ch_s = np.arange(s.n_channels, dtype=float)
    ch_b = np.arange(b.n_channels, dtype=float)
    return ({"counts": np.asarray(s.counts, dtype=np.int64),
             "e_of_ch": np.asarray(s.channel_to_energy(ch_s), dtype=float),
             "live_s": float(s.live_time), "real_s": float(s.real_time),
             "start": str(s.start_datetime)},
            {"counts": np.asarray(b.counts, dtype=np.int64),
             "e_of_ch": np.asarray(b.channel_to_energy(ch_b), dtype=float),
             "live_s": float(b.live_time), "real_s": float(b.real_time)})


def run_method2(library, sums, resp, e, ch_edges, keys):
    """Урезанный, но физически тот же run_method2, что в export_data.py
    (F_B-депопуляция, F_B-нормировка сумм-пиков) -- без канальной
    раскладки и без диагностики peak_area_with_shelf (не нужны лёгкой
    странице)."""
    shape_total = np.zeros_like(e)
    by_nuc_w = {k: np.zeros_like(e) for k in keys}

    def add(nuc_key, weight, shp):
        shape_total[:] += weight * shp
        by_nuc_w[nuc_key] += weight * shp

    # ИСПРАВЛЕНО 09.08.2026 (аудит Б2, коммит df5d178 -- та же находка, что и
    # в export_data.py.run_method2, здесь отдельная, НЕ синхронизированная
    # копия): эффективность ПАРТНЁРА каскада в депопуляции должна быть
    # ПОЛНОЙ (eps_total = shape.sum(), вероятность зарегистрировать хоть
    # что-то от кванта где угодно в спектре), не пиковой -- см. подробное
    # обоснование и цитату (Chehade 2007, IUP Bremen, ур. 2.2-2.3) в
    # export_data.py.run_method2. У сумм-пика ниже (строки 101-105) обе
    # эффективности остаются пиковыми -- там другая физика (полное
    # поглощение ОБОИХ квантов), Б2 её не касается.
    depl = {}
    for E1s, E2s, nuc_keys, I1s, I2s, _note_s, fb_pct_s in sums:
        shp1s, _, _ = resp(E1s)
        shp2s, _, _ = resp(E2s)
        eps1s_tot = float(shp1s.sum())
        eps2s_tot = float(shp2s.sum())
        fb_frac_s = fb_pct_s / 100.0
        k1 = (nuc_keys, round(E1s, 3))
        k2 = (nuc_keys, round(E2s, 3))
        depl[k1] = depl.get(k1, 0.0) + (I1s / 100.0) * (I2s / 100.0) * eps2s_tot / fb_frac_s
        depl[k2] = depl.get(k2, 0.0) + (I2s / 100.0) * (I1s / 100.0) * eps1s_tot / fb_frac_s

    lines_out = []
    for E, I_pct, nuc_key, note in library:
        shp, chans, eps = resp(E)
        w = I_pct / 100.0
        w_depl = depl.get((nuc_key, round(E, 3)), 0.0)
        depl_pct = 0.0
        if w_depl > 0:
            depl_pct = 100.0 * w_depl / max(w, 1e-30)
            w = max(0.0, w - w_depl)
        add(nuc_key, w, shp)
        lines_out.append({"E_keV": E, "nuclide": nuc_key, "I_pct": I_pct,
                          "note": note, "kind": "line",
                          "depleted_pct": depl_pct})

    n_sum_used = 0
    for E1, E2, nuc_key, I1_pct, I2_pct, note, fb_pct in sums:
        Esum = E1 + E2
        if Esum > ed.E_FIT_HI:
            continue
        _, _, eps1 = resp(E1)
        _, _, eps2 = resp(E2)
        shp, chans, eps_sum_node = resp(Esum)
        w = ((I1_pct / 100.0) * (I2_pct / 100.0) * eps1 * eps2
             / max(eps_sum_node, 1e-30) / (fb_pct / 100.0))
        add(nuc_key, w, shp)
        lines_out.append({"E_keV": Esum, "nuclide": nuc_key, "I_pct": None,
                          "note": note, "kind": "sum",
                          "E1_keV": E1, "E2_keV": E2})
        n_sum_used += 1

    return shape_total, by_nuc_w, lines_out, n_sum_used


def main():
    meas, bg = read_pair()
    e = meas["e_of_ch"]
    T = meas["live_s"]
    bg_on_meas = np.interp(e, bg["e_of_ch"], bg["counts"].astype(float),
                           left=0.0, right=0.0)
    bg_scale_time = T / bg["live_s"]
    bg_scaled = bg_on_meas * bg_scale_time

    ch_edges = np.concatenate((
        [e[0] - 0.5 * (e[1] - e[0])],
        0.5 * (e[:-1] + e[1:]),
        [e[-1] + 0.5 * (e[-1] - e[-2])],
    ))

    fwhm_cal = ed.fit_fwhm_calibration(meas["counts"], e)
    ed.FWHM_LAW.update({"kind": "power", "k": fwhm_cal["k"], "p": fwhm_cal["p"]})

    eps_peak = ed.make_eps_peak_interp(os.path.join(ed.BUILD, "grid"))
    resp = ed.make_full_response(os.path.join(ed.BUILD, "grid"), ch_edges,
                                 True, eps_peak)

    sel = (e >= ed.E_FIT_LO) & (e <= ed.E_FIT_HI)
    y_sel = meas["counts"][sel].astype(float)
    bgm = bg_scaled[sel]

    keys = [n[0] for n in ed.NUCS]
    lib_full, _ = ed.load_full_library(nuc_keys=set(keys))

    variants = {}
    for tag, library in (("sel", ed.GAMMA_LIBRARY), ("full", lib_full)):
        shape_total, by_nuc_w, lines_out, n_sum = run_method2(
            library, ed.SUM_PEAKS, resp, e, ch_edges, keys)
        coef, dcoef, chi2, ndof, _ = ed.fit_amplitudes(
            y_sel, [shape_total[sel] * T, bgm], ed.SYS_FLOOR)
        A_Bq, dA_Bq, bg_amp = float(coef[0]), float(dcoef[0]), float(coef[1])
        for ln in lines_out:
            ln["predicted_net"] = ln.get("I_pct")  # placeholder, filled below if needed
        p = ed._CFG["passport"]
        decay_f = ed.decay_factor_years(p["half_life_years"], p["days_pass_to_meas"])
        A_pass = p["bq_per_kg"] * (p["mass_g"] / 1000.0) * decay_f
        dA_pass = A_pass * p["unc_pct"] / 100.0
        variants[tag] = {
            "A_Bq": A_Bq, "dA_Bq": dA_Bq, "bg_amplitude": bg_amp,
            "chi2": chi2, "ndof": ndof, "chi2_ndof": chi2 / ndof,
            "n_lines": len(library), "n_sum_peaks": n_sum,
            "n_sum_peaks_total": len(ed.SUM_PEAKS),
            "ratio_to_passport": A_Bq / A_pass, "d_ratio": dA_Bq / A_pass,
            "lines": lines_out,
            "stack": {k: (by_nuc_w[k] * A_Bq * T).tolist() for k in keys},
        }

    p = ed._CFG["passport"]
    decay_f = ed.decay_factor_years(p["half_life_years"], p["days_pass_to_meas"])
    A_pass = p["bq_per_kg"] * (p["mass_g"] / 1000.0) * decay_f
    dA_pass = A_pass * p["unc_pct"] / 100.0

    reference_lines = [[ln["E_keV"], ln["nuclide"]]
                       for ln in variants["sel"]["lines"] if ln["kind"] == "line"]

    # ── проверка на утечку радона (замечание оператора №3, 09.08.2026) ──
    # Rn-222 (T1/2=3,82 сут) стоит МЕЖДУ Ra-226 и Pb-214/Bi-214 в цепочке.
    # Метод 2 выше считает ОДНУ амплитуду на всю цепочку (вековое
    # равновесие) и утечку в принципе не увидит. ПОПЫТКА независимой
    # проверки: площадь линии-родителя 186,211 кэВ (ДО радона) против
    # чистой линии дочернего 351,932 кэВ (ПОСЛЕ). НЕ ПОЛУЧИЛОСЬ:
    # ПШПВ на 186 кэВ (~22,6 кэВ) настолько широка, что "плечо" окна
    # фона по методу peak_area_with_shelf дотягивается до соседней
    # СИЛЬНОЙ линии 241,995 кэВ Pb-214 (всего 56 кэВ между линиями) --
    # оценка фона получается ЗАВЫШЕННОЙ настолько, что чистая площадь
    # выходит ОТРИЦАТЕЛЬНОЙ при любом разумном наборе окон (roi/shelf
    # от 0,5 до 1,0 x ПШПВ, проверено). Это не признак утечки, это
    # метод, сломанный контаминацией соседним пиком -- число НЕ
    # публикуется как результат (запрет выдумки: сломанное измерение
    # остаётся сломанным измерением, не подгоняется под ответ).
    #
    # Честный вывод: метод 2 в текущем виде утечку радона проверить не
    # может (одна амплитуда на цепочку) и наивная альтернатива тоже не
    # работает на этой паре линий. Нужен метод 1 (независимые МК-
    # амплитуды по нуклидам, ещё не прогнан) либо деконволюция
    # мультиплета 186/242 вместо линейного плеча.
    radon_check = {
        "attempted": True, "reliable": False,
        "reason": "ПШПВ на 186,211 кэВ (~22,6 кэВ) настолько широка, что "
                  "окно фона peak_area_with_shelf дотягивается до соседней "
                  "сильной линии 241,995 кэВ Pb-214 (56 кэВ между линиями) "
                  "-- оценка фона завышена, чистая площадь отрицательна при "
                  "любом разумном наборе окон (roi/shelf 0,5-1,0 x ПШПВ). "
                  "Метод 2 в текущем виде (одна амплитуда на всю цепочку) "
                  "утечку в принципе не видит. Нужен метод 1 (независимые "
                  "МК-амплитуды, ещё не прогнан) или деконволюция "
                  "мультиплета 186/242 вместо окна с линейным плечом.",
    }
    print("проверка утечки радона: НЕ дала надёжного результата (см. radon_check.reason)")

    palette = {n["key"]: n["color"] for n in ed._CFG["nuclides"]}
    label_ru = {n["key"]: n["label_ru"] for n in ed._CFG["nuclides"]}

    data = {
        "meta": {
            "detector": "Гамма-1С (УДС-ГЦ-63х63)",
            "vessel": "Маринелли 1 л, ОИСН-16 ρ=1,6 г/см³ (та же геометрия,"
                     " что и Th-232 -- сетка отклика переиспользована)",
            "live_s": meas["live_s"], "real_s": meas["real_s"],
            "bg_live_s": bg["live_s"], "bg_real_s": bg["real_s"],
            "bg_scale_time": bg_scale_time,
            "start_time": meas["start"],
            "fwhm662_keV": ed.FWHM662,
            "e_fit_lo": ed.E_FIT_LO, "e_fit_hi": ed.E_FIT_HI,
            "level_note": "Библиотека и сумм-пики -- IAEA Live Chart of "
                          "Nuclides (decay_rads), быстрый проход без "
                          "перекрёстной проверки LNHB, в отличие от "
                          "библиотеки Th-232. Метод 1 не строился -- "
                          "МК-шаблоны по нуклидам ветви Ra-226 ещё не "
                          "прогнаны.",
        },
        "fwhm_cal": {"k": fwhm_cal["k"], "p": fwhm_cal["p"],
                    "rms_dev_pct": fwhm_cal["rms_dev_pct"],
                    "n_used": fwhm_cal.get("n_used")},
        "passport": {"A_Bq": A_pass, "dA_Bq": dA_pass,
                    "Bq_per_kg": p["bq_per_kg"], "unc_pct": p["unc_pct"],
                    "mass_g": p["mass_g"], "date_certified": p["passport_date"],
                    "date_measured": p["measured_date"],
                    "decay_factor": decay_f},
        "nuclides": [{"key": k, "label_ru": label_ru[k], "color": palette[k]}
                    for k in keys],
        "spectrum": {
            "e_of_ch": e.tolist(),
            "counts": meas["counts"].tolist(),
            "bg_counts": bg_scaled.tolist(),
        },
        "method2_sel": variants["sel"],
        "method2_full": variants["full"],
        "reference_lines": reference_lines,
        "radon_check": radon_check,
    }

    out = os.path.join(HERE, "g1s_ra226_data.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, separators=(",", ":"))
    print("написано: %s (%d КБ)" % (out, os.path.getsize(out) // 1024))
    for tag in ("sel", "full"):
        v = variants[tag]
        print("  %-4s A=%.0f+-%.0f Бк  ratio=%.3f  chi2/ndof=%.2f  линий=%d сумм=%d"
              % (tag, v["A_Bq"], v["dA_Bq"], v["ratio_to_passport"],
                 v["chi2_ndof"], v["n_lines"], v["n_sum_peaks"]))


if __name__ == "__main__":
    main()
