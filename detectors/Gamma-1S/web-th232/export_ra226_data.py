# -*- coding: utf-8 -*-
"""Экспорт данных для страницы Ra-226 -- метод 2 (двух вариантов библиотеки,
отобранная I>=2% / полная цепочка, сумм-пики те же в обоих) и метод 1
(МК-шаблоны по нуклидам, добавлен 10.08.2026 -- прогон
macros/decay_ra226_isotopes.mac, добро оператора).

Схема JSON остаётся СВОЕЙ, короче g1s_th232_data.json (нет варианта "cs"
по отдельной калибровке цезия, нет масок достоверности МК-статистики R66)
-- под фронтенд ra226.js, не g1s-th232.js. run_method1() -- копия
build_templates()/run_method1() из export_data.py (те объявлены ВНУТРИ
main(), не переиспользуемы напрямую), логика воспроизведена дословно,
включая выделение К-рентгена дочерних отдельной сущностью.

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
             "start": str(s.start_datetime),
             "coefs": [float(c) for c in s.energy_cal],
             "n_channels": s.n_channels},
            {"counts": np.asarray(b.counts, dtype=np.int64),
             "e_of_ch": np.asarray(b.channel_to_energy(ch_b), dtype=float),
             "live_s": float(b.live_time), "real_s": float(b.real_time),
             "coefs": [float(c) for c in b.energy_cal],
             "n_channels": b.n_channels})


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


def load_col(path, name):
    """{E_keV: counts} по именованной колонке и число распадов прогона --
    копия одноимённого вложенного помощника export_data.py.main() (та же
    логика; export_data.py объявляет его ВНУТРИ main(), не переиспользуем
    напрямую -- см. докстринг run_method2 выше про урезанные копии в этом
    лёгком конвейере)."""
    hist, N = {}, None
    with open(path, encoding="utf-8", errors="replace") as fh:
        cols = None
        for ln in fh:
            if ln.startswith("#"):
                if "N_primaries" in ln:
                    N = float(ln.split("=")[1])
                continue
            p = ln.rstrip("\n").split(",")
            if cols is None:
                cols = p
                if name not in cols:
                    raise SystemExit(
                        "в %s нет колонки %s -- файл посчитан сборкой до "
                        "разделения по происхождению кванта"
                        % (os.path.basename(path), name))
                continue
            hist[float(p[0])] = float(p[cols.index(name)])
    if not N:
        raise SystemExit("в %s нет N_primaries" % path)
    return hist, N


def run_method1(e, ch_edges, T, y_sel, bgm, sel, NUCS):
    """Метод 1 для Ra-226 -- прогон 10.08.2026 (оператор дал добро после
    вопроса про ХРИ дочерних, которого в лёгкой странице не было).

    Физика та же, что build_templates()/run_method1() в export_data.py, для
    Th-232 давно проверена и не переизобретается: `chain_<id>.csv` (полный
    физический транспорт всей цепочки за один прогон, уже есть в build/,
    прежняя работа) задаёт АМПЛИТУДУ и полную форму `templ_total`; отдельные
    `iso_<Nuc>.csv` (свежий прогон, macros/decay_ra226_isotopes.mac) дают
    только НОРМИРОВАННУЮ долю каждого нуклида в каждом канале -- суммировать
    их напрямую с амплитудой нельзя (систематика nucleusLimits: одиночный
    нуклид регистрирует энергию отдачи и вторичные, которых в цепочечном
    прогоне нет, см. комментарий в export_data.py). К-рентген дочерних
    (`iso_<Nuc>_shield.csv`, колонка `src_xray`) выделяется вычитанием
    ТОЧНОГО подмножества по признаку рождения кванта в самом Geant4
    (model_RDM_AtomicRelaxation), не энергетическим окном -- то самое,
    что оператор просил учесть в шаблонах нуклидов, не только упомянуть
    как ограничение.

    Копия, не переиспользование: build_templates/run_method1 в
    export_data.py объявлены ВНУТРИ main() (замыкание на локальные
    переменные), импортировать напрямую нельзя -- тот же компромисс, что
    уже принят для run_method2 в этом файле.
    """
    hist_chain, N_chain = ed.load_hist(ed.TEMPLATE_CSV)

    hist_iso = {}
    missing = []
    for key, ru, en, col, br, note in NUCS:
        p = os.path.join(ed.BUILD, "iso_%s.csv" % key)
        if not os.path.isfile(p):
            missing.append(key)
            continue
        hist_iso[key] = ed.load_hist(p)
    if missing:
        raise SystemExit(
            "Нет МК-шаблонов индивидуальных нуклидов Ra-226: %s\n"
            "Запустить: cd %s && ./g1s.exe decay_ra226_isotopes.mac vessel 1.60 OISN16"
            % (missing, ed.BUILD))

    xray_frac_of_branch = {}
    xray_dep = {}
    xray_emit = {}
    for key, ru, en, col, br, note in NUCS:
        _, N_iso = hist_iso[key]
        pe = os.path.join(ed.BUILD, "iso_%s_emitx.csv" % key)
        if os.path.isfile(pe):
            hist_x, N_x = load_col(pe, "x_atomic")
        else:
            hist_x, N_x = {}, N_iso     # нуклид не эмитирует рентген вовсе -- легитимный ноль
        tot = 0.0
        for E0, c in hist_x.items():
            if c <= 0:
                continue
            tot += c
            xray_emit[float(E0)] = xray_emit.get(float(E0), 0.0) + (c / N_x) * br
        xray_frac_of_branch[key] = (tot / N_x) * br

        ps = os.path.join(ed.BUILD, "iso_%s_shield.csv" % key)
        if os.path.isfile(ps):
            hist_d, N_d = load_col(ps, "src_xray")
        else:
            hist_d, N_d = {}, N_iso
        xray_dep[key] = ({E0: c for E0, c in hist_d.items() if c > 0}, N_d)
    XRAY_TOTAL_PER_BRANCH = sum(xray_frac_of_branch.values())

    keys1 = [k for k, _, _, _, _, _ in NUCS] + ["XRAY"]

    templ_total = ed.broaden_and_rebin(hist_chain, N_chain, ch_edges, True)
    by_nuc_raw = {}
    for key, ru, en, col, br, note in NUCS:
        hist, N = hist_iso[key]
        by_nuc_raw[key] = ed.broaden_and_rebin(hist, N, ch_edges, True) * br
    xray_raw = {}
    for key, ru, en, col, br, note in NUCS:
        hist_d, N_d = xray_dep[key]
        xray_raw[key] = ed.broaden_and_rebin(hist_d, N_d, ch_edges, True) * br
        by_nuc_raw[key] = by_nuc_raw[key] - xray_raw[key]
        bad = float(by_nuc_raw[key].min())
        if bad < -1e-9 * float(np.max(np.abs(by_nuc_raw[key])) + 1e-30):
            raise SystemExit(
                "%s: рентген больше самого шаблона (%.3e) -- iso_%s.csv и "
                "iso_%s_shield.csv из разных прогонов" % (key, bad, key, key))
        by_nuc_raw[key] = np.maximum(by_nuc_raw[key], 0.0)
    by_nuc_raw["XRAY"] = sum(xray_raw.values())
    iso_sum = sum(by_nuc_raw.values())

    by_nuc = {}
    for k in by_nuc_raw:
        with np.errstate(divide="ignore", invalid="ignore"):
            share = np.where(iso_sum > 0, by_nuc_raw[k] / iso_sum, 0.0)
        by_nuc[k] = templ_total * share

    neg = {k: float(v.min()) for k, v in by_nuc.items() if float(v.min()) < -1e-12}
    if neg:
        raise SystemExit("отрицательные значения в разложении метода 1: %s" % neg)

    resid = float(np.max(np.abs(sum(by_nuc.values()) - templ_total)))
    if resid > 1e-6 * float(np.max(templ_total)):
        raise SystemExit(
            "XRAY: баланс Σ by_nuc == templ_total нарушен, невязка %.3e" % resid)

    coef, dcoef, chi2, ndof, _ = ed.fit_amplitudes(y_sel, [templ_total[sel] * T, bgm])
    A_branch, dA_branch, bg_amp = float(coef[0]), float(dcoef[0]), float(coef[1])

    stack = {k: (by_nuc[k] * A_branch * T).tolist() for k in keys1}
    lines_out = {"template_decays": [{"nuclide": ru, "n": hist_iso[key][1]}
                                     for key, ru, en, col, br, note in NUCS],
                "chain_decays": N_chain}
    return {
        "A_Bq": A_branch, "dA_Bq": dA_branch,
        "bg_amplitude": bg_amp, "d_bg_amplitude": float(dcoef[1]),
        "chi2": chi2, "ndof": ndof, "chi2_ndof": chi2 / ndof,
        "xray_total_per_branch_pct": 100.0 * XRAY_TOTAL_PER_BRANCH,
        "n_channels_fit": int(sel.sum()),
    }, stack, lines_out


RA226_FWHM_ANCHORS = (186.211, 351.932, 609.321, 1764.491, 2204.100)


def main():
    meas, bg = read_pair()
    e = meas["e_of_ch"]
    T = meas["live_s"]

    # ── перекалибровка энергии по чистым линиям цепочки Ra-226 (замечание
    # оператора 10.08.2026: сторонняя программа СпектраЛайн подсветила
    # что-то похожее на линию Tl-208 2614,5 кэВ, что для Ra-226 физически
    # невозможно -- «2614 на спектре нет»). Прямая проверка на самом .spe
    # (см. сессию) нашла калибровку прибора РАБОТОСПОСОБНОЙ, но с
    # остаточным смещением до -2,6 кэВ на чистых линиях. Кроме того,
    # `ed.FWHM_ANCHORS` по умолчанию -- энергии ветви Th-232 (238,6; 583,2;
    # 2614,5 и т.д.): для Ra-226 реальной линии там нет, деконволюция на
    # этих точках подгоняла шум/континуум (не физику) -- худшая точка
    # прежней таблицы ПШПВ, отклонение -35% на 860,6 кэВ, была именно этим
    # артефактом, не свойством детектора. Свои чистые изолированные якоря
    # (без пометки "+суммирование" в библиотеке): 186,211/351,932/609,321/
    # 1764,491/2204,100 кэВ.
    #
    # ПРОБНЫЙ проход fit_fwhm_calibration на СЫРОЙ (заводской) калибровке
    # даёт истинные центроиды этих линий; линейная поправка E=a+b·E_сырое
    # по невязке (центроида минус номинал) применяется ко всему спектру
    # ДО чего-либо ещё. Коэффициенты .spe (meta.cal_sample) остаются как
    # есть -- это то, что реально хранит прибор; поправка -- отдельное
    # поле meta.energy_correction, прозрачно, не подмешана молча.
    probe = ed.fit_fwhm_calibration(meas["counts"], e, anchors=RA226_FWHM_ANCHORS)
    used_anchors = [q for q in probe["points"] if q["used"]]
    if len(used_anchors) < 3:
        raise SystemExit(
            "перекалибровка энергии: годных чистых якорей меньше трёх (%d "
            "из %d) -- %s" % (len(used_anchors), len(probe["points"]),
                             "; ".join("%.1f: %s" % (q["E_nominal"], q["reject"])
                                       for q in probe["points"] if not q["used"])))
    xm = np.array([q["E_centroid"] for q in used_anchors])
    xt = np.array([q["E_nominal"] for q in used_anchors])
    A_lin = np.vstack([np.ones_like(xm), xm]).T
    corr_coef, _, _, _ = np.linalg.lstsq(A_lin, xt, rcond=None)
    corr_a, corr_b = float(corr_coef[0]), float(corr_coef[1])
    resid = xt - (corr_a + corr_b * xm)
    energy_correction = {
        "a": corr_a, "b": corr_b,
        "n_anchors_used": len(used_anchors), "n_anchors_total": len(probe["points"]),
        "rms_residual_keV": float(np.sqrt(np.mean(resid ** 2))),
        "shift_at_609_keV": corr_a + corr_b * 609.321 - 609.321,
        "shift_at_2204_keV": corr_a + corr_b * 2204.100 - 2204.100,
        "anchors": [{"E_nominal": q["E_nominal"], "E_centroid_raw": q["E_centroid"]}
                   for q in used_anchors],
    }
    print("перекалибровка энергии по %d чистым линиям: a=%.4f b=%.6f "
          "(сдвиг на 609/2204 кэВ = %.3f/%.3f кэВ), СКО невязки %.3f кэВ"
          % (len(used_anchors), corr_a, corr_b,
             energy_correction["shift_at_609_keV"],
             energy_correction["shift_at_2204_keV"],
             energy_correction["rms_residual_keV"]))
    e = corr_a + corr_b * e

    bg_on_meas = np.interp(e, bg["e_of_ch"], bg["counts"].astype(float),
                           left=0.0, right=0.0)
    bg_scale_time = T / bg["live_s"]
    bg_scaled = bg_on_meas * bg_scale_time

    ch_edges = np.concatenate((
        [e[0] - 0.5 * (e[1] - e[0])],
        0.5 * (e[:-1] + e[1:]),
        [e[-1] + 0.5 * (e[-1] - e[-2])],
    ))

    fwhm_cal = ed.fit_fwhm_calibration(meas["counts"], e, anchors=RA226_FWHM_ANCHORS)
    ed.FWHM_LAW.update({"kind": "power", "k": fwhm_cal["k"], "p": fwhm_cal["p"]})

    eps_peak = ed.make_eps_peak_interp(os.path.join(ed.BUILD, "grid"))
    resp = ed.make_full_response(os.path.join(ed.BUILD, "grid"), ch_edges,
                                 True, eps_peak)

    sel = (e >= ed.E_FIT_LO) & (e <= ed.E_FIT_HI)
    y_sel = meas["counts"][sel].astype(float)
    bgm = bg_scaled[sel]

    # ── метод 1: МК-шаблоны по нуклидам (прогон 10.08.2026, добро оператора)
    m1_result, m1_stack, m1_meta = run_method1(e, ch_edges, T, y_sel, bgm,
                                               sel, ed.NUCS)
    print("метод 1: A=%.0f+-%.0f Бк  ratio=%.3f  chi2/ndof=%.2f  "
          "рентген=%.3f%% на распад ветви"
          % (m1_result["A_Bq"], m1_result["dA_Bq"],
             m1_result["A_Bq"] / (ed._CFG["passport"]["bq_per_kg"]
                                  * ed._CFG["passport"]["mass_g"] / 1000.0
                                  * ed.decay_factor_years(
                                        ed._CFG["passport"]["half_life_years"],
                                        ed._CFG["passport"]["days_pass_to_meas"])),
             m1_result["chi2_ndof"], m1_result["xray_total_per_branch_pct"]))

    keys = [n[0] for n in ed.NUCS]
    lib_full, _ = ed.load_full_library(nuc_keys=set(keys))

    variants = {}
    by_nuc_w_sel = None
    for tag, library in (("sel", ed.GAMMA_LIBRARY), ("full", lib_full)):
        shape_total, by_nuc_w, lines_out, n_sum = run_method2(
            library, ed.SUM_PEAKS, resp, e, ch_edges, keys)
        if tag == "sel":
            by_nuc_w_sel = by_nuc_w
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

    # ── проверка на утечку радона (замечание оператора №3 09.08.2026,
    # доведено до числа 10.08.2026) ──────────────────────────────────
    # Rn-222 (T1/2=3,82 сут) стоит МЕЖДУ Ra-226 и Pb-214/Bi-214 в цепочке.
    # Одноамплитудный метод 2 выше (variants["sel"]/["full"]) считает
    # ОДНУ амплитуду на всю цепочку (вековое равновесие встроено в саму
    # модель) и утечку в принципе не увидит.
    #
    # ПЕРВАЯ ПОПЫТКА (09.08.2026, здесь удалена) была локальной: площадь
    # линии-родителя 186,211 кэВ методом peak_area_with_shelf (линейное
    # плечо фона в узком окне) против чистой линии 351,932 кэВ. Провалилась
    # содержательно: ПШПВ на 186 кэВ (~22,6 кэВ) настолько широка, что
    # плечо окна дотягивается до соседней сильной линии 241,995 кэВ
    # Pb-214 (56 кэВ между линиями) -- оценка фона завышена, чистая площадь
    # отрицательна при любом разумном наборе окон.
    #
    # НАСТОЯЩАЯ ПРОВЕРКА: та же модель метода 2 (та же матрица отклика,
    # та же библиотека), но с ДВУМЯ независимыми амплитудами вместо одной
    # -- родитель (Ra-226, единственная линия 186,211 кэВ выше порога
    # библиотеки) отдельно от дочерних (Rn-222/Po-218/Pb-214/Bi-214/
    # Po-214). Подгонка по ВСЕМУ диапазону {{e_fit_lo}}-{{e_fit_hi}} кэВ
    # (как метод 2 в целом), не по узкому окну -- контаминация соседней
    # линией 241,995 корректно разделяется весами по всему спектру, а не
    # локальным плечом. Базисные вектора -- ТЕ ЖЕ by_nuc_w, что уже
    # взвешены F_B-депопуляцией и сумм-пиками внутри run_method2(), просто
    # сгруппированы по родитель/дочерние вместо суммирования в одну
    # shape_total.
    PARENT_KEYS = {"Ra226"}
    daughter_keys = [k for k in keys if k not in PARENT_KEYS]
    shape_parent = np.zeros_like(e)
    shape_daughter = np.zeros_like(e)
    for k in keys:
        if k in PARENT_KEYS:
            shape_parent += by_nuc_w_sel[k]
        else:
            shape_daughter += by_nuc_w_sel[k]
    cols_rn = [shape_parent[sel] * T, shape_daughter[sel] * T, bgm]
    coef_rn, dcoef_rn, chi2_rn, ndof_rn, model_rn = ed.fit_amplitudes(
        y_sel, cols_rn, ed.SYS_FLOOR)
    A_par, dA_par = float(coef_rn[0]), float(dcoef_rn[0])
    A_dtr, dA_dtr = float(coef_rn[1]), float(dcoef_rn[1])
    if A_par > 0 and A_dtr > 0:
        ratio_rn = A_dtr / A_par
        d_ratio_rn = ratio_rn * float(np.sqrt((dA_par / A_par) ** 2
                                              + (dA_dtr / A_dtr) ** 2))
    else:
        ratio_rn, d_ratio_rn = float("nan"), float("nan")
    # Число обусловленности взвешенной нормальной матрицы -- диагностика
    # вырожденности (ridge/regularization здесь НЕ применяется: он бы
    # маскировал вырожденность красивым числом, а не показывал её честно).
    A_rn = np.stack(cols_rn, axis=1)
    sig_rn = np.sqrt(np.maximum(y_sel, 1.0) + (ed.SYS_FLOOR * model_rn) ** 2)
    cond_rn = float(np.linalg.cond((A_rn.T * (1.0 / sig_rn ** 2)) @ A_rn))
    # Проверка на СТАБИЛЬНОСТЬ отношения при разных верхних границах окна
    # подгонки (100 -- испытанные 400/700/1200/2300 кэВ дали 0,47-0,57,
    # не порядковый разброс) -- не пересчитывается на каждый прогон (не
    # тот случай, где нужна автоматизация), проверено вручную 10.08.2026,
    # см. журнал сессии.
    A_par_over_pass = A_par / A_pass
    radon_check = {
        "attempted": True, "method": "method2_split2amp",
        "parent_nuclide": "Ra226", "daughter_nuclides": daughter_keys,
        "A_parent_Bq": A_par, "dA_parent_Bq": dA_par,
        "A_daughter_Bq": A_dtr, "dA_daughter_Bq": dA_dtr,
        "ratio_daughter_to_parent": ratio_rn, "d_ratio": d_ratio_rn,
        "chi2": chi2_rn, "ndof": ndof_rn, "chi2_ndof": chi2_rn / ndof_rn,
        "cond_number": cond_rn, "A_parent_over_passport": A_par_over_pass,
        "reliable": False,
        "caveat": "Метод -- деконволюция ПО ВСЕМУ спектру двумя независимыми "
                  "амплитудами (родитель Ra-226 против дочерних Rn-222+Pb-214"
                  "+Bi-214+Po-218/214), теми же базисами, что и single-"
                  "амплитудный метод 2 выше; не локальное окно "
                  "peak_area_with_shelf (та попытка 09.08.2026 провалилась "
                  "содержательно -- см. git-историю этого файла). Число "
                  "устойчиво к границе окна подгонки (проверено вручную на "
                  "400/700/1200/2300 кэВ: отношение 0,47-0,57, не порядковый "
                  "разброс) -- это НЕ численная случайность обусловленности "
                  "(число обусловленности взвешенной нормальной матрицы "
                  "порядка 1e7, см. cond_number). Но САМА величина "
                  "физически неправдоподобна для утечки радона: "
                  "A(родитель) получается почти вдвое БОЛЬШЕ паспортной "
                  "активности источника (A_parent_over_passport), а не "
                  "долей процента ниже паспорта, как ожидалось бы при "
                  "частичной потере Rn-222. Родительская сторона "
                  "определяется ЕДИНСТВЕННОЙ линией 186,211 кэВ (I=3,565%) "
                  "-- избыток в этой области правдоподобнее объясняется "
                  "несмоделированным вкладом (кандидат, НЕ подтверждён: "
                  "линия U-235 185,715 кэВ, I=57,2%, неразличима с Ra-226 "
                  "186,211 при разрешении NaI ~22,6 кэВ -- сертификаты "
                  "радиевых источников иногда несут следовый уран), чем "
                  "реальным дисбалансом цепочки. Число посчитано и "
                  "приводится, но НЕ читается как подтверждённая утечка "
                  "радона -- reliable=false.",
    }
    print("проверка утечки радона: A(родитель)=%.0f+-%.0f Бк  "
          "A(дочерние)=%.0f+-%.0f Бк  отношение=%.3f+-%.3f  chi2/ndof=%.2f  "
          "cond=%.2e  A_par/паспорт=%.2f -- ЧИСЛО НЕ НАДЁЖНО, см. caveat"
          % (A_par, dA_par, A_dtr, dA_dtr, ratio_rn, d_ratio_rn,
             chi2_rn / ndof_rn, cond_rn, A_par_over_pass))

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
            "cal_sample": {"coefs": meas["coefs"],
                          "order": len(meas["coefs"]) - 1,
                          "n_channels": meas["n_channels"]},
            "cal_bg": {"coefs": bg["coefs"], "order": len(bg["coefs"]) - 1,
                      "n_channels": bg["n_channels"]},
            "energy_correction": energy_correction,
            "level_note": "Библиотека и сумм-пики -- IAEA Live Chart of "
                          "Nuclides (decay_rads), быстрый проход без "
                          "перекрёстной проверки LNHB, в отличие от "
                          "библиотеки Th-232. Метод 1 (МК-шаблоны по "
                          "нуклидам, прогон 10.08.2026) построен, включая "
                          "К-рентген дочерних отдельной сущностью; "
                          "статистика слабых звеньев (R66) на порядки ниже, "
                          "чем у Th-232 (тот прогон занял недели ручной "
                          "донастройки, этот -- первый заход).",
        },
        # Полный словарь fit_fwhm_calibration -- те же поля, что несёт
        # g1s_th232_data.json (points/n_anchors/fwhm662_law/fwhm662_cs/
        # res662_pct), калибровочная вкладка одна и та же на обеих
        # страницах (JS-код общий, см. buildCal/drawFwhm в ra226.js).
        "fwhm_cal": fwhm_cal,
        "passport": {"A_Bq": A_pass, "dA_Bq": dA_pass,
                    "Bq_per_kg": p["bq_per_kg"], "unc_pct": p["unc_pct"],
                    "mass_g": p["mass_g"], "date_certified": p["passport_date"],
                    "date_measured": p["measured_date"],
                    "decay_factor": decay_f},
        "nuclides": [{"key": k, "label_ru": label_ru[k], "color": palette[k]}
                    for k in keys]
                   # XRAY -- та же сущность и тот же цвет, что на Th-232
                   # (export_data.py: {"key":"XRAY","color":"#6b5f4a"}).
                   + [{"key": "XRAY", "label_ru": "K-рентген",
                       "color": "#6b5f4a"}],
        "spectrum": {
            "e_of_ch": e.tolist(),
            "counts": meas["counts"].tolist(),
            "bg_counts": bg_scaled.tolist(),
        },
        "method1": m1_result,
        "method1_meta": m1_meta,
        "method2_sel": variants["sel"],
        "method2_full": variants["full"],
        "reference_lines": reference_lines,
        "radon_check": radon_check,
    }
    data["spectrum"]["stack1"] = m1_stack

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
