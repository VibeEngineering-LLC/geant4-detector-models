# -*- coding: utf-8 -*-
"""Экспорт данных страницы Ra-226 (docs/gamma-1s-ra226): метод 1 (МК-шаблоны) и метод 2 (линии × отклик сетки, два
варианта библиотеки — отобранная I≥2 % и полная цепочка; сумм-пики те же в обоих). Схема JSON — под src/scripts/ra226.js.

Пересчёт 16.09.2026 на действующей модели Гамма-1С (CFG1-2026-09-16-ra226-page.md, «да» оператора):
  - шаблоны g1s_npsm (option4, deex=deex, npsm=on, порог 0,05 мм, защита geo1_2026_09_15, облицовка Cu 3 / Cd 1,5 мм);
    формат «ключ,значение» читает адаптер read_g1s_npsm поверх mix_unfold_g1s.read_template;
  - метод 1 (М1-Б): сумма k24g1_Ra226 + k24g1_Rn222chain с одной амплитудой, доли — r226g1_<нуклид>;
    К-рентген дочерних отдельно не выделяется (Р-А, ключ XRAY нулевой);
  - метод 2 (М2-Б): 24-узловая сетка r226g1_grid/gridon_E*.csv, узлы — light_grid.build_nodes;
  - #CAL-0: шкала света по реперам самого файла (data/light_scale_ra226_2016.json), ПШПВ — трёхчленный закон по точкам
    комплекта 2016 (data/fwhm_law3_g1s_2016.*); критерии A2 / A2V (D-020);
  - разрезы метода 2 — справочно, утечка радона по линии 186 кэВ не выводится (У-Б).
Прежняя версия (10.08.2026: chain_Ra226.csv + iso_*.csv, гауссово уширение, заводская ПШПВ, якорь 77 кэВ) — в истории git.

Запуск:
    SPECTRAVIBE_ROOT=<...> G4MODELS_RA226_BG_SPE=<фон .spe> python export_ra226_data.py
Выход: g1s_ra226_data.json и src/data-ra226.js (window.RA226 = …).
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

SPECTRAVIBE_ROOT = os.environ.get("SPECTRAVIBE_ROOT")
if not SPECTRAVIBE_ROOT:
    raise SystemExit("Задайте SPECTRAVIBE_ROOT -- путь к рабочему каталогу "
                      "gamma-spectrum-analysis (см. README.md).")
sys.path.insert(0, os.path.join(SPECTRAVIBE_ROOT, "scripts"))
from gamma.io.lsrm_spe import read_lsrm_spe  # noqa: E402

# ── r226g1 (16.09.2026): пересчёт на действующей модели Гамма-1С по CFG1-2026-09-16-ra226-page.md («да» оператора:
# В1 М1-Б, В2 М2-Б, В3 Р-А, В4 Ш-Б, В5 У-Б, В6 100–2300 кэВ, В7 #CAL-0). Шаблоны g1s_npsm: option4, deex=deex, npsm=on,
# порог 0,05 мм, защита geo1_2026_09_15. Разбор — core.unfold (левый хвост T, свёртка в каналах, шкала света файла).
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "analysis"))
import yaml                      # noqa: E402
import mix_unfold_g1s as g1s     # noqa: E402  read_template, read_lsrm_spe, make_fwhm, broaden_ch, LSRM_CH_OFFSET
import mix_unfold_core as core   # noqa: E402  unfold
import light_mode                # noqa: E402  load_scale, light_to_energy
import light_grid                # noqa: E402  build_nodes
import crit_bench_m2 as cbm2     # noqa: E402  fit_A2V
BUILD_OUT = os.path.join(os.environ.get("G4MODELS_BUILD_GAMMA_1S_NPSM",
                                        os.path.join(ed.REPO, "build", "Gamma-1S-npsm-1142")), "out")
GRID_DIR = os.path.join(BUILD_OUT, "r226g1_grid")            # М2-Б, run_r226g1.sh
CHAIN_TPLS = (("Ra226", "k24g1_Ra226_marinelli_OISN06_epoxy.csv"),           # М1-Б: форма и амплитуда —
              ("Rn222chain", "k24g1_Rn222chain_marinelli_OISN06_epoxy.csv"))  # сумма двух, одна амплитуда
ISO_TPL = "r226g1_%s_marinelli_OISN06_epoxy.csv"                # доли нуклидов, склейка post_r226g1.sh
FWHM_POINTS = os.path.join(HERE, "data", "fwhm_points_g1s_2016.csv")
FWHM_LAW_CSV = os.path.join(HERE, "data", "fwhm_law3_g1s_2016.csv")     # prep_r226g1_page.py
FWHM_LAW_JSON = os.path.join(HERE, "data", "fwhm_law3_g1s_2016.json")
LIGHT_JSON = os.path.join(HERE, "data", "light_scale_ra226_2016.json")
TAIL_T = 0.75
BLUR = float(yaml.safe_load(open(os.path.join(HERE, "configs", "amticseu.yaml"), encoding="utf-8"))["fit"]["blur"])
SETUP = {"em_cut_mm": "0.05", "em_deex": "deex", "npsm_enabled": "1", "shield_variant": "geo1_2026_09_15",
         "vessel": "marinelli", "sample_matrix": "OISN06_epoxy", "shield_cu_mm": "3", "shield_cd_mm": "1.5"}


def read_g1s_npsm(path, expect=SETUP):
    """Адаптер формата g1s_npsm (CFG1 §5): гистограмма и число событий — g1s.read_template (импорт, не копия);
    шапка «ключ,значение» до маркера bin_keV сверяется с постановкой #CFG-1 — прогон не той постановки даёт отказ."""
    hist, n, npsm = g1s.read_template(path)
    header = {}
    with open(path, encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("bin_keV"):
                break
            key, _, val = ln.strip().partition(",")
            header.setdefault(key, val)
    bad = {k: header.get(k) for k, v in expect.items() if header.get(k) != v}
    if bad or npsm != 1:
        raise SystemExit("ОТКАЗ адаптера: %s не той постановки: %s, npsm=%s" % (os.path.basename(path), bad, npsm))
    return hist, n, header


def hist_array(hist, n_bins=3200):
    """Плотный массив гистограммы по бинам 1 кэВ — для сверки адаптера с read_template (приёмка 2)."""
    a = np.zeros(n_bins)
    for k, v in hist.items():
        a[int(np.floor(k))] += v
    return a


def light_col(hist, n, e, fwhm, file_e, scale):
    """Столбец шаблона в каналах пробы тем же путём, что core.unfold (light_mode.light_to_energy → g1s.broaden_ch),
    но без отказа на пустом переводе: у Po-218 гамма-квантов нет, его доля — законный ноль."""
    e_of = lambda c: float(np.interp(c, np.arange(len(e), dtype=float), e))
    h, st = light_mode.light_to_energy(hist, scale[0], scale[1], e_of, len(e), c=(scale[2] if len(scale) > 2 else 0.0))
    return g1s.broaden_ch(h, n, e, fwhm, file_e, BLUR) if st["sum_used"] > 0 else np.zeros(len(e))


def run_method1_npsm(spec, bgs, e_bg, scale):
    """Метод 1 (М1-Б): столбцы k24g1_Ra226 и k24g1_Rn222chain — core.unfold (шкала света, хвост, свёртка в каналах);
    форма и амплитуда — их СУММА с одной амплитудой (A2, crit_bench.fit_A2, D-020); доли нуклидов по каналам — прогоны
    r226g1_<нуклид> (логика прежнего run_method1). Двухамплитудный разбор core.unfold печатается справочно."""
    paths = [(g, os.path.join(BUILD_OUT, f)) for g, f in CHAIN_TPLS]
    ns = [read_g1s_npsm(p)[1] for _, p in paths]
    if len(set(ns)) != 1:
        raise SystemExit("ОТКАЗ метод 1: у шаблонов суммы разное число событий %s" % ns)
    r = core.unfold(spec, bgs, paths, FWHM_LAW_CSV, lo=ed.E_FIT_LO, hi=ed.E_FIT_HI, recalibrate=False, tail=TAIL_T,
                    bg_energy_of_ch=e_bg, verbose=False, conv="channel", blur=BLUR, ch_offset=g1s.LSRM_CH_OFFSET,
                    light_scale=scale)
    col, cb = r["cols"].sum(axis=0)[None, :], core.cb
    k_bg, cnt, sel, T = spec.live_time / bgs.live_time, r["counts"], r["sel"], r["live_s"]
    a, sd, _ = cb.fit_A2(col, cnt, r["bg_scaled"], k_bg, ns[:1], sel)
    var = cb.fit_A1(col, cnt, r["bg_scaled"], k_bg, ns[:1], sel)[2]["var"] + col[0, sel] / ns[0] * a[0] ** 2
    chi2, ndof = float(np.sum((col[0, sel] * a[0] - r["net"][sel]) ** 2 / var)), int(sel.sum()) - 1
    a_e1 = cb.fit_E1(col, cnt, r["bg_scaled"], k_bg, ns[:1], sel)[0]
    print("метод 1, справочно два шаблона (core.unfold, A2): %s"
          % ", ".join("%s %.1f Бк" % (g, A) for g, A in zip(r["names"], r["activities"])))
    total = col[0] * a[0]
    file_e = lambda c: float(spec.channel_to_energy(c + g1s.LSRM_CH_OFFSET))
    raw, decays = {}, []
    for key, ru, en, color, br, note in ed.NUCS:
        hist, n, _ = read_g1s_npsm(os.path.join(BUILD_OUT, ISO_TPL % key))
        raw[key] = light_col(hist, n, r["e"], r["fwhm"], file_e, scale) * br
        decays.append({"nuclide": ru, "n": n})
    iso_sum = sum(raw.values())
    with np.errstate(divide="ignore", invalid="ignore"):
        stack = {k: np.where(iso_sum > 0, total * v / iso_sum, 0.0) for k, v in raw.items()}
    lost = float(np.sum(total[sel]) - np.sum(sum(stack.values())[sel]))
    if abs(lost) > 1e-6 * float(np.sum(total[sel])):
        raise SystemExit("ОТКАЗ метод 1: Σ долей нуклидов ≠ модели в окне подгонки, потеря %.3e отсчётов" % lost)
    stack["XRAY"] = np.zeros(len(total))    # В3 = Р-А: К-рентген дочерних отдельно не выделяется
    res = {"A_Bq": float(a[0] / T), "dA_Bq": float(sd[0] / T), "bg_amplitude": 1.0, "d_bg_amplitude": 0.0,
           "chi2": chi2, "ndof": ndof, "chi2_ndof": chi2 / ndof, "xray_total_per_branch_pct": None,
           "n_channels_fit": int(sel.sum())}
    extra = {"A_E1_Bq": float(a_e1[0] / T), "two": dict(zip(r["names"], map(float, r["activities"]))), "r": r}
    return (res, {k: [round(float(x), 4) for x in v] for k, v in stack.items()},
            {"template_decays": decays, "chain_decays": int(ns[0])}, extra)


def grid_resp_shift(e, fwhm, file_e, scale):
    """Отклик метода 2 на сетке r226g1_grid (М2-Б, 24 узла): узлы — light_grid.build_nodes (свет → каналы пробы);
    между узлами — форма ближайшего узла, сдвинутая на E − E_узла и приведённая к пиковой эффективности,
    интерполированной лог-лог (приближение export_data.make_full_response, CFG1 §3 ДОП-6)."""
    cq = scale[2] if len(scale) > 2 else 0.0
    l2e = lambda h, a, b, eo, nc: light_mode.light_to_energy(h, a, b, eo, nc, c=cq)
    nodes, skipped, empty = light_grid.build_nodes(GRID_DIR, "gridon_E", g1s.read_template, l2e, g1s.broaden_ch,
                                                   e, file_e, fwhm, scale[0], scale[1], BLUR)
    names = [p[len("gridon_E"):-len(".csv")] for p in os.listdir(GRID_DIR)]
    n_map = {float(s): read_g1s_npsm(os.path.join(GRID_DIR, "gridon_E%s.csv" % s))[1] for s in names}
    print(light_grid.summary(nodes, empty))
    if len(nodes) != 24 or skipped or set(n_map) != set(nodes):
        raise SystemExit("ОТКАЗ сетки метода 2: узлов %d, пропущено %s" % (len(nodes), skipped))
    Es = np.array(sorted(nodes))
    eps = np.array([nodes[E][1] for E in Es])
    lE, lP = np.log(Es[eps > 0]), np.log(eps[eps > 0])

    def eps_at(E):
        x = np.log(E)
        if lE[0] <= x <= lE[-1]:
            return float(np.exp(np.interp(x, lE, lP)))
        i0, i1 = (0, 1) if x < lE[0] else (-2, -1)
        return float(np.exp(lP[i1] + (lP[i1] - lP[i0]) / (lE[i1] - lE[i0]) * (x - lE[i1])))

    def resp(E):
        j = int(np.argmin(np.abs(Es - E)))
        shape, eps_j = nodes[Es[j]][0], nodes[Es[j]][1]
        s = shape if abs(E - Es[j]) < 1e-9 else np.interp(e - (E - Es[j]), e, shape, left=0.0, right=0.0)
        ea = eps_at(E)
        return s * (ea / eps_j if eps_j > 0 else 0.0), {}, ea
    resp.n_of = lambda E: n_map[Es[int(np.argmin(np.abs(Es - E)))]]
    return resp


def fit_m2(W, V, counts, bg_scaled, k_bg, sel):
    """Критерий A2V с дисперсией узлов сетки — как метод 2 смеси 2016 (export_amticseu_data.py, D-020 п.3).
    Возвращает амплитуды (отсчёты за живое время), их sd, χ² по дисперсии критерия, ν, число обусловленности."""
    a, sd, _ = cbm2.fit_A2V(W, counts, bg_scaled, k_bg, V, sel)
    var = np.maximum(counts[sel] + k_bg * bg_scaled[sel], 1.0) + V[:, sel].T @ a ** 2
    chi2 = float(np.sum(((W.T @ a)[sel] - (counts - bg_scaled)[sel]) ** 2 / var))
    Wn = W[:, sel] / np.sqrt(var)
    return a, sd, chi2, int(sel.sum()) - W.shape[0], float(np.linalg.cond(Wn @ Wn.T))


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
    sf = getattr(s, "stored_fwhm_calibration", None)
    fwhm_coefs = (list(sf.coefficients) if sf is not None and sf.coefficients
                 else None)
    fwhm_model = sf.model if sf is not None else None
    return ({"counts": np.asarray(s.counts, dtype=np.int64),
             "e_of_ch": np.asarray(s.channel_to_energy(ch_s), dtype=float),
             "live_s": float(s.live_time), "real_s": float(s.real_time),
             "start": str(s.start_datetime),
             "coefs": [float(c) for c in s.energy_cal],
             "n_channels": s.n_channels,
             "fwhm_coefs": fwhm_coefs, "fwhm_model": fwhm_model},
            {"counts": np.asarray(b.counts, dtype=np.int64),
             "e_of_ch": np.asarray(b.channel_to_energy(ch_b), dtype=float),
             "live_s": float(b.live_time), "real_s": float(b.real_time),
             "coefs": [float(c) for c in b.energy_cal],
             "n_channels": b.n_channels})


def factory_fwhm_keV(coefs, model, E_keV):
    """ПШПВ(E) по ЗАВОДСКОЙ калибровке прибора (полином из шапки .spe),
    без какой-либо собственной деконволюции/подгонки по измеренному
    спектру -- директива оператора 10.08.2026 "сам не калибруй".

    Формат `lsrm_fwhm_polynomial_in_E` (LSRM .spe) вопреки названию
    берёт аргументом НЕ E, а z=sqrt(E) -- задокументированная особенность
    формата (BUG-22 в SpectraVibe, LSRM «Алгоритмические основы» §8.3,
    `gamma/io/lsrm_spe.py:46-58`; прямая подстановка E даёт отрицательные
    ПШПВ). ПШПВ(E) = sum_k c_k * sqrt(E)^k.
    """
    if not coefs or model != "lsrm_fwhm_polynomial_in_E" or E_keV <= 0:
        return None
    z = float(E_keV) ** 0.5
    val = sum(float(ck) * (z ** k) for k, ck in enumerate(coefs))
    return float(val) if val > 0 else None


def fit_power_law_to_factory_fwhm(coefs, model, e_lo=50.0, e_hi=3000.0, n=60):
    """Аппроксимация заводской ПШПВ(E) степенным законом k*E^p -- НЕ
    калибровка (данные не измеренные, а сама заводская функция,
    посчитанная в n точках), а пересчёт формы под существующий
    интерфейс `ed.FWHM_LAW`/`ed.fwhm_kev()` (используется в окнах
    деконволюции К-рентгена, response-грид), который принимает только
    степенной закон. МНК по логарифмам -- невязка (rms_fit_pct)
    сохраняется в JSON как честная мера точности АППРОКСИМАЦИИ, не
    точности самой заводской калибровки.
    """
    Es = np.geomspace(e_lo, e_hi, n)
    vals = [factory_fwhm_keV(coefs, model, E) for E in Es]
    ws = np.array([v if v is not None else float("nan") for v in vals])
    ok = np.isfinite(ws) & (ws > 0)
    Es, ws = Es[ok], ws[ok]
    x = np.log(Es); y = np.log(ws)
    A = np.vstack([np.ones_like(x), x]).T
    coef, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    k = float(np.exp(coef[0])); p = float(coef[1])
    model_ws = k * Es ** p
    rms_fit_pct = float(np.sqrt(np.mean((model_ws / ws - 1.0) ** 2)) * 100.0)
    return k, p, rms_fit_pct




def run_method2(library, sums, resp, e, ch_edges, keys, var_acc=None, n_of=None):
    """Урезанный, но физически тот же run_method2, что в export_data.py
    (F_B-депопуляция, F_B-нормировка сумм-пиков) -- без канальной
    раскладки и без диагностики peak_area_with_shelf (не нужны лёгкой
    странице).
    var_acc/n_of (13.09.2026, D-020 для метода 2, scripts/_spec_crit_bench_m2.md П1): при словаре
    var_acc накапливается дисперсия столбца от шума узлов сетки w**2 * shape / n_of(E узла).
    Без них поведение прежнее."""
    if var_acc is not None and n_of is None:
        raise SystemExit("ОТКАЗ run_method2: var_acc передан без n_of (число событий узла)")
    shape_total = np.zeros_like(e)
    by_nuc_w = {k: np.zeros_like(e) for k in keys}

    def add(nuc_key, weight, shp, E_node):
        shape_total[:] += weight * shp
        by_nuc_w[nuc_key] += weight * shp
        if var_acc is not None:
            if nuc_key not in var_acc:
                var_acc[nuc_key] = np.zeros_like(e)
            var_acc[nuc_key] += weight ** 2 * shp / float(n_of(E_node))

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
        add(nuc_key, w, shp, E)
        lines_out.append({"E_keV": E, "nuclide": nuc_key, "I_pct": I_pct,
                          "note": note, "kind": "line",
                          "depleted_pct": depl_pct,
                          "eps_peak": eps, "weight_per_branch": w * eps})

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
        add(nuc_key, w, shp, Esum)
        lines_out.append({"E_keV": Esum, "nuclide": nuc_key, "I_pct": None,
                          "note": note, "kind": "sum",
                          "E1_keV": E1, "E2_keV": E2,
                          "eps_peak": eps_sum_node, "weight_per_branch": w * eps_sum_node})
        n_sum_used += 1

    return shape_total, by_nuc_w, lines_out, n_sum_used






VESSEL_NOTE = ("Маринелли 1 л, ОИСН-06 ρ=0,60 г/см³ (насыпная эпоксидка, источник -18/2016); шаблоны Geant4 16.09.2026: "
               "g1s_npsm, option4, deex=deex, npsm=on, порог 0,05 мм, защита по замерам оператора (geo1_2026_09_15), "
               "облицовка Cu 3 / Cd 1,5 мм")
LEVEL_NOTE = ("Библиотека и сумм-пики — IAEA Live Chart of Nuclides (decay_rads). Метод 1 — сумма шаблонов распада Ra-226 "
              "и Rn-цепочки (по 10⁷) с одной амплитудой, доли нуклидов — отдельные прогоны по 10⁷ распадов; К-рентген "
              "дочерних отдельно не выделяется (входит в шаблоны своих нуклидов). Метод 2 — 24-узловая моноэнергетическая "
              "сетка по 2·10⁶ квантов, форма между узлами — сдвиг ближайшего узла. Критерии A2/A2V (D-020), шкала света "
              "по реперам самого файла (#CAL-0), ПШПВ — трёхчленный закон по точкам комплекта поверки 2016.")
CAVEAT_SPLIT2 = ("Разрез метода 2 на родителя (Ra-226, линия 186,211 кэВ) и дочерние — справочно. Родитель здесь и группа "
                 "Ra-226 в radon_dpr_vs_ra_check — один фит в разной группировке меток (у Rn-222/Po-218/Po-214 линий в "
                 "библиотеке нет). Амплитуда линии 186,211 кэВ спектром не определяется: линия лежит на пике обратного "
                 "рассеяния Bi-214 609 кэВ, где форма модели не согласована с измерением (METHOD.md, дыра 26; "
                 "RESULT-2026-09-14-kit2024-fit.md §4 п.1).")
CAVEAT_SPLIT3 = ("Три свободные амплитуды (остальные звенья / Pb-214 / Bi-214) — справочно; равновесие цепочки разорвано "
                 "намеренно. Расхождение групп может быть систематикой формы модели и обусловленности совместного фита, "
                 "а не физикой источника.")
CAVEAT_DPR = ("Пересчёт 16.09.2026: прежняя «принятая версия» утечки радона (10.08.2026) снята. Активность Ra-226 "
              "определяется по Rn-цепочке (ДПР); амплитуда группы Ra-226 по линии 186,211 кэВ спектром не определяется "
              "(пик обратного рассеяния Bi-214 609 кэВ, METHOD.md дыра 26), поэтому доля утечки не выводится "
              "(leak_fraction = null). На комплекте поверки 2024 того же источника принудительное равновесие Ra = Rn "
              "ухудшает χ² лишь на 6–9 % (RESULT-2026-09-14-kit2024-fit.md §4 п.1).")


def fin(x):
    """float, если число конечно; иначе None (NaN/inf в JSON не пускать)."""
    return float(x) if x is not None and np.isfinite(x) else None


def main():
    """Пересчёт r226g1: #CAL-0 (шкала файла и шкала света по его реперам), ПШПВ — трёхчленный закон 2016, метод 1 (М1-Б),
    метод 2 (М2-Б), разрезы метода 2 справочно (утечка радона по линии 186 кэВ не выводится, В5 = У-Б)."""
    cfg, p = ed._CFG, ed._CFG["passport"]
    spe_path = os.path.join(str(ed.KIT), *cfg["source"]["measured_sample_spe_rel"].split("/"))
    bg_path = os.environ.get(cfg["source"]["measured_background_spe_env"])
    if not bg_path:
        raise SystemExit("Задайте %s — путь к фоновому .spe." % cfg["source"]["measured_background_spe_env"])
    spec, bgs = g1s.read_lsrm_spe(spe_path), g1s.read_lsrm_spe(bg_path)
    light = json.load(open(LIGHT_JSON, encoding="utf-8"))
    law = json.load(open(FWHM_LAW_JSON, encoding="utf-8"))
    if spec.sample_id != "Ra226_420-7-18" or light["files"] != [os.path.basename(spe_path)] \
            or light["coefs"] != [float(c) for c in spec.energy_cal]:
        raise SystemExit("ОТКАЗ #CAL-0: образец %r или шкала света построена не по этому файлу" % spec.sample_id)
    g1s.TAIL_T, scale = TAIL_T, light_mode.load_scale(LIGHT_JSON)
    print("#CAL-0 %s: коэффициенты файла %s; фона %s" % (spec.sample_id, list(spec.energy_cal), list(bgs.energy_cal)))
    print("  шкала света a=%.6g b=%.6g c=%.6g, реперов %d, СКО невязок %.3f кан"
          % (light["a"], light["b"], light["c"], light["n_used"], light["std_light"]))
    for rf in light["refs"]:
        print("  E=%8.3f c=%8.3f невязка %+.3f кан" % (rf["E"], rf["c"], rf["resid_light"]))
    print("ПШПВ: sqrt(%.6g + %.6g·E + %.6g·E²), точек %d, СКО закона %.2f %%"
          % (*law["params"], law["n_points"], law["rms_rel_pct"]))
    T, k_bg = float(spec.live_time), float(spec.live_time) / float(bgs.live_time)
    e_bg = np.array([bgs.channel_to_energy(i) for i in range(len(bgs.counts))])   # нумерация как в core.unfold
    m1_result, m1_stack, m1_meta, m1x = run_method1_npsm(spec, bgs, e_bg, scale)
    r = m1x["r"]
    e, sel, cnt, bg_scaled, fwhm = r["e"], r["sel"], r["counts"], r["bg_scaled"], r["fwhm"]
    A_pass = p["bq_per_kg"] * (p["mass_g"] / 1000.0) * ed.decay_factor_years(p["half_life_years"], p["days_pass_to_meas"])
    dA_pass = A_pass * p["unc_pct"] / 100.0
    print("метод 1: A = %.1f ± %.1f Бк, к паспорту %.3f, χ²/ν = %.2f; вторая мера E1 %.1f Бк"
          % (m1_result["A_Bq"], m1_result["dA_Bq"], m1_result["A_Bq"] / A_pass, m1_result["chi2_ndof"], m1x["A_E1_Bq"]))
    file_e = lambda c: float(spec.channel_to_energy(c + g1s.LSRM_CH_OFFSET))
    resp = grid_resp_shift(e, fwhm, file_e, scale)
    keys = [n[0] for n in ed.NUCS]
    lib_full, _ = ed.load_full_library(nuc_keys=set(keys))
    variants, var_sel, by_nuc_w_sel = {}, None, None
    for tag, library in (("sel", ed.GAMMA_LIBRARY), ("full", lib_full)):
        var = {}
        shape_total, by_nuc_w, lines_out, n_sum = run_method2(library, ed.SUM_PEAKS, resp, e, r["ch_edges"], keys,
                                                              var_acc=var, n_of=resp.n_of)
        a, sd, chi2, ndof, _ = fit_m2(np.array([shape_total]), np.array([sum(var.values())]), cnt, bg_scaled, k_bg, sel)
        A_Bq, dA_Bq = float(a[0] / T), float(sd[0] / T)
        for ln in lines_out:
            ln["predicted_net"] = ln.get("weight_per_branch", 0.0) * A_Bq * T
        variants[tag] = {"A_Bq": A_Bq, "dA_Bq": dA_Bq, "bg_amplitude": 1.0, "chi2": chi2, "ndof": ndof,
                         "chi2_ndof": chi2 / ndof, "n_lines": len(library), "n_sum_peaks": n_sum,
                         "n_sum_peaks_total": len(ed.SUM_PEAKS), "ratio_to_passport": A_Bq / A_pass,
                         "d_ratio": dA_Bq / A_pass, "lines": lines_out,
                         "stack": {k: [round(float(x), 4) for x in by_nuc_w[k] * a[0]] for k in keys}}
        if tag == "sel":
            var_sel, by_nuc_w_sel = var, by_nuc_w
        print("метод 2 (%s): A = %.1f ± %.1f Бк, к паспорту %.3f, χ²/ν = %.2f, линий %d, сумм-пиков %d"
              % (tag, A_Bq, dA_Bq, A_Bq / A_pass, chi2 / ndof, len(library), n_sum))
    zero = np.zeros_like(e)

    def split(groups):
        W = np.array([sum((by_nuc_w_sel[k] for k in g), zero) for g in groups])
        V = np.array([sum((var_sel.get(k, zero) for k in g), zero) for g in groups])
        a, sd, chi2, ndof, cond = fit_m2(W, V, cnt, bg_scaled, k_bg, sel)
        return a / T, sd / T, chi2 / ndof, ndof, chi2, cond

    daughter_keys = [k for k in keys if k != "Ra226"]
    A2, dA2, c2n, nd2, ch2, cond2 = split([["Ra226"], daughter_keys])
    ratio_rn = fin(A2[1] / A2[0]) if A2[0] > 0 and A2[1] > 0 else None
    d_ratio_rn = fin(ratio_rn * np.hypot(dA2[0] / A2[0], dA2[1] / A2[1])) if ratio_rn else None
    radon_check = {"attempted": True, "method": "method2_split2amp", "parent_nuclide": "Ra226",
                   "daughter_nuclides": daughter_keys, "A_parent_Bq": fin(A2[0]), "dA_parent_Bq": fin(dA2[0]),
                   "A_daughter_Bq": fin(A2[1]), "dA_daughter_Bq": fin(dA2[1]),
                   "ratio_daughter_to_parent": ratio_rn, "d_ratio": d_ratio_rn, "chi2": ch2, "ndof": nd2,
                   "chi2_ndof": c2n, "cond_number": cond2, "A_parent_over_passport": fin(A2[0] / A_pass),
                   "reliable": False, "caveat": CAVEAT_SPLIT2}
    A3, dA3, c3n, _, _, _ = split([["Ra226", "Rn222", "Po218", "Po214"], ["Pb214"], ["Bi214"]])
    split3_check = {"attempted": True, "method": "method2_split3amp",
                    "groups": {"rest_Ra226_Rn222_Po218_Po214": {"A_Bq": fin(A3[0]), "dA_Bq": fin(dA3[0])},
                               "Pb214": {"A_Bq": fin(A3[1]), "dA_Bq": fin(dA3[1]), "A_over_passport": fin(A3[1] / A_pass)},
                               "Bi214": {"A_Bq": fin(A3[2]), "dA_Bq": fin(dA3[2]), "A_over_passport": fin(A3[2] / A_pass)}},
                    "chi2_ndof": c3n, "reliable": False, "caveat": CAVEAT_SPLIT3}
    Ad, dAd, cdn, _, _, _ = split([["Pb214", "Bi214"], ["Ra226", "Rn222", "Po218", "Po214"]])
    radon_dpr_vs_ra_check = {"attempted": True, "method": "dpr_vs_ra226_both_free",
                             "A_DPR_Bq": fin(Ad[0]), "dA_DPR_Bq": fin(dAd[0]), "A_DPR_over_passport": fin(Ad[0] / A_pass),
                             "A_Ra226_Bq": fin(Ad[1]), "dA_Ra226_Bq": fin(dAd[1]),
                             "A_Ra226_over_passport": fin(Ad[1] / A_pass),
                             "ratio_Ra226_to_DPR": fin(Ad[1] / Ad[0]) if Ad[0] > 0 else None,
                             "leak_fraction": None, "d_leak_fraction": None, "leak_fraction_window_range": None,
                             "chi2_ndof": cdn, "accepted": False, "caveat": CAVEAT_DPR}
    print("разрезы метода 2 (справочно): родитель %.1f / дочерние %.1f Бк; Pb214 %.3f, Bi214 %.3f паспорта; "
          "ДПР %.3f, группа Ra-226 %.3f паспорта — утечка не выводится (В5 = У-Б)"
          % (A2[0], A2[1], A3[1] / A_pass, A3[2] / A_pass, Ad[0] / A_pass, Ad[1] / A_pass))
    pts = sorted(tuple(map(float, ln.split(",")[:2])) for ln in open(FWHM_POINTS, encoding="utf-8") if ln[:1].isdigit())
    Eg = np.geomspace(50.0, 3000.0, 60)
    wg = np.array([fwhm(E) for E in Eg])
    p_pl, lnk = np.polyfit(np.log(Eg), np.log(wg), 1)
    k_pl, p_pl = float(np.exp(lnk)), float(p_pl)
    w_cs = next(w for E, w in pts if abs(E - 661.657) < 0.01)
    fwhm_cal = {"source": "измеренные точки комплекта поверки 2016 (%d линий), трёхчленный закон" % len(pts),
                "coefs": law["params"], "model": "sqrt(a + b*E + c*E^2)", "k": k_pl, "p": p_pl,
                "fit_rms_pct": float(np.sqrt(np.mean((k_pl * Eg ** p_pl / wg - 1.0) ** 2)) * 100.0),
                "fwhm662_law": k_pl * 661.657 ** p_pl, "fwhm662_cs": w_cs,
                "res662_pct": 100.0 * k_pl * 661.657 ** p_pl / 661.657,
                "reference_points": [{"E_keV": E, "fwhm_factory_keV": w, "fwhm_power_law_keV": k_pl * E ** p_pl}
                                     for E, w in pts]}
    energy_correction = {"applied": False, "method": "cal0_light_scale_file_refs", "anchor_kev": None,
                         "sample_anchor_channels": [rf["c"] for rf in light["refs"]],
                         "sample_anchor_kev": [rf["E"] for rf in light["refs"]], "blend_lo_ch": None, "blend_hi_ch": None,
                         "note": "#CAL-0: шкала файла проверена по реперам его таблицы пиков; шаблоны npsm=on переведены "
                                 "шкалой «канал = a + b·свет + c·свет²» по тем же реперам. Прежняя локальная поправка "
                                 "по якорю 77 кэВ снята."}
    palette = {n["key"]: n["color"] for n in cfg["nuclides"]}
    label_ru = {n["key"]: n["label_ru"] for n in cfg["nuclides"]}
    data = {
        "meta": {"detector": "Гамма-1С (УДС-ГЦ-63х63)", "vessel": VESSEL_NOTE,
                 "live_s": T, "real_s": float(spec.real_time), "bg_live_s": float(bgs.live_time),
                 "bg_real_s": float(bgs.real_time), "bg_scale_time": k_bg, "start_time": str(spec.start_datetime),
                 "fwhm662_keV": w_cs, "e_fit_lo": ed.E_FIT_LO, "e_fit_hi": ed.E_FIT_HI,
                 "cal_sample": {"coefs": [float(c) for c in spec.energy_cal], "order": len(spec.energy_cal) - 1,
                                "n_channels": int(spec.n_channels)},
                 "cal_bg": {"coefs": [float(c) for c in bgs.energy_cal], "order": len(bgs.energy_cal) - 1,
                            "n_channels": len(bgs.counts)},
                 "energy_correction": energy_correction, "level_note": LEVEL_NOTE},
        "fwhm_cal": fwhm_cal,
        "passport": {"A_Bq": A_pass, "dA_Bq": dA_pass, "Bq_per_kg": p["bq_per_kg"], "unc_pct": p["unc_pct"],
                     "mass_g": p["mass_g"], "date_certified": p["passport_date"], "date_measured": p["measured_date"],
                     "decay_factor": ed.decay_factor_years(p["half_life_years"], p["days_pass_to_meas"])},
        "nuclides": [{"key": k, "label_ru": label_ru[k], "color": palette[k]} for k in keys]
                    + [{"key": "XRAY", "label_ru": "K-рентген", "color": "#6b5f4a"}],
        # Ось энергии канала: полином файла в нумерации таблицы пиков (индекс + LSRM_CH_OFFSET, #CH-1), проверен реперами.
        "spectrum": {"e_of_ch": [round(float(spec.channel_to_energy(i + g1s.LSRM_CH_OFFSET)), 4)
                                 for i in range(spec.n_channels)],
                     "counts": [int(c) for c in spec.counts],
                     "bg_counts": [round(float(x), 4) for x in bg_scaled], "stack1": m1_stack},
        "method1": m1_result, "method1_meta": m1_meta,
        "method2_sel": variants["sel"], "method2_full": variants["full"],
        "reference_lines": [[ln["E_keV"], ln["nuclide"]] for ln in variants["sel"]["lines"] if ln["kind"] == "line"],
        "radon_check": radon_check, "split3_check": split3_check, "radon_dpr_vs_ra_check": radon_dpr_vs_ra_check,
    }
    js = json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    with open(os.path.join(HERE, "g1s_ra226_data.json"), "w", encoding="utf-8") as f:
        f.write(js)
    with open(os.path.join(HERE, "src", "data-ra226.js"), "w", encoding="utf-8") as f:
        f.write("window.RA226 = " + js + ";")
    print("написано: g1s_ra226_data.json и src/data-ra226.js (%d КБ)" % (len(js.encode("utf-8")) // 1024))
    notes = []
    if abs(m1x["A_E1_Bq"] / m1_result["A_Bq"] - 1.0) > 0.05:
        notes.append("метод 1: вторая мера E1 %.1f Бк против A2 %.1f Бк" % (m1x["A_E1_Bq"], m1_result["A_Bq"]))
    for tag in ("sel", "full"):
        if abs(variants[tag]["A_Bq"] / m1_result["A_Bq"] - 1.0) > 0.05:
            notes.append("метод 2 (%s) и метод 1 расходятся: %.3f против %.3f паспорта"
                         % (tag, variants[tag]["ratio_to_passport"], m1_result["A_Bq"] / A_pass))
        if variants[tag]["chi2_ndof"] > 5.0:
            notes.append("метод 2 (%s): χ²/ν = %.2f > 5" % (tag, variants[tag]["chi2_ndof"]))
    two = m1x["two"]
    notes.append("метод 1 двумя шаблонами (справочно): Ra226 %.3f, Rn222chain %.3f паспорта — группа Ra226 спектром "
                 "не определяется (дыра 26)" % (two["Ra226"] / A_pass, two["Rn222chain"] / A_pass))
    notes.append("bg_amplitude = 1: в критериях A2/A2V фон приведён по живому времени и не подгоняется")
    notes.append("xray_total_per_branch_pct = null (В3 = Р-А): панель метода 1 выводит его как «0,000 %»")
    print("\nТРЕБУЕТ ТОЛКОВАНИЯ:\n" + "\n".join(notes))


if __name__ == "__main__":
    main()
