"""Выгрузка данных страницы «Разложение спектра смеси Am-241/Ti-44/Cs-137/Eu-152»
(src/amticseu.html) по разбору 11.09.2026: полная геометрия с кюветой, шаблоны
npsm=off, шкала пробы по реперам файла, шкала фона — рефит по 11 якорям (#CAL-0),
измеренная кривая ПШПВ, левый хвост T = 0,75; путь 1 — полный спектр, путь 2 —
нетто-площади. Сгенерировано qwen3-coder:30b по _spec_export_amticseu.md, принято
с правками. Запуск: SPECTRAVIBE_ROOT=<...> python export_amticseu_data.py
"""
import datetime, json, math, os, sys
import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
GAMMA1S = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(GAMMA1S))
sys.path.insert(0, os.path.join(GAMMA1S, "analysis"))
if not os.environ.get("SPECTRAVIBE_ROOT"):
    raise SystemExit("ОТКАЗ: задайте SPECTRAVIBE_ROOT — каталог gamma-spectrum-analysis (ридер .spe)")
import mix_unfold_g1s as g1s
import mix_unfold_core as core
import bg_seven_line_anchor_check as bga
import glob
from scipy.optimize import nnls
os.environ.setdefault("G4MODELS_SOURCE_CONFIG", os.path.join(HERE, "configs", "amticseu.yaml"))
sys.path.insert(0, HERE)
import export_data as ed          # библиотека линий и сумм-пики с F_B — из конфига, тем же кодом, что прежний метод 2
import export_ra226_data as erd   # run_method2: F_B-депопуляция линий-партнёров и сумм-пики

KIT = os.path.join(GAMMA1S, "reference", "lsrm", "raw_lsrm", "Work", "BG",
                   "Gamma-1S", "Spe - поверки", "Поверка 2016")
SPE_SAMPLE   = os.path.join(KIT, "Маринелли", "Смесь_AmTiCsEu_Маринелли.spe")
SPE_BG       = os.path.join(KIT, "Фон вода", "фон вода_13.spe")
SPE_POINT_AM = os.path.join(KIT, "Точка 5см", "Am-241 42.13_Точечная-5см_5cm.spe")
SPE_POINT_CS = os.path.join(KIT, "Точка 5см", "Cs-137 #SRC-07_Точечная-5см_5cm.spe")
SPE_PETRI    = os.path.join(KIT, "Чашка Петри 60мл", "РИСН №SRC-04_Am-Ti-Eu-Cs_Петри-60.spe")
BUILD_OUT = os.path.join(os.environ.get("G4MODELS_BUILD_GAMMA_1S_NPSM",
                         os.path.join(REPO, "build", "Gamma-1S-npsm-1142")), "out")
FWHM_POINTS = os.path.join(HERE, "data", "fwhm_points_g1s_2016.csv")
CONFIG = os.path.join(HERE, "configs", "amticseu.yaml")
TAIL_T = 0.75
KEYS = ["Am241", "Ti44chain", "Cs137chain", "Eu152"]
NET_MAIN = {"Am241": 59.541, "Ti44chain": 1157.022, "Cs137chain": 661.657, "Eu152": 344.279}
NET_EXTRA = {"Ti44chain": [511.0], "Eu152": [121.7817, 1408.013]}
# W-077 (11.09.2026): прогон цепочки Ti-44→Sc-44 одним событием суммировал γ Ti-44
# с Sc-44 (T½ 4 ч; временного окна в g1s_npsm нет). Столбец Ti44chain —
# сумма раздельных прогонов Ti44only + Sc44only (scripts/run_mix_ti_sc_split.sh).
TEMPLATE_FILE = {"Ti44chain": "mix_Ti44split_npsmoff.csv"}


def tpl(k, suffix="npsmoff"):
    """Шаблон нуклида. suffix: npsmoff — прогон без непропорциональности (гистограмма отложенной
    энергии), npsmon — с ней (гистограмма СВЕТА, её переводит core.unfold(light_scale=...))."""
    name = TEMPLATE_FILE.get(k, "mix_%s_npsmoff.csv" % k).replace("npsmoff", suffix)
    return os.path.join(BUILD_OUT, name)
NET_WIN, NET_SIDE, BLEND_K = 1.25, (1.6, 3.0), 1.5

def fin(x):
    """float, если число конечно; иначе None (NaN/inf в JSON не пускать)."""
    if x is None:
        return None
    v = float(x)
    return v if math.isfinite(v) else None


def zjson(z):
    """Выжимка core.zone без numpy-массивов."""
    return {"ratio": fin(z["ratio"]), "sum_meas": fin(z["sum_meas"]), "sum_model": fin(z["sum_model"]),
            "model_shares_pct": rl(z["model_shares_pct"], 2), "peak_meas_keV": fin(z["peak_meas_keV"]),
            "peak_model_keV": fin(z["peak_model_keV"]), "names": list(z["names"])}

def rl(a, d=4):
    return [round(float(v), d) for v in a]

def decay_factor(half_life_years, days):
    return 0.5 ** (days / (half_life_years * 365.25))

def days_between(d_ref, d_meas):
    return (d_meas - d_ref).days

def model_grid():
    return np.arange(20.0, 3000.0, 1.0)

def model_col(template_path, fwhm):
    hist, n, _ = g1s.read_template(template_path)
    e = model_grid()
    col = g1s.broaden(hist, n, core.channel_edges(e), fwhm)
    return (e, col, n)

def eff_model(template_path, fwhm, E0):
    e, col, _ = model_col(template_path, fwhm)
    return core.net_area(e, col, E0, fwhm(E0), NET_WIN, NET_SIDE)

def device_row(spec, E0, fwhm):
    """Строка таблицы пиков прибора у E0 (допуск max(6 кэВ, 0,6·ПШПВ кривой))."""
    table = spec.extras.get("lsrm_peaks_table") or []
    return core.peak_row(table, E0, max(6.0, 0.6 * fwhm(E0)))


def dev(row):
    """Поля строки прибора, которые идут на страницу."""
    if row is None:
        return None
    return {"E_keV": fin(row["energy_keV"]), "fwhm_keV": fin(row["fwhm_keV"]),
            "area": fin(row["area"]), "d_area": fin(row["d_area"])}

def blends(E0, fwhm, lib_lines, sum_peaks, own_key):
    # Ближе 0,5 кэВ — это сама линия (в библиотеке 344,2785 против 344,279),
    # а не соседка; сумма, совпавшая с прямой линией, не дублируется.
    res, seen = [], set()
    cand = [(line["e_kev"], line["nuclide"], "линия") for line in lib_lines]
    cand += [(p["e1_kev"] + p["e2_kev"], p["nuclide"], "сумма") for p in sum_peaks]
    for E, nuc, kind in cand:
        if 0.5 < abs(E - E0) < BLEND_K * fwhm(E0) and (round(E), nuc) not in seen:
            seen.add((round(E), nuc))
            res.append({"E_keV": E, "nuclide": nuc, "kind": kind})
    return res

def grid_response(ch_edges, fwhm, light_scale=None, e=None, file_e=None, blur=1.0):
    """resp(E) метода 2 на ПРЯМЫХ прогонах сетки grid_mar_E*.csv (run_mix_grid.sh):
    без интерполяции и сдвига формы — энергия обязана совпасть с прогоном до 0,01 кэВ.
    Форма — свёртка измеренной кривой ПШПВ с левым хвостом (g1s.broaden), на квант.
    Эффективность пика — строго полное поглощение: бин 1 кэВ с E0 и соседний снизу
    (аннигиляционный квант оседает в 510,5); окно ±2,5 ПШПВ старого движка захватывало
    пик вылета иода и ближний комптон.
    При light_scale узлы строятся на прогонах gridon_E*.csv (npsm=on): их гистограмма — СВЕТ,
    ключи переводятся в энергию шкалы пробы, свёртка идёт в каналах (analysis/light_grid.py).
    Пиковая эффективность там берётся по НАЙДЕННОМУ положению пика, а не по номиналу кванта:
    после перевода полное поглощение стоит там, куда его ставит шкала."""
    if light_scale is not None:
        import light_grid, light_mode
        nodes_l, skipped, empty_l = light_grid.build_nodes(BUILD_OUT, "gridon_E", g1s.read_template,
                                                           light_mode.light_to_energy, g1s.broaden_ch,
                                                           e, file_e, fwhm, light_scale[0], light_scale[1], blur)
        print(light_grid.summary(nodes_l, empty_l))
        if skipped:
            print("  пропущено файлов (имя не разбирается в энергию): %s" % ", ".join(skipped))
        resp_l = light_grid.make_resp(nodes_l)
        # D-020 для метода 2 (13.09.2026, _spec_crit_bench_m2.md П2): узлы и число событий узла для
        # дисперсии столбцов и closure-стенда. Нечисловые имена (gridon_eplusg_* под регистронезависимым
        # глобом Windows) отбрасываются так же, как в light_grid.build_nodes.
        resp_l.nodes = {E: (v[0], v[1]) for E, v in nodes_l.items()}
        resp_l.n_map = {}
        for p in glob.glob(os.path.join(BUILD_OUT, "gridon_E*.csv")):
            try:
                E_n = float(os.path.basename(p)[len("gridon_E"):-len(".csv")])
            except ValueError:
                continue
            resp_l.n_map[E_n] = g1s.read_template(p)[1]
        return resp_l, len(nodes_l)
    nodes, n_map = {}, {}
    for p in glob.glob(os.path.join(BUILD_OUT, "grid_mar_E*.csv")):
        E0 = float(os.path.basename(p)[len("grid_mar_E"):-len(".csv")])
        hist, n, _ = g1s.read_template(p)
        n_map[E0] = n
        # Бин 0,5 кэВ — квант не отложил энергии (90–94 % событий сетки): это не регистрация.
        # С ним shape.sum() ≈ 1, и F_B-депопуляция run_method2 (eps_total партнёра) выедала
        # каскадные линии на 90–97 % (Ti 67,9/78,3, Eu 121,8/244,7/...; 11.09.2026).
        hist = {k: v for k, v in hist.items() if k >= 1.0}
        if sum(hist.values()) / n > 0.6:
            raise SystemExit("ОТКАЗ: полная эффективность %.3f на квант %.3f кэВ — бин нулевого отложения?"
                             % (sum(hist.values()) / n, E0))
        b = math.floor(E0) + 0.5
        nodes[E0] = (g1s.broaden(hist, n, ch_edges, fwhm), (hist.get(b, 0.0) + hist.get(b - 1.0, 0.0)) / n)
    if not nodes:
        raise SystemExit("ОТКАЗ: нет прогонов сетки grid_mar_E*.csv в %s" % BUILD_OUT)
    Es = np.array(sorted(nodes))

    def resp(E):
        j = int(np.argmin(np.abs(Es - E)))
        if abs(Es[j] - E) > 0.01:
            raise SystemExit("ОТКАЗ: нет прямого прогона сетки для %.3f кэВ (ближайший %.3f)" % (E, Es[j]))
        shape, eps = nodes[Es[j]]
        return shape, {}, eps
    resp.nodes, resp.n_map = nodes, n_map   # D-020 метод 2, _spec_crit_bench_m2.md П2
    return resp, len(Es)


# fit_cols (NNLS метода 2 весами 1/√max(y,1), прежний критерий) удалён 13.09.2026: после переноса метода 2
# на A2 (D-020 п.3, scripts/_spec_m2_transfer.md) вызовов не осталось; дублирующая подгонка — дефект (§33).


def main():
    for path in [SPE_SAMPLE, SPE_BG, SPE_POINT_AM, SPE_POINT_CS, SPE_PETRI, FWHM_POINTS, CONFIG]:
        if not os.path.exists(path):
            raise SystemExit("ОТКАЗ: отсутствует файл %s" % path)
    for key in KEYS:
        template = tpl(key)
        if not os.path.exists(template):
            raise SystemExit("ОТКАЗ: отсутствует шаблон %s" % template)

    cfg = yaml.safe_load(open(CONFIG, encoding="utf-8"))
    g1s.TAIL_T = TAIL_T
    fwhm = g1s.make_fwhm(FWHM_POINTS)
    passport = cfg["passport"]
    nuclides = cfg["nuclides"]
    library = cfg["library"]
    sum_peaks = cfg["sum_peaks"]
    fit = cfg["fit"]
    # #CH-1 (14.09.2026): сдвиг нумерации LSRM — единый источник в ядре; ключ конфига обязан совпадать, иначе отказ
    # (прежние чтения fit.get("lsrm_ch_offset", 0.0) при потере ключа молча ушли бы на 0).
    if "lsrm_ch_offset" not in fit or float(fit["lsrm_ch_offset"]) != g1s.LSRM_CH_OFFSET:
        raise SystemExit("ОТКАЗ: fit.lsrm_ch_offset=%r не равен ядру LSRM_CH_OFFSET=%r (#CH-1)"
                         % (fit.get("lsrm_ch_offset"), g1s.LSRM_CH_OFFSET))
    # Внешние реперы шкалы пробы ниже 59,5 кэВ (#AMT-3, analysis/low_scale_donors.py).
    low_refs = [(float(x["ch"]), float(x["e_kev"])) for x in (fit.get("low_refs") or [])] \
        if fit.get("use_low_refs") else []
    comp = {c["key"]: c for c in passport["components"]}

    # Паспорт маринелли
    d_ref = datetime.date.fromisoformat(passport["passport_date"])
    d_meas = datetime.date.fromisoformat(passport["measured_date"])
    days = days_between(d_ref, d_meas)
    A_pass = {}
    dA_pass, passport_out = {}, {}
    for key in KEYS:
        A_Bq = comp[key]["bq_per_kg"] * passport["mass_g"]/1000 * decay_factor(
            comp[key]["half_life_years"], days)
        dA_Bq = A_Bq * comp[key]["unc_pct"] / 100
        A_pass[key] = A_Bq
        dA_pass[key] = dA_Bq
        passport_out[key] = {"A_Bq": A_Bq, "dA_Bq": dA_Bq, "Bq_per_kg": comp[key]["bq_per_kg"],
                             "unc_pct": comp[key]["unc_pct"], "half_life_years": comp[key]["half_life_years"],
                             "decay_factor": decay_factor(comp[key]["half_life_years"], days),
                             "ref_date": passport["passport_date"], "meas_date": passport["measured_date"]}

    # Шкала фона
    bgs = g1s.read_lsrm_spe(SPE_BG)
    coefs = list(bgs.energy_cal)
    if len(coefs) != 2:
        raise SystemExit("ОТКАЗ: шкала фона не линейная")
    c0, c1 = coefs
    ch_a, E_a, w_a = bga.measure_anchors(np.asarray(bgs.counts, float), c0, c1,
                                        fwhm, verbose=False)
    W = np.sqrt(w_a)
    M = np.vstack([np.ones_like(ch_a), ch_a]).T
    cb, *_ = np.linalg.lstsq(M * W[:, None], E_a * W, rcond=None)
    resid = M @ cb - E_a
    rms = np.sqrt(np.mean(resid**2))
    E_bg = cb[0] + cb[1] * np.arange(len(bgs.counts))

    # Подгонка
    # Шкала в СВЕТЕ (12.09.2026). При fit.light_mode шаблоны берутся с суффиксом npsmon: их
    # гистограммы — свет в своих единицах, и core.unfold переводит ключи прямой «канал = a + b·свет»
    # (analysis/light_mode.py), а свёртка идёт в каналах с множителем blur.
    light_on = bool(fit.get("light_mode", False))
    suffix = str(fit.get("template_suffix", "npsmoff"))
    conv = str(fit.get("conv", "energy"))
    blur = float(fit.get("blur", 1.0))
    light_scale = None
    if light_on:
        if conv != "channel":
            raise SystemExit("ОТКАЗ: light_mode: true требует conv: channel, в конфиге %r" % conv)
        if suffix == "npsmoff":
            raise SystemExit("ОТКАЗ: light_mode: true с шаблонами npsmoff — в них света нет")
        import light_mode
        light_scale = light_mode.load_scale(os.path.join(
            os.path.dirname(os.path.abspath(__file__)), str(fit["light_scale_json"])))

    templates = [(k, tpl(k, suffix)) for k in KEYS]
    r = core.unfold(SPE_SAMPLE, SPE_BG, templates, FWHM_POINTS, lo=fit["e_lo_kev"],
         hi=fit["e_hi_kev"], recalibrate=True, extra_refs=low_refs, ch_offset=float(fit.get("lsrm_ch_offset", 0.0)), tail=TAIL_T, bg_energy_of_ch=E_bg, verbose=False,
         conv=conv, blur=blur, light_scale=light_scale)
    r0 = core.unfold(SPE_SAMPLE, SPE_BG, templates, FWHM_POINTS, lo=fit["e_lo_kev"],
         hi=fit["e_hi_kev"], recalibrate=True, extra_refs=low_refs, ch_offset=float(fit.get("lsrm_ch_offset", 0.0)), tail=TAIL_T, bg_energy_of_ch=None, verbose=False,
         conv=conv, blur=blur, light_scale=light_scale)

    if r["names"] != KEYS:
        raise SystemExit("ОТКАЗ: не совпадают ключи в результатах")

    # Погрешности даёт САМ критерий (D-020): A2 считает их с учётом дисперсии шаблонов.
    # Прежний amplitude_errors выводил их из матрицы плана в весах 1/√max(y,1) — весах другого
    # критерия; после переноса подгонки на A2 это было бы рассогласовано.
    err = r["sd"]
    # Множитель Бирге — от ЕДИНОЙ метрики (дисперсия A1), а не от χ² самого критерия: у A2 в
    # знаменателе стоит ещё и дисперсия шаблонов, и √(χ²/ν) там означает не то же самое.
    birge = np.sqrt(r["chi2_ref"] / r["ndof"])
    # Вторая оценка (E1) публикуется РЯДОМ, а не вместо: расхождение критериев по Am-241 —
    # факт о модели, а не повод выбрать удобное число (D-020, основание 3).
    groups_e1 = {}
    for i, key in enumerate(KEYS):
        A_e1 = r["e1"]["activities"][i]
        groups_e1[key] = {"A_Bq": fin(A_e1), "dA_Bq": fin(r["e1"]["sd"][i] / r["live_s"]),
                          "A_over_passport": fin(A_e1 / A_pass[key])}

    groups = {}
    stack1 = {}
    for i, key in enumerate(KEYS):
        A = r["activities"][i]
        dA_stat = err[i] / r["live_s"]
        dA = dA_stat * birge
        groups[key] = {"A_Bq": A, "dA_Bq": fin(dA), "dA_stat_Bq": fin(dA_stat),
                       "A_over_passport": A / A_pass[key]}
        stack1[key] = rl(r["cols"][i] * r["coef"][i], 3)

    # B-7: полосы χ² — от нижней границы окна подгонки (с 40 кэВ терялось 15,4 % χ²).
    bands = core.band_shares(r, bands=((fit["e_lo_kev"], 40), (40, 90), (90, 200), (200, 500),
                                       (500, 1000), (1000, fit["e_hi_kev"])))
    zn = core.zone(r, 50.0, 70.0)
    zn40 = core.zone(r, 40.0, 90.0)

    # Метод 2: линии библиотеки конфига × ПРЯМОЙ отклик сетки, сумм-пики с F_B
    file_e_sample = lambda c: float(r["spec"].channel_to_energy(c + float(fit.get("lsrm_ch_offset", 0.0))))
    resp, n_nodes = grid_response(r["ch_edges"], fwhm, light_scale=light_scale,
                                  e=r["e"], file_e=file_e_sample, blur=blur)
    # β⁺ Sc-44 (11.09.2026): пара 511 и γ 1157 одного распада — прогон primary=eplus_gamma
    # (scripts/run_mix_eplus.sh). Линия 511 (2·β⁺) → составная с долей β⁺; 1157 → остаток без β⁺.
    if light_scale is not None:
        # Пара β⁺ на световом прогоне: тот же перевод, что у сетки (иначе 511 и 1157 остались бы
        # на старой шкале и метод 2 оказался бы переведён наполовину).
        import light_grid as _lg
        nodes_p, _skipped_p, _empty_p = _lg.build_nodes(BUILD_OUT, "gridon_eplusg_E", g1s.read_template,
                                              light_mode.light_to_energy, g1s.broaden_ch,
                                              r["e"], file_e_sample, fwhm,
                                              light_scale[0], light_scale[1], blur)
        shp_p, eps_p, _sh, _nu, _ni = nodes_p[1157.022]
        pair = (shp_p, {}, eps_p)
        beta = next(I for E, I, nk, _ in ed.GAMMA_LIBRARY if nk == "Ti44chain" and abs(E - 511.0) < 0.01) / 2.0
    else:
        pg = os.path.join(BUILD_OUT, "grid_eplusg_E1157.022.csv")
        if not os.path.exists(pg):
            raise SystemExit("ОТКАЗ: нет прогона β⁺ Sc-44 %s" % pg)
        hp, npair, _ = g1s.read_template(pg)
        hp = {k: v for k, v in hp.items() if k >= 1.0}   # без событий без отложения (см. grid_response)
        pair = (g1s.broaden(hp, npair, r["ch_edges"], fwhm), {}, (hp.get(510.5, 0.0) + hp.get(511.5, 0.0)) / npair)
        beta = next(I for E, I, nk, _ in ed.GAMMA_LIBRARY if nk == "Ti44chain" and abs(E - 511.0) < 0.01) / 2.0
    lib_m2 = []
    for E, I, nk, note in ed.GAMMA_LIBRARY:
        if nk == "Ti44chain" and abs(E - 511.0) < 0.01:
            lib_m2.append((E, beta, nk, "β⁺ Sc-44: пара 511 кэВ и γ 1157 одного распада (прогон eplus_gamma)"))
        elif nk == "Ti44chain" and abs(E - 1157.022) < 0.01:
            lib_m2.append((E, I - beta, nk, "γ 1157 после захвата электрона (без β⁺)"))
        else:
            lib_m2.append((E, I, nk, note))
    resp_m2 = lambda E: pair if abs(E - 511.0) < 0.01 else resp(E)
    # D-020 п.3 (13.09.2026, scripts/_spec_m2_transfer.md, правка А): дисперсия столбцов метода 2 от шума
    # узлов сетки копится тем же вызовом. Принято стендом analysis/crit_bench_m2.py (closure с шумом узлов:
    # покрытие 68 % у A1 0,33–0,41 против A2V 0,61–0,70).
    pair_file = "gridon_eplusg_E1157.022.csv" if light_scale is not None else "grid_eplusg_E1157.022.csv"
    pair_n = g1s.read_template(os.path.join(BUILD_OUT, pair_file))[1]
    def n_of(E):
        Ek = min(resp.n_map, key=lambda x: abs(x - E))
        if abs(E - 511.0) >= 0.01 and abs(Ek - E) > 0.01:
            raise SystemExit("ОТКАЗ: нет числа событий узла %.3f кэВ" % E)
        return pair_n if abs(E - 511.0) < 0.01 else resp.n_map[Ek]
    var2 = {}
    _, by_nuc_w, lines_m2, n_sum = erd.run_method2(lib_m2, ed.SUM_PEAKS, resp_m2, r["e"], r["ch_edges"],
                                                   KEYS, var_acc=var2, n_of=n_of)
    cols2 = [by_nuc_w[k] * r["live_s"] for k in KEYS]
    dump_path = os.environ.get("G1S_M2_DUMP")
    if dump_path:  # D-020 метод 2 (_spec_crit_bench_m2.md П3; правка Д: by_nuc_w, var2, n_of, pair_n — из правки А)
        import pickle
        w_ref, acc = by_nuc_w, var2
        with open(dump_path, "wb") as fh:
            pickle.dump({"nodes": resp.nodes, "n_map": resp.n_map, "pair": (pair[0], pair[2]), "pair_n": pair_n,
                         "lib_m2": lib_m2, "sum_peaks": ed.SUM_PEAKS, "e": r["e"], "ch_edges": r["ch_edges"],
                         "sel": r["sel"], "keys": KEYS, "live_s": r["live_s"], "passport": A_pass,
                         "counts": np.asarray(r["spec"].counts, dtype=float), "bg_scaled": r["bg_scaled"],
                         "k": r["spec"].live_time / g1s.read_lsrm_spe(SPE_BG).live_time,
                         "cols2_ref": [w_ref[k] for k in KEYS], "V2_ref": [acc[k] for k in KEYS]}, fh)
        print("дамп метода 2: %s" % dump_path)
    # D-020 п.3 (правка Б): A2 с дисперсией узлов — импортом из принятого стенда; E1 — вторая мера.
    # Псевдоним crb, НЕ cb: в main уже есть локальная cb — коэффициенты перекалибровки ("coefs_refit").
    # Импорт под именем cb перезаписал её модулем и уронил выгрузку (W-101, 13.09.2026).
    import crit_bench as crb
    import crit_bench_m2 as cbm2
    k_bg = r["spec"].live_time / g1s.read_lsrm_spe(SPE_BG).live_time
    W2, V2 = np.array([by_nuc_w[k] for k in KEYS]), np.array([var2[k] for k in KEYS])
    cnt2, sel2, live = np.asarray(r["spec"].counts, dtype=float), r["sel"], r["live_s"]
    a2, sd2, _ = cbm2.fit_A2V(W2, cnt2, r["bg_scaled"], k_bg, V2, sel2)
    a2e, sd2e, ex2e = crb.fit_E1(W2, cnt2, r["bg_scaled"], k_bg, np.ones(len(KEYS)), sel2)
    met2 = crb.metrics_all(a2, W2, cnt2, r["bg_scaled"], k_bg, sel2, r["e"], KEYS, A_pass, live)
    met2e = crb.metrics_all(a2e, W2, cnt2, r["bg_scaled"], k_bg, sel2, r["e"], KEYS, A_pass, live)
    model2, net2 = W2.T @ a2, cnt2 - r["bg_scaled"]
    var_crit2 = np.maximum(cnt2[sel2] + k_bg * r["bg_scaled"][sel2], 1.0) + V2[:, sel2].T @ a2 ** 2
    chi2_2 = float(np.sum((model2[sel2] - net2[sel2]) ** 2 / var_crit2))
    ndof2 = int(sel2.sum()) - len(KEYS)
    birge2 = math.sqrt(met2["chi2_ref"] / ndof2)
    c2 = a2 / live
    m2_groups = {k: {"A_Bq": fin(c2[i]), "dA_Bq": fin(sd2[i] / live * birge2), "dA_stat_Bq": fin(sd2[i] / live),
                     "A_over_passport": fin(met2["A_over_passport"][k])} for i, k in enumerate(KEYS)}
    m2e_groups = {k: {"A_Bq": fin(a2e[i] / live), "dA_Bq": fin(sd2e[i] / live),
                      "A_over_passport": fin(met2e["A_over_passport"][k])} for i, k in enumerate(KEYS)}
    m2_stack = {k: rl(W2[i] * a2[i], 3) for i, k in enumerate(KEYS)}
    zmask = (r["e"] >= 50.0) & (r["e"] <= 70.0)
    m2_zone = fin(np.sum(model2[zmask]) / np.sum(r["y"][zmask]))
    for ln in lines_m2:
        ln["predicted_net"] = fin(ln["weight_per_branch"] * c2[KEYS.index(ln["nuclide"])] * r["live_s"])
        for f in ("eps_peak", "weight_per_branch", "depleted_pct"):
            if f in ln:
                ln[f] = fin(ln[f])

    groups0 = {}
    for i, key in enumerate(KEYS):
        A = r0["activities"][i]
        groups0[key] = {"A_Bq": fin(A), "A_over_passport": fin(A / A_pass[key])}
    chi2_ndof = r["chi2"] / r["ndof"]
    zone_ratio = core.zone(r0, 50.0, 70.0)["ratio"]

    # Реперы шкалы пробы
    spec = r["spec"]
    table = spec.extras["lsrm_peaks_table"]
    refs = []
    for E in g1s.RECAL_REFS:
        row = core.peak_row(table, E, max(6.0, 0.04 * E))
        if row is not None:
            refs.append({
                "E_true": E,
                "ch": row["position_ch"],
                "E_file": row["energy_keV"],
                "shift_keV": E - row["energy_keV"],
                "shift_fwhm": (E - row["energy_keV"]) / fwhm(E)
            })

    # Путь 2
    netarea_lines = []
    for key in KEYS:
        for E0 in [NET_MAIN[key]] + NET_EXTRA.get(key, []):
            row = device_row(spec, E0, fwhm)
            eff_meas = None
            if row is not None:
                eff_meas = row["area"] / (A_pass[key] * r["live_s"])
            eff_mod = eff_model(tpl(key), fwhm, E0)
            ratio = None
            d_ratio_stat = None
            if eff_meas is not None and eff_mod > 0:
                ratio = eff_meas / eff_mod
                d_ratio_stat = ratio * row["d_area"] / row["area"]
            netarea_lines.append({
                "nuclide": key,
                "E_keV": E0,
                "main": E0 == NET_MAIN[key],
                "device": dev(row),
                "fwhm_curve_keV": fin(fwhm(row["energy_keV"]) if row else fwhm(E0)),
                "fwhm_ratio": row["fwhm_keV"] / fwhm(row["energy_keV"]) if row else None,
                "eff_meas_pct": 100 * eff_meas if eff_meas is not None else None,
                "eff_model_pct": 100 * eff_mod,
                "ratio": ratio,
                "d_ratio_stat": d_ratio_stat,
                "blends": blends(E0, fwhm, library["lines"], sum_peaks, key)
            })

    # Сверка по трём постановкам
    cross_geometry = []
    for record, path, nuclide, E0, fname in [
        ("точечный 5 см", SPE_POINT_AM, "Am241", 59.541, "pt_Am241.csv"),
        ("точечный 5 см", SPE_POINT_CS, "Cs137chain", 661.657, "pt_Cs137.csv"),
        ("Петри-60",      SPE_PETRI,    "Am241", 59.541, "petri_Am241.csv"),
        ("Петри-60",      SPE_PETRI,    "Cs137chain", 661.657, "petri_Cs137chain.csv"),
        ("маринелли",     SPE_SAMPLE,   "Am241", 59.541, "mix_Am241_npsmoff.csv"),
        ("маринелли",     SPE_SAMPLE,   "Cs137chain", 661.657, "mix_Cs137chain_npsmoff.csv")
    ]:
        s = g1s.read_lsrm_spe(path)
        row = device_row(s, E0, fwhm)
        if row is None:
            raise SystemExit("ОТКАЗ: не найдена линия %s в %s" % (E0, path))
        t_meas = s.start_datetime.date()
        T_half = comp[nuclide]["half_life_years"]
        if record == "маринелли":
            A = A_pass[nuclide]
        elif record == "Петри-60":
            A = comp[nuclide]["bq_per_kg"] * s.sample_mass_kg * decay_factor(
                T_half, days_between(datetime.date.fromisoformat(passport["passport_date"]), t_meas))
        else:
            p = (s.extras.get("lsrm_passport") or [None])[0]
            want = {"Am241": "Am-241", "Cs137chain": "Cs-137"}[nuclide]
            if not p or not p.get("value") or not str(p.get("nuclide", "")).startswith(want):
                raise SystemExit("ОТКАЗ: паспорт %s не найден в %s" % (want, os.path.basename(path)))
            A = p["value"] * decay_factor(
                T_half, days_between(datetime.date.fromisoformat(p["reference_date"]), t_meas))
        eff_meas = row["area"] / (A * s.live_time)
        eff_mod = eff_model(os.path.join(BUILD_OUT, fname), fwhm, E0)
        cross_geometry.append({
            "record": record,
            "nuclide": nuclide,
            "E_keV": E0,
            "A_Bq": A,
            "device": dev(row),
            "fwhm_curve_keV": fin(fwhm(row["energy_keV"]) if row else fwhm(E0)),
            "fwhm_ratio": row["fwhm_keV"] / fwhm(row["energy_keV"]),
            "eff_meas_pct": 100 * eff_meas,
            "eff_model_pct": 100 * eff_mod,
            "ratio": eff_meas / eff_mod,
            "template": fname
        })

    # ПШПВ
    points = []
    with open(FWHM_POINTS, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("E_keV"):
                continue
            parts = line.strip().split(",")
            points.append({
                "E_keV": float(parts[0]),
                "fwhm_keV": float(parts[1]),
                "d_fwhm_keV": float(parts[2]) if len(parts) > 2 else None,
                "source": parts[3] if len(parts) > 3 else ""
            })
    curve = [[E, fwhm(E)] for E in range(40, 3001, 10)]
    f662 = [p["fwhm_keV"] for p in points if abs(p["E_keV"] - 661.657) < 0.1][0]
    sqrt_law = [[E, f662 * np.sqrt(E/661.657)] for E in range(40, 3001, 10)]
    device_peaks = []
    for row in table:
        device_peaks.append({
            "E_keV": row["energy_keV"],
            "fwhm_keV": row["fwhm_keV"],
            "fwhm_curve_keV": fwhm(row["energy_keV"]),
            "ratio": row["fwhm_keV"] / fwhm(row["energy_keV"])
        })

    # Сборка данных
    data = {
        "meta": {
            "detector": "Гамма-1С (УДС-ГЦ-63х63)",
            "vessel": "маринелли 1 л по чертежу, проба risn379, ρ = 1,0 г/см³, 1000 см³, источник по объёму пробы",
            "model": {
                "exe": "g1s_npsm",
                "physics": "G4EmStandardPhysics_option4, deex=deex, порог 0,05 мм",
                "npsm": "off",
                "templates": [os.path.basename(t) for _, t in templates],
                "decays": r["n_events"]
            },
            "live_s": r["live_s"],
            "real_s": float(spec.real_time),
            "bg_live_s": bgs.live_time,
            "bg_real_s": bgs.real_time,
            "bg_scale_time": spec.live_time / bgs.live_time,
            "start_time": str(spec.start_datetime),
            "e_fit_lo": fit["e_lo_kev"],
            "e_fit_hi": fit["e_hi_kev"],
            "tail_T": TAIL_T,
        "cal_sample": {
            "coefs_file": list(spec.energy_cal),
            "refs": refs
        },
        "cal_bg": {
            "coefs_file": [c0, c1],
            "coefs_refit": [cb[0], cb[1]],
            "rms_keV": rms,
            "anchors": [
                {"E_keV": fin(E_a[i]), "ch": fin(ch_a[i]), "netsum": fin(w_a[i]), "resid_keV": fin(resid[i])}
                for i in range(len(ch_a))
            ]
        },
        "generated": datetime.date.today().isoformat(),
        },
        "passport": passport_out,
        "nuclides": [{"key": n["key"], "label_ru": n["label_ru"], "color": n["color"]} for n in nuclides],
        "spectrum": {
            "e_of_ch": rl(r["e"], 3),
            "counts": rl(r["counts"], 1),
            "bg_counts": rl(r["bg_scaled"], 3),
            "stack1": stack1
        },
        "method1": {
            "groups": groups,
            "crit": r["crit"],
            "crit_note": "A2: нетто без обрезки, вес 1/√(p + k·b + Σaₖ²·varₖ), 6 итераций (D-020)",
            "chi2": r["chi2"],
            "ndof": r["ndof"],
            "chi2_ndof": chi2_ndof,
            "chi2_ref": r["chi2_ref"],
            "chi2_ref_ndof": r["chi2_ref"] / r["ndof"],
            "chi2_note": "chi2 — по дисперсии самого критерия (с дисперсией шаблонов); chi2_ref — "
                         "единая метрика сравнения критериев (дисперсия A1). Меньшее значение "
                         "первого НЕ означает лучшего согласия: у него шире знаменатель",
            "birge": birge,
            "n_channels_fit": int(r["sel"].sum()),
            "bands": bands,
            "zone_50_70": zjson(zn),
            "zone_40_90": zjson(zn40)
        },
        "method1_e1": {
            "groups": groups_e1,
            "crit": "E1",
            "crit_note": "E1 (вторая мера D-020): минимум полной вариации нормированных форм при "
                         "Σ модель = Σ нетто, точное ЛП; погрешности — параметрический бутстрап",
            "tv": r["e1"]["tv"],
            "chi2_ref": r["e1"]["chi2_ref"],
            "chi2_ref_ndof": r["e1"]["chi2_ref"] / r["ndof"]
        },
        "method1_meta": {
            "template_decays": [{"nuclide": {n["key"]: n["label_ru"] for n in nuclides}[k], "n": int(ne)} for k, ne in zip(KEYS, r["n_events"])]
        },
        "method1_bg_by_channel": {
            "groups": groups0,
            "chi2_ndof": r0["chi2"] / r0["ndof"],
            "zone_ratio": zone_ratio
        },
        "method2": {"groups": m2_groups, "crit": "A2",
                    "crit_note": "A2 (D-020 п.3): нетто без обрезки, вес 1/√(p + k·b + Σaₖ²·Vₖ), "
                                 "Vₖ = Σ w²·s/n по узлам сетки, 6 итераций",
                    "chi2": chi2_2, "ndof": ndof2, "chi2_ndof": chi2_2 / ndof2,
                    "chi2_ref": met2["chi2_ref"], "chi2_ref_ndof": met2["chi2_ref"] / ndof2,
                    "chi2_note": "chi2 — по дисперсии самого критерия (с дисперсией узлов сетки); chi2_ref — "
                                 "единая метрика (дисперсия A1), сопоставима с методом 1",
                    "birge": birge2,
                    "n_lines": len(ed.GAMMA_LIBRARY), "n_sum_peaks": n_sum, "n_nodes": n_nodes,
                    "zone_50_70_ratio": m2_zone, "stack": m2_stack, "lines": lines_m2},
        "method2_e1": {"groups": m2e_groups, "crit": "E1",
                       "crit_note": "E1 (вторая мера D-020): минимум полной вариации нормированных форм при "
                                    "Σ модель = Σ нетто, точное ЛП; погрешности — параметрический бутстрап",
                       "tv": ex2e["tv"], "chi2_ref": met2e["chi2_ref"],
                       "chi2_ref_ndof": met2e["chi2_ref"] / ndof2},
        "netarea": {
            "definition": {
                "win_fwhm": NET_WIN,
                "side_fwhm": list(NET_SIDE),
                "grid_keV": 1.0,
                "blend_fwhm": BLEND_K
            },
            "lines": netarea_lines
        },
        "cross_geometry": cross_geometry,
        "fwhm_cal": {
            "points": points,
            "curve": curve,
            "sqrt_law": sqrt_law,
            "f662_keV": f662,
            "device_peaks": device_peaks
        },
        "reference_lines": [[line["e_kev"], line["nuclide"]] for line in library["lines"] if line["nuclide"] in KEYS]
    }

    # Запись файлов
    with open(os.path.join(HERE, "g1s_amticseu_data.json"), "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    with open(os.path.join(HERE, "src/data-amticseu.js"), "w", encoding="utf-8") as f:
        f.write("window.AMTICSEU = " + json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False) + ";\n")

    # Печать
    print("#CAL-0")
    print("коэффициенты файла: %s" % list(spec.energy_cal))
    for ref in refs:
        print("E_true=%.3f E_file=%.3f shift_keV=%.3f shift_fwhm=%.3f" % (
            ref["E_true"], ref["E_file"], ref["shift_keV"], ref["shift_fwhm"]))
    print()
    print("#CAL-0 фон")
    print("коэффициенты файла: [%s, %s]" % (c0, c1))
    print("коэффициенты рефита: [%s, %s]" % (cb[0], cb[1]))
    print("RMS = %.3f кэВ" % rms)
    for anchor in data["meta"]["cal_bg"]["anchors"]:
        print("E=%.3f ch=%d resid=%.3f" % (anchor["E_keV"], anchor["ch"], anchor["resid_keV"]))
    print()
    print("путь 1")
    for key in KEYS:
        print("%s: A=%.1f ± %s Бк, к паспорту %.3f" % (
            key, groups[key]["A_Bq"], groups[key]["dA_Bq"], groups[key]["A_over_passport"]))
    print("χ²/ν = %.3f (дисперсия критерия A2); единая метрика χ²_ref/ν = %.3f"
          % (chi2_ndof, r["chi2_ref"] / r["ndof"]))
    print("вторая мера E1: полная вариация %.6f, χ²_ref/ν = %.3f; к паспорту %s"
          % (r["e1"]["tv"], r["e1"]["chi2_ref"] / r["ndof"],
             ", ".join("%s %.3f" % (k, groups_e1[k]["A_over_passport"]) for k in KEYS)))
    print("зона 50–70: модель/измерение = %.3f" % zn["ratio"])
    print("максимумы: измерение %.1f кэВ, модель %.1f кэВ" % (zn["peak_meas_keV"], zn["peak_model_keV"]))
    print("зона 50–70 (фон канал в канал): отношение = %.3f" % zone_ratio)
    print()
    print("метод 2 (линии × прямой отклик, узлов %d, сумм-пиков %d): χ²/ν = %.1f, зона 50–70 = %.3f"
          % (n_nodes, n_sum, chi2_2 / ndof2, m2_zone))
    for k in KEYS:
        print("  %s: A=%.1f ± %s Бк, к паспорту %.3f" % (k, m2_groups[k]["A_Bq"], m2_groups[k]["dA_Bq"],
                                                       m2_groups[k]["A_over_passport"]))
    print()
    print("путь 2")
    for line in netarea_lines:
        if line["device"]:
            print("%s %.3f: eff_meas=%.2f%% eff_model=%.2f%% ratio=%.3f fwhm_ratio=%.3f" % (
                line["nuclide"], line["E_keV"],
                line["eff_meas_pct"], line["eff_model_pct"], line["ratio"], line["fwhm_ratio"]))
        else:
            print("%s %.3f: device=None" % (line["nuclide"], line["E_keV"]))
    print()
    print("сверка по трём постановкам")
    for item in cross_geometry:
        print("%s %s %.3f: ratio=%.3f" % (
            item["record"], item["nuclide"], item["E_keV"], item["ratio"]))
    print()

    # Требует толкования
    issues = []
    if r["chi2_ref"] / r["ndof"] > 5:
        issues.append("χ²_ref/ν = %.1f > 5: погрешности пути 1 даны с множителем √(χ²_ref/ν)"
                      % (r["chi2_ref"] / r["ndof"]))
    issues.append("две величины χ²/ν: %.1f по дисперсии критерия A2 и %.1f по единой метрике — "
                  "разница от дисперсии шаблонов в знаменателе, не от согласия модели"
                  % (chi2_ndof, r["chi2_ref"] / r["ndof"]))
    issues.append("метод 2 и метод 1 считаются одним критерием A2 (D-020 п.3); дисперсия столбцов метода 2 — "
                  "от шума узлов сетки (Σ w²·s/n), проверено closure-стендом analysis/crit_bench_m2.py")
    issues.append("при шуме узлов критерий смещает амплитуды метода 2 до 0,1 %% (Ti −0,05…−0,06 %%, "
                  "Eu +0,08…+0,09 %%) — в 2–5 раз меньше статистической погрешности; механизм не установлен"
                  % ())
    for key in KEYS:
        x2, z2 = m2_groups[key]["A_over_passport"], m2e_groups[key]["A_over_passport"]
        if x2 and z2 is not None and abs(z2 - x2) / x2 > 0.05:
            issues.append("метод 2, %s: вторая мера E1 %.3f против A2 %.3f — расхождение %.0f %%"
                          % (key, z2, x2, 100 * abs(z2 - x2) / x2))
    for line in netarea_lines:
        if line["blends"]:
            issues.append("%s %.1f кэВ — бленд: %s" % (line["nuclide"], line["E_keV"], ", ".join("%.1f %s" % (b["E_keV"], b["nuclide"]) for b in line["blends"])))
    for key in KEYS:
        ratio1 = groups[key]["A_over_passport"]
        ln = [x for x in netarea_lines if x["nuclide"] == key and x["main"]][0]
        if ln["ratio"] is not None and abs(ln["ratio"] - ratio1) / ratio1 > 0.10:
            issues.append("%s: пути расходятся на %.0f %% (путь 1 %.3f, путь 2 %.3f)" % (
                key, 100 * abs(ln["ratio"] - ratio1) / ratio1, ratio1, ln["ratio"]))
        r2m = m2_groups[key]["A_over_passport"]
        if r2m is not None and abs(r2m - ratio1) / ratio1 > 0.10:
            issues.append("%s: метод 2 и путь 1 расходятся на %.0f %% (%.3f против %.3f)" % (
                key, 100 * abs(r2m - ratio1) / ratio1, r2m, ratio1))
    for item in netarea_lines + cross_geometry:
        if item["fwhm_ratio"] is not None and (item["fwhm_ratio"] < 0.9 or item["fwhm_ratio"] > 1.1):
            issues.append("ширина прибора / кривая = %.3f: %s %s %.1f кэВ" % (item["fwhm_ratio"], item.get("record", "маринелли"), item["nuclide"], item["E_keV"]))
    dmax = max(abs(groups[k]["A_over_passport"] - groups0[k]["A_over_passport"]) for k in KEYS)
    if dmax > 0.005:
        issues.append("приведение фона по его шкале меняет отношение к паспорту до %.3f" % dmax)
    issues = list(dict.fromkeys(issues))  # маринелли-Am есть и в пути 2, и в сверке постановок
    print("ТРЕБУЕТ ТОЛКОВАНИЯ: %s" % (", ".join(issues) if issues else "пусто"))
    print()
    size1 = os.path.getsize(os.path.join(HERE, "g1s_amticseu_data.json"))
    size2 = os.path.getsize(os.path.join(HERE, "src/data-amticseu.js"))
    print("написано: g1s_amticseu_data.json (%d КБ)" % (size1 // 1024))
    print("написано: src/data-amticseu.js (%d КБ)" % (size2 // 1024))

if __name__ == "__main__":
    main()
