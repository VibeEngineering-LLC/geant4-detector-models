Напиши ОДИН модуль Python 3 целиком. Верни ТОЛЬКО код, без объяснений и без
markdown-заборов. Комментарии и docstring — на РУССКОМ языке. Никаких голых
`except`, никаких `except: pass`. Любая ошибка входа — громкий отказ
`raise SystemExit("ОТКАЗ: ...")` с русским текстом. В JSON не должно попасть
ни одного абсолютного пути к файлу (только базовые имена файлов).

ИМЯ ФАЙЛА: export_amticseu_data.py (папка detectors/Gamma-1S/web-th232).

DOCSTRING МОДУЛЯ: выгрузка данных веб-страницы «Разложение спектра смеси
Am-241/Ti-44/Cs-137/Eu-152» (src/amticseu.html) по разбору 11.09.2026:
полная геометрия Гамма-1С с кюветой маринелли, источник по объёму пробы,
шаблоны npsm=off по 10⁷ распадов, шкала пробы — кусочно-линейная по реперам
таблицы пиков самого файла, шкала фона — рефит по 11 якорям фона (#CAL-0),
ПШПВ — измеренная кривая комплекта (16 точек), отклик «гаусс + левый хвост»
T = 0,75. Два пути: (1) подгонка полного нетто-спектра суммой шаблонов,
(2) нетто-площади пиков прибора против чистых площадей пиков модели. Метод 2
прежней страницы (библиотека линий × сетка отклика) на новой модели НЕ
пересчитан — сетка в полной геометрии ещё не посчитана; в выгрузку не входит.
Запуск: SPECTRAVIBE_ROOT=<...> python export_amticseu_data.py

ИМПОРТЫ И ПУТИ (ровно так):
    import datetime, json, math, os, sys
    import numpy as np
    import yaml
    HERE = os.path.dirname(os.path.abspath(__file__))
    GAMMA1S = os.path.dirname(HERE)
    REPO = os.path.dirname(os.path.dirname(GAMMA1S))
    sys.path.insert(0, os.path.join(GAMMA1S, "analysis"))
    if not os.environ.get("SPECTRAVIBE_ROOT"):
        raise SystemExit("ОТКАЗ: задайте SPECTRAVIBE_ROOT — каталог gamma-spectrum-analysis (ридер .spe)")
    import mix_unfold_g1s as g1s            # read_lsrm_spe, make_fwhm, RECAL_REFS
    import mix_unfold_core as core          # unfold, amplitude_errors, band_shares, zone, net_area, peak_row
    import bg_seven_line_anchor_check as bga  # measure_anchors

КОНСТАНТЫ:
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
    NET_WIN, NET_SIDE, BLEND_K = 1.25, (1.6, 3.0), 1.5
    Проверить существование каждого файла (SPE_*, FWHM_POINTS, CONFIG) и каждого
    шаблона; иначе SystemExit с именем отсутствующего.

ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ:
  fin(x): float(x), если конечно, иначе None (NaN/inf в JSON не пускать).
  rl(a, d=4): [round(float(v), d) for v in a].
  decay_factor(half_life_years, days): 0.5 ** (days / (half_life_years * 365.25)).
  days_between(d_ref, d_meas): (d_meas - d_ref).days, аргументы — datetime.date.
  model_grid(): np.arange(20.0, 3000.0, 1.0) — центры каналов равномерной сетки 1 кэВ
      для чистых площадей МОДЕЛИ (площадь не должна зависеть от шкалы записи).
  model_col(template_path, fwhm): hist, n, _ = g1s.read_template(template_path);
      e = model_grid(); col = g1s.broaden(hist, n, core.channel_edges(e), fwhm);
      вернуть (e, col, n)   (col — отсчёты на распад по каналам сетки).
  eff_model(template_path, fwhm, E0): e, col, _ = model_col(...);
      вернуть core.net_area(e, col, E0, fwhm(E0), NET_WIN, NET_SIDE)   (на распад).
  device_row(spec, E0): таблица spec.extras["lsrm_peaks_table"]; строка
      core.peak_row(table, E0, max(6.0, 0.6 * ПШПВ_кривой(E0))); None, если нет.
  blends(E0, fwhm, lib_lines, sum_peaks, own_key):
      все линии библиотеки конфига (любого нуклида) и суммы сумм-пиков конфига
      (e1_kev + e2_kev), у которых 0 < |E − E0| < BLEND_K * fwhm(E0); вернуть
      список dict {"E_keV", "nuclide", "kind": "линия" | "сумма"}.

ПОРЯДОК РАБОТЫ main():
 1. cfg = yaml.safe_load(open(CONFIG, encoding="utf-8")). Из него:
    passport (mass_g, passport_date, measured_date, components[key, bq_per_kg,
    unc_pct, half_life_years]), nuclides[key, label_ru, color], library.lines
    [e_kev, i_pct, nuclide], sum_peaks[e1_kev, e2_kev, nuclide], fit[e_lo_kev, e_hi_kev].
    g1s.TAIL_T = TAIL_T. fwhm = g1s.make_fwhm(FWHM_POINTS).
 2. Паспорт маринелли на дату измерения: days = days_between(date(passport_date),
    date(measured_date)); для каждого ключа A_Bq = bq_per_kg * mass_g/1000 *
    decay_factor(T½, days), dA_Bq = A_Bq * unc_pct / 100. Сохранить также
    Bq_per_kg, unc_pct, half_life_years, decay_factor, ref_date, meas_date.
 3. Шкала фона (#CAL-0): bgs = g1s.read_lsrm_spe(SPE_BG); coefs = list(bgs.energy_cal);
    если len(coefs) != 2 — SystemExit («шкала фона не линейная»). c0, c1 = coefs.
    ch_a, E_a, w_a = bga.measure_anchors(np.asarray(bgs.counts, float), c0, c1,
    fwhm, verbose=False). Взвешенный МНК: W = np.sqrt(w_a);
    M = np.vstack([np.ones_like(ch_a), ch_a]).T;
    cb, *_ = np.linalg.lstsq(M * W[:, None], E_a * W, rcond=None);
    resid = M @ cb − E_a; rms = sqrt(mean(resid²)).
    E_bg = cb[0] + cb[1] * np.arange(len(bgs.counts)).
 4. templates = [(k, os.path.join(BUILD_OUT, "mix_%s_npsmoff.csv" % k)) for k in KEYS].
    r  = core.unfold(SPE_SAMPLE, SPE_BG, templates, FWHM_POINTS, lo=fit.e_lo_kev,
         hi=fit.e_hi_kev, recalibrate=True, tail=TAIL_T, bg_energy_of_ch=E_bg, verbose=False)
    r0 = то же, но bg_energy_of_ch=None (фон канал в канал — так считал отчёт 11.09;
         нужен только для сравнения).
    Проверить r["names"] == KEYS, иначе SystemExit.
 5. Путь 1 (по r): err = core.amplitude_errors(r["A"], r["coef"]);
    birge = sqrt(r["chi2"] / r["ndof"]); для каждого i, key:
      A = r["activities"][i]; dA_stat = err[i] / r["live_s"]; dA = dA_stat * birge;
      groups[key] = {"A_Bq", "dA_Bq": fin(dA), "dA_stat_Bq": fin(dA_stat),
                     "A_over_passport": A / A_pass}.
    stack1[key] = rl(r["cols"][i] * r["coef"][i], 3).
    bands = core.band_shares(r); zn = core.zone(r, 50.0, 70.0); zn40 = core.zone(r, 40.0, 90.0).
    То же для r0 — только {"groups": {key: {"A_Bq", "A_over_passport"}},
    "chi2_ndof", "zone_ratio": core.zone(r0,50,70)["ratio"]}.
 6. Реперы шкалы пробы: spec = r["spec"]; table = spec.extras["lsrm_peaks_table"];
    для каждого E из g1s.RECAL_REFS: row = core.peak_row(table, E, max(6.0, 0.04 * E));
    если есть — {"E_true": E, "ch": row["position_ch"], "E_file": row["energy_keV"],
    "shift_keV": E − row["energy_keV"], "shift_fwhm": (E − row["energy_keV"]) / fwhm(E)}.
 7. Путь 2 (нетто-площади), запись маринелли: для каждого key и каждой линии
    E0 из [NET_MAIN[key]] + NET_EXTRA.get(key, []):
      row = device_row(spec, E0); если None — строку всё равно выдать с device=None;
      eff_meas = row["area"] / (A_pass[key] * r["live_s"]);
      eff_mod = eff_model(шаблон key, fwhm, E0);
      ratio = eff_meas / eff_mod; d_ratio_stat = ratio * row["d_area"] / row["area"];
      fwhm_ratio = row["fwhm_keV"] / fwhm(row["energy_keV"]);
      записать {"nuclide": key, "E_keV": E0, "main": E0 == NET_MAIN[key],
        "device": {"E_keV", "fwhm_keV", "area", "d_area"}, "fwhm_curve_keV": fwhm(row E),
        "fwhm_ratio", "eff_meas_pct": 100*eff_meas, "eff_model_pct": 100*eff_mod,
        "ratio", "d_ratio_stat", "blends": blends(E0, ...)}.
 8. Сверка по трём постановкам (разгадка Am-241): для каждой записи
      ("точечный 5 см", SPE_POINT_AM, "Am241", 59.541, "pt_Am241.csv"),
      ("точечный 5 см", SPE_POINT_CS, "Cs137chain", 661.657, "pt_Cs137.csv"),
      ("Петри-60",      SPE_PETRI,    "Am241", 59.541, "petri_Am241.csv"),
      ("Петри-60",      SPE_PETRI,    "Cs137chain", 661.657, "petri_Cs137chain.csv"),
      ("маринелли",     SPE_SAMPLE,   "Am241", 59.541, "mix_Am241_npsmoff.csv"),
      ("маринелли",     SPE_SAMPLE,   "Cs137chain", 661.657, "mix_Cs137chain_npsmoff.csv"):
      s = g1s.read_lsrm_spe(путь); row = device_row(s, E0) (None — SystemExit);
      активность на дату измерения s.start_datetime.date():
        точечные — из s.extras["lsrm_passport"][0]: value (Бк), reference_date
          ("ГГГГ-ММ-ДД"), uncertainty_pct; нуклид в строке паспорта обязан
          начинаться с "Am-241" или "Cs-137" в соответствии с key, иначе SystemExit;
        Петри — bq_per_kg из конфига для key × s.sample_mass_kg, дата отсчёта —
          passport_date конфига; маринелли — A_pass[key] из шага 2;
        T½ — из конфига для key.
      eff_meas = row["area"] / (A * s.live_time); eff_mod = eff_model(BUILD_OUT/шаблон, fwhm, E0);
      записать {"record", "nuclide", "E_keV", "A_Bq", "device": {...}, "fwhm_curve_keV",
        "fwhm_ratio", "eff_meas_pct", "eff_model_pct", "ratio": eff_meas/eff_mod,
        "template": базовое имя шаблона}.
 9. ПШПВ: points из FWHM_POINTS (E_keV, fwhm_keV, d_fwhm_keV, source);
    curve = [[E, fwhm(E)] for E in range(40, 3001, 10)];
    f662 = fwhm_keV точки с E_keV = 661.657; sqrt_law = [[E, f662*sqrt(E/661.657)] ...];
    device_peaks = для каждой строки таблицы пиков пробы: {"E_keV", "fwhm_keV",
    "fwhm_curve_keV": fwhm(E), "ratio"}.
10. Собрать data (все числа — fin/rl, никаких numpy-типов):
    "meta": detector "Гамма-1С (УДС-ГЦ-63х63)"; vessel "маринелли 1 л по чертежу,
      проба risn379, ρ = 1,0 г/см³, 1000 см³, источник по объёму пробы";
      model {"exe": "g1s_npsm", "physics": "G4EmStandardPhysics_option4, deex=deex,
      порог 0,05 мм", "npsm": "off", "templates": базовые имена, "decays": r["n_events"]};
      live_s, real_s (spec), bg_live_s, bg_real_s (bgs), bg_scale_time
      (spec.live_time/bgs.live_time), start_time str(spec.start_datetime),
      e_fit_lo, e_fit_hi, tail_T; "cal_sample": {"coefs_file": list(spec.energy_cal),
      "refs": шаг 6}; "cal_bg": {"coefs_file": [c0, c1], "coefs_refit": [cb0, cb1],
      "rms_keV": rms, "anchors": [{"E_keV", "ch", "netsum", "resid_keV"}]};
      "generated": datetime.date.today().isoformat().
    "passport", "nuclides" [{key, label_ru, color}] из конфига,
    "spectrum": {"e_of_ch": rl(r["e"], 3), "counts": rl(r["counts"], 1),
      "bg_counts": rl(r["bg_scaled"], 3), "stack1": stack1},
    "method1": {"groups", "chi2", "ndof", "chi2_ndof", "birge", "n_channels_fit":
      int(r["sel"].sum()), "bands": bands, "zone_50_70": zn, "zone_40_90": zn40}
      (из zone — только числа и списки: ratio, sum_meas, sum_model, model_shares_pct,
      peak_meas_keV, peak_model_keV, names),
    "method1_meta": {"template_decays": [{"nuclide": label_ru, "n": n}]},
    "method1_bg_by_channel": шаг 5 для r0,
    "netarea": {"definition": {"win_fwhm": NET_WIN, "side_fwhm": list(NET_SIDE),
      "grid_keV": 1.0, "blend_fwhm": BLEND_K}, "lines": шаг 7},
    "cross_geometry": шаг 8,
    "fwhm_cal": {"points", "curve", "sqrt_law", "f662_keV", "device_peaks"},
    "reference_lines": [[e_kev, nuclide] из library.lines, где nuclide in KEYS].
11. Записать HERE/g1s_amticseu_data.json (json.dump(..., ensure_ascii=False,
    separators=(",", ":"), allow_nan=False)) и HERE/src/data-amticseu.js ровно так:
    "window.AMTICSEU = " + тот же JSON-текст + ";\n" (encoding utf-8).
12. ПЕЧАТЬ (#CAL-0 и #SA-4), в таком порядке:
    - шкала пробы: коэффициенты файла; по каждому реперу E_true, E_file,
      сдвиг в кэВ и в долях ПШПВ;
    - шкала фона: коэффициенты файла, рефита, RMS и невязка каждого якоря;
    - путь 1: по ключу A, ±dA, отношение к паспорту; χ²/ν; зона 50–70: отношение
      модель/измерение и максимумы; то же отношение для варианта «фон канал в канал»;
    - путь 2: таблица линий; сверка по трём постановкам;
    - блок «ТРЕБУЕТ ТОЛКОВАНИЯ:» — строки для каждого сработавшего условия:
      (а) χ²/ν > 5; (б) у линии пути 2 есть бленды — перечислить;
      (в) |отношение пути 2 (главная линия) − отношение пути 1| / отношение пути 1 > 0,10;
      (г) fwhm_ratio вне [0,9; 1,1] у любой строки путей 2 и 8;
      (д) max |отношение(фон по шкале) − отношение(фон канал в канал)| > 0,005.
      Если не сработало ничего — строка «ТРЕБУЕТ ТОЛКОВАНИЯ: пусто».
    - «написано: <базовое имя> (N КБ)» для обоих файлов.

if __name__ == "__main__": main()
