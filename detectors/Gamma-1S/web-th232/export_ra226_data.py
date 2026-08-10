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
Переменные окружения: G4MODELS_BUILD_GAMMA_1S (сетка отклика), G4MODELS_
RA226_BG_SPE (приватный фон, в репозитории его нет).

Источник переведён 10.08.2026 на -18 (поверка 2016, оператор: "-19" в
чёрный список, дефект): СВОЯ матрица/плотность, ОТДЕЛЬНАЯ от Th-232 --
ОИСН-06 (насыпная эпоксидка, ро=0,60) вместо ОИСН-16 (ро=1,60). Сетка
метода 2 и МК-шаблоны метода 1 поэтому СВОИ, не общие с Th-232 (см.
GRID_TAG ниже) -- geometry/G1SDetector.cc, матрица "OISN06_epoxy",
ra226-remarks.md §14/15.
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

# Тег сетки метода 2 -- см. drivers/run_grid.py: tag = "<матрица><ро>",
# т.е. для матрицы OISN06_epoxy при ро=0,60 файлы grid/OISN06_epoxy0.60_
# E*.csv. НЕ "rho1.60" (умолчание модуля export_data.py) -- та сетка
# посчитана под ОИСН-16 источника -19, для -18 не годится (другая
# плотность и состав матрицы, см. докстринг модуля).
GRID_MATRIX = "OISN06_epoxy"
GRID_RHO = 0.60
GRID_PATTERN = "%s%.2f_E*.csv" % (GRID_MATRIX, GRID_RHO)

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


# ── НИЗКОЭНЕРГЕТИЧНАЯ ПОПРАВКА ШКАЛЫ (10.08.2026, оператор) ─────────────────
# Отдельная история от energy_correction ниже (та была ошибкой и снята) --
# ЭТА поправка живая, локальная, по одному внешне проверенному якорю.
#
# Находка: реальный пик К-рентгена Pb/Bi (канал ~28-30 у образца, ~29-31 у
# фона) по заводскому полиному попадал на 65-69 кэВ, а физическая энергия
# линии -- 74,8-90,9 кэВ (K-рентген Pb-214->Bi и Bi-214->Po, дублет+ветвь,
# ✅ ПРЯМОЙ запрос IAEA NDS 10.08.2026: nds.iaea.org/relnsd/v1/data?
# nuclides=214pb|214bi&rad_types=x -- цифры см. ra226-remarks.md §16).
# Интенсивностно-взвешенное среднее всей нерасщепимой на NaI группы --
# 79,08 кэВ (образец, смесь Pb214+Bi214); у фона (чистый Pb защиты, БЕЗ
# вклада Bi214 -- фон снят без источника) взвешенное среднее чище: 76,74
# кэВ. Оператор 10.08.2026: «при нашем разрешении фит не имеет смысла,
# просто прими 77» -- округление обоих чисел до общей опорной точки,
# обоснованное (76,74 и 79,08 оба близки к 77, расхождение <2,4 кэВ --
# меньше половины ПШПВ на этой энергии).
#
# ПОЧЕМУ НЕ ГЛОБАЛЬНЫЙ РЕФИТ. Прямая попытка перерешить квадратик ТОЛЬКО
# по (29,77)+(186/242/295/352) даёт чистые невязки на этих точках, но на
# 2204 кэВ уводит результат на +313 кэВ -- квадратика с одним новым нижним
# якорем не хватает степеней свободы одновременно сесть внизу и остаться
# верной наверху (проверено численно). Добавление точки 609,3 кэВ (канал
# 211, статистика надёжная) в тот же рефит делает квадратик снова хорошим
# ВПЛОТЬ ДО ~ канала 250 (невязки <1,3 кэВ на всех 6 якорях), но при
# экстраполяции дальше (туда, где статистика уже шумная -- 1764/2204 кэВ
# дают 306/73 отсчёта, ARGMAX там ненадёжен) старая и новая кривая
# расходятся: +80 кэВ на канале 598, +141 на 743, +306 на 1023. Верхнюю
# область НЕЗАЧЕМ и НЕЛЬЗЯ трогать без надёжного якоря там -- оставлена
# ЗАВОДСКОЙ, локальная поправка сшивается плавно и уходит в ноль к
# каналу 250.
_XRAY_ANCHOR_KEV = 77.0
_SAMPLE_XRAY_CH = 29        # реальный максимум пика (net-counts), не курсор
_SAMPLE_ANCHORS_CH = np.array([29, 67, 87, 105, 124, 211], dtype=float)
_SAMPLE_ANCHORS_KEV = np.array([77.0, 186.211, 241.995, 295.224,
                                351.932, 609.321], dtype=float)
_BLEND_LO, _BLEND_HI = 150.0, 250.0   # ниже BLEND_LO -- новая кривая,
                                     # выше BLEND_HI -- заводская, между -- плавный переход


def apply_xray_anchor_correction_sample(e_of_ch, ch):
    """Локальная поправка шкалы образца по одному внешнему якорю (см.
    комментарий выше). Возвращает НОВЫЙ массив энергий, не трогая ch."""
    A = np.vstack([_SAMPLE_ANCHORS_CH ** k for k in range(3)]).T
    coef, *_ = np.linalg.lstsq(A, _SAMPLE_ANCHORS_KEV, rcond=None)
    e_new = coef[0] + coef[1] * ch + coef[2] * ch ** 2

    # плавная сшивка (smoothstep, кубический Эрмит 0->1) между новой и
    # заводской кривой на [BLEND_LO, BLEND_HI]; вне интервала -- 0/1 ровно.
    t = np.clip((ch - _BLEND_LO) / (_BLEND_HI - _BLEND_LO), 0.0, 1.0)
    w_old = t * t * (3 - 2 * t)     # 0 при ch<=LO (чистая новая), 1 при ch>=HI (чистая старая)
    return (1.0 - w_old) * e_new + w_old * e_of_ch


def apply_xray_anchor_correction_bg(e_of_ch_bg, bg_energy_cal_coefs, xray_ch_bg):
    """Поправка шкалы фона: линейная калибровка (2 коэффициента) не даёт
    вторую степень свободы под один якорь -- сдвигаем ТОЛЬКО ноль (c0),
    наклон (кэВ/канал) заводской оставляем как есть. Физически это
    соответствует смещению нуля АЦП/усилителя, а не ошибке шкалы
    коэффициента усиления -- самый экономный вариант объяснения одной
    точкой. Сдвиг КОНСТАНТНЫЙ по всему диапазону (в отличие от квадратика
    образца, экстраполяция линии сама по себе не расходится)."""
    if len(bg_energy_cal_coefs) != 2:
        return e_of_ch_bg  # не линейная калибровка -- поправка не определена, не трогаем
    c0, c1 = bg_energy_cal_coefs
    e_old_at_anchor = c0 + c1 * xray_ch_bg
    shift = _XRAY_ANCHOR_KEV - e_old_at_anchor
    return e_of_ch_bg + shift


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
        add(nuc_key, w, shp)
        lines_out.append({"E_keV": Esum, "nuclide": nuc_key, "I_pct": None,
                          "note": note, "kind": "sum",
                          "E1_keV": E1, "E2_keV": E2,
                          "eps_peak": eps_sum_node, "weight_per_branch": w * eps_sum_node})
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
            "Запустить: cd %s && ./g1s.exe decay_ra226_isotopes.mac vessel %.2f %s"
            % (missing, ed.BUILD, GRID_RHO, GRID_MATRIX))

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


def main():
    meas, bg = read_pair()
    T = meas["live_s"]

    # ── энергетическая шкала образца: ЗАВОДСКАЯ + ЛОКАЛЬНАЯ поправка низа.
    #
    # ИСТОРИЯ ГЛОБАЛЬНОЙ поправки (важна, чтобы не наступить второй раз --
    # касается файла -19, ρ=1,60, снят с конвейера 10.08.2026). Тогда сюда
    # была добавлена линейная поправка E = a + b·E_заводское, снятая по
    # центроидам пяти чистых линий из ed.fit_peak_multiplet. Она оказалась
    # ОШИБКОЙ и удалена в тот же день: поправка растягивала шкалу на 0,74%
    # (+4,9 кэВ на 609 кэВ, +16,7 кэВ на 2204 кэВ) и уводила пики с их
    # мест на графике -- замечание оператора «этот пик не на месте»,
    # «уехал влево почти на 200 кэВ». Причина: центроиды fit_peak_multiplet
    # систематически занижены (модель окна -- гауссиана + ПРЯМАЯ подложка,
    # а под пиком реально комптоновская СТУПЕНЬ). Артефакт подгонки был
    # принят за смещение шкалы. На файле -19 заводская шкала была тогда
    # проверена тремя независимыми способами (СпектраЛайн, SpectraVibe
    # seven_line_check, гладкость полинома) и оказалась верна ЦЕЛИКОМ, без
    # поправок -- см. git-историю (коммит aae647b) для деталей на -19.
    #
    # Для ТЕКУЩЕГО файла -18 (2016) описанная выше ГЛОБАЛЬНАЯ трёхкратная
    # проверка НЕ ПОВТОРЯЛАСЬ -- сравнение реперных линий 186-609 кэВ по
    # argmax сошлось в пределах 1-2 каналов (заводская верна), а вот
    # К-рентген-пик (канал ~29, физика 74,8-90,9 кэВ) заводская давала на
    # 65-69 кэВ -- разрыв ~10-14 кэВ, устойчивый к методу (не артефакт
    # подгонки, проверялось: гаусс+линия, гаусс+экспонента, центр масс
    # разными окнами -- всё в одном месте, и ВАЖНО: тот же центр-масс метод
    # на заведомо верных 186-2204 кэВ точках сам давал ложный уход
    # 3-15 кэВ, поэтому решающим стал argmax + прямая проверка IAEA NDS,
    # не подгонка). Оператор 10.08.2026: «на заводе не всегда идеально
    # делают... поправь калибровку с учётом этого пика... просто прими 77»
    # -- см. apply_xray_anchor_correction_sample/_bg выше, ra226-remarks.md
    # §16 (полный разбор, численные пробы, обоснование локальности).
    e = apply_xray_anchor_correction_sample(meas["e_of_ch"],
                                            np.arange(meas["n_channels"],
                                                     dtype=float))
    bg["e_of_ch"] = apply_xray_anchor_correction_bg(
        bg["e_of_ch"], bg["coefs"], 30)  # канал 30 -- argmax пика Pb-ХРИ фона
    energy_correction = None

    bg_on_meas = np.interp(e, bg["e_of_ch"], bg["counts"].astype(float),
                           left=0.0, right=0.0)
    bg_scale_time = T / bg["live_s"]
    bg_scaled = bg_on_meas * bg_scale_time

    ch_edges = np.concatenate((
        [e[0] - 0.5 * (e[1] - e[0])],
        0.5 * (e[:-1] + e[1:]),
        [e[-1] + 0.5 * (e[-1] - e[-2])],
    ))

    # ── ПШПВ(E): ЗАВОДСКАЯ калибровка прибора, не собственная деконволюция
    # (директива оператора 10.08.2026 «сам не калибруй» -- та же логика,
    # что и для энергетической шкалы выше: fit_peak_multiplet на широких
    # NaI-пиках систематически смещает результат, см. историю оси E).
    # `factory_fwhm_keV`/`fit_power_law_to_factory_fwhm` берут ГОТОВЫЙ
    # полином из шапки .spe (LSRM «Алгоритмические основы» §8.3), никакой
    # подгонки по измеренному спектру. Перевод в степенной закон k·E^p --
    # чисто техническая аппроксимация ФОРМЫ под существующий интерфейс
    # ed.FWHM_LAW, не калибровка: fit_power_law_rms_pct -- невязка этой
    # аппроксимации к самой заводской функции (обычно <1%), НЕ невязка
    # заводской калибровки к реальности.
    if not meas.get("fwhm_coefs"):
        raise SystemExit(
            "в шапке .spe нет заводской ПШПВ-калибровки "
            "(stored_fwhm_calibration пуст) -- посчитать нечем.")
    fwhm_k, fwhm_p, fwhm_fit_rms_pct = fit_power_law_to_factory_fwhm(
        meas["fwhm_coefs"], meas["fwhm_model"])
    ed.FWHM_LAW.update({"kind": "power", "k": fwhm_k, "p": fwhm_p})
    _ref_pts = (186.211, 351.932, 609.321, 661.657, 1764.491, 2204.100, 2447.7)
    fwhm_cal = {
        "source": "заводская (LSRM, полином в шапке .spe, z=sqrt(E))",
        "coefs": meas["fwhm_coefs"], "model": meas["fwhm_model"],
        "k": fwhm_k, "p": fwhm_p, "fit_rms_pct": fwhm_fit_rms_pct,
        "fwhm662_law": fwhm_k * 661.657 ** fwhm_p,
        "fwhm662_cs": ed.FWHM662,
        "res662_pct": 100.0 * fwhm_k * 661.657 ** fwhm_p / 661.657,
        "reference_points": [
            {"E_keV": E,
             "fwhm_factory_keV": factory_fwhm_keV(meas["fwhm_coefs"],
                                                  meas["fwhm_model"], E),
             "fwhm_power_law_keV": fwhm_k * E ** fwhm_p}
            for E in _ref_pts],
    }
    print("ПШПВ: заводская калибровка, аппроксимация степенным законом "
          "k=%.4f p=%.4f (невязка аппроксимации %.2f%%), ПШПВ(662)=%.1f кэВ"
          % (fwhm_k, fwhm_p, fwhm_fit_rms_pct, fwhm_cal["fwhm662_law"]))

    eps_peak = ed.make_eps_peak_interp(os.path.join(ed.BUILD, "grid"),
                                       GRID_PATTERN)
    resp = ed.make_full_response(os.path.join(ed.BUILD, "grid"), ch_edges,
                                 True, eps_peak, GRID_PATTERN)

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
            # ИСПРАВЛЕНО 10.08.2026 (замечание оператора «сортировка по
            # вкладу»): было ln.get("I_pct") -- копия интенсивности без
            # смысла (не отсчёты, не Бк), заглушка так и осталась не
            # заполненной. Настоящий предсказанный вклад строки в счёт
            # спектра -- weight_per_branch (площадь отклика на распад
            # ветви, из run_method2) x амплитуда x живое время.
            ln["predicted_net"] = ln.get("weight_per_branch", 0.0) * A_Bq * T
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
        # ПЕРЕСЧИТАНО 10.08.2026 вечером (источник -18, после локальной
        # поправки низкоэнергетичной шкалы, см. §16 remarks) -- ПРЕЖНИЙ
        # текст ("почти вдвое") относился к файлу -19 и после смены
        # источника стал ФАКТИЧЕСКИ НЕВЕРЕН (число здесь -- живое,
        # A_par_over_pass -- НЕ хардкод; но подстрока-объяснение ниже
        # НАПИСАНА руками под текущий прогон -- если A_par_over_pass ещё
        # раз всерьёз изменится, переписать текст заново, не оставлять
        # рассинхрон, как уже было один раз).
        "caveat": "Метод -- деконволюция ПО ВСЕМУ спектру двумя независимыми "
                  "амплитудами (родитель Ra-226 против дочерних Rn-222+Pb-214"
                  "+Bi-214+Po-218/214), теми же базисами, что и single-"
                  "амплитудный метод 2 выше. Проверка устойчивости к границе "
                  "окна подгонки (E_FIT_HI, вручную 10.08.2026 вечером, "
                  "источник -18): 400 кэВ -> A_par/паспорт=1,05, отношение="
                  "0,92; 700 -> 1,33 / 0,77; 1200 -> 1,33 / 0,77; 2300 (рабочее"
                  " значение) -> 1,30 / 0,81. Отношение дочерние/родитель "
                  "ВСЮДУ <1 (направление СОГЛАСУЕТСЯ с гипотезой утечки), но "
                  "величина ЗАВИСИТ от ширины окна (1,05...1,33) -- сама эта "
                  "зависимость означает остаточную систематику модели формы, "
                  "не чистое измерение. Мягче прежнего вывода на файле -19 "
                  "(там было 'почти вдвое', здесь 1,05-1,33) -- возможно, "
                  "частично снято правкой низкоэнергетичной шкалы (родитель "
                  "определяется в т.ч. континуумом К-рентгена Ra-226 в "
                  "низкоэнергетичной части, которая перекалибрована §16). "
                  "Гипотеза U-235 185,715 кэВ как источник избытка на "
                  "родительской линии ОТКЛОНЕНА оператором 10.08.2026: "
                  "«в нашем КИ этого нет» (паспорт контрольного источника "
                  "U-235 не декларирует) -- причина избытка на 186,211 кэВ "
                  "ОТКРЫТА, не заменена новой недоказанной гипотезой. "
                  "Оператор отдельно указал СТАНДАРТНУЮ практику: активность "
                  "Ra-226 по линии 186,211 кэВ (I=3,565%, самая слабая в "
                  "модели, самопоглощение/интерференции) НЕ ОПРЕДЕЛЯЮТ -- т.е. "
                  "эта проверка изначально построена на линии, которую "
                  "профессионально считают непригодной для количественной "
                  "оценки; отклонение теста именно на ней ожидаемо по этой "
                  "причине, не обязательно физическая утечка радона. Число "
                  "посчитано и приводится, но НЕ "
                  "читается как подтверждённая утечка радона -- "
                  "reliable=false, нужна более узкая модель (3+ амплитуды, "
                  "разделяющая Pb214/Bi214 отдельно от Rn222/Po218/Po214) "
                  "прежде чем публиковать как факт.",
    }
    print("проверка утечки радона: A(родитель)=%.0f+-%.0f Бк  "
          "A(дочерние)=%.0f+-%.0f Бк  отношение=%.3f+-%.3f  chi2/ndof=%.2f  "
          "cond=%.2e  A_par/паспорт=%.2f -- ЧИСЛО НЕ НАДЁЖНО, см. caveat"
          % (A_par, dA_par, A_dtr, dA_dtr, ratio_rn, d_ratio_rn,
             chi2_rn / ndof_rn, cond_rn, A_par_over_pass))

    # ── 3-амплитудный разрез (10.08.2026, оператор -- визуально Bi-214
    # выше данных, Pb-214 ниже на графике метода 1; "может Ra226 нужно
    # добавить" -- отдельная группа, не смешивать с Pb214/Bi214) ──────────
    # ГРУППЫ: {Ra226+Rn222+Po218+Po214} (хвост -- 0,04-2,2% модели каждый,
    # см. таблицу метода 1 -- фактически несёт только линию 186,211 кэВ)
    # | Pb214 (отдельно) | Bi214 (отдельно). Прямая проверка: РЕАЛЬНО ли
    # Pb214 занижена, а Bi214 завышена относительно друг друга, а не
    # просто визуальное впечатление от локального шума.
    G3 = {"rest": ["Ra226", "Rn222", "Po218", "Po214"],
         "Pb214": ["Pb214"], "Bi214": ["Bi214"]}
    shapes3 = {}
    for gname, gkeys in G3.items():
        s = np.zeros_like(e)
        for k in gkeys:
            s += by_nuc_w_sel[k]
        shapes3[gname] = s
    cols3 = [shapes3["rest"][sel] * T, shapes3["Pb214"][sel] * T,
            shapes3["Bi214"][sel] * T, bgm]
    coef3, dcoef3, chi2_3, ndof_3, model_3 = ed.fit_amplitudes(
        y_sel, cols3, ed.SYS_FLOOR)
    A_rest3, A_pb3, A_bi3 = (float(coef3[0]), float(coef3[1]),
                             float(coef3[2]))
    dA_rest3, dA_pb3, dA_bi3 = (float(dcoef3[0]), float(dcoef3[1]),
                                float(dcoef3[2]))
    print("3-амплитудный разрез: rest(Ra226+..)=%.0f+-%.0f  Pb214=%.0f+-%.0f "
          " Bi214=%.0f+-%.0f Бк  Pb214/паспорт=%.2f  Bi214/паспорт=%.2f  "
          "chi2/ndof=%.2f -- СПРАВОЧНО, та же оговорка о неоднозначности, "
          "что у 2-амплитудного разреза выше"
          % (A_rest3, dA_rest3, A_pb3, dA_pb3, A_bi3, dA_bi3,
             A_pb3 / A_pass, A_bi3 / A_pass, chi2_3 / ndof_3))
    split3_check = {
        "attempted": True, "method": "method2_split3amp",
        "groups": {"rest_Ra226_Rn222_Po218_Po214": {
                       "A_Bq": A_rest3, "dA_Bq": dA_rest3},
                  "Pb214": {"A_Bq": A_pb3, "dA_Bq": dA_pb3,
                            "A_over_passport": A_pb3 / A_pass},
                  "Bi214": {"A_Bq": A_bi3, "dA_Bq": dA_bi3,
                            "A_over_passport": A_bi3 / A_pass}},
        "chi2_ndof": chi2_3 / ndof_3,
        "reliable": False,
        "caveat": "3-амплитудный разрез (10.08.2026 вечером) -- прямая "
                  "проверка визуального замечания оператора (Bi-214 чуть "
                  "выше данных, Pb-214 чуть ниже на графике метода 1). "
                  "Та же методологическая оговорка, что у 2-амплитудного "
                  "разреза выше: секулярное равновесие цепочки ПРЕДПОЛАГАЕТСЯ "
                  "одноамплитудной моделью метода 1/2 по умолчанию, здесь "
                  "оно НАМЕРЕННО разорвано на три независимых числа, чтобы "
                  "увидеть расхождение -- само расхождение (если есть) "
                  "может быть как реальной физикой (разная самопоглощение "
                  "линий разных энергий, недостаточная статистика МК-шаблона "
                  "Pb214 -- 200 тыс. распадов против 2 млн у Bi214), так и "
                  "артефактом регуляризации/обусловленности совместного "
                  "3-амплитудного фита. НЕ читать как подтверждённый вывод "
                  "без дополнительной проверки устойчивости к окну "
                  "(аналогично 400/700/1200/2300 кэВ выше).",
    }

    # ── ДПР и Ra-226 -- ДВЕ НЕЗАВИСИМЫЕ СВОБОДНЫЕ амплитуды по форме
    # спектра, СРАВНИВАЮТСЯ друг с другом и с паспортом (10.08.2026
    # вечером, три захода оператора до верной формулировки: "подгоняй
    # радий под сумму" -> "не так, паспорт получен ПО ДПР, не по радию,
    # закреплять надо ДПР" -> "подгоняй отдельно активность ДПР и
    # активность радия по форме спектра и сравнивай"). Ключевой факт из
    # sci-search (§19 remarks): активность Ra-226 в гамма-спектрометрии
    # НЕ измеряется по его собственной линии 186,211 кэВ вообще (та
    # линия неразделимо накладывается на U-235 185,715 кэВ, разница
    # 0,5 кэВ, не берёт даже HPGe) -- значит и ПАСПОРТНАЯ активность
    # источника традиционно получена ПО ДПР (Pb-214/Bi-214) или прямым
    # альфа-счётом, но НЕ по линии Ra-226. Поэтому предыдущая версия
    # теста (фиксация "родителя" на паспорте) была КОНЦЕПТУАЛЬНО
    # НЕВЕРНОЙ -- паспорт логичнее сопоставлять с ДПР-стороной, а Ra-226
    # мерить НЕЗАВИСИМО и смотреть, совпадает ли он с обоими.
    #
    # Здесь -- НИЧЕГО не фиксируется. Две группы, обе свободны:
    #   ДПР   = Pb214 + Bi214 (те самые линии, что реально используют)
    #   Ra226 = Ra226 (+ Rn222/Po218/Po214 -- гамма практически нет,
    #           0,00-0,04% модели, физически неотделимы от Ra226 по
    #           спектру, включены для полноты баланса, не меняют число)
    shape_dpr = by_nuc_w_sel["Pb214"] + by_nuc_w_sel["Bi214"]
    shape_ra = (by_nuc_w_sel["Ra226"] + by_nuc_w_sel["Rn222"]
               + by_nuc_w_sel["Po218"] + by_nuc_w_sel["Po214"])
    cols_2 = [shape_dpr[sel] * T, shape_ra[sel] * T, bgm]
    coef_2, dcoef_2, chi2_2, ndof_2, model_2 = ed.fit_amplitudes(
        y_sel, cols_2, ed.SYS_FLOOR)
    A_dpr, dA_dpr = float(coef_2[0]), float(dcoef_2[0])
    A_ra, dA_ra = float(coef_2[1]), float(dcoef_2[1])
    ratio_ra_dpr = A_ra / A_dpr if A_dpr > 0 else float("nan")
    print("ДПР против Ra-226 (обе свободно, по форме): ДПР=%.0f+-%.0f Бк "
          "(%.2fx паспорта)  Ra-226=%.0f+-%.0f Бк (%.2fx паспорта)  "
          "Ra226/ДПР=%.3f  chi2/ndof=%.2f -- СПРАВОЧНО"
          % (A_dpr, dA_dpr, A_dpr / A_pass, A_ra, dA_ra, A_ra / A_pass,
             ratio_ra_dpr, chi2_2 / ndof_2))
    radon_dpr_vs_ra_check = {
        "attempted": True, "method": "dpr_vs_ra226_both_free",
        "A_DPR_Bq": A_dpr, "dA_DPR_Bq": dA_dpr,
        "A_DPR_over_passport": A_dpr / A_pass,
        "A_Ra226_Bq": A_ra, "dA_Ra226_Bq": dA_ra,
        "A_Ra226_over_passport": A_ra / A_pass,
        "ratio_Ra226_to_DPR": ratio_ra_dpr,
        "chi2_ndof": chi2_2 / ndof_2,
        "reliable": False,
        "caveat": "Обе амплитуды СВОБОДНЫ (по форме спектра, метод 2), "
                  "НИЧЕГО не зафиксировано на паспорте -- ДПР (Pb-214+"
                  "Bi-214, линии, реально используемые для активности "
                  "Ra-226 по стандартной практике, см. §19 remarks) и "
                  "Ra-226 (+ Rn222/Po218/Po214, гамма практически нет) "
                  "измерены НЕЗАВИСИМО друг от друга. Результат: ДПР=%.2f "
                  "паспорта, Ra-226=%.2f паспорта -- ДПР близко к 1 "
                  "(в пределах ожидаемого для метода), Ra-226 заметно "
                  "выше. Поскольку паспортная активность источника "
                  "традиционно получена ПО ДПР или альфа-счётом, а не по "
                  "линии Ra-226 (та профессионально не используется, §19),"
                  " согласие ДПР~паспорт и рассогласование Ra-226~паспорт "
                  "СОГЛАСУЕТСЯ с гипотезой (б): проблема именно в "
                  "измерении Ra-226 через 186,211 кэВ в этой модели, не в "
                  "балансе цепочки. Утечка радона (гипотеза а) не имеет "
                  "предсказательной силы для ЭТОГО расхождения: утечка "
                  "занизила бы ДПР ОТНОСИТЕЛЬНО Ra-226, а здесь "
                  "расходится ровно наоборот -- Ra-226 читается высоко, "
                  "не ДПР низко. Не проверено независимым способом "
                  "(нужна проверка чувствительности эффективности "
                  "186,211 кэВ к неопределённости состава/плотности новой "
                  "матрицы) -- reliable=false."
                  % (A_dpr / A_pass, A_ra / A_pass),
    }

    palette = {n["key"]: n["color"] for n in ed._CFG["nuclides"]}
    label_ru = {n["key"]: n["label_ru"] for n in ed._CFG["nuclides"]}

    data = {
        "meta": {
            "detector": "Гамма-1С (УДС-ГЦ-63х63)",
            "vessel": "Маринелли 1 л, ОИСН-06 ρ=0,60 г/см³ (насыпная"
                     " эпоксидка, источник -18/2016; СВОЙ Geant4-прогон"
                     " 10.08.2026, НЕ переиспользован от Th-232/ОИСН-16)",
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
        "split3_check": split3_check,
        "radon_dpr_vs_ra_check": radon_dpr_vs_ra_check,
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
