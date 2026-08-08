"""Выгрузка данных для страницы «Разложение спектра Th-232 в Маринелли».

Разложение — по индивидуальным нуклидам ветви Th-232 (директива оператора
08.08.2026): Th-232, Ac-228, Th-228, Ra-224, Rn-220, Pb-212, Bi-212, Tl-208.
Каждый нуклид — отдельный МК-прогон полного распада в той же геометрии
(macros/decay_th232_isotopes.mac). В вековом равновесии активность каждого
звена равна активности родителя, поэтому сумма восьми шаблонов даёт полный
цепочечный шаблон.

Выход `g1s_th232_data.json`:
- measurement / background:  каналы, отсчёты, времена, калибровка
- template.total:            сумма индивидуальных, на 1 распад ветви
- template.by_nuc:            каждый нуклид отдельно (stacked area)
- method1:                    МК-шаблоны по нуклидам — NNLS подгонка
                              полного спектра (амплитуда ветви + фон)
- method2:                    функция пика полного поглощения (сетка
                              grid/rho1.60_E*.csv) × библиотека выходов
                              γ-линий эмиссии + сумм-пики каскадного
                              совпадения — МНК по ансамблю измеренных
                              площадей (директива оператора 08.08.2026,
                              документ wt20-methods-compare.md)
- passport, meta:             как раньше
"""
import glob
import json
import math
import os
import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

import numpy as np
from scipy.optimize import curve_fit, nnls

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
REPO = os.path.abspath(os.path.join(ROOT, "..", ".."))

sys.path.insert(0, os.path.join(REPO, "common", "py"))
import paths  # noqa: E402

sys.path.insert(0, os.path.join(ROOT, "analysis"))
from detector_params import FWHM662, ESCAPE_KEV  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))
KIT = str(paths.ref("Gamma-1S"))

XML_MEAS = os.path.join(
    KIT, "reference_spectra", "reference_kits_becqmoni",
    "Marinelli_1L", "Th-232", "Th232_420-7-17_Маринелли_0cm.xml")

# Полный шаблон цепочки Th-232 в маринелли — единый прогон (150 000
# распадов Th-232 с прохождением всей цепочки через nucleusLimits
# 208..232 81..90). Нормирован на 1 распад РОДИТЕЛЯ ветви со всеми
# ветвлениями внутри — им пользуются оба метода. Индивидуальные iso_
# файлы дают только ДОЛИ нуклидов в отклике (визуально), в амплитуду не
# входят: суммирование iso с учётом branching совпадает с chain_Th232
# по интегралу в 2× (систематика Geant4 с nucleusLimits для одиночных
# нуклидов — регистрация энергии отдачи ядра и вторичных, зависит от
# `applyToAllProcesses`; в chain_Th232 её нет, а физический ответ там).
TEMPLATE_CSV = os.path.join(BUILD, "chain_Th232.csv")

# ─── Ветвь Th-232: восемь звеньев с γ-эмиссией + сама голова ветви ─────────
# ключ · подпись RU · подпись EN · цвет (для stacked area) · комментарий ·
#   BR — ветвление от родителя ветви Th-232 (на 1 распад ветви — сколько
#   распадов ЭТОГО нуклида в равновесии).
# Для Th-232…Bi-212 линия распада без развилок: BR = 1.
# Bi-212 → Tl-208 (α, 35,94 %) или Po-212 (β⁻, 64,06 %). Ниже — только
# Tl-208 распадается γ (Po-212 α, γ нет).
# Числа: ENSDF/DDEP для 212Bi (2020).
BR_TL208 = 0.3594
# Восемь легко различимых цветов; проверено на соседних парах в стеке.
# Основа — палитра с максимальным попарным контрастом (не ColorBrewer Set2
# и не Category10 — там пары дублируются на 8 категориях).
NUCS = [
    ("Th232", "Th-232", "Th-232", "#4a3e6b", 1.0,
     "α, γ практически отсутствуют — тонкая полоска у нуля"),
    ("Ac228", "Ac-228", "Ac-228", "#c4479e", 1.0,
     "β⁻; γ 338, 911, 964, 969 кэВ — доминирует в средней полосе"),
    ("Th228", "Th-228", "Th-228", "#b04e2c", 1.0,
     "α; γ малой интенсивности 84 и 216 кэВ"),
    ("Ra224", "Ra-224", "Ra-224", "#1b8a8f", 1.0,
     "α; γ 241 кэВ"),
    ("Rn220", "Rn-220", "Rn-220", "#7a5cb0", 1.0,
     "α; γ 549 кэВ малой интенсивности"),
    ("Pb212", "Pb-212", "Pb-212", "#e8801c", 1.0,
     "β⁻; γ 238,6 / 300 / 415 кэВ"),
    ("Bi212", "Bi-212", "Bi-212", "#2f7a3c", 1.0,
     "β⁻; γ 727, 785, 1620 кэВ; α-ветка 36 % идёт в Tl-208"),
    ("Tl208", "Tl-208", "Tl-208", "#d2b93a", BR_TL208,
     "β⁻; γ 2614,5 / 583 / 860 / 510 / 277 кэВ — весь жёсткий край; "
     "образуется только в 35,94 % распадов Bi-212 (α-ветка)"),
]

# ─── Паспорт партии 420-7-17 ────────────────────────────────────────────────
PASSPORT_BQ_KG   = 1940.0
PASSPORT_UNC_PCT = 6.0
PASSPORT_MASS_G  = 1600.0
PASSPORT_DATE    = "2007-09-17"
MEAS_DATE        = "2024-10-24"
T12_TH232_YEARS  = 1.405e10
DAYS_PASS_TO_MEAS = 6247

def decay_factor_years(t12, dt_days):
    return math.exp(-math.log(2.0) * (dt_days / 365.25) / t12)

DECAY_FACTOR = decay_factor_years(T12_TH232_YEARS, DAYS_PASS_TO_MEAS)
PASSPORT_A_BQ  = PASSPORT_BQ_KG * (PASSPORT_MASS_G / 1000.0) * DECAY_FACTOR
PASSPORT_dA_BQ = PASSPORT_A_BQ * PASSPORT_UNC_PCT / 100.0

# ─── Метод 2: библиотека γ-линий эмиссии ветви + сумм-пики ─────────────────
# Собрана и перепроверена 08.08.2026 двумя независимыми оценёнными базами
# (ENSDF через IAEA Live Chart API, LNHB/DDEP через PDF-таблицы) —
# расхождение везде <5 %, кроме отмеченных. Значения ниже — ENSDF.
# Порог отбора (директива оператора, R46, 08.08.2026): I_γ ≥ 2 % на распад
# НУКЛИДА (не ветви). Заменил прежний порог 0,5 % — список короче вдвое.
# Th-228 из фиксированной библиотеки ВЫПАЛ ЦЕЛИКОМ: обе его линии (84,373 —
# 1,188 %; 215,985 — 0,2469 %) ниже 2 %. Прежнее ручное исключение из общего
# порога для 215,985 кэВ снято — новый порог применён без исключений.
# Линия Tl-208 1093,900 кэВ (прямая компонента, 0,43 %) тоже ниже порога и
# выпала, но сумм-пик 510,77+583,19→1093,9 (SUM_PEAKS ниже) остаётся: это
# независимая физика (совпадение каскада, а не branching ratio прямого
# перехода) и по интенсивности как раз доминирует над выпавшей прямой линией.
# ИСКЛЮЧЕНЫ отдельно, по данным (нет подтверждения независимым источником,
# порог I≥2 % тут ни при чём): Ac-228 214,850 кэВ (0,76 %, отсутствует в
# LNHB, нет мультипольности в ENSDF) и Ac-228 674,750 кэВ (2,1 % — прошла бы
# порог, но отсутствует в LNHB, погрешность ENSDF ~33 %).
# Полная (нефильтрованная) библиотека — load_full_library() ниже, отдельно,
# без порога вообще: её единственный смысл — сравнение с этой куцей.
# Порог как именованная константа — не только для фильтра ниже, но и чтобы
# число в прозе страницы шло через {{m2_ithresh}}, а не было набрано цифрой
# в разметке (сторож build_page.py такое ловит и останавливает сборку).
I_THRESHOLD_PCT = 2.0
# (E_keV, I_gamma_percent, nuclide_key, note)
GAMMA_LIBRARY = [
    # Ac-228 (BR=1.0), ENSDF cutoff 2012 (K. Abusaleem)
    (129.065, 2.42,  "Ac228", ""),
    (209.253, 3.89,  "Ac228", ""),
    (270.245, 3.46,  "Ac228", ""),
    (328.000, 2.95,  "Ac228", ""),
    (338.320, 11.27, "Ac228", "главная линия Ac-228"),
    (463.004, 4.40,  "Ac228", ""),
    (794.947, 4.25,  "Ac228", ""),
    (911.204, 25.8,  "Ac228", "главная линия Ac-228"),
    (964.766, 4.99,  "Ac228", "+ суммирование 835,710+129,065"),
    (968.971, 15.8,  "Ac228", ""),
    (1588.200, 3.22, "Ac228", ""),
    # Ra-224 (BR=1.0), ENSDF cutoff 2010
    (240.986, 4.10, "Ra224", ""),
    # Pb-212 (BR=1.0), ENSDF cutoff 2020
    (238.632, 43.6,  "Pb212", "главная линия Pb-212"),
    (300.087, 3.301, "Pb212", ""),
    # Bi-212 (BR=1.0), ENSDF cutoff 2007
    (727.330, 6.67, "Bi212", "главная линия Bi-212"),
    # Tl-208 (BR=0,3594), ENSDF cutoff 2007
    (277.371, 6.6,   "Tl208", ""),
    (510.770, 22.6,  "Tl208", "ядерный переход M1+E2, не аннигиляционная линия"),
    (583.187, 85.0,  "Tl208", "главная линия ветви (аналитическая, purity=1,000)"),
    (860.557, 12.5,  "Tl208", "+ суммирование 277,371+583,187"),
    (2614.511, 99.754, "Tl208", "реперная линия ветви"),
]

# Подтверждённые каскады (true coincidence summing) — доказаны через общий
# уровень схемы распада (start/end levels ENSDF) + существование прямой
# crossover-линии, независимо в ENSDF и LNHB. Проверены и ОТКЛОНЕНЫ как
# физически невозможные: Ac-228 338+911, 911+969, 338+964 и Bi-212 727+1620
# — это ветвления ОДНОГО уровня на разные конечные (общий верхний либо
# общий нижний уровень), а не последовательные переходы одного каскада;
# оба кванта пары никогда не испускаются в одном акте распада.
# (E1_keV, E2_keV, nuclide_key, I1_percent, I2_percent, note)
SUM_PEAKS = [
    (510.770, 583.187, "Tl208", 22.6, 85.0,
     "уровень 3197,717 кэВ; crossover 1093,9 подтверждён (0,43-0,44 %)"),
    (277.371, 583.187, "Tl208", 6.6, 85.0,
     "уровень 3197,717 кэВ; суммирование в прямую линию 860,557"),
    (463.004, 911.204, "Ac228", 4.40, 25.8,
     "уровень 968,972 кэВ; crossover 1374,19 подтверждён (~0,014 %)"),
    (794.947, 270.245, "Ac228", 4.25, 3.46,
     "уровень 328,006 кэВ; один из трёх путей на ~1065,18"),
    (726.863, 338.320, "Ac228", 0.62, 11.27,
     "уровень 396,083 кэВ; один из трёх путей на ~1065,18"),
    (153.977, 911.204, "Ac228", 0.722, 25.8,
     "уровень 968,972 кэВ; один из трёх путей на ~1065,19"),
    (835.710, 129.065, "Ac228", 1.61, 2.42,
     "уровень 1022,531 кэВ; суммирование в прямую линию 964,766"),
    # Tl-208 583,187+2614,511=3197,7 — самый значимый канал во всей
    # цепочке (I1·I2 = 8479 %²), но сумма ВНЕ диапазона подгонки
    # E_FIT_HI=2900 кэВ — не включён в активную модель.
]

# ─── Полная библиотека: ВСЕ известные γ-линии ветви, без порога ────────────
# Выгрузка ENSDF через IAEA Live Chart of Nuclides (REST `decay_rads`,
# rad_types=g), дата извлечения в самом файле; провенанс каждой строки —
# оценщик и дата оценки ENSDF — в колонке `source`.
#
# Две особенности выгрузки, обе проверены сверкой с отдельным запросом
# rad_types=x и отражены в файле:
#   * запрос «g» отдаёт вместе с ядерными переходами и характеристический
#     рентген дочернего атома (внутренняя конверсия) — 54 строки из 417;
#     они помечены `line_type=xray` и в γ-библиотеку НЕ идут: рентген в
#     методе 2 учитывается отдельной сущностью, из эмиссионных спектров
#     Geant4, и попал бы в счёт дважды;
#   * часть линий известна по энергии, но интенсивность в ENSDF не
#     оценена (16 строк, Pb-212 и Bi-212) — поле пустое, не ноль. Такие
#     строки пропускаются: подставить им интенсивность нечем.
FULL_LIBRARY_CSV = os.path.join(HERE, "data", "ensdf_th232_chain_lines.csv")


def load_full_library(path=FULL_LIBRARY_CSV, nuc_keys=None):
    """(E_keV, I_percent, nuclide_key, note) по всем линиям с оценённой I."""
    import csv
    if not os.path.isfile(path):
        raise SystemExit("Нет полной библиотеки линий: " + path)
    out, skipped_noI, skipped_xray, skipped_nuc = [], 0, 0, 0
    with io_open(path) as f:
        for row in csv.DictReader(f):
            if row.get("line_type") != "gamma":
                skipped_xray += 1
                continue
            if not (row.get("I_percent") or "").strip():
                skipped_noI += 1
                continue
            key = row["nuclide"]
            if nuc_keys is not None and key not in nuc_keys:
                skipped_nuc += 1
                continue
            unc = (row.get("unc_I_percent") or "").strip()
            note = "ENSDF"
            if unc:
                note += ", σ(I) = " + unc + " %"
            lvl = (row.get("level") or "").strip()
            if lvl:
                note += ", уровень " + lvl + " кэВ"
            out.append((float(row["E_keV"]), float(row["I_percent"]),
                        key, note))
    out.sort(key=lambda r: r[0])
    return out, {"xray": skipped_xray, "no_intensity": skipped_noI,
                 "other_nuclide": skipped_nuc}


def io_open(path):
    return open(path, encoding="utf-8", newline="")


# Диапазон подгонки методов 1 и 2.
E_FIT_LO, E_FIT_HI = 100.0, 2900.0
SYS_FLOOR = 0.03
FINE_STEP_KEV = 0.25


# ─── читатель BecqMoni XML ─────────────────────────────────────────────────

def _text_of(node, tag):
    n = node.find(tag)
    return None if n is None or n.text is None else n.text.strip()


def read_becqmoni_pair(xml_path):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    rd = root.find("./ResultDataList/ResultData")
    out = {}
    for tag, key in [("EnergySpectrum", "sample"),
                     ("BackgroundEnergySpectrum", "background_inline")]:
        node = rd.find(tag)
        if node is None:
            continue
        n_ch = int(_text_of(node, "NumberOfChannels"))
        live = float(_text_of(node, "LiveTime"))
        real = float(_text_of(node, "MeasurementTime"))
        coefs = [float((c.text or "0").strip())
                 for c in node.findall(
                     "./EnergyCalibration/Coefficients/Coefficient")]
        dps = node.findall("./Spectrum/DataPoint")
        counts = np.array([int((dp.text or "0").strip()) for dp in dps],
                          dtype=np.int64)
        if len(counts) != n_ch:
            raise SystemExit("XML: длина Spectrum != NumberOfChannels: %s vs %s"
                             % (len(counts), n_ch))
        ch = np.arange(n_ch, dtype=float)
        e_of_ch = np.polynomial.polynomial.polyval(ch, coefs)
        out[key] = {"counts": counts, "e_of_ch": e_of_ch,
                    "live_s": live, "real_s": real,
                    "start": _text_of(rd, "StartTime"), "coefs": coefs}
    return out


# ─── читатель CSV-выхода Geant4 ────────────────────────────────────────────

def load_hist(path):
    hist = {}
    N = None
    for line in open(path, encoding="utf-8"):
        if line.startswith("#"):
            if "N_primaries" in line:
                N = int(line.split("=")[1])
            continue
        s = line.strip()
        if not s or not s[0].isdigit():
            continue
        e, c = s.split(",")
        hist[float(e)] = int(c)
    return hist, N


def broaden_and_rebin(hist, N_prim, ch_edges, broaden=True):
    """Уширяет депозитный спектр приборным разрешением и интегрирует по
    границам каналов записи. Возвращает счёт в канале на ОДИН распад.

    `broaden=False` — уширение НЕ применяется: депозитные линии остаются
    такими, какими их выдал Geant4 (ширина в один шаг гистограммы), и
    только перекладываются по каналам записи. Это второй режим страницы
    («без учёта аппаратной ПШПВ»): интеграл каждой линии сохраняется, но
    вся площадь сидит в одном-двух каналах вместо реального пика.
    """
    egrid = np.arange(0.0, 3300.0, FINE_STEP_KEV)
    out = np.zeros_like(egrid)
    if not broaden:
        idx = np.clip(np.searchsorted(egrid, list(hist.keys())),
                      0, len(egrid) - 1)
        np.add.at(out, idx, np.asarray(list(hist.values()), dtype=float))
    else:
        for E0, c in hist.items():
            if c == 0 or E0 <= 0:
                continue
            sig = fwhm_kev(max(E0, 8.0)) / 2.3548
            lo = np.searchsorted(egrid, E0 - 5 * sig)
            hi = np.searchsorted(egrid, E0 + 5 * sig)
            if hi <= lo:
                continue
            g = np.exp(-0.5 * ((egrid[lo:hi] - E0) / sig) ** 2)
            s = g.sum()
            if s > 0:
                out[lo:hi] += c * g / s
    cum = np.cumsum(out)
    v = np.interp(ch_edges, egrid, cum, left=cum[0], right=cum[-1])
    return np.diff(v) / max(N_prim, 1)


# ─── площадь пика с полками ─────────────────────────────────────────────────

# Действующий закон ширины. До калибровки (см. fit_fwhm_calibration) —
# одноточечный корневой закон по записи цезия комплекта: ПШПВ(E) =
# FWHM662·√(E/661,657). После калибровки по самому́ измеренному спектру
# Th-232 сюда кладётся степенной закон ПШПВ(E) = k·E^p, и им пользуются
# ВСЕ потребители ширины: уширение шаблонов, окна масок, окна съёма
# площадей и окно эффективности ПП на моно-сетке.
FWHM_LAW = {"kind": "sqrt", "k": FWHM662 / math.sqrt(661.657), "p": 0.5}


def fwhm_kev(E):
    """ПШПВ прибора на энергии E, кэВ — по действующему закону FWHM_LAW."""
    return FWHM_LAW["k"] * float(E) ** FWHM_LAW["p"]


# ─── калибровка по ПШПВ: ширины пиков снимаются с самого спектра ───────────
# Опорные линии — те же шесть, по которым разрешение снималось штатным ПО
# прибора (detector_params.FWHM_MEASURED); здесь они переснимаются НАШИМ
# кодом с этой же записи, и совпадение служит независимой проверкой.
#
# Одиночной гауссианой с линейной подложкой обойтись нельзя, и это не
# теоретическое опасение: при такой попытке линия 300,087 кэВ, сидящая на
# спаде дублета 238,6/241,0, дала ПШПВ 344 кэВ вместо ожидаемых 18 —
# подгонка ушла в подложку; линия 860,557 кэВ, зажатая между 835,7/840,4 и
# 904,2/911,2, дала 69,5 вместо 54. Поэтому окно разбирается ДЕКОНВОЛЮЦИЕЙ:
# сумма гауссиан на позициях ВСЕХ библиотечных линий окна (амплитуды
# свободные, неотрицательные), общая ширина по действующему закону, общий
# сдвиг шкалы — два нелинейных параметра, амплитуды при них линейны.
FWHM_ANCHORS = (238.632, 300.087, 583.187, 727.330, 860.557, 2614.511)


def _multiplet_chi2(theta, x, y, lines, E0, w0, p_law):
    """χ² окна при (σ на E0, сдвиг шкалы); амплитуды — NNLS при них.

    Подложка задана тремя неотрицательными столбцами (константа и наклон
    в обе стороны): у NNLS нет отрицательных коэффициентов, а наклон
    континуума под пиком бывает любого знака.
    """
    sig0, shift = float(theta[0]), float(theta[1])
    if sig0 <= 1e-3 or sig0 > 5.0 * w0 or abs(shift) > w0:
        return 1e30, None
    cols = []
    for Ei in lines:
        si = sig0 * (Ei / E0) ** p_law
        cols.append(np.exp(-0.5 * ((x - (Ei + shift)) / si) ** 2))
    cols.append(np.ones_like(x))
    cols.append((x - E0) / w0)
    cols.append(-(x - E0) / w0)
    A = np.stack(cols, axis=1)
    s = np.sqrt(np.maximum(y, 1.0))
    coef, _ = nnls(A / s[:, None], y / s)
    r = (y - A @ coef) / s
    return float((r * r).sum()), coef


def fit_peak_multiplet(counts, e_of_ch, E0, library_E, half_win_fwhm=2.2):
    """Ширина линии E0 из деконволюции её окна. dict или None."""
    from scipy.optimize import minimize

    w0 = fwhm_kev(E0)
    p_law = FWHM_LAW["p"]
    lo, hi = E0 - half_win_fwhm * w0, E0 + half_win_fwhm * w0
    m = (e_of_ch >= lo) & (e_of_ch <= hi)
    if m.sum() < 9:
        return None
    x = e_of_ch[m]
    y = counts[m].astype(float)
    # Соседи ЗА краем окна тоже заводятся в модель: их хвост входит внутрь
    # и иначе был бы приписан подложке (а через неё — ширине якоря).
    lines = sorted(set(E for E in library_E if lo - w0 <= E <= hi + w0)
                   | {float(E0)})
    args = (x, y, lines, float(E0), w0, p_law)
    best = minimize(lambda t: _multiplet_chi2(t, *args)[0],
                    x0=[w0 / 2.3548, 0.0], method="Nelder-Mead",
                    options={"xatol": 1e-3, "fatol": 1e-3, "maxiter": 800})
    sig0, shift = float(best.x[0]), float(best.x[1])
    chi2, coef = _multiplet_chi2(best.x, *args)
    if coef is None or not np.isfinite(chi2):
        return None
    # Погрешность ширины — по кривизне χ²(σ) при найденном сдвиге:
    # Δσ = √(2/H), H — вторая производная (χ² растёт на 1 при 1σ).
    h = max(0.01 * sig0, 1e-3)
    cp = _multiplet_chi2([sig0 + h, shift], *args)[0]
    cm = _multiplet_chi2([sig0 - h, shift], *args)[0]
    H = (cp - 2.0 * chi2 + cm) / h ** 2
    d_sig = math.sqrt(2.0 / H) if H > 0 and np.isfinite(H) else float("nan")
    amp_anchor = float(coef[lines.index(float(E0))])
    return {"fwhm_keV": 2.3548 * sig0,
            "d_fwhm_keV": 2.3548 * d_sig if np.isfinite(d_sig) else float("nan"),
            "shift_keV": shift, "E_centroid": float(E0) + shift,
            "amp_anchor": amp_anchor, "n_lines_window": len(lines),
            "chi2_ndof": chi2 / max(1, m.sum() - len(lines) - 3)}


def fit_fwhm_calibration(counts, e_of_ch, anchors=FWHM_ANCHORS, n_pass=2):
    """Степенной закон ПШПВ(E) = k·E^p по снятым ширинам опорных линий.

    Форма не произвольна: `sqrt(a + b·E)`, стандартная для полупроводников,
    на этих точках даёт отрицательное `a` и обращается в ноль внутри
    рабочего диапазона (разобрано в detector_params). Степенной закон
    остаётся конечным и монотонным вниз по энергии.

    Два прохода: окна и относительные ширины соседей в деконволюции
    строятся по действующему закону, поэтому после первой калибровки они
    пересчитываются заново уже по ней.

    МНК по логарифмам, веса — обратные квадраты относительной погрешности
    ширины: жёсткая точка снимается втрое точнее мягких, и равные веса
    отдали бы ей подгонку целиком.
    """
    library_E = sorted(set(float(E) for E, _, _, _ in GAMMA_LIBRARY))
    law_saved = dict(FWHM_LAW)
    out = None
    try:
        for _ in range(max(1, n_pass)):
            pts = []
            for E0 in anchors:
                expect = fwhm_kev(E0)
                r = fit_peak_multiplet(counts, e_of_ch, E0, library_E)
                if r is None:
                    pts.append({"E_nominal": float(E0), "used": False,
                                "reject": "окно вырождено, подгонка не строится"})
                    continue
                w, dw = r["fwhm_keV"], r["d_fwhm_keV"]
                q = {"E_nominal": float(E0), "E_centroid": r["E_centroid"],
                     "fwhm_keV": w, "d_fwhm_keV": dw,
                     "res_pct": 100.0 * w / r["E_centroid"],
                     "shift_keV": r["shift_keV"],
                     "n_lines_window": r["n_lines_window"],
                     "chi2_ndof": r["chi2_ndof"], "used": True, "reject": ""}
                # Отбраковка — по признакам, что подогналось не то: якорь
                # обнулён деконволюцией, ширина вне физически возможного
                # коридора, ширина не определена данными.
                if r["amp_anchor"] <= 0:
                    q["used"] = False
                    q["reject"] = "деконволюция обнулила амплитуду линии"
                elif not (0.4 * expect <= w <= 2.5 * expect):
                    q["used"] = False
                    q["reject"] = ("ширина вне коридора 0,4…2,5 от ожидаемой "
                                   "(%.1f кэВ)" % expect)
                elif not np.isfinite(dw) or dw > 0.20 * w:
                    q["used"] = False
                    q["reject"] = "ширина не определена данными (погрешность >20 %)"
                pts.append(q)

            used = [q for q in pts if q["used"]]
            if len(used) < 3:
                raise SystemExit(
                    "калибровка ПШПВ: годных опорных линий меньше трёх (%d из "
                    "%d) — степенной закон по двум точкам не проверяем.\n%s"
                    % (len(used), len(pts),
                       "\n".join("  %.1f: %s" % (q["E_nominal"], q["reject"])
                                 for q in pts if not q["used"])))
            x = np.log(np.array([q["E_centroid"] for q in used]))
            yl = np.log(np.array([q["fwhm_keV"] for q in used]))
            rel = np.array([max(q["d_fwhm_keV"], 1e-6) / q["fwhm_keV"]
                            for q in used])
            wgt = 1.0 / rel ** 2
            Sw = wgt.sum()
            Sx = (wgt * x).sum(); Sy = (wgt * yl).sum()
            Sxx = (wgt * x * x).sum(); Sxy = (wgt * x * yl).sum()
            p = (Sw * Sxy - Sx * Sy) / (Sw * Sxx - Sx * Sx)
            k = math.exp((Sy - p * Sx) / Sw)
            for q in pts:
                if not q["used"]:
                    continue
                q["fwhm_model_keV"] = k * q["E_centroid"] ** p
                q["dev_pct"] = 100.0 * (q["fwhm_model_keV"] / q["fwhm_keV"] - 1.0)
            rms = math.sqrt(sum(q["dev_pct"] ** 2 for q in used) / len(used))
            out = {"k": float(k), "p": float(p), "points": pts,
                   "n_used": len(used), "n_anchors": len(pts),
                   "rms_dev_pct": rms,
                   "fwhm662_law": float(k * 661.657 ** p),
                   "fwhm662_cs": FWHM662,
                   "res662_pct": float(100.0 * k * 661.657 ** p / 661.657)}
            FWHM_LAW.update({"kind": "power", "k": out["k"], "p": out["p"]})
    finally:
        FWHM_LAW.clear()
        FWHM_LAW.update(law_saved)
    return out


def channel_of(E, e_of_ch):
    return float(np.interp(E, e_of_ch, np.arange(len(e_of_ch), dtype=float)))


def peak_area_with_shelf(counts, e_of_ch, E, roi=1.0, shelf=1.0):
    w = fwhm_kev(E)
    lo, hi = E - roi * w, E + roi * w
    lo_s1, lo_s2 = lo - shelf * w, lo
    hi_s1, hi_s2 = hi, hi + shelf * w
    def integ(a, b):
        i0 = max(0, int(round(channel_of(a, e_of_ch))))
        i1 = min(len(counts), int(round(channel_of(b, e_of_ch))) + 1)
        return float(counts[i0:i1].sum()), i0, i1
    gross, i0, i1 = integ(lo, hi)
    n_roi = i1 - i0
    left_sum, _, _  = integ(lo_s1, lo_s2)
    right_sum, _, _ = integ(hi_s1, hi_s2)
    n_left = max(1, int(round(channel_of(lo_s2, e_of_ch)
                              - channel_of(lo_s1, e_of_ch))))
    n_right = max(1, int(round(channel_of(hi_s2, e_of_ch)
                               - channel_of(hi_s1, e_of_ch))))
    background = 0.5 * (left_sum / n_left + right_sum / n_right) * n_roi
    net = gross - background
    dnet = math.sqrt(max(gross, 0)
                     + (n_roi / (2.0 * n_left)) ** 2 * left_sum
                     + (n_roi / (2.0 * n_right)) ** 2 * right_sum)
    return gross, background, net, dnet


# ─── эффективность пика полного поглощения — интерполяция по сетке ────────
# grid/rho1.60_E*.csv — 24 моноэнергетических прогона того же объёмного
# источника в маринелли (ρ=1,6 г/см³ — тот же раствор, что и образец),
# 400 000 историй каждый. Геометрия подтверждена шапкой файлов
# (`mode = vessel`) и разбором режима в main.cc: `vessel` без двоеточия —
# это `VesselGeom::Preset("marinelli")`. Раньше здесь стояла сетка
# `denta1.60` — снятая в ДЕНТЕ 120 мл (`mode = vessel:denta`), то есть в
# чужой геометрии: ε_ПП там выше маринелли на 5 % на 2614,5 кэВ и на 14 %
# на 238,6 кэВ, и расхождение энергозависимо — смещало не только
# амплитуду метода 2, но и наклон его модели по энергии. ε_ПП(E) = счёт в окне ±2,5 ПШПВ вокруг узла,
# делённый на число историй. Между узлами — кусочно-линейная интерполяция
# в координатах (ln E, ln ε): кривая эффективности детектора гладкая и
# без изломов, лог-лог спрямляет степенной участок падения после ~200 кэВ.
# За краями сетки (45,3…3552,5 кэВ) — линейная экстраполяция тем же
# наклоном, что на крайнем интервале (библиотека линий из неё не выходит).

def _grid_main_csvs(grid_dir, pattern):
    """glob(pattern), но БЕЗ файлов-спутников _chan.csv/_emit.csv.

    R45 (08.08.2026): main.cc стал писать rho1.60_E00661.7_chan.csv рядом
    с основным rho1.60_E00661.7.csv — новый файл текстуально подходит под
    старый шаблон "rho1.60_E*.csv" (glob не видит границу токена), и
    load_eps_peak_grid пытался распарсить 12-колоночный файл разложения по
    каналам как 2-колоночный спектр — ValueError на первой же строке.
    Найдено этим же прогоном на первом запуске после правки."""
    return [f for f in glob.glob(os.path.join(grid_dir, pattern))
            if not (f.endswith("_chan.csv") or f.endswith("_emit.csv"))]


def load_eps_peak_grid(grid_dir):
    files = sorted(_grid_main_csvs(grid_dir, "rho1.60_E*.csv"))
    if not files:
        raise SystemExit("Нет моно-сетки эффективности: " + grid_dir)
    Es, Eps = [], []
    for f in files:
        hist = {}
        N = None
        E0 = None
        for line in open(f, encoding="utf-8"):
            if line.startswith("#"):
                if "N_primaries" in line:
                    N = int(line.split("=")[1].strip())
                if "E_prim_keV" in line:
                    E0 = float(line.split("=")[1].strip())
                continue
            s = line.strip()
            if not s or not s[0].isdigit():
                continue
            e_s, c_s = s.split(",")
            hist[float(e_s)] = int(c_s)
        fw = fwhm_kev(E0)
        lo, hi = E0 - 2.5 * fw, E0 + 2.5 * fw
        peak = sum(c for e0, c in hist.items() if lo <= e0 <= hi)
        Es.append(E0)
        Eps.append(peak / N)
    order = np.argsort(Es)
    return np.asarray(Es)[order], np.asarray(Eps)[order]


def make_eps_peak_interp(grid_dir):
    Es, Eps = load_eps_peak_grid(grid_dir)
    logE, logEps = np.log(Es), np.log(Eps)

    def eps_peak(E):
        scalar_in = np.ndim(E) == 0
        Ea = np.atleast_1d(np.asarray(E, dtype=float))
        out = np.exp(np.interp(np.log(Ea), logE, logEps))
        below = Ea < Es[0]
        above = Ea > Es[-1]
        if below.any():
            k = (logEps[1] - logEps[0]) / (logE[1] - logE[0])
            out[below] = np.exp(logEps[0] + k * (np.log(Ea[below]) - logE[0]))
        if above.any():
            k = (logEps[-1] - logEps[-2]) / (logE[-1] - logE[-2])
            out[above] = np.exp(logEps[-1] + k * (np.log(Ea[above]) - logE[-1]))
        return float(out[0]) if scalar_in else out

    return eps_peak


# ─── R45: полная матрица отклика метода 2 ──────────────────────────────────
# Замена узких гауссиан (R19/R33/R37, включая отдельную сущность SECOND
# для вылета аннигиляции) на методику из wt20-methods-compare.md:
#     r_i = Σ_линий y_l · R(E_l -> i)
# R(E->i) — ПОЛНЫЙ депонированный спектр монохроматического кванта энергии E
# (пик + комптоновский континуум + все вторичные процессы), а не только его
# площадь под пиком. Вторичные пики (одиночный/двойной вылет аннигиляции)
# перестают быть отдельной сущностью — они физически часть отклика ЛЮБОЙ
# линии выше 1022 кэВ, канал pair_esc1/pair_esc2 (см. ниже).

# Разложение отклика по каналам взаимодействия — те же 11 каналов, что
# main.cc пишет в *_chan.csv (см. enum Chan там же); порядок и имена
# ОБЯЗАНЫ совпадать построчно, иначе zip() молча перепутает столбцы.
CHAN_NAMES = ["photo", "compt_full", "compt_esc1", "compt_escN", "xray_esc",
              "brems_esc", "pair_full", "pair_esc1", "pair_esc2", "external",
              "other"]
CHAN_LABEL_RU = {
    "photo":      "фотоэффект",
    "compt_full": "комптон, поглощён",
    "compt_esc1": "комптон, квант ушёл (1×)",
    "compt_escN": "комптон, квант ушёл (N×)",
    "xray_esc":   "вылет характеристического рентгена",
    "brems_esc":  "вылет тормозного",
    "pair_full":  "пары, оба кванта поглощены",
    "pair_esc1":  "пары, один квант 511 кэВ ушёл",
    "pair_esc2":  "пары, оба кванта 511 кэВ ушли",
    "external":   "вторичные извне (защита/MgO/корпус/резина)",
    "other":      "остаточный канал (сторож)",
}
# Палитра — визуально различима на соседних парах (тот же приём, что и у
# NUCS выше): проверялась на живом рендере стека, не подбором на глаз.
CHAN_COLOR = {
    "photo":      "#2f6fb0",
    "compt_full": "#5aa8d6",
    "compt_esc1": "#1b8a8f",
    "compt_escN": "#7ec9a8",
    "xray_esc":   "#c4479e",
    "brems_esc":  "#e07bb0",
    "pair_full":  "#b04e2c",
    "pair_esc1":  "#e8801c",
    "pair_esc2":  "#f0b93a",
    "external":   "#6b5f4a",
    "other":      "#9a9284",
}


def load_chan_hist(path):
    """*_chan.csv -> ({имя_канала: {E_keV: count}}, N_primaries)."""
    N = None
    out = {c: {} for c in CHAN_NAMES}
    header = None
    with io_open(path) as f:
        for line in f:
            if line.startswith("#"):
                if "N_primaries" in line:
                    N = int(line.split("=")[1].strip())
                continue
            s = line.strip()
            if not s:
                continue
            if header is None:
                header = s.split(",")
                assert header[1:] == CHAN_NAMES, (
                    "заголовок %s разошёлся с CHAN_NAMES: %s" % (path, header))
                continue
            parts = s.split(",")
            E0 = float(parts[0])
            for name, val in zip(header[1:], parts[1:]):
                c = int(val)
                if c:
                    out[name][E0] = c
    if N is None:
        raise SystemExit("Нет N_primaries в шапке: " + path)
    return out, N


def load_grid_nodes(grid_dir, pattern="rho1.60_E*.csv"):
    """Все узлы моно-сетки по возрастанию энергии:
    [(E0, hist_полный, N, {канал: hist_канала}), ...]."""
    files = sorted(_grid_main_csvs(grid_dir, pattern))
    if not files:
        raise SystemExit("Нет моно-сетки: " + grid_dir)
    nodes = []
    for f in files:
        hist, N = load_hist(f)
        E0 = None
        for line in open(f, encoding="utf-8"):
            if line.startswith("#") and "E_prim_keV" in line:
                E0 = float(line.split("=")[1].strip())
                break
        if E0 is None:
            raise SystemExit("Нет E_prim_keV в шапке: " + f)
        cf = f[:-4] + "_chan.csv"
        if not os.path.isfile(cf):
            raise SystemExit(
                "Нет разложения по каналам: " + cf + " — сетку нужно "
                "пересчитать инструментированным exe (R45)")
        chan, Nc = load_chan_hist(cf)
        if Nc != N:
            raise SystemExit("N_primaries спектра и каналов разошлись: %s "
                              "(%d vs %d)" % (f, N, Nc))
        # Та же проверка баланса, что при записи в main.cc, но здесь ещё
        # раз — ловит и повреждение файла при копировании, не только
        # дефект самого прогона.
        sum_chan = sum(sum(d.values()) for d in chan.values())
        sum_hist = sum(hist.values())
        if sum_chan != sum_hist:
            raise SystemExit("Узел %s: сумма каналов %d != спектру %d"
                              % (f, sum_chan, sum_hist))
        nodes.append((E0, hist, N, chan))
    nodes.sort(key=lambda t: t[0])
    return nodes


def make_full_response(grid_dir, ch_edges, broaden, eps_peak_interp):
    """Полная матрица отклика R(E->i): для произвольной энергии линии E —
    ГОТОВАЯ форма на сетке каналов записи (пик + континуум + все вторичные
    процессы, не окно вокруг пика) и её разложение по 11 каналам
    взаимодействия (см. CHAN_NAMES).

    Форма между 24 узлами сетки НЕ интерполируется по форме (задача,
    которую wt20-methods-compare.md не специфицирует) — берётся форма
    БЛИЖАЙШЕГО по энергии узла, сдвинутая на (E-E_узла) и промасштабированная
    так, чтобы её СОБСТВЕННАЯ площадь пика совпала с ИНТЕРПОЛИРОВАННОЙ
    (лог-лог, eps_peak_interp — тот же интерполятор, что и раньше) на E.
    Точна на самих узлах (сдвиг <1 кэВ для большинства линий отобранной
    библиотеки); худший случай в ней — сдвиг ~176 кэВ (Ac-228 1588,2 кэВ к
    узлу 1764,5). Приближение первого порядка: континуум детектора меняется
    с энергией плавно, без изломов на масштабе одного шага сетки — тот же
    класс допущения, что уже стоит за eps_peak между узлами."""
    nodes = load_grid_nodes(grid_dir)
    node_Es = np.array([n[0] for n in nodes])

    # ОДИН broaden_and_rebin на узел (не на линию библиотеки: 20-347 линий
    # против 24 узлов) — экономия на порядки при той же физике. Сдвиг ГОТОВОЙ
    # формы коммутирует со свёрткой ПШПВ с точностью до медленной
    # зависимости ПШПВ(E) от энергии — второй порядок на шаге сетки.
    shape_node, chan_shape_node, peak_node = [], [], []
    for E0, hist, N, chan in nodes:
        shape_node.append(broaden_and_rebin(hist, N, ch_edges, broaden))
        chan_shape_node.append({
            c: broaden_and_rebin(chan[c], N, ch_edges, broaden)
            for c in CHAN_NAMES})
        fw = fwhm_kev(E0)
        lo, hi = E0 - 2.5 * fw, E0 + 2.5 * fw
        peak_node.append(sum(cnt for e0, cnt in hist.items()
                              if lo <= e0 <= hi) / N)

    e_centers = 0.5 * (ch_edges[:-1] + ch_edges[1:])

    def _shift(arr, shift_kev):
        if shift_kev == 0.0:
            return arr
        return np.interp(e_centers - shift_kev, e_centers, arr,
                          left=0.0, right=0.0)

    def response(E_target):
        """-> (shape[len(e)], {канал: shape[len(e)]}, eps_peak(E_target))."""
        j = int(np.argmin(np.abs(node_Es - E_target)))
        shift = E_target - node_Es[j]
        scale = float(eps_peak_interp(E_target)) / max(peak_node[j], 1e-30)
        shape = _shift(shape_node[j], shift) * scale
        chans = {c: _shift(chan_shape_node[j][c], shift) * scale
                 for c in CHAN_NAMES}
        return shape, chans, float(eps_peak_interp(E_target))

    return response


# ─── общая процедура подгонки амплитуд ─────────────────────────────────────

def fit_amplitudes(y, cols, sys_floor=SYS_FLOOR):
    """Двухпроходный NNLS: сперва пуассоновские веса, затем с систематическим
    полом по модели первого прохода.

    `cols` — столбцы матрицы плана, каждый уже в отсчётах измерения.
    Возвращает (coef, d_coef, chi2, ndof, model).
    """
    A = np.stack(cols, axis=1)
    sig = np.sqrt(np.maximum(y, 1.0))
    coef, _ = nnls(A / sig[:, None], y / sig)
    model = A @ coef
    sig = np.sqrt(np.maximum(y, 1.0) + (sys_floor * model) ** 2)
    coef, _ = nnls(A / sig[:, None], y / sig)
    model = A @ coef
    chi2 = float((((y - model) / sig) ** 2).sum())
    ndof = max(1, int(len(y) - A.shape[1]))
    W = 1.0 / sig ** 2
    try:
        cov = np.linalg.inv((A.T * W) @ A)
        d_coef = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    except np.linalg.LinAlgError:
        d_coef = np.full(A.shape[1], float("nan"))
    return coef, d_coef, chi2, ndof, model


# ─── главная процедура ─────────────────────────────────────────────────────

def main():
    if not os.path.isfile(XML_MEAS):
        raise SystemExit("Нет XML измерения: %s" % XML_MEAS)

    pair = read_becqmoni_pair(XML_MEAS)
    meas = pair["sample"]
    bg   = pair["background_inline"]

    # фон в шкалу измерения по энергии, масштаб по живому времени
    bg_on_meas = np.interp(meas["e_of_ch"], bg["e_of_ch"],
                           bg["counts"].astype(float), left=0.0, right=0.0)
    bg_scale_time = meas["live_s"] / bg["live_s"]
    bg_scaled = bg_on_meas * bg_scale_time

    # границы каналов измерения
    e = meas["e_of_ch"]
    ch_edges = np.concatenate((
        [e[0] - 0.5 * (e[1] - e[0])],
        0.5 * (e[:-1] + e[1:]),
        [e[-1] + 0.5 * (e[-1] - e[-2])],
    ))
    T = meas["live_s"]

    # ── калибровка по ПШПВ ──────────────────────────────────────────────
    # Снимается с САМОГО разбираемого спектра и заменяет одноточечный
    # корневой закон по записи цезия. Это не косметика: закон по цезию
    # (FWHM662 = 49,9 кэВ) и ширины линий тория расходятся на ~10 % —
    # расхождение отмечено в detector_params и до сих пор не разобрано.
    # Уширение шаблонов входит в обе подгонки, поэтому ошибка в законе
    # ширины смещает обе оценки активности одинаково и незаметно.
    fwhm_cal = fit_fwhm_calibration(meas["counts"], e)
    LAW_LINES = {"kind": "power", "k": fwhm_cal["k"], "p": fwhm_cal["p"]}
    # Закон по записи цезия комплекта — тот же одноточечный корень, которым
    # FWHM_LAW инициализирован при определении модуля (строка ~364), снятый
    # заново явным словарём: после калибровки модуль-глобал будет перезаписан
    # степенным законом, и без этой копии исходное значение было бы потеряно.
    LAW_CS = {"kind": "sqrt", "k": FWHM662 / math.sqrt(661.657), "p": 0.5}
    FWHM_LAW.update(LAW_LINES)

    # ── сырые гистограммы Geant4 читаются ОДИН раз ─────────────────────
    if not os.path.isfile(TEMPLATE_CSV):
        raise SystemExit("Нет chain_Th232.csv: " + TEMPLATE_CSV)
    hist_chain, N_chain = load_hist(TEMPLATE_CSV)

    hist_iso = {}
    missing = []
    for key, ru, en, col, br, note in NUCS:
        p = os.path.join(BUILD, "iso_%s.csv" % key)
        if not os.path.isfile(p):
            missing.append(key)
            continue
        hist_iso[key] = load_hist(p)
    if missing:
        raise SystemExit(
            "Нет МК-шаблонов индивидуальных нуклидов: %s\n"
            "Запустить: cd %s && ./g1s.exe decay_th232_isotopes.mac"
            % (missing, BUILD))

    # Доля K-рентгена каждого нуклида-источника в потоке ветви — из
    # ЭМИССИОННЫХ спектров (энергия кванта до детектора). Депозитный
    # спектр этого не хранит: 60–110 кэВ там набирается и рассеянными
    # квантами жёстких линий, отделить рентген постфактум нельзя.
    K_XRAY_LO, K_XRAY_HI = 60.0, 110.0
    XRAY_GRID_CSV = os.path.join(BUILD, "grid", "rho1.60_E00088.0.csv")
    xray_frac_of_branch = {}
    # Здесь же собирается СПЕКТР эмиссии рентгена (не только его интеграл):
    # методу 2 нужны отдельные энергии, чтобы взять ε_ПП на каждой, — K-серия
    # Z = 80…83 разнесена от 72 до 91 кэВ, и эффективность на её краях
    # отличается заметно.
    xray_emit = defaultdict(float)
    for key, ru, en, col, br, note in NUCS:
        pe = os.path.join(BUILD, "iso_%s_emit.csv" % key)
        if not os.path.isfile(pe):
            continue
        hist_e, N_e = load_hist(pe)
        win = 0
        for E0, c in hist_e.items():
            if K_XRAY_LO <= E0 <= K_XRAY_HI:
                win += c
                xray_emit[float(E0)] += (c / N_e) * br
        xray_frac_of_branch[key] = (win / N_e) * br
    XRAY_TOTAL_PER_BRANCH = sum(xray_frac_of_branch.values())
    xray_share = {k: v / XRAY_TOTAL_PER_BRANCH
                  for k, v in xray_frac_of_branch.items()}
    if not os.path.isfile(XRAY_GRID_CSV):
        raise SystemExit("Нет моно-отклика для XRAY: " + XRAY_GRID_CSV)
    hist_grid_xray, N_grid_xray = load_hist(XRAY_GRID_CSV)

    eps_peak = make_eps_peak_interp(os.path.join(BUILD, "grid"))
    BR_of = {k: br for k, _, _, _, br, _ in NUCS}
    keys = [k for k, _, _, _, _, _ in NUCS] + ["XRAY"]

    sel = (e >= E_FIT_LO) & (e <= E_FIT_HI)
    y_sel = meas["counts"][sel].astype(float)
    bgm = bg_scaled[sel]

    # ── шаблоны и разложение при заданном режиме уширения ───────────────

    def build_templates(broaden):
        """templ_total и by_nuc (включая XRAY) для режима уширения."""
        templ_total = broaden_and_rebin(hist_chain, N_chain, ch_edges, broaden)
        by_nuc_raw = {}
        for key, ru, en, col, br, note in NUCS:
            hist, N = hist_iso[key]
            by_nuc_raw[key] = broaden_and_rebin(hist, N, ch_edges, broaden) * br
        iso_sum = sum(by_nuc_raw.values())
        # Амплитуду даёт chain_Th232 (полный физический транспорт), iso —
        # только НОРМИРОВАННУЮ долю нуклида в канале: суммирование iso с
        # branching расходится с chain по интегралу вдвое (систематика
        # nucleusLimits для одиночного нуклида — регистрируется энергия
        # отдачи ядра и вторичных, которых в цепочечном прогоне нет).
        by_nuc = {}
        for k in by_nuc_raw:
            with np.errstate(divide="ignore", invalid="ignore"):
                share = np.where(iso_sum > 0, by_nuc_raw[k] / iso_sum, 0.0)
            by_nuc[k] = templ_total * share

        # Псевдо-нуклид XRAY: форма — моно-отклик ~88 кэВ (K-серия Z=80..83
        # лежит в узкой полосе 72–91 кэВ), масштаб — выход рентгена на распад
        # ветви. Вычитается из нуклидов-источников пропорционально их доле в
        # эмиссии, иначе рентген был бы учтён дважды.
        xray_shape = broaden_and_rebin(hist_grid_xray, N_grid_xray,
                                       ch_edges, broaden)
        xray_template = xray_shape * XRAY_TOTAL_PER_BRANCH
        for k, sh in xray_share.items():
            if sh > 0 and k in by_nuc:
                by_nuc[k] = by_nuc[k] - xray_template * sh
        by_nuc["XRAY"] = xray_template
        # Приближение одной моноэнергией не обязано точно повторять форму
        # рентгеновской доли каждого источника: на краях полосы вычет уходит
        # в минус. Недостачу возвращаем в XRAY, откуда она пришла, — тогда
        # сумма по компонентам остаётся равна templ_total ТОЧНО в каждом
        # канале, независимо от точности приближения формы.
        for k, sh in xray_share.items():
            if sh <= 0 or k not in by_nuc:
                continue
            deficit = np.clip(-by_nuc[k], 0.0, None)
            if deficit.any():
                by_nuc[k] = by_nuc[k] + deficit
                by_nuc["XRAY"] = by_nuc["XRAY"] - deficit
        neg = {k: float(v.min()) for k, v in by_nuc.items()
               if float(v.min()) < -1e-12}
        if neg:
            raise SystemExit(
                "XRAY: после компенсации всё ещё отрицательно (ошибка в "
                "схеме перераспределения): %s" % neg)
        resid = float(np.max(np.abs(sum(by_nuc.values()) - templ_total)))
        if resid > 1e-9 * float(np.max(templ_total)):
            raise SystemExit(
                "XRAY: баланс Σ by_nuc == templ_total нарушен, невязка %.3e"
                % resid)
        return templ_total, by_nuc

    # ── метод 1: МК-шаблоны по нуклидам ────────────────────────────────
    # ОДНА амплитуда на всю ветвь: в вековом равновесии активности звеньев
    # тождественно равны, а разложение с восемью независимыми амплитудами
    # вырождено на коллинеарных шаблонах соседних звеньев.

    def run_method1(templ_total, by_nuc):
        coef, dcoef, chi2, ndof, _ = fit_amplitudes(
            y_sel, [templ_total[sel] * T, bgm])
        A_branch, bg_amp = float(coef[0]), float(coef[1])
        per_nuc = {}
        for k in keys:
            share = float(by_nuc[k][sel].sum()
                          / max(templ_total[sel].sum(), 1e-30))
            per_nuc[k] = {"A_Bq": A_branch, "dA_Bq": float(dcoef[0]),
                          "share": share}
        return {
            "A_Bq": A_branch, "dA_Bq": float(dcoef[0]),
            "per_nuclide": per_nuc,
            "bg_amplitude": bg_amp, "d_bg_amplitude": float(dcoef[1]),
            "chi2": chi2, "ndof": ndof, "chi2_ndof": chi2 / ndof,
            "E_fit_lo": E_FIT_LO, "E_fit_hi": E_FIT_HI,
            "sys_floor": SYS_FLOOR,
            "ratio_to_passport": A_branch / PASSPORT_A_BQ,
            "d_ratio": float(dcoef[0]) / PASSPORT_A_BQ,
        }, (templ_total * A_branch * T + bg_scaled * bg_amp)

    # ── метод 2 (R45): библиотека линий свёрнута с ПОЛНОЙ матрицей ──────
    # отклика R(E->i) (make_full_response выше) вместо узких гауссиан на
    # позициях линий (R19-R37). Континуум и все вторичные процессы (в т.ч.
    # вылет аннигиляции — ранее отдельная сущность SECOND) входят в отклик
    # КАЖДОЙ линии естественно, отдельно моделировать их больше не нужно.
    # Директива оператора 08.08.2026: "Метод 2 это и есть функция отклика
    # по каналам взаимодействия" — по образцу AtomSpectra Nano 16 PRO.

    def run_method2(library, sums, resp, with_diag=True):
        shape_total = np.zeros_like(e)
        by_nuc_w = {k: np.zeros_like(e) for k in keys}
        by_chan_w = {c: np.zeros_like(e) for c in CHAN_NAMES}
        photon_lines = []

        def add(nuc_key, weight, shp, chans):
            shape_total[:] += weight * shp
            by_nuc_w[nuc_key] += weight * shp
            for c in CHAN_NAMES:
                by_chan_w[c] += weight * chans[c]

        for E, I_pct, nuc_key, note in library:
            shp, chans, eps = resp(E)
            w = BR_of[nuc_key] * (I_pct / 100.0)
            add(nuc_key, w, shp, chans)
            photon_lines.append({
                "E_keV": E, "nuclide": nuc_key, "I_gamma_pct": I_pct,
                "note": note, "eps_peak": eps, "weight_per_branch": w * eps,
                "kind": "line"})

        # ── K-рентген как отдельная сущность метода 2 ──────────────────
        # Выход берётся из ЭМИССИОННЫХ спектров Geant4 (энергия кванта до
        # детектора) в окне 60-110 кэВ — там же, где он выделяется для
        # метода 1. Отличие в том, что метод 2 применяет к каждой энергии
        # K-серии свою полную форму отклика, а не общий моно-отклик 88 кэВ.
        xray_w_total = 0.0
        for E0, yld in xray_emit.items():
            shp, chans, eps = resp(E0)
            add("XRAY", yld, shp, chans)
            xray_w_total += yld * eps
        if xray_w_total > 0:
            photon_lines.append({
                "E_keV": sum(E0 * y for E0, y in xray_emit.items())
                         / max(sum(xray_emit.values()), 1e-30),
                "nuclide": "XRAY", "I_gamma_pct": 100.0 * XRAY_TOTAL_PER_BRANCH,
                "note": "K-серия дочерних атомов (Bi, Pb, Tl, Po, Ra, Th), "
                        "%d энергий в окне %.0f-%.0f кэВ; выход на распад "
                        "ветви — из эмиссионных спектров Geant4"
                        % (len(xray_emit), K_XRAY_LO, K_XRAY_HI),
                "eps_peak": xray_w_total / max(sum(xray_emit.values()), 1e-30),
                "weight_per_branch": xray_w_total, "kind": "xray"})

        n_sum_used = 0
        for E1, E2, nuc_key, I1_pct, I2_pct, note in sums:
            Esum = E1 + E2
            if Esum > E_FIT_HI:
                continue
            _, _, eps1 = resp(E1)
            _, _, eps2 = resp(E2)
            # Форма — с узла, ближайшего к СУММАРНОЙ энергии Esum: континуум
            # суммарного пика физически размазан похоже на континуум
            # одиночного кванта той же полной энергии (то же приближение
            # первого порядка, что и у обычных линий). Абсолютная величина —
            # честное произведение eps1*eps2 (эффективность КАЖДОГО из двух
            # квантов на своей энергии), НЕ eps_peak(Esum).
            shp, chans, eps_sum_node = resp(Esum)
            w = (BR_of[nuc_key] * (I1_pct / 100.0) * (I2_pct / 100.0)
                 * eps1 * eps2 / max(eps_sum_node, 1e-30))
            add(nuc_key, w, shp, chans)
            photon_lines.append({
                "E_keV": Esum, "nuclide": nuc_key, "I_gamma_pct": None,
                "note": note, "eps_peak": eps1 * eps2,
                "weight_per_branch": BR_of[nuc_key] * (I1_pct / 100.0)
                                      * (I2_pct / 100.0) * eps1 * eps2,
                "kind": "sum", "E1_keV": E1, "E2_keV": E2,
                "I1_pct": I1_pct, "I2_pct": I2_pct})
            n_sum_used += 1

        # Подгонка по ВСЕМУ диапазону, как метод 1: модель теперь физически
        # полная (континуум + все вторичные процессы), окно ±2 ПШПВ вокруг
        # пиков было нужно только peak-only модели (R19-R37) и с R45 снято.
        coef, dcoef, chi2, ndof, _ = fit_amplitudes(
            y_sel, [shape_total[sel] * T, bgm])
        A_ph, bg_amp = float(coef[0]), float(coef[1])

        for ph in photon_lines:
            ph["predicted_net"] = float(ph["weight_per_branch"] * A_ph * T)
        if with_diag:
            # Справочно, ПОСЛЕ подгонки и вне её: полковая оценка net — та
            # самая, что на мультиплете неверна; имя поля это и говорит.
            for ph in photon_lines:
                gross, bgroi, net, dnet = peak_area_with_shelf(
                    meas["counts"], meas["e_of_ch"], ph["E_keV"])
                ph["gross"] = int(gross)
                ph["shelf_bg"] = float(bgroi)
                ph["net_shelf_unreliable"] = float(net)
                ph["dnet_shelf_unreliable"] = float(dnet)

        stack2 = {k: (by_nuc_w[k] * A_ph * T).tolist() for k in keys}
        stack2_chan = {c: (by_chan_w[c] * A_ph * T).tolist()
                       for c in CHAN_NAMES}

        return {
            "A_Bq": A_ph, "dA_Bq": float(dcoef[0]),
            "lines": photon_lines,
            "n_lines": len(library),
            "n_sum_peaks": n_sum_used, "n_sum_peaks_total": len(sums),
            "n_channels_fit": int(sel.sum()),
            "xray_weight_per_branch": xray_w_total,
            "n_xray_energies": len(xray_emit),
            "bg_amplitude": bg_amp,
            "chi2": chi2, "ndof": ndof, "chi2_ndof": chi2 / ndof,
            "E_fit_lo": E_FIT_LO, "E_fit_hi": E_FIT_HI,
            "sys_floor": SYS_FLOOR,
            "ratio_to_passport": A_ph / PASSPORT_A_BQ,
            "d_ratio": float(dcoef[0]) / PASSPORT_A_BQ,
            "chan_names": CHAN_NAMES,
        }, (shape_total * A_ph * T + bg_scaled * bg_amp), stack2, stack2_chan, by_nuc_w

    # ── четыре расчёта: {закон ширины: по линиям спектра, по цезию} ×
    #    {отобранная библиотека, полная} ─────────────────────────────────
    # Уширение включено ВСЕГДА (реальная гауссиана приборного разрешения);
    # варианты отличаются ЗАКОНОМ ширины, не наличием свёртки. «По линиям» —
    # степенной закон, снятый деконволюцией этого же спектра (действующий,
    # используется по умолчанию everywhere). «По цезию» — одноточечный
    # корневой закон FWHM662·√(E/661,657) по записи цезия комплекта; две
    # оценки расходятся на 662 кэВ ~10 % (detector_params.py), сравнение
    # показывает, сколько в ответе держится на выборе закона ширины.
    # Метод 1 от библиотеки не зависит — считается по одному разу на режим.
    lib_full, lib_skip = load_full_library(nuc_keys=set(BR_of))
    variants = {}
    for tag, law in (("lines", LAW_LINES), ("cs", LAW_CS)):
        FWHM_LAW.clear()
        FWHM_LAW.update(law)
        templ_total, by_nuc = build_templates(True)
        m1, model1 = run_method1(templ_total, by_nuc)
        resp = make_full_response(os.path.join(BUILD, "grid"), ch_edges,
                                   True, eps_peak)
        m2, model2, stack2, stack2_chan, by_nuc_w2 = run_method2(
            GAMMA_LIBRARY, SUM_PEAKS, resp, with_diag=True)
        m2f, model2f, stack2f, stack2f_chan, by_nuc_w2f = run_method2(
            lib_full, SUM_PEAKS, resp, with_diag=False)
        stack = {k: (by_nuc[k] * m1["A_Bq"] * T).tolist() for k in keys}
        variants[tag] = {
            "method1": m1, "method2": m2, "method2_full": m2f,
            "stack": stack, "stack2": stack2, "stack2_full": stack2f,
            "stack2_chan": stack2_chan, "stack2_chan_full": stack2f_chan,
            "model_counts": model1.tolist(),
            "model2_counts": model2.tolist(),
            "model2_full_counts": model2f.tolist(),
            "templ_total": templ_total, "by_nuc": by_nuc,
            "by_nuc_w2": by_nuc_w2, "by_nuc_w2f": by_nuc_w2f,
        }
    # Восстановить действующий закон на «по линиям» — вариант по умолчанию,
    # на нём же и дальнейший код модуля (ничего ниже fwhm_kev не вызывает,
    # но пусть глобал не остаётся в последнем состоянии цикла — ловушка).
    FWHM_LAW.clear()
    FWHM_LAW.update(LAW_LINES)

    V = variants["lines"]
    templ_total, by_nuc = V["templ_total"], V["by_nuc"]
    method1_mc, method2_photon = V["method1"], V["method2"]
    A_branch = method1_mc["A_Bq"]

    # ── подгонка с гипотезой утечки торона (Rn-220 — газ, T½ = 55,6 с) ──
    # Ветвь разделена по узлу Rn-220: утечка отрубает всё, что образуется
    # ИЗ ушедшего газа, — группу «after» пропорционально её доле. Две
    # амплитуды NNLS дают A_before и A_after независимо; η = 1 − A_after/A_before.
    # Применяется к ОБОИМ методам (уточнение оператора, R56) — метод 1
    # (МК-шаблоны по нуклидам) и метод 2 (библиотека линий) дважды,
    # отобранная и полная библиотека: физика деления по группам одна и та
    # же, входные шаблоны — разные.
    KEYS_BEFORE = ["Th232", "Ac228", "Th228", "Ra224"]
    KEYS_AFTER  = ["Rn220", "Pb212", "Bi212", "Tl208"]
    # XRAY распределён между группами по факту происхождения, иначе
    # templ_before + templ_after перестаёт замыкаться на templ_total.
    xray_share_before = sum(xray_share.get(k, 0.0) for k in KEYS_BEFORE)
    xray_share_after  = sum(xray_share.get(k, 0.0) for k in KEYS_AFTER)

    def _leak_refit(by_nuc_src, m_single):
        """Подгонка с гипотезой утечки на ЛЮБОМ источнике по-нуклидных
        форм (by_nuc метода 1 или by_nuc_w2/by_nuc_w2f метода 2) — та же
        схема групп KEYS_BEFORE/KEYS_AFTER, тот же взвешенный XRAY.
        m_single — однoамплитудный результат ЭТОГО ЖЕ метода (источник
        χ² для сравнения). Возвращает (fit-словарь, stack по нуклидам)."""
        tb = (sum(by_nuc_src[k] for k in KEYS_BEFORE)
              + xray_share_before * by_nuc_src["XRAY"])
        ta = (sum(by_nuc_src[k] for k in KEYS_AFTER)
              + xray_share_after * by_nuc_src["XRAY"])
        c, d, chi, ndof_l, _ = fit_amplitudes(
            y_sel, [tb[sel] * T, ta[sel] * T, bgm])
        Ab, Aa = float(c[0]), float(c[1])
        dAb, dAa = float(d[0]), float(d[1])
        if Ab > 0:
            eta = 1.0 - Aa / Ab
            rel = math.sqrt((dAa / max(Aa, 1e-9)) ** 2
                            + (dAb / max(Ab, 1e-9)) ** 2)
            d_eta = abs(Aa / Ab) * rel
        else:
            eta = float("nan"); d_eta = float("nan")
        xray_amp = xray_share_before * Ab + xray_share_after * Aa
        xray_damp = xray_share_before * dAb + xray_share_after * dAa
        fit = {
            "A_before_Bq": Ab, "dA_before_Bq": dAb,
            "A_after_Bq":  Aa, "dA_after_Bq":  dAa,
            "eta_leak": eta, "d_eta": d_eta,
            "bg_amplitude": float(c[2]),
            "chi2": chi, "chi2_ndof": chi / ndof_l, "ndof": ndof_l,
            # Δχ² однoамплитудной модели ЭТОГО ЖЕ метода минус эта подгонка
            # (одним параметром больше) — знак: положительное = утечка лучше.
            "delta_chi2_vs_single": m_single["chi2"] - chi,
            "keys_before": KEYS_BEFORE, "keys_after": KEYS_AFTER,
            # XRAY не принадлежит целиком ни одной группе (шаблон — уже
            # сумма вкладов всех нуклидов-источников) — взвешенная по факту
            # происхождения амплитуда; погрешность — линейная комбинация той
            # же весовой смесью (консервативно, без учёта корреляции
            # A_before/A_after — вклад XRAY в спектр ~0,1 %, точнее не нужно).
            "xray_amp_leak": xray_amp, "xray_damp_leak": xray_damp,
        }
        stack_leak = {}
        for k in keys:
            if k == "XRAY":
                amp = xray_amp
            elif k in KEYS_BEFORE:
                amp = Ab
            else:
                amp = Aa
            stack_leak[k] = (by_nuc_src[k] * amp * T).tolist()
        return fit, stack_leak

    # Посчитано ОДИН раз, на законе ширины «по линиям» — при переключении
    # на закон цезия страница показывает этот же результат (упрощение,
    # объявлено явно, тот же класс, что у остальных полей V=variants["lines"]).
    leak_fit, stack_leak = _leak_refit(by_nuc, method1_mc)
    leak_fit2, stack2_leak = _leak_refit(V["by_nuc_w2"], method2_photon)
    leak_fit2_full, stack2_leak_full = _leak_refit(
        V["by_nuc_w2f"], V["method2_full"])

    # ── упаковка JSON ─────────────────────────────────────────────────────
    CS = variants["cs"]
    data = {
        "meta": {
            "detector": "Гамма-1С (УДС-ГЦ-63х63)",
            "vessel": "Маринелли 1 л, ОИСН-16 ρ=1,6 г/см³",
            "source_batch": "420-7-17",
            "start_time": meas["start"],
            "live_s": meas["live_s"], "real_s": meas["real_s"],
            "bg_source": "BackgroundEnergySpectrum того же XML",
            "bg_live_s": bg["live_s"], "bg_real_s": bg["real_s"],
            "bg_scale_time": bg_scale_time,
            "fwhm662_keV": FWHM662, "escape_keV": ESCAPE_KEV,
            "k_xray_lo_keV": K_XRAY_LO, "k_xray_hi_keV": K_XRAY_HI,
            "template_source": "iso_*.csv (200 000 распадов на нуклид)",
            "sys_floor_pct": SYS_FLOOR * 100.0,
            # Калибровка энергии из XML: полином E(канал) = Σ c_i · канал^i.
            # Порядок и коэффициенты свои у образца и фона (разные записи).
            "cal_sample": {
                "coefs": meas["coefs"], "order": len(meas["coefs"]) - 1,
                "n_channels": int(len(meas["counts"])),
            },
            "cal_bg": {
                "coefs": bg["coefs"], "order": len(bg["coefs"]) - 1,
                "n_channels": int(len(bg["counts"])),
            },
        },
        "fwhm_cal": fwhm_cal,
        "passport": {
            "batch": "420-7-17", "Bq_per_kg": PASSPORT_BQ_KG,
            "unc_pct": PASSPORT_UNC_PCT, "mass_g": PASSPORT_MASS_G,
            "date_certified": PASSPORT_DATE, "date_measured": MEAS_DATE,
            "decay_factor": DECAY_FACTOR,
            "A_Bq": PASSPORT_A_BQ, "dA_Bq": PASSPORT_dA_BQ,
        },
        "nuclides": [
            {"key": k, "label_ru": ru, "label_en": en, "color": col,
             "note": note, "branching": br}
            for k, ru, en, col, br, note in NUCS
        ] + [
            {"key": "XRAY", "label_ru": "K-рентген", "label_en": "K X-ray",
             "color": "#6b5f4a", "branching": XRAY_TOTAL_PER_BRANCH,
             "note": "характеристический K-рентген дочерних атомов после "
                     "распада (Bi, Th, Ra, Pb, Po); выделен из шаблонов-"
                     "источников свёрткой их доли эмиссии с моно-откликом "
                     "~88 кэВ — без этого рентген учитывался бы дважды"},
        ],
        # Разложение отклика метода 2 по каналам взаимодействия (R45) — та
        # же классификация, что в main.cc (enum Chan): фотоэффект, комптон
        # (с вылетом / без), рождение пар (с вылетом одного или двух
        # аннигиляционных квантов), вылет рентгена/тормозного, вторичные
        # извне. Заменяет прежнюю отдельную сущность SECOND (R33/R37) —
        # вылет аннигиляции теперь виден напрямую как pair_esc1/pair_esc2.
        "channels": [
            {"key": c, "label_ru": CHAN_LABEL_RU[c], "color": CHAN_COLOR[c]}
            for c in CHAN_NAMES
        ],
        "spectrum": {
            "e_of_ch": e.tolist(),
            "counts": meas["counts"].tolist(),
            "bg_counts": bg_scaled.tolist(),
            "model_counts": V["model_counts"],
            # Модель метода 2 (R45) — полная матрица отклика: пик каждой
            # библиотечной линии со своим комптоновским континуумом и
            # всеми вторичными процессами, не узкая гауссиана.
            "model2_counts": V["model2_counts"],
            "model2_full_counts": V["model2_full_counts"],
            "stack": V["stack"],
            "stack_leak": stack_leak,
            "stack2": V["stack2"],
            "stack2_full": V["stack2_full"],
            "stack2_leak": stack2_leak,
            "stack2_leak_full": stack2_leak_full,
            "stack2_chan": V["stack2_chan"],
            "stack2_chan_full": V["stack2_chan_full"],
        },
        # Тот же расчёт с ДРУГИМ законом ширины: одноточечный корень по
        # записи цезия комплекта (FWHM662·√(E/661,657)) вместо степенного
        # закона, снятого деконволюцией самого́ разбираемого спектра тория.
        # Уширение применяется в обоих случаях — отличие только в законе.
        "cs": {
            "method1": CS["method1"], "method2": CS["method2"],
            "method2_full": CS["method2_full"],
            "spectrum": {"model_counts": CS["model_counts"],
                         "model2_counts": CS["model2_counts"],
                         "model2_full_counts": CS["model2_full_counts"],
                         "stack": CS["stack"],
                         "stack2": CS["stack2"],
                         "stack2_full": CS["stack2_full"],
                         "stack2_chan": CS["stack2_chan"],
                         "stack2_chan_full": CS["stack2_chan_full"]},
        },
        "method1": method1_mc,
        "method2": method2_photon,
        "method2_full": V["method2_full"],
        "library": {
            "fixed_n": len(GAMMA_LIBRARY), "full_n": len(lib_full),
            "i_threshold_pct": I_THRESHOLD_PCT,
            "skipped_xray": lib_skip["xray"],
            "skipped_no_intensity": lib_skip["no_intensity"],
            "source": "ENSDF через IAEA Live Chart of Nuclides "
                      "(decay_rads, rad_types=g)",
            "csv": os.path.basename(FULL_LIBRARY_CSV),
        },
        "leak_fit": leak_fit,
        "leak_fit2": leak_fit2,
        "leak_fit2_full": leak_fit2_full,
        "reference_lines": [
            # (E_keV, нуклид, короткая метка) — для реперов в калибровке
            (238.632, "Pb-212", "Pb-212 238"),
            (241.997, "Ra-224", "Ra-224 242"),
            (277.351, "Tl-208", "Tl-208 277"),
            (300.087, "Pb-212", "Pb-212 300"),
            (338.320, "Ac-228", "Ac-228 338"),
            (463.004, "Ac-228", "Ac-228 463"),
            (510.770, "Tl-208", "Tl-208 511"),
            (583.187, "Tl-208", "Tl-208 583"),
            (727.330, "Bi-212", "Bi-212 727"),
            (785.370, "Bi-212", "Bi-212 785"),
            (860.564, "Tl-208", "Tl-208 861"),
            (911.204, "Ac-228", "Ac-228 911"),
            (964.766, "Ac-228", "Ac-228 965"),
            (968.971, "Ac-228", "Ac-228 969"),
            (1620.500, "Bi-212", "Bi-212 1621"),
            (2614.511, "Tl-208", "Tl-208 2615"),
        ],
    }

    # Спектральных массивов теперь три десятка (два режима ПШПВ × два
    # варианта библиотеки × разложения по нуклидам). Двойная точность в
    # JSON раздувает выгрузку втрое, не добавляя ни одной значащей цифры:
    # это отсчёты в канале, у них своя пуассоновская погрешность.
    def _round_arrays(node, nd=2):
        if isinstance(node, list):
            if node and all(isinstance(v, float) for v in node):
                return [round(v, nd) for v in node]
            return [_round_arrays(v, nd) for v in node]
        if isinstance(node, dict):
            return {k: _round_arrays(v, nd) for k, v in node.items()}
        return node

    for sec in (data["spectrum"], data["cs"]["spectrum"]):
        for k, v in list(sec.items()):
            sec[k] = _round_arrays(v, 4 if k == "e_of_ch" else 2)

    out = os.path.join(HERE, "g1s_th232_data.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, separators=(",", ":"))
    print("записано:", out, "(%.0f КБ)" % (os.path.getsize(out) / 1024))
    print("калибровка ПШПВ:  %.3f·E^%.4f  (СКО %.1f %%, на 662 кэВ %.1f кэВ "
          "против %.1f по цезию)"
          % (fwhm_cal["k"], fwhm_cal["p"], fwhm_cal["rms_dev_pct"],
             fwhm_cal["fwhm662_law"], fwhm_cal["fwhm662_cs"]))
    for q in fwhm_cal["points"]:
        print("   %8.1f кэВ  ПШПВ %6.2f ± %.2f кэВ  (%.2f %%)  модель %6.2f  "
              "откл. %+5.1f %%"
              % (q["E_nominal"], q["fwhm_keV"], q["d_fwhm_keV"],
                 q["res_pct"], q["fwhm_model_keV"], q["dev_pct"]))
    print("паспорт            A = %.0f ± %.0f Бк"
          % (PASSPORT_A_BQ, PASSPORT_dA_BQ))
    print("библиотека: отобранная %d линий, полная %d (пропущено: рентген %d, "
          "без оценённой I %d)"
          % (len(GAMMA_LIBRARY), len(lib_full), lib_skip["xray"],
             lib_skip["no_intensity"]))
    for tag, title in (("lines", "ПШПВ линий"), ("cs", "ПШПВ цезия")):
        m1 = variants[tag]["method1"]
        print("%s метод 1          A = %6.0f ± %3.0f Бк, ratio %.3f, "
              "χ²/ν = %8.2f"
              % (title, m1["A_Bq"], m1["dA_Bq"],
                 m1["ratio_to_passport"], m1["chi2_ndof"]))
        for lk, lt in (("method2", "метод 2 (отобр.) "),
                       ("method2_full", "метод 2 (полная) ")):
            m2 = variants[tag][lk]
            print("%s %sA = %6.0f ± %3.0f Бк, ratio %.3f, χ²/ν = %8.2f, "
                  "линий %d, каналов %d"
                  % (title, lt, m2["A_Bq"], m2["dA_Bq"],
                     m2["ratio_to_passport"], m2["chi2_ndof"],
                     m2["n_lines"], m2["n_channels_fit"]))
    print("по нуклидам (метод 1, с ПШПВ):")
    for k in keys:
        v = method1_mc["per_nuclide"][k]
        print("  %-6s  доля %5.1f %%   A = %8.1f ± %6.1f  Бк"
              % (k, 100 * v["share"], v["A_Bq"], v["dA_Bq"]))


if __name__ == "__main__":
    main()
