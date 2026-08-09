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
import yaml
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

# ─── Конфиг источника (этап 2 обобщения конвейера, задача #175/#177) ──
# Всё, что специфично ИМЕННО этому источнику/цепочке (нуклиды, паспорт,
# библиотека линий, пары суммирования, диапазон подгонки), лежит в
# configs/<id>.yaml, а не в теле скрипта. Ниже — только раскрытие в те же
# имена переменных, что использует остальной (общий) движок файла, чтобы
# нижняя часть модуля не менялась вовсе. Какой конфиг использовать —
# переменная G4MODELS_SOURCE_CONFIG, по умолчанию th232.yaml (этот же
# источник, что и раньше — рефакторинг не должен ничего сдвинуть).
CONFIG_PATH = os.environ.get(
    "G4MODELS_SOURCE_CONFIG", os.path.join(HERE, "configs", "th232.yaml"))
with open(CONFIG_PATH, encoding="utf-8") as _f:
    _CFG = yaml.safe_load(_f)

XML_MEAS = os.path.join(KIT, *_CFG["source"]["measured_xml_rel"].split("/"))

# Полный шаблон цепочки Th-232 в маринелли — единый прогон (150 000
# распадов Th-232 с прохождением всей цепочки через nucleusLimits
# 208..232 81..90). Нормирован на 1 распад РОДИТЕЛЯ ветви со всеми
# ветвлениями внутри — им пользуются оба метода. Индивидуальные iso_
# файлы дают только ДОЛИ нуклидов в отклике, в амплитуду не входят.
#
# Прежде здесь стояло объяснение расхождения iso и chain систематикой
# Geant4 с nucleusLimits — регистрацией энергии отдачи ядра и вторичных.
# Объяснение было НЕВЕРНЫМ и держалось потому, что звучало правдоподобно
# и никто не сверил его с числом. Настоящая причина (R61, 08.08.2026):
# восемь iso_-прогонов шли в режиме геометрии без сосуда, где тома
# Sample не существует, поэтому /gps/pos/confine Sample был снят самим
# Geant4 с предупреждением в консоли, и 12,5 % распадов рождались внутри
# кристалла NaI. Доля объёма розыгрыша, попадающая в кристалл, считается
# из чертежа: 12,517 %; наблюдалось 12,494…12,562 %.
# Файлы перегнаны в режиме vessel; в шапке каждого прогона теперь есть
# поле src_in_crystal, и оно обязано быть нулём (проверяется ниже).
TEMPLATE_CSV = os.path.join(BUILD, *_CFG["source"]["chain_template_csv"].split("/"))

NUCS = [(n["key"], n["label_ru"], n["label_en"], n["color"], n["br"], n["note_ru"])
        for n in _CFG["nuclides"]]

# ─── Паспорт источника ──────────────────────────────────────────────
# Номер партии — НЕ публикуется (директива оператора 08.08.2026): в
# отличие от каталожного номера набора ОИСН, партия конкретного
# аттестационного источника этой же схемой лишь условно защищена от
# идентификации владельца комплекта поверки, риск невыгодно мал против
# цены ошибки. Сам номер остаётся только в приватном имени входного XML
# (XML_MEAS выше) — он не публикуется, туда путь не уходит.
_P = _CFG["passport"]
PASSPORT_BQ_KG   = _P["bq_per_kg"]
PASSPORT_UNC_PCT = _P["unc_pct"]
PASSPORT_MASS_G  = _P["mass_g"]
PASSPORT_DATE    = _P["passport_date"]
MEAS_DATE        = _P["measured_date"]
T12_TH232_YEARS  = _P["half_life_years"]
DAYS_PASS_TO_MEAS = _P["days_pass_to_meas"]

def decay_factor_years(t12, dt_days):
    return math.exp(-math.log(2.0) * (dt_days / 365.25) / t12)

DECAY_FACTOR = decay_factor_years(T12_TH232_YEARS, DAYS_PASS_TO_MEAS)
PASSPORT_A_BQ  = PASSPORT_BQ_KG * (PASSPORT_MASS_G / 1000.0) * DECAY_FACTOR
PASSPORT_dA_BQ = PASSPORT_A_BQ * PASSPORT_UNC_PCT / 100.0

# ─── Метод 2: библиотека γ-линий эмиссии ветви + сумм-пики ───────────
# Собрана и перепроверена 08.08.2026 двумя независимыми оценёнными
# базами (ENSDF через IAEA Live Chart API, LNHB/DDEP через PDF-таблицы)
# — расхождение везде <5 %, кроме отмеченных. Значения ниже — ENSDF.
# Порог отбора (директива оператора, R46, 08.08.2026): I_γ ≥ 2 % на
# распад НУКЛИДА (не ветви). Полная (нефильтрованная) библиотека —
# load_full_library() ниже, отдельно, без порога вообще: её единственный
# смысл — сравнение с этой куцей. Какие нуклиды/линии вошли или были
# исключены вручную (например Th-228 целиком из-за порога, Ac-228
# 214,850/674,750 без подтверждения LNHB) — см. configs/th232.yaml и
# историю коммитов этого файла до этапа 2.
# Порог как именованная константа — не только для фильтра ниже, но и
# чтобы число в прозе страницы шло через {{m2_ithresh}}, а не было
# набрано цифрой в разметке (сторож build_page.py такое ловит).
I_THRESHOLD_PCT = _CFG["library"]["intensity_threshold_pct"]
# (E_keV, I_gamma_percent, nuclide_key, note)
GAMMA_LIBRARY = [(l["e_kev"], l["i_pct"], l["nuclide"], l["note"])
                 for l in _CFG["library"]["lines"]]

# Подтверждённые каскады (true coincidence summing) — доказаны через
# общий уровень схемы распада (start/end levels ENSDF) + существование
# прямой crossover-линии, независимо в ENSDF и LNHB. Проверены и
# ОТКЛОНЕНЫ как физически невозможные: Ac-228 338+911, 911+969, 338+964
# и Bi-212 727+1620 — это ветвления ОДНОГО уровня на разные конечные
# (общий верхний либо общий нижний уровень), а не последовательные
# переходы одного каскада; оба кванта пары никогда не испускаются в
# одном акте распада. Сами пары — configs/th232.yaml (ручной отбор,
# категория C инвентаризации, задача #176); F_B считает движок ниже.
# (E1_keV, E2_keV, nuclide_key, I1_percent, I2_percent, note)
SUM_PEAKS = [(s["e1_kev"], s["e2_kev"], s["nuclide"], s["i1_pct"], s["i2_pct"], s["note"])
             for s in _CFG["sum_peaks"]]

def _sum_peaks_with_fb(pairs, path=None):
    """Добавляет седьмое поле fb_pct — F_B, суммарную ДЕПОПУЛЯЦИЮ (=ЗАСЕЛЕНИЕ,
    по закону сохранения) общего уровня каскада, гамма-канал + внутренняя
    конверсия.

    Зачем. P(оба кванта пары полностью поглощены за один распад) —
    НЕ произведение маргинальных выходов I1·I2 (это двойной счёт
    заселения уровня), а I1·I2/F_B, где F_B — доля распадов, в которых
    уровень вообще заселяется (см. #166/R78, находка внутреннего аудита
    09.08.2026, P1).

    ИСПРАВЛЕНО 09.08.2026 (внешний аудит, находка Б1, коммит df5d178):
    F_B считался как сумма I_percent ТОЛЬКО гамма-квантов, исходящих из
    уровня — это депопуляция ОДНИМ каналом (радиационным), а не полная.
    Для низкоэнергетичных переходов в тяжёлых ядрах внутренняя конверсия
    может на порядок превышать гамма-канал: уровень 186,827 кэВ Ac-228
    депопулируется ЕДИНСТВЕННЫМ гамма-переходом 129,065 кэВ (I=2,42 %),
    но входящих в тот же уровень (каскадные гамма сверху + прямое
    бета-питание) — 11,01-11,47 % → старая формула занижала F_B в
    4,5-4,7 раза для этого уровня (проверено ДВУМЯ независимыми путями:
    population = Σ прямого бета-питания уровня, IAEA Live Chart
    decay_rads rad_types=bm, + Σ каскадных гамма сверху; и depopulation =
    Σ I_gamma·(1+CC) по коэффициентам полной внутренней конверсии, тот же
    decay_rads rad_types=g, поле conversion_coeff — сошлись в пределах
    4 %). Для уровней с преимущественно высокоэнергетичными переходами
    (609,318 кэВ Bi-214, CC≈0,02; большинство уровней Ac-228 выше 300 кэВ)
    поправка на порядки меньше — единицы процентов (до ~4,3 % на уровне
    396,083 кэВ Ac-228), в пределах обычного разброса самосогласованности
    ENSDF (баланс население/деполяция не обязан сходиться день в день,
    известное свойство оценённых данных).

    Верная формула: F_B = Σ I_gamma·(1+CC) по ВСЕМ гамма-переходам,
    исходящим из уровня (CC — полный коэффициент внутренней конверсии
    перехода). Источник CC — data/conversion_coeff_sum_peak_levels.csv
    (tools/fetch_conversion_coeff.py, тот же decay_rads запрос, что и
    основная библиотека, просто с сохранённой колонкой conversion_coeff,
    которую прежний экспорт не забирал). Список нуклидов там намеренно
    узкий — только те, что фигурируют в SUM_PEAKS (Ac228, Tl208, Bi214);
    для перехода без найденного CC подставляется 0 (чистый радиационный
    канал, поведение как в старом коде) — молча, но со счётчиком в конце
    функции, чтобы пропуск не остался незамеченным.

    Уровень для каждой пары находится не по тексту `note` (там уже был
    минимум один опечатанный номер — 1022,531 вместо 186,827 для пары
    835,710+129,065, найдено и исправлено при том же аудите), а
    программно: конец перехода E1 обязан совпасть с началом перехода
    E2 (или наоборот) в колонке `level` полной выгрузки.
    """
    import csv
    # FULL_LIBRARY_CSV определяется НИЖЕ по файлу (после этого блока) —
    # здесь нельзя ссылаться на неё как на значение по умолчанию
    # (NameError при загрузке модуля), поэтому путь собирается заново
    # из HERE и source.id конфига (не хардкод "th232" — найдено при
    # Пилоте 1, задача #182: до этой правки F_B всегда считался по
    # библиотеке Th-232 независимо от загруженного источника).
    csv_path = path or os.path.join(
        HERE, "data", "ensdf_%s_chain_lines.csv" % _CFG["source"]["id"])
    cc_path = os.path.join(HERE, "data", "conversion_coeff_sum_peak_levels.csv")
    # (nuc, round(E_keV,3)) -> CC полного коэффициента внутренней конверсии;
    # см. tools/fetch_conversion_coeff.py и докстринг выше (аудит Б1).
    cc_of = {}
    if os.path.isfile(cc_path):
        with open(cc_path, encoding="utf-8", newline="") as f:
            for row in csv.DictReader(
                    ln for ln in f if not ln.startswith("#")):
                cc_raw = (row.get("conversion_coeff") or "").strip()
                if not cc_raw:
                    continue
                try:
                    cc_of[(row["nuclide"], round(float(row["E_keV"]), 3))] = \
                        float(cc_raw)
                except ValueError:
                    continue

    # Нуклиды, реально фигурирующие в SUM_PEAKS этого конфига — диагностика
    # "нет CC" имеет смысл только для них (только их уровни идут в F_B);
    # для остальных нуклидов ветви (Bi212, Pb212, ... в th232) отсутствие
    # CC ожидаемо (tools/fetch_conversion_coeff.py их не выгружал) и не
    # влияет на F_B, шуметь диагностикой о них незачем.
    pairs_nuclides = set(p[2] for p in pairs)

    by_energy = {}   # (nuc, E_keV округлённая) -> (start, end)
    by_start = []     # (nuc, start, I_percent*(1+CC)) — для депопуляции уровня (F_B)
    no_cc_rows = []   # (nuc, E_keV, I_percent) — только для pairs_nuclides
    with open(csv_path, encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            if row.get("line_type") != "gamma":
                continue
            lvl = (row.get("level") or "").strip()
            if not lvl or "->" not in lvl:
                continue
            try:
                start_s, end_s = lvl.split("->")
                start, end = float(start_s), float(end_s)
                E_row = float(row["E_keV"])
                ip = float(row["I_percent"])
            except ValueError:
                continue
            nuc = row["nuclide"]
            by_energy[(nuc, round(E_row, 3))] = (start, end)
            cc = cc_of.get((nuc, round(E_row, 3)))
            if cc is None:
                if nuc in pairs_nuclides:
                    no_cc_rows.append((nuc, E_row, ip))
                cc = 0.0  # нет данных о конверсии — как в старом коде (только γ)
            by_start.append((nuc, start, ip * (1.0 + cc)))
    if no_cc_rows:
        sys.stderr.write(
            "_sum_peaks_with_fb: %d гамма-строк нуклидов %s без CC в "
            "data/conversion_coeff_sum_peak_levels.csv — учтён только "
            "радиационный канал для НИХ, F_B уровней с их участием может "
            "быть занижен (см. докстринг функции): %s\n"
            % (len(no_cc_rows), sorted(pairs_nuclides),
               ["%s %.3f (I=%.3g%%)" % r for r in no_cc_rows]))

    def transition(nuc, E, tol=0.05):
        # Матч по КОЛОНКЕ E_keV самой строки (то же число, что в SUM_PEAKS),
        # а не пересчёт из start-end: разница уровней и заявленная E_keV
        # округлены в ENSDF независимо и не обязаны совпасть день в день
        # (пример: 510,77 кэВ vs 3708,41-3197,717=510,693 — 0,08 кэВ мимо).
        best = None
        for (n, e_row), se in by_energy.items():
            if n != nuc:
                continue
            d = abs(e_row - E)
            if d < tol and (best is None or d < best[0]):
                best = (d, se)
        return best[1] if best else None

    def depopulation(nuc, level, tol=0.02):
        return sum(ip for n, s, ip in by_start
                   if n == nuc and abs(s - level) < tol)

    out = []
    for E1, E2, nuc, I1, I2, note in pairs:
        t1, t2 = transition(nuc, E1), transition(nuc, E2)
        if t1 is None or t2 is None:
            raise SystemExit(
                "SUM_PEAKS: переход %.3f или %.3f нуклида %s не найден в %s"
                % (E1, E2, nuc, csv_path))
        s1, e1 = t1
        s2, e2 = t2
        if abs(e1 - s2) < 0.05:
            level = s2
        elif abs(e2 - s1) < 0.05:
            level = s1
        else:
            raise SystemExit(
                "SUM_PEAKS: %.3f (%.3f→%.3f) и %.3f (%.3f→%.3f) "
                "нуклида %s не стыкуются общим уровнем"
                % (E1, s1, e1, E2, s2, e2, nuc))
        fb = depopulation(nuc, level)
        if fb <= 0 or fb < max(I1, I2) - 0.5:
            raise SystemExit(
                "SUM_PEAKS: депопуляция уровня %.3f кэВ (%s) = %.3f %% — "
                "меньше входящей в неё линии, не может быть верно"
                % (level, nuc, fb))
        out.append((E1, E2, nuc, I1, I2, note, fb))
    return out


SUM_PEAKS = _sum_peaks_with_fb(SUM_PEAKS)

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
FULL_LIBRARY_CSV = os.path.join(
    HERE, "data", "ensdf_%s_chain_lines.csv" % _CFG["source"]["id"])


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


# Диапазон подгонки методов 1 и 2 — configs/<id>.yaml, поле fit.
E_FIT_LO, E_FIT_HI = _CFG["fit"]["e_lo_kev"], _CFG["fit"]["e_hi_kev"]
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

def load_hist(path, require_vessel=True):
    """Гистограмма энерговыделения из выхода Geant4, с проверкой шапки.

    `require_vessel` — прогон обязан быть сделан с построенным сосудом и
    без единой первичной вершины внутри кристалла. Проверка не косметика:
    ровно этих двух полей не хватало, чтобы поймать R61 — восемь шаблонов
    ветви были посчитаны в режиме без сосуда, /gps/pos/confine молча снят
    самим Geant4, и 12,5 % распадов разыграно внутри NaI. Дефект прожил
    недели и всплыл случайно, по «лишнему» горбу на графике страницы.
    Файл без поля `src_in_crystal` считается прогнанным до появления
    сторожа в main.cc и потому непроверенным — тоже отказ.
    """
    hist = {}
    N = None
    mode = None
    src_in_crystal = None
    for line in open(path, encoding="utf-8"):
        if line.startswith("#"):
            if "N_primaries" in line:
                N = int(line.split("=")[1])
            elif line.startswith("# mode ="):
                mode = line.split("=", 1)[1].strip()
            elif line.startswith("# src_in_crystal ="):
                src_in_crystal = int(line.split("=", 1)[1])
            continue
        s = line.strip()
        if not s or not s[0].isdigit():
            continue
        e, c = s.split(",")
        hist[float(e)] = int(c)
    if require_vessel:
        name = os.path.basename(path)
        if mode is None or not mode.startswith("vessel"):
            raise SystemExit(
                "%s: прогон в режиме '%s', а нужен vessel. В режиме без\n"
                "сосуда тома 'Sample' не существует, /gps/pos/confine Sample\n"
                "снимается самим Geant4 с предупреждением, и источник\n"
                "разыгрывается по всему цилиндру, включая кристалл (R61)."
                % (name, mode))
        if src_in_crystal is None:
            raise SystemExit(
                "%s: в шапке нет поля src_in_crystal — файл посчитан сборкой\n"
                "без сторожа 'источник в детекторе'. Перегнать текущим exe."
                % name)
        if src_in_crystal != 0:
            raise SystemExit(
                "%s: src_in_crystal = %d, первичные вершины попали внутрь\n"
                "кристалла. Шаблон негоден." % (name, src_in_crystal))
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


def tcs_report(templ_mc, shape_lib, e_of_ch, ch_edges):
    """Измерить, чего методу 2 недостаёт по каскадному совпадению (R78).

    Метод 1 строит шаблон полным МК-транспортом ЦЕПОЧКИ: два кванта одного
    распада попадают в кристалл одним событием, и суммирование совпадений
    учтено само собой — и уход из одиночных пиков, и приход в сумм-пик.

    Метод 2 складывает отклики отдельных линий с их библиотечными выходами.
    Приход учтён (список SUM_PEAKS), УХОД — нет: каждая линия входит с полным
    выходом I, как если бы её квант всегда регистрировался в одиночку. Модель
    поэтому завышает площади одиночных пиков, а подгонка одной амплитудой ко
    всему спектру занижает активность.

    Обе кривые здесь — НА РАСПАД РОДИТЕЛЯ, без подгонки, поэтому отношение
    сравнимо напрямую. Приводится по полосам энергии: где отношение ниже
    единицы, метод 2 «нарисовал» больше отсчётов, чем даёт полный транспорт.
    """
    e = np.asarray(e_of_ch, dtype=float)
    mc = np.asarray(templ_mc, dtype=float)
    lb = np.asarray(shape_lib, dtype=float)
    bands = [(50, 150), (150, 350), (350, 600), (600, 1000),
             (1000, 1800), (1800, 2900)]
    print("   каскадное суммирование (R78): МК-цепочка / свёртка библиотеки, "
          "на распад родителя", flush=True)
    for lo, hi in bands:
        m = (e >= lo) & (e < hi)
        a, b = float(mc[m].sum()), float(lb[m].sum())
        if b <= 0:
            continue
        print("      %4d-%4d кэВ: %6.3f" % (lo, hi, a / b), flush=True)
    m = (e >= 50) & (e < 2900)
    a, b = float(mc[m].sum()), float(lb[m].sum())
    if b > 0:
        print("      весь диапазон: %6.3f  (ниже 1 — метод 2 завышает "
              "модель и потому занижает активность)" % (a / b), flush=True)


def resolve_mask(mask, e_of_ch):
    """Огрубить булеву маску до разрешения прибора (сторож R66).

    Поканальное сравнение n_eff с порогом даёт рваную маску: n_eff — сумма
    гауссиан от отдельных МК-отсчётов, и между ними она проседает. На графике
    это выглядело как дыры, пробитые сквозь плотный слой (замечание
    оператора 09.08.2026 по Tl-208 около 2500 кэВ) — при том что дыра уже
    ПШПВ, то есть уже того, что прибор в принципе способен разрешить.

    Поэтому провал или островок короче ПШПВ на этой энергии не считается
    признаком: сначала закрываются короткие провалы, затем убираются короткие
    островки. Оба порога — одна и та же ширина, так что операция не сдвигает
    границы протяжённых участков, а только гасит колебание около порога.
    """
    m = np.asarray(mask, dtype=bool).copy()
    n = m.size
    e = np.asarray(e_of_ch, dtype=float)
    # ширина ПШПВ в каналах на каждом канале
    dE = np.gradient(e)
    wch = np.maximum(1.0, np.array([fwhm_kev(max(x, 1.0)) for x in e]) / dE)

    def runs_of(arr):
        out, i = [], 0
        while i < n:
            j = i
            while j + 1 < n and arr[j + 1] == arr[i]:
                j += 1
            out.append((i, j, bool(arr[i])))
            i = j + 1
        return out

    for want in (False, True):          # сперва закрыть провалы, потом островки
        for a, b, val in runs_of(m):
            if val is not want:
                continue
            if a == 0 or b == n - 1:    # краевой участок не трогаем: он не
                continue                # окружён противоположным значением
            if (b - a + 1) < wch[(a + b) // 2]:
                m[a:b + 1] = not want
    return m.tolist()


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
    """glob(pattern), но БЕЗ файлов-спутников _chan/_emit/_emitx/_shield.csv.

    R45 (08.08.2026): main.cc стал писать rho1.60_E00661.7_chan.csv рядом
    с основным rho1.60_E00661.7.csv — новый файл текстуально подходит под
    старый шаблон "rho1.60_E*.csv" (glob не видит границу токена), и
    load_eps_peak_grid пытался распарсить 12-колоночный файл разложения по
    каналам как 2-колоночный спектр — ValueError на первой же строке.
    Найдено этим же прогоном на первом запуске после правки.

    R69 (09.08.2026), тот же класс: main.cc стал писать rho1.60_E00088.0_
    shield.csv (4-колоночный, src_xray) — glob снова его подобрал, снова
    ValueError, снова на первом прогоне после правки."""
    SUFFIXES = ("_chan.csv", "_emit.csv", "_emitx.csv", "_shield.csv")
    return [f for f in glob.glob(os.path.join(grid_dir, pattern))
            if not f.endswith(SUFFIXES)]


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

    # ── K-рентген: отбор ПО ПРОЦЕССУ РОЖДЕНИЯ, не по энергетическому окну ──
    #
    # Было (до 09.08.2026): доля рентгена считалась как интеграл спектра
    # эмиссии в полосе 60–110 кэВ, а форма приближалась одним моно-откликом
    # 88 кэВ. Окно не различает, ЧЕМ рождён квант, и забирало ядерные линии
    # Th-228 84,373 и Th-232 63,8 — то есть лишало Th-228 единственной
    # наблюдаемой линии, на которой держится проверка возраста ряда (R69).
    # Приближение формы одной энергией вдобавок перебирало в максимуме и
    # недобирало по краям, отчего вычет уходил в минус, обрезался по нулю, и
    # слои нуклидов ложились в РОВНЫЙ ноль на 73–105 кэВ.
    #
    # Стало: разделение делает сам Geant4 по модели, породившей трек
    # (main.cc, Tracking: model_RDM_AtomicRelaxation против model_RDM_IT).
    # Модель пишет два файла на каждый нуклид:
    #   iso_<X>_emitx.csv  — испущенное, колонки x_atomic и g_nuclear;
    #   iso_<X>_shield.csv — ОТКЛИК, колонка src_xray: срабатывания, энергию
    #                        в которые принёс рентген атомной релаксации.
    # Второй — точное ПОДМНОЖЕСТВО шаблона нуклида, поэтому вычитание не
    # может дать отрицательного и обрезки по нулю не требует вовсе.
    def load_col(path, name):
        """{E_keV: counts} по ИМЕНОВАННОЙ колонке и число распадов прогона."""
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
                            "в %s нет колонки %s — файл посчитан сборкой до "
                            "разделения по происхождению кванта (R69)"
                            % (os.path.basename(path), name))
                    continue
                hist[float(p[0])] = float(p[cols.index(name)])
        if not N:
            raise SystemExit("в %s нет N_primaries" % path)
        return hist, N

    xray_frac_of_branch = {}
    xray_dep = {}          # депозитный спектр рентгена нуклида, на распад
    # Спектр ЭМИССИИ рентгена (не только интеграл): методу 2 нужны отдельные
    # энергии, чтобы взять ε_ПП на каждой, — K-серия Z = 80…83 разнесена от
    # 72 до 91 кэВ, и эффективность на её краях отличается заметно.
    xray_emit = defaultdict(float)
    for key, ru, en, col, br, note in NUCS:
        # N_primaries того же прогона — из уже проверенного load_hist(iso_X.csv)
        # (mode=vessel*, src_in_crystal=0, см. hist_iso чуть выше). emitx.csv и
        # shield.csv — файлы ТОГО ЖЕ вызова EndOfRunAction, отдельно их не
        # штампуют: доверие к главному файлу распространяется на спутников.
        _, N_iso = hist_iso[key]

        pe = os.path.join(BUILD, "iso_%s_emitx.csv" % key)
        pemit = os.path.join(BUILD, "iso_%s_emit.csv" % key)
        if os.path.isfile(pe):
            hist_x, N_x = load_col(pe, "x_atomic")
        elif os.path.isfile(pemit):
            # _emit.csv есть, _emitx.csv нет: main.cc пишет их одним блоком
            # `if (emitted > 0)`, разойтись при текущем exe они не могут —
            # значит спектр посчитан сборкой ДО разделения по происхождению.
            raise SystemExit(
                "%s: есть iso_%s_emit.csv, но нет iso_%s_emitx.csv — "
                "посчитано сборкой до разделения эмиссии по происхождению "
                "(R69). Перепрогнать rerun_th232_R61.mac текущим exe."
                % (key, key, key))
        else:
            hist_x, N_x = {}, N_iso     # нуклид не эмитирует вовсе — легитимный ноль
        tot = 0.0
        for E0, c in hist_x.items():
            if c <= 0:
                continue
            tot += c
            xray_emit[float(E0)] += (c / N_x) * br
        xray_frac_of_branch[key] = (tot / N_x) * br

        ps = os.path.join(BUILD, "iso_%s_shield.csv" % key)
        if os.path.isfile(ps):
            hist_d, N_d = load_col(ps, "src_xray")
        else:
            # Файл пишется, только если хоть одно срабатывание несёт признак
            # свинца/защиты/рентгена (main.cc: `if (fPbXHits||fShXHits||
            # fSrcXHits)`). Его отсутствие для нуклида с редкой эмиссией —
            # легитимный ноль, а не признак старого прогона: пары в
            # load_col уже отказывают на файле без колонки src_xray, и
            # отдельно на «есть emit, нет emitx» выше.
            hist_d, N_d = {}, N_iso
        xray_dep[key] = ({E0: c for E0, c in hist_d.items() if c > 0}, N_d)
    XRAY_TOTAL_PER_BRANCH = sum(xray_frac_of_branch.values())
    # Диапазон энергий, реально классифицированных как рентген атомной
    # релаксации, — для отчёта на странице. Не окно отбора (его больше нет,
    # разбор R69): это то, что фактически вернул признак Geant4.
    XRAY_SPAN_LO = min(xray_emit) if xray_emit else 0.0
    XRAY_SPAN_HI = max(xray_emit) if xray_emit else 0.0

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
        # Рентген атомной релаксации выделяется ДО нормировки, вычитанием
        # подмножества: в шаблоне нуклида он уже есть, и, если его не вынуть,
        # сущность XRAY учла бы его вторично.
        xray_raw = {}
        for key, ru, en, col, br, note in NUCS:
            hist_d, N_d = xray_dep[key]
            xray_raw[key] = broaden_and_rebin(hist_d, N_d, ch_edges,
                                              broaden) * br
            by_nuc_raw[key] = by_nuc_raw[key] - xray_raw[key]
            # Подмножество не может превысить целое; отрицательное здесь
            # означало бы рассогласование файлов, а не неточность формы.
            bad = float(by_nuc_raw[key].min())
            if bad < -1e-9 * float(np.max(np.abs(by_nuc_raw[key])) + 1e-30):
                raise SystemExit(
                    "%s: рентген больше самого шаблона (%.3e) — iso_%s.csv и "
                    "iso_%s_shield.csv из разных прогонов" % (key, bad, key, key))
            by_nuc_raw[key] = np.maximum(by_nuc_raw[key], 0.0)
        by_nuc_raw["XRAY"] = sum(xray_raw.values())
        iso_sum = sum(by_nuc_raw.values())

        # ── сторож статистики на долю нуклида (R66/R76) ──────────────────
        # Доля нуклида в канале считается как by_nuc_raw[k]/iso_sum и потом
        # умножается на templ_total. Там, где шаблон нуклида набран единицами
        # отсчётов МК, эта доля — не физика, а пуассоновский шум; умноженная
        # на большой полный отклик, она даёт слой, ПОВТОРЯЮЩИЙ ФОРМУ отклика
        # и читающийся на лог-шкале как настоящие пики. Так вылезли «фантомы»
        # слабых нуклидов выше 2000 кэВ и структура сущности рентгена на
        # 1600-2600 кэВ, где рентгена быть не может (разбор R76).
        #
        # Мера доверия — ожидаемое число МК-отсчётов, попавших в канал: тот же
        # broaden_and_rebin, но БЕЗ деления на число розыгрышей. Уширение
        # размазывает отсчёт по ~ПШПВ (на 2614 кэВ это ~38 каналов), поэтому
        # одиночный отсчёт даёт n_eff ~ 0,03, а населённая линия — сотни.
        # Относительная пуассоновская погрешность доли ~ 1/sqrt(n_eff).
        #
        # ОГРУБЛЕНИЕ, объявленное явно: отсчёт, размазанный по нескольким
        # каналам, вносит в них КОРРЕЛИРОВАННЫЕ доли, и точная дисперсия
        # ниже пуассоновской по n_eff. Для порога «шум или не шум» этого
        # достаточно; полноценная ковариация шаблонов здесь не считается.
        n_eff = {}
        for key, ru, en, col, br, note in NUCS:
            hist, _ = hist_iso[key]
            n_eff[key] = broaden_and_rebin(hist, 1.0, ch_edges, broaden)
        xr_hist = {}
        for key, ru, en, col, br, note in NUCS:
            for E0, c in xray_dep[key][0].items():
                xr_hist[E0] = xr_hist.get(E0, 0.0) + c
        n_eff["XRAY"] = broaden_and_rebin(xr_hist, 1.0, ch_edges, broaden)
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

        # XRAY вошёл в by_nuc_raw наравне с нуклидами и потому уже
        # отнормирован вместе с ними: обрезки по нулю и возврата недостачи
        # больше нет — они и создавали ровные нули в слоях (R69).
        neg = {k: float(v.min()) for k, v in by_nuc.items()
               if float(v.min()) < -1e-12}
        if neg:
            raise SystemExit("отрицательные значения в разложении: %s" % neg)
        # Сторож R69 п.4: вычет рентгена не должен обнулять слой там, где до
        # вычета вклад был. Проверяется именно ЭТО, а не «нули вообще»: ноль
        # там, где у нуклида нет линий (Th-232 выше 100 кэВ), законен и
        # физичен, а ноль на месте бывшего вклада — след обрезки.
        for k, v in by_nuc.items():
            if k == "XRAY" or k not in xray_raw:
                continue
            full = by_nuc_raw[k] + xray_raw[k]
            eaten = (full > 0) & (v <= 0)
            if int(eaten.sum()) > 3:
                lo = float(np.asarray(e)[eaten][0])
                raise SystemExit(
                    "%s: вычет рентгена обнулил слой в %d каналах, где вклад "
                    "был (около %.1f кэВ) — R69" % (k, int(eaten.sum()), lo))
        resid = float(np.max(np.abs(sum(by_nuc.values()) - templ_total)))
        if resid > 1e-9 * float(np.max(templ_total)):
            raise SystemExit(
                "XRAY: баланс Σ by_nuc == templ_total нарушен, невязка %.3e"
                % resid)
        return templ_total, by_nuc, n_eff

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

        # ── истощение одиночных линий совпадением (R78/R85, #166) ───────
        # Каждая подтверждённая пара (сумм-пик) уносит часть вероятности
        # ИЗ ОБОИХ одиночных пиков E1 и E2: в том распаде, где оба кванта
        # каскада захвачены одновременно, событие идёт в сумм-пик и не
        # может ещё раз войти в одиночный пик ни одной из двух линий.
        # Раньше цикл ниже брал I_pct библиотеки НЕИЗМЕНЁННЫМ, то есть
        # эту потерю не учитывал — метод 2 предсказывал БОЛЬШЕ отсчётов
        # на активность, чем есть на самом деле, и подгонка той же
        # измеренной площади занижала активность (это и есть суть R78).
        # Стандартная формула поправки (Debertin/Helmer, TCS correction),
        # ИСПРАВЛЕНО 09.08.2026 (внутренний аудит перед внешним, находка P1,
        # #174): наивное произведение маргинальных выходов I1·I2 дважды
        # учитывает заселение общего уровня, когда тот заселяется НЕ в
        # 100% распадов (F_B<100%). Верно — I1·I2/F_B, где F_B — суммарная
        # депопуляция уровня (см. SUM_PEAKS, седьмое поле fb_pct,
        # _sum_peaks_with_fb выше). При F_B≈100% (583+2614 Tl-208)
        # поправка на F_B ничтожна; при F_B=42% (пары 911,204 Ac-228) —
        # растёт разложение в 2,4 раза; проверено также обратным счётом
        # по независимому прогону внутреннего аудита — совпало день в день.
        #
        # ИСПРАВЛЕНО 09.08.2026, вечер (внешний аудит, находка Б2, коммит
        # df5d178): эффективность ПАРТНЁРА каскада в депопуляции должна
        # быть ПОЛНОЙ (εT — вероятность зарегистрировать ХОТЬ ЧТО-ТО от
        # этого кванта где угодно в спектре), а не пиковой. Из совпадения
        # выпадает событие уже при ЛЮБОМ энерговыделении партнёра в то же
        # окно разрешения — не только при его полном поглощении. Источник:
        # Chehade (2007, IUP Bremen, магистерская диссертация, ур. 2.2-2.3,
        # метод Debertin & Schötzig 1990) — ✅ прочитано по содержанию,
        # дословно: «the probability of γ2 being detected and appearing
        # anywhere in the spectrum (εT2). Therefore all the coincidences
        # whether giving rise to a sum peak count or not, should be taken
        # into consideration and hence the final term uses εT2, the total
        # efficiency for the detection of γ2» (стр. 31, ур. 2.2: n'1 =
        # Af1ε1 − Af1ε1εT2). До правки код использовал eps_peak партнёра
        # (третий элемент resp()) — недоучёт истощения в 1,1-7,5 раза в
        # зависимости от энергии (проверено на сетке: E=129 кэВ ratio=1,06,
        # E=2614,511 ratio=7,53, E=583,187 ratio=3,01, E=835,71 ratio=3,63,
        # E=609,321 ratio=3,11 — отношение εT/ε_peak для этой геометрии,
        # независимо перепроверено адверсариальным аудитом субагента).
        # ВАЖНО: у сумм-пика (Esum-вес ниже) партнёр эффективности остаётся
        # ПИКОВОЙ у ОБОИХ квантов — сумм-пик требует ПОЛНОГО поглощения
        # обоих, это другая величина и другая физика, Б2 её не касается.
        # eps_total(E) = shape.sum() — интеграл ГОТОВОЙ формы response()
        # по ВСЕЙ сетке каналов (пик+континуум+все вторичные процессы).
        # Баланс sum(shape)==sum(chans) на САМИХ УЗЛАХ сетки (до сдвига/
        # масштаба/уширения под конкретную E) проверен инвариантом
        # load_grid_nodes; ИТОГОВЫЙ инвариант ПОСЛЕ response() (то, что
        # реально используется здесь) отдельно инвариантом кода не
        # проверяется — проверен ЧИСЛЕННО при аудите (5 энергий, машинный
        # ноль расхождения), но это разовая проверка, не постоянный guard.
        depl = {}
        for E1s, E2s, nuc_keys, I1s, I2s, _note_s, fb_pct_s in sums:
            shp1s, _, _ = resp(E1s)
            shp2s, _, _ = resp(E2s)
            eps1s_tot = float(shp1s.sum())
            eps2s_tot = float(shp2s.sum())
            fb_frac_s = fb_pct_s / 100.0
            k1 = (nuc_keys, round(E1s, 3))
            k2 = (nuc_keys, round(E2s, 3))
            depl[k1] = depl.get(k1, 0.0) + (
                BR_of[nuc_keys] * (I1s / 100.0) * (I2s / 100.0)
                * eps2s_tot / fb_frac_s)
            depl[k2] = depl.get(k2, 0.0) + (
                BR_of[nuc_keys] * (I2s / 100.0) * (I1s / 100.0)
                * eps1s_tot / fb_frac_s)

        for E, I_pct, nuc_key, note in library:
            shp, chans, eps = resp(E)
            w = BR_of[nuc_key] * (I_pct / 100.0)
            w_depl = depl.get((nuc_key, round(E, 3)), 0.0)
            if w_depl > 0:
                # Доля веса ЭТОЙ строки библиотеки, не всего нуклида —
                # правка формулировки по той же находке аудита (было
                # ошибочно подписано «от нуклида»).
                note = (note + "; " if note else "") + (
                    "истощена совпадением на %.3f %% от линии (R78/R85)"
                    % (100.0 * w_depl / max(w, 1e-30)))
                w = max(0.0, w - w_depl)
            add(nuc_key, w, shp, chans)
            photon_lines.append({
                "E_keV": E, "nuclide": nuc_key, "I_gamma_pct": I_pct,
                # Выход линии дан НА РАСПАД СВОЕГО НУКЛИДА — так он и стоит в
                # ENSDF. Ветвление от родителя ряда — отдельное число, и в
                # модель оно входит множителем (BR_of выше). Складывать их в
                # одно поле нельзя: читатель сверяет выход с библиотекой.
                "branch": BR_of[nuc_key],
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
                        "%d энергий, %.1f-%.1f кэВ; отбор по модели-родителю "
                        "трека Geant4 (не по энергетическому окну, R69)"
                        % (len(xray_emit), XRAY_SPAN_LO, XRAY_SPAN_HI),
                "eps_peak": xray_w_total / max(sum(xray_emit.values()), 1e-30),
                "weight_per_branch": xray_w_total, "kind": "xray"})

        n_sum_used = 0
        for E1, E2, nuc_key, I1_pct, I2_pct, note, fb_pct in sums:
            Esum = E1 + E2
            if Esum > E_FIT_HI:
                continue
            _, _, eps1 = resp(E1)
            _, _, eps2 = resp(E2)
            fb_frac = fb_pct / 100.0
            # Форма — с узла, ближайшего к СУММАРНОЙ энергии Esum: континуум
            # суммарного пика физически размазан похоже на континуум
            # одиночного кванта той же полной энергии (то же приближение
            # первого порядка, что и у обычных линий). Абсолютная величина —
            # честное произведение eps1*eps2 (эффективность КАЖДОГО из двух
            # квантов на своей энергии), НЕ eps_peak(Esum), делённое на F_B
            # уровня — та же поправка на неполное заселение, что и в
            # истощении одиночных линий выше (одна и та же физическая
            # величина с двух сторон, см. комментарий там).
            shp, chans, eps_sum_node = resp(Esum)
            w = (BR_of[nuc_key] * (I1_pct / 100.0) * (I2_pct / 100.0)
                 * eps1 * eps2 / max(eps_sum_node, 1e-30) / fb_frac)
            add(nuc_key, w, shp, chans)
            photon_lines.append({
                "E_keV": Esum, "nuclide": nuc_key, "I_gamma_pct": None,
                "branch": BR_of[nuc_key],
                "note": note, "eps_peak": eps1 * eps2,
                "weight_per_branch": BR_of[nuc_key] * (I1_pct / 100.0)
                                      * (I2_pct / 100.0) * eps1 * eps2
                                      / fb_frac,
                "kind": "sum", "E1_keV": E1, "E2_keV": E2,
                "I1_pct": I1_pct, "I2_pct": I2_pct, "fb_pct": fb_pct})
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
        }, (shape_total * A_ph * T + bg_scaled * bg_amp), stack2, stack2_chan, \
            by_nuc_w, shape_total

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
        templ_total, by_nuc, n_eff = build_templates(True)
        m1, model1 = run_method1(templ_total, by_nuc)
        resp = make_full_response(os.path.join(BUILD, "grid"), ch_edges,
                                   True, eps_peak)
        m2, model2, stack2, stack2_chan, by_nuc_w2, shape2 = run_method2(
            GAMMA_LIBRARY, SUM_PEAKS, resp, with_diag=True)
        m2f, model2f, stack2f, stack2f_chan, by_nuc_w2f, shape2f = run_method2(
            lib_full, SUM_PEAKS, resp, with_diag=False)
        if tag == "lines":
            tcs_report(templ_total, shape2f, e, ch_edges)
        stack = {k: (by_nuc[k] * m1["A_Bq"] * T).tolist() for k in keys}
        # Сторож статистики (R66/R76): порог доверия к доле нуклида в канале.
        # N_EFF_MIN = 4 отсчёта МК -> относительная пуассоновская погрешность
        # доли 50 %. Ниже порога слой рисуется пунктиром без заливки: значение
        # остаётся в данных (баланс Σ слоёв = сумме не рвётся), но читателю
        # видно, что это шум шаблона, а не структура отклика.
        N_EFF_MIN = 4.0
        trusted = {k: resolve_mask(n_eff[k] >= N_EFF_MIN, e) for k in keys}
        # Мера «сколько слоя ненадёжно» — доля ИНТЕГРАЛА ниже порога, а не
        # доля каналов. Считать каналы нельзя: у нуклида с короткой шкалой
        # (Pb-212 обрывается на 479 кэВ) почти все каналы пусты по физике, и
        # счёт каналов объявил бы ненадёжными 89 % шаблона, надёжного на
        # 99,7 % по вкладу. Ноль — не недостоверность, а отсутствие вклада.
        noise_frac = {}
        for k in keys:
            v = by_nuc[k]
            m = ~np.asarray(trusted[k], dtype=bool)   # ПОСЛЕ огрубления до
            noise_frac[k] = float(v[m].sum()          # разрешения прибора,
                                  / max(v.sum(), 1e-30))  # иначе метка в
                                                      # легенде разойдётся с
                                                      # тем, что на графике
            if noise_frac[k] > 0.01:
                print("   сторож статистики: %-6s %.1f %% интеграла слоя "
                      "ниже %g отсчётов МК на канал" % (k, 100 * noise_frac[k],
                                                        N_EFF_MIN), flush=True)
        variants[tag] = {
            "trusted": trusted, "n_eff_min": N_EFF_MIN,
            "noise_frac": noise_frac,
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

    # ── состав матрицы пробы: из ПОСТРОЕННОЙ геометрии ────────────────────
    # Имя «ОИСН-16» состав не определяет: под ним в комплекте ходят две
    # разные рецептуры (G1SDetector::MakeMatrix), и ни одна не подтверждена
    # первичным документом. Значит на страницу идут не переписанные руками
    # доли, а те, которыми прогон ДЕЙСТВИТЕЛЬНО посчитан — и снимаются они
    # теми же аргументами запуска, что стоят в шапке шаблонов.
    # Читается ФАЙЛ выгрузки, а не запускается exe: экспорт страницы не должен
    # требовать окружения Geant4 (без него запуск падает с 0xC0000135 — DLL не
    # найдены, и падение выглядит как ошибка данных, хотя дело в PATH).
    def _run_args(path, what):
        with open(path, encoding="utf-8") as fh:
            for ln in fh:
                if not ln.startswith("#"):
                    break
                if "run_args" in ln:
                    return ln.split("=", 1)[1].split(";")[0].split()
        raise SystemExit("в шапке %s нет run_args — %s" % (path, what))

    tpl_args = _run_args(os.path.join(BUILD, "iso_Tl208.csv"),
                         "перепрогнать шаблоны текущей сборкой")
    dump_path = os.path.join(BUILD, "geom_dump.csv")
    if not os.path.isfile(dump_path):
        raise SystemExit(
            "нет %s — состав матрицы взять неоткуда. Снять выгрузку тем же "
            "запуском, что и шаблоны:\n"
            "  $env:G1S_DUMP_GEOM='geom_dump.csv'; .\\g1s.exe vtest.mac %s"
            % (dump_path, " ".join(tpl_args[1:])))
    dump_args = _run_args(dump_path, "снять выгрузку текущей сборкой")
    # Сверка ИМЕННО режима и материала: имя макроса у выгрузки своё (она
    # снимается пустым прогоном), а всё остальное обязано совпасть — иначе
    # страница назовёт состав, которым шаблоны не считались.
    if dump_args[1:] != tpl_args[1:]:
        raise SystemExit(
            "выгрузка геометрии снята другими аргументами, чем шаблоны:\n"
            "  шаблоны: %s\n  выгрузка: %s"
            % (" ".join(tpl_args[1:]), " ".join(dump_args[1:])))
    matrix_name = tpl_args[3] if len(tpl_args) > 3 else "?"
    matrix_rho, matrix_els = None, []
    with open(dump_path, encoding="utf-8") as fh:
        for ln in fh:
            p = ln.rstrip("\n").split(",")
            if len(p) >= 5 and p[0] == "MAT" and p[1] == "Sample":
                matrix_rho = float(p[2])
                matrix_els.append((p[3], float(p[4])))
    if matrix_rho is None:
        raise SystemExit("в %s нет состава материала Sample" % dump_path)
    matrix_els.sort(key=lambda t: -t[1])

    # ── упаковка JSON ─────────────────────────────────────────────────────
    CS = variants["cs"]
    data = {
        "meta": {
            "detector": "Гамма-1С (УДС-ГЦ-63х63)",
            "vessel": "Маринелли 1 л, ОИСН-16 ρ=1,6 г/см³",
            "start_time": meas["start"],
            "live_s": meas["live_s"], "real_s": meas["real_s"],
            "bg_source": "BackgroundEnergySpectrum того же XML",
            "bg_live_s": bg["live_s"], "bg_real_s": bg["real_s"],
            "bg_scale_time": bg_scale_time,
            "fwhm662_keV": FWHM662, "escape_keV": ESCAPE_KEV,
            # Больше не окно отбора (снято в R69) — фактический диапазон
            # энергий, которые модель Geant4 разобрала как атомную релаксацию.
            "xray_span_lo_keV": XRAY_SPAN_LO, "xray_span_hi_keV": XRAY_SPAN_HI,
            # Число розыгрышей у звеньев РАЗНОЕ: слабым добавлено статистики
            # (decay_th232_isotopes_hi.mac, R66), и писать общее «200 000»
            # значило бы называть чужое число. Берётся из шапок самих файлов.
            "template_source": "iso_*.csv, распадов на нуклид: "
                               + ", ".join("%s %d" % (n[0], hist_iso[n[0]][1])
                                           for n in NUCS),
            # То же самое структурированно — для прозы страницы (build_page
            # группирует по числу розыгрышей и подставляет русское имя,
            # текст руками не пишет).
            "template_decays": [{"nuclide": n[1], "n": hist_iso[n[0]][1]}
                                for n in NUCS],
            # Перечисление нуклидов ветви для прозы (t-method) — раньше
            # набиралось текстом руками мимо системы токенов и молча
            # расходилось бы с NUCS при смене состава (найдено этапом 5
            # обобщения конвейера, задача #175/#176, агент C).
            "nuclide_list_ru": ", ".join(n[1] for n in NUCS),
            # Состав матрицы — из ПОСТРОЕННОЙ геометрии, а не из текста
            # шаблона: под именем «ОИСН-16» в комплекте ходят две разные
            # рецептуры, и страница обязана называть ту, которой посчитано.
            "matrix_name": matrix_name,
            "matrix_density_g_cm3": matrix_rho,
            "matrix_composition": [{"element": s, "mass_fraction": w}
                                   for s, w in matrix_els],
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
            "Bq_per_kg": PASSPORT_BQ_KG,
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
            # Сторож статистики (R66/R76): по каналу на нуклид — набран ли
            # его МК-шаблон достаточно, чтобы доле в этом канале верить.
            # false — доля определяется пуассоновским шумом шаблона, и слой
            # там рисуется без заливки. Значение в stack при этом сохранено:
            # выбрасывать его нельзя, иначе сумма слоёв перестанет сходиться
            # с полным откликом.
            "trusted": V["trusted"],
            "n_eff_min": V["n_eff_min"],
            # Доля ИНТЕГРАЛА слоя, попавшая в недостоверную область. Именно
            # она идёт в легенду: доля КАНАЛОВ обманывает — у нуклида с
            # короткой шкалой почти все каналы пусты по физике.
            "noise_frac": V["noise_frac"],
            "stack2": V["stack2"],
            "stack2_full": V["stack2_full"],
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
                         # Своя маска сторожа, а не заимствованная у закона
                         # «по линиям»: n_eff считается ПОСЛЕ свёртки, и при
                         # другой ширине линии тот же набор отсчётов МК
                         # размазывается по другому числу каналов — граница
                         # доверия смещается вместе с законом.
                         "trusted": CS["trusted"],
                         "n_eff_min": CS["n_eff_min"],
                         "noise_frac": CS["noise_frac"],
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
