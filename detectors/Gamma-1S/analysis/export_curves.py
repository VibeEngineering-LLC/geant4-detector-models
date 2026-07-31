"""Выгрузка расчётных кривых ППП в таблицы, пригодные к повторному использованию.

Расчётные спектры (сотни файлов, десятки мегабайт) в репозиторий не
коммитятся: они воспроизводятся драйверами. Но САМИ КРИВЫЕ должны лежать
готовыми — иначе, чтобы узнать эффективность на 662 кэВ, надо собрать Geant4
и посчитать сутки. Этот скрипт превращает сетки прогонов в таблицы в
detectors/Gamma-1S/results/.

    python detectors/Gamma-1S/analysis/export_curves.py

Что получается:

  efficiency_curves.csv   все кривые одной длинной таблицей (для обработки)
  eff_<метка>.csv         по кривой на геометрию (для чтения и построения)
  runs_manifest.csv       чем и с какими доводами посчитана каждая сетка
  summing_C.csv           поправки на каскадное суммирование по линиям

ЧТО ТАКОЕ eps В ЭТИХ ТАБЛИЦАХ. Абсолютная эффективность регистрации в пике
полного поглощения: доля испущенных в пробе квантов данной энергии, давших
отсчёт в ППП. На активность и на выход линии НЕ поделено — это чистая
характеристика «проба+сосуд+детектор+защита».

ДВА СТОЛБЦА ЭФФЕКТИВНОСТИ, И ЭТО НАМЕРЕННО.
  eps_net   — площадь пика за вычетом левой полки континуума правилом
              common/py/peakwin (окно ±6 кэВ, полка [E−25; E−10], счёт в
              каналах). Тем же правилом строит кривую compare_lsrm.py.
  eps_gross — та же площадь без вычета полки. Так считает compare_point.py.
Разница между ними — систематика способа взятия площади, а не разброс
расчёта. Публиковать одно число, умолчав о втором, значило бы выдать выбор
обработки за свойство детектора.

СТОЛБЕЦ in_range. Паспортный диапазон регистрируемых энергий прибора —
50…3000 кэВ (п. 2.2), а сетка выходит за него в обе стороны: 45,3 кэВ снизу,
3304,8 и 3552,5 сверху. За границами прибор не аттестован, поэтому согласие
или расхождение расчёта с измерением там не довод ни за модель, ни против
неё. Узлы НЕ удаляются — они годятся для сверки с другим кодом (EffCalcMC), —
но помечаются: in_range = 1 внутри диапазона, 0 вне его. Границы берутся из
detector_params.ATTESTED_RANGE_KEV, то есть оттуда же, откуда их читает любой
другой скрипт.

ТОЧЕЧНЫЕ ГЕОМЕТРИИ СЧИТАНЫ В КОНУС. Изотропный источник на 25 см тратит
99,4 % событий впустую, поэтому кванты разыгрывались в конус вокруг
направления на детектор, а эффективность приведена к полному телесному углу
умножением на его долю (файл grid/<метка>_solidangle.txt). Для ППП это
законно: в пик идут практически только прямые кванты. Полная (не пиковая)
эффективность при таком розыгрыше была бы занижена.
"""
import glob
import math
import os
import re
import sys

# Корни путей — из переменных окружения (common/py/paths.py), чтобы в коде не
# было ни одного пути, привязанного к конкретной машине.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
import paths  # noqa: E402
import stamp  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import detector_params as dp  # noqa: E402
import peakwin  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))
OUT = str(paths.results("Gamma-1S"))

# Окно ППП и полка континуума — ЕДИНСТВЕННАЯ реализация common/py/peakwin
# (полка [E−25; E−10], счёт в каналах). Прежде здесь жили ДВЕ собственные
# копии правила с полкой E−30, захватывавшей пик вылета иода E−28,6, — тот
# дефект, который §4 отчёта объявлял устранённым; найдено внутренним аудитом
# 31.07.2026 (сверх ТЗ п.1-2). Заодно лгало объявление наблюдаемой в
# summing_C.csv («в каналах; [E-25;E-10]» при коде «в кэВ; E-30»). Констант
# WIN/BG0/BG1 здесь больше нет; объявление — из peakwin.declare().

# Описание сеток: метка -> (сосуд, матрица, плотность г/см3, объём мл,
#                           режим защиты, чем задана, пояснение)
# Значения обязаны совпадать с drivers/run_grid.py и drivers/run_all_grids.py.
#
# ЗАПЯТЫХ ВНУТРИ ЗНАЧЕНИЙ ЗДЕСЬ БОЛЬШЕ НЕТ, И ЭТО НАМЕРЕННО. Запись идёт через
# csv.writer, который такое поле закавычит, — но потребитель читает выгрузку
# СВОИМ кодом, и кавычки его не спасут, если он разбирает строку простым
# split(","). Правильное экранирование обязательно, а отсутствие запятой в
# значении — второй рубеж, который работает и для наивного читателя. Было
# «точечный, 25 см» и «рабочая кривая, сверялась с .efr ЛСРМ».
GRIDS = [
    ("rho1.00", "Маринелли 1 л", "ОИСН-16", 1.00, 1000.0, "закрыта",
     "run_grid.py", "нижняя плотность для подгонки d_eff"),
    ("rho1.60", "Маринелли 1 л", "ОИСН-16", 1.60, 1000.0, "закрыта",
     "run_grid.py", "рабочая кривая; сверялась с .efr ЛСРМ"),
    ("water1.00", "Маринелли 1 л", "вода", 1.00, 1000.0, "закрыта",
     "run_all_grids.py", "проверка влияния состава матрицы"),
    ("denta0.60", "Дента 120 мл", "ОИСН-16", 0.60, 120.0, "закрыта",
     "run_all_grids.py", ""),
    ("denta1.60", "Дента 120 мл", "ОИСН-16", 1.60, 120.0, "закрыта",
     "run_all_grids.py", ""),
    ("petri0.60", "Петри 60 мл", "ОИСН-16", 0.60, 60.0, "закрыта",
     "run_all_grids.py", ""),
    ("petri1.60", "Петри 60 мл", "ОИСН-16", 1.60, 60.0, "закрыта",
     "run_all_grids.py", ""),
    ("p5cm", "точечный 5 см", "—", None, None, "закрыта",
     "run_all_grids.py", "источник внутри полости защиты"),
    ("p25cm", "точечный 25 см", "—", None, None, "ОТКРЫТА",
     "run_all_grids.py", "источник над защитой; с закрытой крышкой "
     "50 мм свинца дают почти нулевой счёт"),
]


def read_run(path):
    """(E0, чистая площадь, погрешность, N первичных, полная площадь).

    Площадь — общей реализацией `peakwin.area`; здесь остаётся только чтение
    файла и погрешность. Дисперсия чистой площади: пуассон окна плюс
    перенесённый шум полки — bg = side·(n/nside), D(bg) = (n/nside)²·side =
    (n/nside)·bg (вывод см. common/py/becqmoni.py).
    """
    N = E0 = None
    hist = {}
    for line in open(path, encoding="utf-8"):
        if line.startswith("#"):
            if "N_primaries" in line:
                N = int(line.split("=")[1])
            elif "E_prim_keV" in line:
                E0 = float(line.split("=")[1])
            continue
        if line and line[0].isdigit():
            e, c = line.split(",")
            hist[float(e)] = int(c)
    if E0 is None or N is None:
        return None
    det = {}
    net = peakwin.area(hist, E0, detail=det)
    gross, side = det["gross"], det["side"]
    n, nside = det["n_peak"], det["n_side"]
    bg = side / nside * n if nside else 0.0
    var = gross + (n / nside) * bg if nside else gross
    return E0, net, math.sqrt(max(var, 1.0)), N, gross


def solid_angle(tag):
    """Доля телесного угла для точечных прогонов; 1.0 для объёмных."""
    p = os.path.join(BUILD, "grid", "%s_solidangle.txt" % tag)
    if os.path.exists(p):
        return float(open(p).read().strip())
    return 1.0


# Объявление наблюдаемой для выгружаемых кривых. Одно на все eff_*.csv и на
# сводную efficiency_curves.csv: обе таблицы содержат ОДНУ величину, снятую
# ОДНИМ правилом. Строки окна/полки/размытия СОБИРАЮТСЯ из констант peakwin
# (peakwin.declare()) — набранное руками объявление уже лгало один раз
# (внутренний аудит 31.07.2026, сверх ТЗ п.2).
OBS = dict(
    {
        "quantity": "абсолютная эффективность регистрации в ППП — доля"
                    " испущенных квантов данной энергии; давших отсчёт"
                    " в пике",
        "area": "eps_net — площадь пика за вычетом полки континуума"
                " (правило common/py/peakwin); eps_gross — та же площадь"
                " без вычета (столбцы намеренно оба)",
        "solid_angle": "точечные сетки разыграны в конус; eps приведена к"
                       " 4pi множителем из grid/<метка>_solidangle.txt",
        "in_range": "1 внутри паспортного диапазона 50…3000 кэВ (п. 2.2);"
                    " 0 вне его — там прибор не аттестован",
    },
    **peakwin.declare())


def curve(tag, dropped=None, used=None):
    """Кривая по сетке. dropped — сюда складываются ПРОПУЩЕННЫЕ узлы.

    Пропуск обязан быть слышен. Узел с нулевой площадью молча выпадал из
    кривой, и так потерялись две верхние точки сетки (3304,8 и 3552,5 кэВ):
    гистограмма модели была обрезана на 3200 кэВ, пик полного поглощения
    уезжал в канал переполнения. Ни прогон, ни выгрузка не сказали ни слова,
    а кривая просто оказалась короче, чем задумано.
    """
    frac = solid_angle(tag)
    rows = []
    for p in sorted(glob.glob(os.path.join(BUILD, "grid", tag + "_E*.csv"))):
        r = read_run(p)
        if r is None:
            if dropped is not None:
                dropped.append((tag, os.path.basename(p), "нет заголовка"))
            continue
        E0, net, dnet, N, gross = r
        if net <= 0:
            if dropped is not None:
                dropped.append((tag, "%.1f кэВ" % E0,
                                "чистая площадь ≤ 0 — пик вне гистограммы?"))
            continue
        rows.append({
            "E_keV": round(E0, 3),
            "eps_net": net / N * frac,
            "d_eps": dnet / N * frac,
            "eps_gross": gross / N * frac,
            "N_primaries": N,
            "net_counts": round(net, 1),
            "solid_angle_fraction": frac,
            # 1/0, а не текст: столбец читается фильтром потребителя. Признак
            # ВЫВОДИТСЯ из паспортных границ в detector_params, а не вписан
            # списком энергий — вписанный разъедется с границами при первой же
            # правке, и это ровно тот класс дефекта, который в отчёте назван
            # главным уроком линии.
            "in_range": 1 if dp.in_attested_range(E0) else 0,
        })
        if used is not None:
            used.append(p)
    rows.sort(key=lambda r: r["E_keV"])
    return rows


# Манифест и поправки на суммирование — ДРУГИЕ наблюдаемые, и объявление у
# них своё. Один OBS на три разные таблицы означал бы, что штамп лжёт ровно
# там, где заводится против лжи.
OBS_MANIFEST = {
    "quantity": "перечень сеток прогонов: чем и с какими доводами посчитана"
                " каждая — величина НЕ числовая",
    "area": "не применимо — площади здесь не снимаются",
    "window": "не применимо",
    "shelf": "не применимо",
    "blurred": "не применимо",
}
OBS_SUMMING = {
    "quantity": "C — поправка на каскадное суммирование: отношение"
                " эффективности на моноэнергии к эффективности той же линии"
                " в прогоне полного распада",
    "area": "чистая площадь пика за вычетом полки; ОДНО правило на обе"
            " стороны отношения",
    "window": "+-6 кэВ в каналах; полка [E-25; E-10]",
    "shelf": "односторонняя слева; вычитается одинаково с обеих сторон",
    "blurred": "нет — депозит-спектры как есть",
}


def _stamp(inputs, obs=None):
    """Строки штампа для выгружаемой кривой.

    ЗАЧЕМ ОН ЗДЕСЬ. Таблицы этого экспортёра штампа не несли, и это дорого
    обошлось: сетки в каталоге прогонов были пересчитаны 30.07 на исправленной
    геометрии плоских кювет (коммит 2df1eb2) и на текущем exe, а выгрузка в
    results/ осталась от 28.07. Опубликованные кривые Петри и Дента разошлись
    с расчётом на 15…28 %, Маринелли на 8 %, и НИЧТО об этом не сказало:
    незаштампованная таблица внешне не отличается от свежей. Скрипты анализа
    читают grid/ напрямую, поэтому выводы отчёта на этом не стояли, — но
    потребитель выгрузки получал устаревшие числа. Со штампом такая таблица
    сама объявляет, каким отпечатком посчитаны её входы.
    """
    return stamp.lines(
        "detectors/Gamma-1S/analysis/export_curves.py", obs or OBS,
        inputs=inputs,
        geometry_dir=str(paths.geometry("Gamma-1S")),
        names=stamp.SRC_LISTS["Gamma-1S"], repo_dir=str(paths.REPO))


def write_csv(path, header, rows, inputs=None, obs=None):
    """Запись таблицы — общей реализацией из common/py/csvio.py.

    Ручной ",".join не экранирует ничего, и любое поле с запятой рвёт строку.
    Так и вышло: geometry у точечных сеток было «точечный, 25 см», и все 48
    строк обеих точечных кривых в efficiency_curves.csv получили 14 полей при
    шапке в 13. csv.DictReader при этом не падал, а ТИХО сдвигал всё правее на
    поле: E_keV становился словом «ОТКРЫТА», eps_net — числом энергии.

    Своей копии записи и своего сторожа здесь больше нет: и то и другое живёт
    в csvio, потому что копия сторожа успела разъехаться с копией в
    tools/check_csv.py по обращению с комментариями.
    """
    csvio.write(path, header, rows,
                stamp=_stamp(inputs, obs) if obs is not None
                or inputs is not None else ())


def export_curves():
    os.makedirs(OUT, exist_ok=True)
    long_rows = []
    made = []
    dropped = []
    all_inputs = []
    for (tag, vessel, matrix, rho, vol, lid, driver, note) in GRIDS:
        used = []
        rows = curve(tag, dropped, used)
        if not rows:
            print("нет сетки %s — пропущена" % tag)
            continue
        all_inputs.extend(used)
        head = ["E_keV", "eps_net", "d_eps", "eps_gross", "N_primaries",
                "net_counts", "solid_angle_fraction", "in_range"]
        write_csv(os.path.join(OUT, "eff_%s.csv" % tag), head, rows, used)
        made.append((tag, len(rows)))
        for r in rows:
            long_rows.append(dict(r, geometry=vessel, matrix=matrix,
                                  density_g_cm3="" if rho is None else rho,
                                  fill_ml="" if vol is None else vol,
                                  shield_lid=lid, grid=tag))
    if long_rows:
        head = ["grid", "geometry", "matrix", "density_g_cm3", "fill_ml",
                "shield_lid", "E_keV", "eps_net", "d_eps", "eps_gross",
                "N_primaries", "net_counts", "solid_angle_fraction",
                "in_range"]
        long_rows.sort(key=lambda r: (r["grid"], r["E_keV"]))
        write_csv(os.path.join(OUT, "efficiency_curves.csv"), head, long_rows,
                  all_inputs)
    out = sorted({r["E_keV"] for r in long_rows if not r["in_range"]})
    if out:
        lo, hi = dp.ATTESTED_RANGE_KEV
        print("\nВНЕ ПАСПОРТНОГО ДИАПАЗОНА %.0f…%.0f кэВ (in_range=0): %s"
              % (lo, hi, "; ".join("%.1f" % E for E in out)))
        print("   Там прибор не аттестован — расхождение расчёта с измерением"
              " не довод ни за модель, ни против.")
    if dropped:
        print("\nПРОПУЩЕНЫ УЗЛЫ СЕТКИ (%d) — кривая короче задуманной:"
              % len(dropped))
        for tag, what, why in dropped:
            print("   %-12s %-14s %s" % (tag, what, why))
    return made


def export_manifest(made):
    got = dict(made)
    rows = []
    for (tag, vessel, matrix, rho, vol, lid, driver, note) in GRIDS:
        rows.append({
            "grid": tag,
            "geometry": vessel,
            "matrix": matrix,
            "density_g_cm3": "" if rho is None else rho,
            "fill_ml": "" if vol is None else vol,
            "shield_lid": lid,
            "driver": driver,
            "n_points": got.get(tag, 0),
            "solid_angle_fraction": solid_angle(tag),
            "note": note,
        })
    write_csv(os.path.join(OUT, "runs_manifest.csv"),
              ["grid", "geometry", "matrix", "density_g_cm3", "fill_ml",
               "shield_lid", "driver", "n_points", "solid_angle_fraction",
               "note"], rows, obs=OBS_MANIFEST)
    return rows


# Линии, по которым считалась поправка на каскадное суммирование, и прогон
# распада, из которого берётся и площадь пика, и число ИСПУЩЕННЫХ квантов.
# Cs-137 и K-40 каскада не имеют: для них C обязан выйти 1,00 — это контроль
# метода, а не результат.
SUM_LINES = {
    "decay_Cs137": ("Cs-137", [661.657], "контроль: каскада нет"),
    "decay_K40": ("K-40", [1460.822], "контроль: каскада нет"),
    "decay_Tl208": ("Tl-208", [583.187, 2614.511], ""),
    "decay_Bi214": ("Bi-214", [609.32, 768.36, 1120.294, 1764.491], ""),
}


def emitted(run, E):
    """Сколько квантов этой энергии реально испущено за прогон.

    Берётся из того же расчёта (*_emit.csv), а не из справочника выходов:
    иначе сравнивались бы числа из разных источников.
    """
    p = os.path.join(BUILD, "%s_emit.csv" % run)
    if not os.path.exists(p):
        return None, None
    tot, N = 0, None
    for line in open(p, encoding="utf-8"):
        if line.startswith("#"):
            if "N_primaries" in line:
                N = int(line.split("=")[1])
            continue
        if line and line[0].isdigit():
            e, c = line.split(",")
            if abs(float(e) - E) <= 2.0:
                tot += int(c)
    return tot, N


def mono_eff(E):
    """(эффективность, погрешность) из моноэнергетической сетки Маринелли."""
    for p in glob.glob(os.path.join(BUILD, "grid", "rho1.60_E*.csv")):
        r = read_run(p)
        if r and abs(r[0] - E) < 1.0:
            return r[1] / r[3], r[2] / r[3]
    return None, None


def export_summing():
    rows = []
    for run, (nuc, lines, note) in sorted(SUM_LINES.items()):
        p = os.path.join(BUILD, "%s.csv" % run)
        if not os.path.exists(p):
            continue
        hist, N = {}, None
        for line in open(p, encoding="utf-8"):
            if line.startswith("#"):
                if "N_primaries" in line:
                    N = int(line.split("=")[1])
                continue
            if line and line[0].isdigit():
                e, c = line.split(",")
                hist[float(e)] = int(c)
        for E in lines:
            # Та же единственная реализация правила, что и в read_run: прежде
            # здесь была ВТОРАЯ собственная копия (аудит 31.07, сверх ТЗ п.1).
            det = {}
            net = peakwin.area(hist, E, detail=det)
            bg = (det["side"] / det["n_side"] * det["n_peak"]
                  if det["n_side"] else 0.0)
            dnet = math.sqrt(max(
                det["gross"] + (det["n_peak"] / det["n_side"]) * bg
                if det["n_side"] else det["gross"], 1.0))
            nem, _ = emitted(run, E)
            if not nem or net <= 0:
                continue
            eps_decay = net / nem
            d_decay = dnet / nem
            eps_mono, d_mono = mono_eff(E)
            C = dC = ""
            if eps_mono:
                C = eps_mono / eps_decay
                # погрешности двух независимых прогонов складываются
                dC = C * math.sqrt((d_mono / eps_mono) ** 2
                                   + (d_decay / eps_decay) ** 2)
                C, dC = round(C, 4), round(dC, 4)
            rows.append({
                "nuclide": nuc,
                "E_keV": E,
                "run": run,
                "N_decays": N,
                "N_emitted": nem,
                "eps_decay": round(eps_decay, 6),
                "eps_mono": "" if not eps_mono else round(eps_mono, 6),
                "C_summing": C,
                "d_C": dC,
                "note": note,
            })
    if rows:
        write_csv(os.path.join(OUT, "summing_C.csv"),
                  ["nuclide", "E_keV", "run", "N_decays", "N_emitted",
                   "eps_decay", "eps_mono", "C_summing", "d_C", "note"],
                  rows, obs=OBS_SUMMING)
    return rows


if __name__ == "__main__":
    if not os.path.isdir(os.path.join(BUILD, "grid")):
        raise SystemExit(
            "Нет каталога расчётных спектров: %s\n"
            "Задайте G4MODELS_BUILD на каталог сборки или посчитайте сетки:\n"
            "    python detectors/Gamma-1S/drivers/run_grid.py\n"
            "    python detectors/Gamma-1S/drivers/run_all_grids.py"
            % os.path.join(BUILD, "grid"))
    made = export_curves()
    for tag, n in made:
        print("eff_%s.csv — точек: %d" % (tag, n))
    export_manifest(made)
    print("runs_manifest.csv — сеток: %d" % len(GRIDS))
    srows = export_summing()
    print("summing_C.csv — линий: %d" % len(srows))
    print("всё записано в", OUT)
