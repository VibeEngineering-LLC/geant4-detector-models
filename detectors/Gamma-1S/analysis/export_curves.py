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
  eps_net   — площадь пика за вычетом левой полки континуума (окно ±6 кэВ,
              полка E−30…E−10 кэВ). Так строит кривую compare_lsrm.py, и с
              этим числом сверялись с ЛСРМ.
  eps_gross — та же площадь без вычета полки. Так считает compare_point.py.
Разница между ними — систематика способа взятия площади, а не разброс
расчёта. Публиковать одно число, умолчав о втором, значило бы выдать выбор
обработки за свойство детектора.

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
import paths  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))
OUT = str(paths.results("Gamma-1S"))

# Окно ППП и полка континуума — те же, что в compare_lsrm.py. Расчёт без
# уширения, пик острый; края окна учитывают утечку в соседний канал.
WIN = 6.0
BG0, BG1 = 30.0, 10.0

# Описание сеток: метка -> (сосуд, матрица, плотность г/см3, объём мл,
#                           режим защиты, чем задана, пояснение)
# Значения обязаны совпадать с drivers/run_grid.py и drivers/run_all_grids.py.
GRIDS = [
    ("rho1.00", "Маринелли 1 л", "ОИСН-16", 1.00, 1000.0, "закрыта",
     "run_grid.py", "нижняя плотность для подгонки d_eff"),
    ("rho1.60", "Маринелли 1 л", "ОИСН-16", 1.60, 1000.0, "закрыта",
     "run_grid.py", "рабочая кривая, сверялась с .efr ЛСРМ"),
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
    ("p5cm", "точечный, 5 см", "—", None, None, "закрыта",
     "run_all_grids.py", "источник внутри полости защиты"),
    ("p25cm", "точечный, 25 см", "—", None, None, "ОТКРЫТА",
     "run_all_grids.py", "источник над защитой; с закрытой крышкой "
     "50 мм свинца дают почти нулевой счёт"),
]


def read_run(path):
    """(E0, чистая площадь, погрешность, N первичных, полная площадь)."""
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
    gross = sum(c for e, c in hist.items() if abs(e - E0) <= WIN)
    side = sum(c for e, c in hist.items() if E0 - BG0 <= e <= E0 - BG1)
    nside = BG0 - BG1                      # ширина полки в каналах по 1 кэВ
    bg = side / nside * (2 * WIN + 1)
    var = gross + (bg / nside) * bg        # пуассон пика + шум полки
    return E0, gross - bg, math.sqrt(max(var, 1.0)), N, gross


def solid_angle(tag):
    """Доля телесного угла для точечных прогонов; 1.0 для объёмных."""
    p = os.path.join(BUILD, "grid", "%s_solidangle.txt" % tag)
    if os.path.exists(p):
        return float(open(p).read().strip())
    return 1.0


def curve(tag):
    frac = solid_angle(tag)
    rows = []
    for p in sorted(glob.glob(os.path.join(BUILD, "grid", tag + "_E*.csv"))):
        r = read_run(p)
        if r is None:
            continue
        E0, net, dnet, N, gross = r
        if net <= 0:
            continue
        rows.append({
            "E_keV": round(E0, 3),
            "eps_net": net / N * frac,
            "d_eps": dnet / N * frac,
            "eps_gross": gross / N * frac,
            "N_primaries": N,
            "net_counts": round(net, 1),
            "solid_angle_fraction": frac,
        })
    rows.sort(key=lambda r: r["E_keV"])
    return rows


def write_csv(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(",".join(header) + "\n")
        for r in rows:
            fh.write(",".join(str(r[h]) for h in header) + "\n")


def export_curves():
    os.makedirs(OUT, exist_ok=True)
    long_rows = []
    made = []
    for (tag, vessel, matrix, rho, vol, lid, driver, note) in GRIDS:
        rows = curve(tag)
        if not rows:
            print("нет сетки %s — пропущена" % tag)
            continue
        head = ["E_keV", "eps_net", "d_eps", "eps_gross", "N_primaries",
                "net_counts", "solid_angle_fraction"]
        write_csv(os.path.join(OUT, "eff_%s.csv" % tag), head, rows)
        made.append((tag, len(rows)))
        for r in rows:
            long_rows.append(dict(r, geometry=vessel, matrix=matrix,
                                  density_g_cm3="" if rho is None else rho,
                                  fill_ml="" if vol is None else vol,
                                  shield_lid=lid, grid=tag))
    if long_rows:
        head = ["grid", "geometry", "matrix", "density_g_cm3", "fill_ml",
                "shield_lid", "E_keV", "eps_net", "d_eps", "eps_gross",
                "N_primaries", "net_counts", "solid_angle_fraction"]
        long_rows.sort(key=lambda r: (r["grid"], r["E_keV"]))
        write_csv(os.path.join(OUT, "efficiency_curves.csv"), head, long_rows)
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
               "note"], rows)
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
            peak = sum(c for e, c in hist.items() if abs(e - E) <= WIN)
            side = sum(c for e, c in hist.items() if E - BG0 <= e <= E - BG1)
            nside = BG0 - BG1
            bg = side / nside * (2 * WIN + 1)
            net = peak - bg
            dnet = math.sqrt(max(peak + (bg / nside) * bg, 1.0))
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
                   "eps_decay", "eps_mono", "C_summing", "d_C", "note"], rows)
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
