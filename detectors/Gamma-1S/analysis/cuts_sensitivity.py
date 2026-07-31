"""Чувствительность эффективности к порогу рождения вторичных (задача 100).

ЗАЧЕМ. Порог рождения (`SetDefaultCutValue`, у модели 0,05 мм) задаёт, с какой
энергии Geant4 рождает вторичный электрон или фотон отдельной частицей; ниже
порога энергия выделяется на месте. Величина не измеряется, а НАЗНАЧАЕТСЯ, и
пока её влияние не измерено, она остаётся неучтённой систематикой выбора
расчётчика. Мягкий край — то место, где эффект наиболее вероятен: там пробеги
вторичных электронов сопоставимы с масштабами слоёв входного торца.

ЧТО СРАВНИВАЕТСЯ. Четыре порога — 0,01 / 0,05 / 0,50 / 1,00 мм — на трёх
мягких линиях (59,5 / 88,0 / 122,1 кэВ) и на 661,7 как контроле, где эффекта
быть не должно. Отношения берутся к РАБОЧЕМУ порогу 0,05 мм, поэтому
геометрия, телесный угол и нормировка из сравнения уходят.

ПОЧЕМУ БЕЗ ПЕРЕСБОРКИ. `/run/setCut` — штатная команда Geant4: порог меняется
в работающем исполняемом файле. Правка `main.cc` сменила бы отпечаток
исходников и обесценила все посчитанные сетки (задача 139); здесь этого не
требуется — тот же exe, тот же отпечаток, различаются только прогоны.

ЧТО ЭТОТ ТЕСТ НЕ ГОВОРИТ. Он измеряет ЧУВСТВИТЕЛЬНОСТЬ к параметру, а не
правильность его значения. Малая чувствительность означает, что выбор порога
не входит в бюджет неопределённостей заметным членом; большая означала бы, что
порог надо обосновывать отдельно, а не назначать.
"""
import glob
import math
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
import paths  # noqa: E402
import peakwin  # noqa: E402
import stamp  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))
OUT = str(paths.results("Gamma-1S"))
BASE = "0.05"          # рабочий порог модели, мм

OBS = {
    "quantity": "отношение эффективности ППП при данном пороге рождения к"
                " эффективности при рабочем пороге 0.05 мм",
    "area": "чистая площадь пика за вычетом полки; одно правило на все"
            " пороги (common/py/peakwin.area)",
    "window": "+-6 кэВ в каналах; полка [E-25; E-10] — конвенция peakwin",
    "shelf": "односторонняя слева; вычитается одинаково на всех порогах",
    "blurred": "нет — депозит-спектры как есть",
    "geometry": "точечный источник z = 91 мм; полный 4pi; 400 тыс. первичных"
                " на точку — как в рабочих сетках",
}


def read_hist(path):
    hist, N = {}, None
    for ln in open(path, encoding="utf-8"):
        if ln.startswith("#"):
            if "N_primaries" in ln:
                N = int(ln.split("=")[1])
            continue
        if ln and ln[0].isdigit():
            e, c = ln.split(",")
            hist[float(e)] = float(c)
    return hist, N


def main():
    files = sorted(glob.glob(os.path.join(BUILD, "cut*_E*.csv")))
    if not files:
        raise SystemExit(
            "Не найдены прогоны cut*_E*.csv в %s.\n"
            "Прогон: g1s.exe cuts_scan.mac shield" % BUILD)
    data = {}
    for p in files:
        m = re.match(r"cut([0-9.]+)_E0*([0-9.]+)\.csv",
                     os.path.basename(p))
        if not m:
            continue
        cut, E = m.group(1), float(m.group(2))
        hist, N = read_hist(p)
        a = peakwin.area(hist, E)
        data.setdefault(E, {})[cut] = (a / N, math.sqrt(max(a, 1.0)) / N)

    cuts = sorted({c for v in data.values() for c in v},
                  key=lambda s: float(s))
    print("Чувствительность к порогу рождения вторичных."
          " Отношение к рабочему порогу %s мм.\n" % BASE)
    print("%10s %s" % ("E, кэВ", "  ".join("%9s мм" % c for c in cuts)))
    rows = []
    for E in sorted(data):
        if BASE not in data[E]:
            print("%10.1f  нет прогона на рабочем пороге — пропущена" % E)
            continue
        e0, d0 = data[E][BASE]
        cells = []
        for c in cuts:
            if c not in data[E]:
                cells.append("%12s" % "—")
                continue
            e, d = data[E][c]
            r = e / e0
            dr = r * math.sqrt((d / e) ** 2 + (d0 / e0) ** 2)
            cells.append("%7.4f+-%.4f" % (r, dr))
            rows.append((E, c, e, r, dr))
        print("%10.1f %s" % (E, "  ".join(cells)))

    print("\nВЫВОД.")
    worst_E, worst_dev, worst_sig = None, 0.0, 0.0
    for E, c, e, r, dr in rows:
        if c == BASE:
            continue
        dev = abs(r - 1.0)
        sig = dev / dr if dr > 0 else float("inf")
        if dev > worst_dev:
            worst_E, worst_dev, worst_sig = (E, c), dev, sig
    if worst_E is None:
        print("  Нет данных для сравнения.")
    else:
        print("  Наибольшее отклонение от рабочего порога: %.1f кэВ при"
              " пороге %s мм — %.2f %% (%.1f сигмы)."
              % (worst_E[0], worst_E[1], 100 * worst_dev, worst_sig))
        if worst_sig < 3.0:
            # У отрицательного результата обязан быть предел обнаружения:
            # «влияния нет» без указания; НАСКОЛЬКО малого влияния хватило бы;
            # чтобы его не заметить; — утверждение ни о чём.
            dr_typ = sorted(dr for E, c, e, r, dr in rows if c != BASE)
            lim = 100 * 2.0 * dr_typ[len(dr_typ) // 2]
            print("  Значимого влияния НЕТ во всём проверенном диапазоне"
                  " порогов — от 0;01 до 1;0 мм; то есть")
            print("  в двадцать раз грубее рабочего.")
            print("  ПРЕДЕЛ ОБНАРУЖЕНИЯ ЭТОГО ТЕСТА: %.1f %% (2 сигмы при"
                  " 400 тыс. первичных на точку)." % lim)
            print("  Эффект больше %.1f %% был бы виден; меньше — нет. В"
                  " бюджете неопределённостей порог рождения" % lim)
            print("  учитывается этой границей; а не нулём.")
        else:
            print("  Влияние ЗНАЧИМО. Порог рождения перестаёт быть свободным"
                  " параметром: его значение надо")
            print("  обосновывать физически; а не назначать; и включать"
                  " в бюджет неопределённостей.")

    csvio.write(
        os.path.join(OUT, "cuts_sensitivity.csv"),
        ["E_keV", "cut_mm", "eps", "ratio_to_base", "d_ratio"],
        [("%.1f" % E, c, "%.6e" % e, "%.4f" % r, "%.4f" % dr)
         for E, c, e, r, dr in rows],
        comments=[
            "Порог рождения вторичных задан командой /run/setCut в работающем"
            " exe — пересборки нет; отпечаток исходников не менялся.",
            "Отношение берётся к рабочему порогу 0;05 мм: геометрия и"
            " нормировка из сравнения уходят.",
            "Мягкие линии 59;5 / 88;0 / 122;1 — там эффект наиболее вероятен;"
            " 661;7 — контроль; где его быть не должно.",
            "Тест измеряет ЧУВСТВИТЕЛЬНОСТЬ к параметру; а не правильность"
            " его значения.",
        ],
        stamp=stamp.lines(
            "detectors/Gamma-1S/analysis/cuts_sensitivity.py", OBS,
            inputs=files,
            geometry_dir=str(paths.geometry("Gamma-1S")),
            names=stamp.SRC_LISTS["Gamma-1S"], repo_dir=str(paths.REPO)))
    print("\nтаблица: %s" % os.path.join(OUT, "cuts_sensitivity.csv"))


if __name__ == "__main__":
    main()
