# -*- coding: utf-8 -*-
"""Отклик двух групп кривых на смену геометрии и корреляция откликов (задача 134).

ВОПРОС ОПЕРАТОРА, три части подряд:
  (1) чем различаются группы — аттестованная ЛСРМ (`.efr`) и наша Geant4-модель —
      в каждой геометрии;
  (2) как меняется кривая ВНУТРИ каждой группы при переходе
      точка -> Петри -> Дента -> Маринелли (именно в этой последовательности);
  (3) коррелируют ли эти изменения между группами: совпадают или есть
      антисовпадения.

ПОЧЕМУ ОТКЛИК СЧИТАЕТСЯ ОТНОШЕНИЕМ К ПРЕДЫДУЩЕЙ ГЕОМЕТРИИ, А НЕ ПО АБСОЛЮТУ.
Абсолютное расхождение групп содержит всё сразу: нормировку, конвенцию съёма
площади, плато ~7…8 % и разрыв на жёстком крае. Вопрос (2) — про ФОРМУ отклика
на геометрию, и она выделяется отношением кривых ОДНОЙ группы к ОДНОЙ опорной
геометрии: тогда всё, что не зависит от геометрии, сокращается тождественно, в
обеих группах одинаково. Это то же соображение, по которому §7 требует наклонов
вместо отношений: сравнивать надо ход, а не точку.

Шаги последовательности берутся ПОСЛЕДОВАТЕЛЬНЫМИ (Петри/точка, Дента/Петри,
Маринелли/Дента), как задал оператор. Каждый шаг — физически разное изменение:
первый добавляет объём и самопоглощение, второй увеличивает объём при той же
плоской укладке, третий переводит источник из «перед торцом» в «вокруг
кристалла».

ЕДИНИЦЫ. Обе группы — эффективность НА ИСПУЩЕННЫЙ КВАНТ (не на распад): `.efr`
по своему формату (`E=eff,dpct,...`), наша сетка — моноэнергии, поделённые на
долю телесного угла для конусных розыгрышей. Выход линии не участвует ни там,
ни там, поэтому библиотека нуклида из сравнения исключена.

ПЛОТНОСТЬ — ПРОЧИТАНА, А НЕ ПРИНЯТА. Условия аттестации записаны в заголовке
самих файлов `.efr`: `Density,g/cm3`, `Material`, `Volume,ml`, `Thick,mm`. Петри
— 60 мл, 1,60 г/см³, слой 10 мм; Дента — 120 мл, 1,658 г/см³, слой 33 мм;
Маринелли — 1000 мл, 1,60 г/см³, слой 31 мм; матрица везде ОИСН-16.
Сопоставимая сетка поэтому одна — на 1,60, и она стоит ПЕРВОЙ в списке меток.
Сетки на 0,60 остаются как проба чувствительности к самопоглощению, сравнением
с аттестацией они не являются.

Первая редакция этого расчёта называла плотность неизвестной и заменяла её
вилкой 0,60/1,60 — ошибка того же класса, что разбирает method-rules:
утверждение выводилось по памяти, тогда как ответ читается из файла одной
командой. Уточнение вывод не смягчило: при сопоставимой плотности расхождение
наклонов БОЛЬШЕ, чем при благоприятной.

Остаточное несовпадение: у Денты 1,658 против 1,60 в сетке — проба легче на
3,5 %; отдельной сетки на 1,658 не считалось.

АНТИСОВПАДЕНИЕ. Ключевой признак — не величина корреляции, а ЗНАК наклона
отклика по энергии. Если у одной группы отклик растёт с энергией, а у другой
падает, модель воспроизводит геометрию с обратным ходом, и никакая нормировка
этого не исправит. Поэтому печатаются оба наклона со своими sigma, и лишь потом
коэффициент корреляции.

ОГОВОРКА О ВХОДАХ. Точечная сетка регенерируется; вердикт провенанса печатается
рядом с числами. Отношения внутри одной группы устойчивее абсолютов, но если
вердикт не `ok`, числа предварительные.
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

sys.path.insert(0, str(paths.tools()))
from fetch_efr import parse_efr  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import point_recalc as pr  # noqa: E402

# Последовательность геометрий — как задал оператор. Для каждой: имя кривой
# ЛСРМ и метки наших сеток. ПЕРВАЯ метка — сопоставимая аттестации по прочитанной
# из `.efr` плотности; вторая, если есть, — проба чувствительности, не сравнение.
SEQ = [
    ("точка 5 см", "Точечная-5см", ("p5cm",)),
    ("Петри 60 мл", "Петри", ("petri1.60", "petri0.60")),
    ("Дента 120 мл", "Дента", ("denta1.60", "denta0.60")),
    ("Маринелли 1 л", "Маринелли", ("rho1.60", "rho1.00")),
]

# Плотность и слой аттестованной пробы — ИЗ ЗАГОЛОВКА .efr, для печати рядом с
# числами. Пересчитывается при каждом запуске, вручную здесь ничего не стоит.
EFR_KEYS = ("Volume,ml", "Density,g/cm3", "Thick,mm", "Distance,cm")

MATCH_KEV = 1.0   # допуск сопоставления узла кривой и узла сетки


def efr_conditions(name):
    """Условия аттестации из заголовка `.efr` — читаются, а не принимаются."""
    p = paths.efficiency_curve(name)
    if p is None:
        return {}
    out = {}
    for ln in paths.read_text(p).splitlines()[:20]:
        for k in EFR_KEYS:
            if ln.startswith(k + "="):
                out.setdefault(k, ln.split("=", 1)[1].strip())
    return out


def efr_nodes(name):
    """{E: (eps, sigma)} по аттестованной кривой; sigma из столбца dpct."""
    p = paths.efficiency_curve(name)
    if p is None:
        return {}
    out = {}
    for sec in parse_efr(paths.read_text(p)):
        for row in sec["points"]:
            E, eff = float(row[0]), float(row[1])
            dpct = float(row[2]) if len(row) > 2 else 0.0
            if eff > 0:
                out[round(E, 3)] = (eff, eff * dpct / 100.0)
    return out


def grid_nodes(gtag):
    """{E: (eps, sigma)} по нашей сетке моноэнергий; sigma — пуассон по площади.

    Площадь снимается ТОЙ ЖЕ функцией `peakwin.area`, что в пересчёте: своя
    копия правила здесь была бы сверкой копии с копией.
    """
    out = {}
    for p in sorted(glob.glob(os.path.join(pr.BUILD, "grid",
                                           gtag + "_E*.csv"))):
        m = re.search(r"_E(\d+\.\d)\.csv$", p)
        if not m:
            continue
        hist, N = pr.load_hist(p)
        if not N:
            continue
        frac = pr.solid_angle_frac(p, gtag)
        if frac is None:
            continue
        pr.USED.add(p)
        E = float(m.group(1))
        det = {}
        a = peakwin.area(hist, E, detail=det)
        if a <= 0:
            continue
        # sigma площади: gross + вклад полки, приведённый по числу каналов.
        var = det["gross"] + det["side"] * (det["n_peak"] / max(
            det["n_side"], 1)) ** 2
        out[E] = ((a / N) * frac, (math.sqrt(max(var, 1.0)) / N) * frac)
    return out


def pair(a, b):
    """Общие узлы двух наборов: [(E, va, sa, vb, sb)] с допуском MATCH_KEV."""
    out = []
    for Ea, (va, sa) in sorted(a.items()):
        best = min(b, key=lambda x: abs(x - Ea), default=None)
        if best is None or abs(best - Ea) > MATCH_KEV:
            continue
        vb, sb = b[best]
        out.append((Ea, va, sa, vb, sb))
    return out


def ratio_series(num, den):
    """[(E, отношение, sigma)] по общим узлам — отклик на смену геометрии."""
    out = []
    for E, vn, sn, vd, sd in pair(num, den):
        r = vn / vd
        out.append((E, r, r * math.hypot(sn / vn, sd / vd)))
    return out


def slope_lnE(series):
    """Наклон ln(отношение) по ln E со своей sigma — «ход отклика по энергии».

    Взвешенная линейная регрессия. Наклон, а не отношение концов: отношение
    взрывается арифметически, когда одно из чисел проходит через единицу
    (method-rules §7).
    """
    pts = [(math.log(E), math.log(r), s / r) for E, r, s in series
           if E > 0 and r > 0 and s > 0]
    if len(pts) < 3:
        return None, None, len(pts)
    sw = sum(1 / e ** 2 for _x, _y, e in pts)
    sx = sum(x / e ** 2 for x, _y, e in pts)
    sy = sum(y / e ** 2 for _x, y, e in pts)
    sxx = sum(x * x / e ** 2 for x, _y, e in pts)
    sxy = sum(x * y / e ** 2 for x, y, e in pts)
    den = sw * sxx - sx * sx
    if abs(den) < 1e-30:
        return None, None, len(pts)
    k = (sw * sxy - sx * sy) / den
    return k, math.sqrt(sw / den), len(pts)


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if sxx <= 0 or syy <= 0:
        return None
    return sxy / math.sqrt(sxx * syy)


def main():
    lsrm = {}
    ours = {}
    for label, efr, tags in SEQ:
        lsrm[label] = efr_nodes(efr)
        ours[label] = {t: grid_nodes(t) for t in tags}

    missing = [l for l, _e, _t in SEQ if not lsrm[l]]
    if missing:
        print("нет аттестованных кривых: %s" % ", ".join(missing))
    empty = [(l, t) for l, _e, tags in SEQ for t in tags
             if not ours[l].get(t)]
    for l, t in empty:
        print("нет нашей сетки %s (геометрия %s)" % (t, l))

    rows = []

    print("=" * 78)
    print("0. УСЛОВИЯ АТТЕСТАЦИИ ИЗ ЗАГОЛОВКОВ .efr (прочитаны, не приняты)")
    print("=" * 78)
    for label, efr, tags in SEQ:
        c = efr_conditions(efr)
        print("%-14s %-42s -> сопоставимая сетка %s"
              % (label, "; ".join("%s=%s" % (k, c[k]) for k in EFR_KEYS
                                  if k in c) or "—", tags[0]))
    print()

    print("=" * 78)
    print("1. РАЗЛИЧИЕ ГРУПП В КАЖДОЙ ГЕОМЕТРИИ: наша/ЛСРМ по общим узлам")
    print("=" * 78)
    print("%-14s %-10s %5s %10s %10s %s"
          % ("геометрия", "сетка", "узлов", "медиана", "размах", "наклон"
             " ln(наша/ЛСРМ) по ln E"))
    for label, _efr, tags in SEQ:
        for t in tags:
            if not ours[label].get(t) or not lsrm[label]:
                continue
            s = ratio_series(ours[label][t], lsrm[label])
            if not s:
                continue
            vals = sorted(r for _E, r, _sg in s)
            med = vals[len(vals) // 2]
            k, dk, n = slope_lnE(s)
            print("%-14s %-10s %5d %10.4f %10s %s"
                  % (label, t, len(s), med,
                     "%.3f..%.3f" % (vals[0], vals[-1]),
                     "—" if k is None else "%+.4f +- %.4f" % (k, dk)))
            for E, r, sg in s:
                rows.append(("группы", label, t, "%.3f" % E, "%.5f" % r,
                             "%.5f" % sg))
            del n

    print()
    print("=" * 78)
    print("2. ОТКЛИК НА ГЕОМЕТРИЮ ВНУТРИ ГРУППЫ: шаг к предыдущей геометрии")
    print("=" * 78)
    steps = []
    for i in range(1, len(SEQ)):
        cur_l, _e, cur_tags = SEQ[i]
        prv_l, _e2, prv_tags = SEQ[i - 1]
        name = "%s / %s" % (cur_l, prv_l)
        # ЛСРМ: одна кривая на геометрию, вилки нет.
        s_l = (ratio_series(lsrm[cur_l], lsrm[prv_l])
               if lsrm[cur_l] and lsrm[prv_l] else [])
        # Наша: вилка — рабочая плотность против альтернативной.
        s_o = {}
        for t in cur_tags:
            base = prv_tags[0]
            if ours[cur_l].get(t) and ours[prv_l].get(base):
                s_o[t] = ratio_series(ours[cur_l][t], ours[prv_l][base])
        steps.append((name, s_l, s_o))

    print("%-30s %-10s %5s %10s %s"
          % ("шаг", "набор", "узлов", "медиана", "наклон ln R по ln E"))
    for name, s_l, s_o in steps:
        for tag, s in [("ЛСРМ", s_l)] + sorted(s_o.items()):
            if not s:
                continue
            vals = sorted(r for _E, r, _sg in s)
            k, dk, _n = slope_lnE(s)
            print("%-30s %-10s %5d %10.4f %s"
                  % (name, tag, len(s), vals[len(vals) // 2],
                     "—" if k is None else "%+.4f +- %.4f" % (k, dk)))
            for E, r, sg in s:
                rows.append(("отклик", name, tag, "%.3f" % E, "%.5f" % r,
                             "%.5f" % sg))

    print()
    print("=" * 78)
    print("3. КОРРЕЛЯЦИЯ ОТКЛИКОВ: совпадение или АНТИСОВПАДЕНИЕ")
    print("=" * 78)
    print("Признак антисовпадения — ПРОТИВОПОЛОЖНЫЙ ЗНАК наклона отклика по")
    print("энергии, а не малая корреляция: нормировкой знак не исправляется.\n")
    print("%-30s %-10s %11s %11s %8s %s"
          % ("шаг", "наша сетка", "наклон ЛСРМ", "наклон наш", "r", "вердикт"))
    for name, s_l, s_o in steps:
        if not s_l:
            continue
        kl, dkl, _ = slope_lnE(s_l)
        for tag, s in sorted(s_o.items()):
            if not s:
                continue
            ko, dko, _ = slope_lnE(s)
            if kl is None or ko is None:
                continue
            # Корреляция считается по ОБЩИМ узлам обоих откликов.
            dl = {round(E, 1): r for E, r, _s in s_l}
            do = {round(E, 1): r for E, r, _s in s}
            com = sorted(set(dl) & set(do))
            r = pearson([math.log(dl[E]) for E in com],
                        [math.log(do[E]) for E in com]) if len(com) >= 3 \
                else None
            same = (kl > 0) == (ko > 0)
            sig = abs(kl - ko) / math.hypot(dkl, dko) if (dkl and dko) else 0.0
            if not same:
                verdict = "АНТИСОВПАДЕНИЕ знака"
            elif sig > 3.0:
                verdict = "знак один; величина расходится на %.1f sigma" % sig
            else:
                verdict = "совпадает (%.1f sigma)" % sig
            print("%-30s %-10s %+11.4f %+11.4f %8s %s"
                  % (name, tag, kl, ko,
                     "—" if r is None else "%+.3f" % r, verdict))
            rows.append(("корреляция", name, tag, "",
                         "%+.4f" % ko, "%+.4f" % kl))

    if pr.NO_FRAC:
        print("\nДОЛЯ ТЕЛЕСНОГО УГЛА ПРИНЯТА 1,0 (полный 4pi) для сеток: %s."
              % ", ".join(sorted(pr.NO_FRAC)))
        print("  Ни шапки прогона; ни файла-сателлита — прогон старого exe."
              " Для объёмных сеток это верно по построению драйвера,")
        print("  но остаётся ДОПУЩЕНИЕМ о чужом прогоне до их регенерации.")

    verdict, detail = stamp.check_inputs(
        sorted(pr.USED), str(paths.geometry("Gamma-1S")),
        stamp.SRC_LISTS["Gamma-1S"])
    print("\nПРОВЕНАНС входов: %s (%s)" % (verdict, detail))
    if verdict != "ok":
        print("  Числа ПРЕДВАРИТЕЛЬНЫЕ. Отношения внутри одной группы"
              " устойчивее абсолютов, но вердикт не подтверждён.")

    if rows:
        op = os.path.join(str(paths.results("Gamma-1S")),
                          "geometry_response.csv")
        csvio.write(
            op, ["block", "case", "series", "E_keV", "value", "sigma"], rows,
            comments=[
                "Отклик двух групп кривых на смену геометрии в"
                " последовательности точка -> Петри -> Дента -> Маринелли.",
                "block=группы: наша/ЛСРМ по общим узлам в каждой геометрии.",
                "block=отклик: отношение кривой к ПРЕДЫДУЩЕЙ геометрии внутри"
                " своей группы (форма отклика;",
                "  всё; что от геометрии не зависит; сокращается).",
                "block=корреляция: value — наклон нашего отклика по ln E;"
                " sigma — наклон отклика ЛСРМ.",
                "Обе группы — эффективность НА ИСПУЩЕННЫЙ КВАНТ; выход линии"
                " не участвует.",
                "Плотность объёмных проб в аттестации против нашей сетки —"
                " ДОПУЩЕНИЕ; отклик дан вилкой 0;60/1;60.",
            ],
            stamp=stamp.lines(
                "detectors/Gamma-1S/analysis/geometry_response.py",
                dict(peakwin.declare(),
                     quantity="эффективность на испущенный квант; отношения"
                              " между геометриями и между группами",
                     area="пик полного поглощения в депозит-спектре (наша)"
                          " против узлов аттестованной кривой .efr (ЛСРМ)",
                     cone="конусные сетки приведены делением на долю"
                          " телесного угла из шапки спектра",
                     density_assumption="плотность объёмных проб в аттестации"
                                        " не совпадает с сеткой; дана вилка"
                                        " 0;60 и 1;60 г/см3",
                     sequence="точка 5 см -> Петри 60 мл -> Дента 120 мл ->"
                              " Маринелли 1 л"),
                inputs=sorted(pr.USED),
                geometry_dir=str(paths.geometry("Gamma-1S")),
                names=stamp.SRC_LISTS["Gamma-1S"],
                repo_dir=str(paths.REPO)))
        print("таблица: %s" % op)


if __name__ == "__main__":
    main()
