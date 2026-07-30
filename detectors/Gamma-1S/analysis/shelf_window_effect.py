# -*- coding: utf-8 -*-
"""Две правки окна площади, измеренные ПО ОТДЕЛЬНОСТИ (задачи 129, 132).

ЗАЧЕМ ПО ОТДЕЛЬНОСТИ. В `area_sim` изменились две независимые вещи:
  (а) окно полки `[E−30, E−10]` -> `[E−25, E−10]` — чтобы плоская подложка не
      вычиталась поверх пика вылета K-рентгена иода (`E−28,6`);
  (б) ширина окна пика: множитель `2*win+1` = 13 каналов -> `nchan()` = 12,
      потому что каналы стоят в серединах `(i+0.5)*bin` и правило,
      записанное в кэВ, применялось к каналам.
Считать их вместе нельзя: тогда неизвестно, что именно сдвинуло кривую, а
контрольный признак починки обязан отличаться от починяемого правила ровно на
одну вещь (method-rules §5). Обе правки трогают ТОЛЬКО вычитаемую полку,
поэтому обе поднимают `eps` и обе сильнее на мягком крае.

ЧТО СЧИТАЕТСЯ. Одна и та же сетка моноэнергий, один и тот же файл, четыре
правила съёма площади. Прогонов не делается вовсе — это чистая сверка
конвенций, самая дешёвая из возможных: различие не может быть отнесено ни к
геометрии, ни к статистике, ни к версии exe.

ОГОВОРКА О ВХОДАХ. Штамп провенанса печатается вместе с числами. Если вердикт
не `ok`, отношения ниже всё равно осмысленны (числитель и знаменатель считаны
с ОДНОГО файла), а вот абсолютные `eps` — нет.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
import paths  # noqa: E402
import stamp  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import point_recalc as pr  # noqa: E402

WIN = 6.0
BG1 = 10.0

# Четыре правила: (метка, край полки, считать ли каналы по сетке).
RULES = (("старое", 30.0, False),
         ("окно25", 25.0, False),
         ("каналы", 30.0, True),
         ("оба", 25.0, True))


def area(hist, E, bg0, exact):
    gross = sum(c for e, c in hist.items() if abs(e - E) <= WIN)
    side = sum(c for e, c in hist.items() if E - bg0 <= e <= E - BG1)
    if exact:
        n_peak = pr.nchan(E - WIN, E + WIN)
        n_side = pr.nchan(E - bg0, E - BG1)
    else:
        n_peak = 2 * WIN + 1
        n_side = bg0 - BG1
    return gross - side / n_side * n_peak


def main():
    rows = []
    for gtag in ("p5cm", "p25cm"):
        saf = os.path.join(pr.BUILD, "grid", "%s_solidangle.txt" % gtag)
        if not os.path.exists(saf):
            print("нет сетки %s — пропуск" % gtag)
            continue
        frac = float(open(saf).read().strip())
        import glob
        import re
        for p in sorted(glob.glob(os.path.join(pr.BUILD, "grid",
                                               gtag + "_E*.csv"))):
            m = re.search(r"_E(\d+\.\d)\.csv$", p)
            if not m:
                continue
            E = float(m.group(1))
            hist, N = pr.load_hist(p)
            if not N:
                continue
            pr.USED.add(p)
            eps = {}
            for label, bg0, exact in RULES:
                a = area(hist, E, bg0, exact)
                eps[label] = (a / N) * frac if a > 0 else None
            if None in eps.values():
                continue
            rows.append((gtag, E, eps))

    if not rows:
        raise SystemExit("Нет сеток моноэнергий — нечего сверять.")

    print("Влияние двух правок окна на eps сетки, ПО ОТДЕЛЬНОСТИ.")
    print("окно25 — только край полки 30 -> 25 кэВ;"
          " каналы — только 13 -> 12 каналов пика; оба — вместе.\n")
    print("%-7s %9s %12s %9s %9s %9s"
          % ("сетка", "E, кэВ", "eps старое", "окно25", "каналы", "оба"))
    out = []
    for gtag, E, eps in rows:
        d = {k: 100.0 * (eps[k] / eps["старое"] - 1.0) for k in
             ("окно25", "каналы", "оба")}
        print("%-7s %9.1f %12.5e %+8.2f%% %+8.2f%% %+8.2f%%"
              % (gtag, E, eps["старое"], d["окно25"], d["каналы"], d["оба"]))
        out.append((gtag, "%.1f" % E, "%.6e" % eps["старое"],
                    "%.6e" % eps["оба"], "%+.3f" % d["окно25"],
                    "%+.3f" % d["каналы"], "%+.3f" % d["оба"]))

    # ПЕРЕКРЁСТНЫЙ ЧЛЕН ВЫВЕДЕН, А НЕ УГАДАН. Первая версия проверки ждала
    # произведения поправок `a*b` — так было бы, если правки действовали
    # МУЛЬТИПЛИКАТИВНО. Здесь они обе уменьшают ОДНО И ТО ЖЕ вычитаемое
    # слагаемое, и точный ответ другой. Пусть
    #   net = G − S/n_s · n_p,  X = S30/20·13 (прежнее вычитание),
    #   u = S25/15·13 (только окно),  v = S30/20·12 (только каналы).
    # Тогда a = (X−u)/net, b = (X−v)/net, совместное = (X − u·12/13)/net, и
    #   a + b − совместное = (X − v − u/13)/net = (X/13 − u/13)/net = a/13,
    # поскольку X − v = S30/20 = X/13. То есть ожидаемый перекрёстный член
    # равен a·(1 − 12/13) и на 10,9 % даёт 0,84 п.п. — ровно наблюдаемое.
    # Это тот самый случай «репродукция ≠ верификация»: числа были верны,
    # неверно было моё ожидание, и поймал это только вывод формулы.
    N_OLD, N_NEW = 2 * WIN + 1, 12.0
    print("\nПерекрёстный член (сумма поодиночке минус совместное):")
    worst = (0.0, 0.0, None, 0.0)
    for gtag, E, eps in rows:
        a = 100.0 * (eps["окно25"] / eps["старое"] - 1.0)
        b = 100.0 * (eps["каналы"] / eps["старое"] - 1.0)
        both = 100.0 * (eps["оба"] / eps["старое"] - 1.0)
        d = abs(a + b - both)
        if d > worst[0]:
            worst = (d, E, gtag, a * (1.0 - N_NEW / N_OLD))
    print("  наибольший %.3f п.п. на %s %.1f кэВ; вывод даёт a*(1-12/13) ="
          " %.3f п.п." % (worst[0], worst[2], worst[1], worst[3]))
    ok = abs(worst[0] - worst[3]) < 0.02 + 0.02 * abs(worst[3])
    print("  сошлось -> правки трогают только полку, лишнего нет." if ok else
          "  !! НЕ сошлось — правки задевают не только полку, разобрать.")

    verdict, detail = stamp.check_inputs(
        sorted(pr.USED), str(paths.geometry("Gamma-1S")),
        stamp.SRC_LISTS["Gamma-1S"])
    print("\nПРОВЕНАНС входов: %s (%s)" % (verdict, detail))
    if verdict != "ok":
        print("  Отношения выше осмысленны (числитель и знаменатель с одного"
              " файла); абсолютные eps — нет.")

    csvio.write(
        os.path.join(str(paths.results("Gamma-1S")), "shelf_window_effect.csv"),
        ["grid", "E_keV", "eps_old", "eps_new", "d_window25_pct",
         "d_channels_pct", "d_both_pct"],
        out,
        comments=[
            "Влияние двух правок окна съёма площади на eps сеток моноэнергий;"
            " прогоны НЕ переделывались.",
            "d_window25_pct — только край полки 30 -> 25 кэВ (пик вылета"
            " K-рентгена иода уходит из окна полки).",
            "d_channels_pct — только ширина окна пика: 2*win+1 = 13 каналов"
            " -> nchan() = 12 (сетка смещена на полканала).",
            "Обе правки только УМЕНЬШАЮТ вычитаемую полку; отсюда знак плюс и"
            " рост эффекта к мягкому краю.",
        ],
        stamp=stamp.lines(
            "detectors/Gamma-1S/analysis/shelf_window_effect.py",
            {"quantity": "eps сетки моноэнергий при четырёх правилах съёма"
                         " площади; отношения к прежнему правилу",
             "area": "пик полного поглощения в депозит-спектре",
             "window": "±6;0 кэВ",
             "shelf": "сравниваются [E-30; E-10] и [E-25; E-10]"
                      " односторонние",
             "blurred": "нет (депозит как есть)",
             "cone": "приведено делением на долю телесного угла из"
                     " <сетка>_solidangle.txt"},
            inputs=sorted(pr.USED),
            geometry_dir=str(paths.geometry("Gamma-1S")),
            names=stamp.SRC_LISTS["Gamma-1S"], repo_dir=str(paths.REPO)))
    print("таблица: %s"
          % os.path.join(str(paths.results("Gamma-1S")),
                         "shelf_window_effect.csv"))


if __name__ == "__main__":
    main()
