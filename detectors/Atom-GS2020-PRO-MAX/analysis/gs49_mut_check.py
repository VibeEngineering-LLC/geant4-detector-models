# -*- coding: utf-8 -*-
"""#GS-49 мутационная приёмка DECAY_TALLY (#SA-3): по tally-файлу печатает число сигнатур, долю событий с фотоном 661,657 кэВ
(Cs-137 → Ba-137m; справочно ≈0,851) и долю событий с электронной частью сигнатуры E:. Запуск: python gs49_mut_check.py tally.txt [...]"""
import sys

sys.stdout.reconfigure(encoding="utf-8")
for path in sys.argv[1:]:
    n_ev, tot, g661, withE, nsig = 0, 0, 0, 0, 0
    for ln in open(path, encoding="utf-8"):
        if ln.startswith("n_events,"):
            n_ev = int(ln.split(",")[1])
        if ";" not in ln or ln.startswith("count"):
            continue
        cnt, *parts = ln.strip().split(";")
        cnt = int(cnt)
        nsig += 1
        tot += cnt
        ph = [float(x) for x in next((p[2:].replace(",", " ").split() for p in parts if p.startswith("P:")), []) if x]
        if any(abs(e - 661.657) < 0.5 for e in ph):
            g661 += cnt
        if any(p.startswith("E:") and len(p) > 2 for p in parts):
            withE += cnt
    print("%s: событий %d, сигнатур %d, Σсчёт %d, доля γ661,657 %.4f, доля с E: %.4f" % (path, n_ev, nsig, tot, g661 / n_ev, withE / n_ev))
