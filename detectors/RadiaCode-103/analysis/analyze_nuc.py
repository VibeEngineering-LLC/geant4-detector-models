# -*- coding: utf-8 -*-
"""Разделение гамма- и бета-вклада по прогонам полного распада.

Идея: полный прогон нуклида даёт ВЕСЬ отклик (гамма, бета, конверсионные
электроны, тормозное). Отдельно предсказываем чисто гамма-часть, свернув
известные линии с посчитанными кривыми eps_t(E) и eps_p(E). Разность — вклад
заряженных частиц. Заодно совпадение по площадям пиков проверяет сами кривые:
это два независимых расчёта одной величины.
"""
import csv
import os

import numpy as np

import nucdata
import sys
# Модули прибора лежат в двух каталогах: разбор в analysis/, запуск прогонов
# в drivers/. Импорт через каталог-сосед иначе не находится: python кладёт
# в sys.path только каталог запускаемого файла.
for _d in ("analysis", "drivers"):
    _p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import rcspec

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.abspath(os.path.join(HERE, "..", "results"))
THRESH = 20


def load_curves(cfg):
    """eps_p(E), eps_t(E) для нужной конфигурации из efficiency.csv."""
    E, ep, et = [], [], []
    for r in csv.DictReader(open(rcspec.rdir("efficiency.csv"),
                                 encoding="utf-8")):
        if r["config"] == cfg:
            E.append(float(r["E_keV"]))
            ep.append(float(r["eps_p"]))
            et.append(float(r["eps_t"]))
    o = np.argsort(E)
    return np.array(E)[o], np.array(ep)[o], np.array(et)[o]


def interp_log(x, xs, ys):
    """Интерполяция кривой эффективности в логарифмических координатах."""
    return np.exp(np.interp(np.log(np.clip(x, xs[0], xs[-1])),
                            np.log(xs), np.log(ys)))


def main():
    base = rcspec.rdir("nuclides")
    for cfgdir in sorted(os.listdir(base)):
        d = os.path.join(base, cfgdir)
        if not os.path.isdir(d):
            continue
        cfg = "full_" + cfgdir.replace("_", "_") if not cfgdir.startswith("full_") \
            else cfgdir
        cfg = "full_" + cfgdir if not cfgdir.startswith("full_") else cfgdir
        try:
            Eg, epg, etg = load_curves(cfg)
        except Exception:
            print("нет кривых для", cfg)
            continue

        print("\n=== %s ===" % cfgdir)
        print("%-22s %11s %11s %11s %8s" %
              ("нуклид", "всего/расп", "гамма/расп", "бета/расп", "доля β"))
        for nuc in ("K40", "Cs137", "Ra226", "Th232", "U238"):
            p = os.path.join(d, "nuc_%s.csv" % nuc)
            if not os.path.exists(p):
                continue
            meta, hist = rcspec.read_spec(p)
            n = float(meta["N_primaries"])
            total = hist[THRESH:].sum() / n

            lines = nucdata.LINES[nuc]
            gam = sum(y * interp_log(e, Eg, etg) for e, y in lines)
            beta = total - gam
            print("%-22s %11.4e %11.4e %11.4e %7.0f %%"
                  % (nucdata.TITLE[nuc], total, gam, beta,
                     100 * beta / total if total else 0))

        # проверка кривых по площадям пиков
        print("\n  проверка кривых: площадь пика, расчёт по линии против прогона")
        print("  %-9s %-8s %11s %11s %8s" %
              ("нуклид", "линия", "по кривой", "по прогону", "разн."))
        for nuc in ("K40", "Cs137", "Ra226", "Th232"):
            p = os.path.join(d, "nuc_%s.csv" % nuc)
            if not os.path.exists(p):
                continue
            meta, hist = rcspec.read_spec(p)
            n = float(meta["N_primaries"])
            for E0, y in sorted(nucdata.LINES[nuc], key=lambda t: -t[1])[:2]:
                if E0 < 100:
                    continue
                pred = y * interp_log(E0, Eg, epg)
                lo, hi = int(E0 - 4), int(E0 + 4) + 1
                obs = hist[lo:hi].sum() / n
                print("  %-9s %-8.1f %11.4e %11.4e %+7.0f %%"
                      % (nuc, E0, pred, obs,
                         100 * (obs / pred - 1) if pred > 0 else 0))


if __name__ == "__main__":
    main()
