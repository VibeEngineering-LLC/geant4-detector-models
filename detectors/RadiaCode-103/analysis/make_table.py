# -*- coding: utf-8 -*-
"""Рабочая таблица эффективности по линиям ЕРН и техногенных нуклидов."""
import csv
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
import rcspec
RESULTS = rcspec.RESULTS

ROWS = [(238.6, "Pb-212", "Th-232"), (295.2, "Pb-214", "Ra-226"),
        (351.9, "Pb-214", "Ra-226"), (583.2, "Tl-208", "Th-232"),
        (609.3, "Bi-214", "Ra-226"), (661.7, "Cs-137", "Cs-137"),
        (911.2, "Ac-228", "Th-232"), (968.9, "Ac-228", "Th-232"),
        (1120.3, "Bi-214", "Ra-226"), (1460.8, "K-40", "K-40"),
        (1764.5, "Bi-214", "Ra-226"), (2614.5, "Tl-208", "Th-232")]

COLS = [("full_air_0.00", "предел"), ("full_organic_0.50", "0,50"),
        ("full_soil_0.80", "0,80"), ("full_water_1.00", "1,00"),
        ("full_soil_1.20", "1,20"), ("full_soil_1.60", "1,60")]


def main():
    data = {}
    for r in csv.DictReader(open(rcspec.rdir("efficiency.csv"),
                                 encoding="utf-8")):
        data.setdefault(r["config"], {})[round(float(r["E_keV"]), 1)] = (
            float(r["eps_p"]), float(r["d_eps_p"]))

    out = rcspec.rdir("efficiency_table.md")
    with open(out, "w", encoding="utf-8") as f:
        f.write("| Линия, кэВ | Нуклид | Определяет | "
                + " | ".join("ρ = %s" % c[1] for c in COLS) + " |\n")
        f.write("|---|---|---|" + "---|" * len(COLS) + "\n")
        for E, nuc, det in ROWS:
            cells = []
            for cfg, _ in COLS:
                v = data.get(cfg, {}).get(E)
                cells.append("%.2e" % v[0] if v else "—")
            f.write("| %.1f | %s | %s | %s |\n" % (E, nuc, det, " | ".join(cells)))
    print(open(out, encoding="utf-8").read())
    print("таблица:", out)


if __name__ == "__main__":
    main()
