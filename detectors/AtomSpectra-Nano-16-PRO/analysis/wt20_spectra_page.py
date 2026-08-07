# -*- coding: utf-8 -*-
"""Страница проверки калибровки: образец и фон на ОБЩЕЙ энергетической сетке.

Директива оператора 07.08.2026: «выведи откалиброванные спектры фона и
образца с наложением и возможностью вычитания для проверки». Оба спектра
приводятся к сетке 1 кэВ ЧЕРЕЗ СВОИ поправки калибровки (§3 методики), фон
масштабируется к времени образца; страница рисует их с наложением,
разностью, зумом и маркерами опорных линий.

    python analysis/wt20_spectra_page.py <спектр.xml> <каталог с calibration_fitted.csv> [вывод.html]
"""
import csv
import io
import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from wt20_unfold import (read_atomspectra_xml, read_correction,   # noqa: E402
                         rebin_to_grid)

E_TOP = 3380.0   # ниже канала переполнения: последний канал шкалы (~3394)
                 # копит оверфлоу (974 отсчёта одним выбросом) и пиком не является
STEP = 1.0


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    src, caldir = sys.argv[1], sys.argv[2]
    out_html = (sys.argv[3] if len(sys.argv) > 3
                else os.path.join(caldir, "wt20_spectra.html"))

    spec = read_atomspectra_xml(src)
    if isinstance(spec, (list, tuple)):
        spec = spec[0]
    bg = getattr(spec, "background_embedded", None)
    if bg is None:
        raise SystemExit("во входном XML нет встроенного фона")
    t_smp, t_bg = float(spec.real_time), float(bg.real_time)

    corr = read_correction(os.path.join(caldir, "calibration_fitted.csv"))
    if not corr:
        raise SystemExit("нет calibration_fitted.csv в %s — сперва прогнать "
                         "wt20_calibration.py" % caldir)

    edges = np.arange(0.0, E_TOP + STEP, STEP)
    y_smp = rebin_to_grid(np.asarray(spec.counts, float),
                          list(spec.energy_cal), corr.get("sample"), edges)
    y_bg = rebin_to_grid(np.asarray(bg.counts, float),
                         list(bg.energy_cal), corr.get("background"), edges)
    y_bgs = y_bg * (t_smp / t_bg)

    # якоря и невязки — из calibration_check.csv той же ревизии, если есть
    anchors = {"sample": [], "background": []}
    pc = os.path.join(caldir, "calibration_check.csv")
    if os.path.exists(pc):
        for row in csv.reader(io.open(pc, encoding="utf-8")):
            if not row or row[0].startswith(("спектр", "#")):
                continue
            tgt = ("sample" if row[0].startswith("ОБРАЗЕЦ")
                   else "background")
            anchors[tgt].append({"nuc": row[1], "E": float(row[2]),
                                 "dE": float(row[6])})

    # зона K-серии вольфрама — окно XW из каталога конструктора: рисуется
    # полосой, чтобы было видно, что локального максимума там нет
    import roi_lines as R
    xw = next((r for r in R.parse_xml(R.DEFAULT_XML)
               if r["key"] == "XW" and r.get("window")), None)

    # сумм-пики каскада Tl-208 — энергии из библиотеки МАГАТЭ
    import wt20_calibration as C
    e26, _ = C.lib_energy("208tl", 2614.51)
    sums = [{"E": round(e26 + C.lib_energy("208tl", q)[0], 2),
             "lab": "Σ 2615+%d" % round(q)}
            for q in (277.37, 510.77, 583.19)]

    data = {
        "step": STEP,
        "t_smp": t_smp, "t_bg": t_bg,
        "xw_zone": list(xw["window"]) if xw else None,
        "sums": sums,
        "smp": [round(float(v), 3) for v in y_smp],
        "bg": [round(float(v), 3) for v in y_bgs],
        "anchors": anchors,
        "corr": {k: list(map(float, v)) for k, v in corr.items()},
        "n_smp": int(y_smp.sum()), "n_bg": int(y_bgs.sum()),
    }

    tpl = io.open(os.path.join(_HERE, "wt20_spectra_template.html"),
                  encoding="utf-8").read()
    html = tpl.replace("@DATA@", json.dumps(data, ensure_ascii=False))
    with io.open(out_html, "w", encoding="utf-8", newline="\n") as f:
        f.write(html)
    print("записано: %s (%d КБ)" % (out_html, os.path.getsize(out_html)
                                    // 1024))
    return 0


if __name__ == "__main__":
    sys.exit(main())
