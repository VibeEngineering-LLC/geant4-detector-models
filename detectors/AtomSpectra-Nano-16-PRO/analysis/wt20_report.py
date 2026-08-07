# -*- coding: utf-8 -*-
"""Сборка интерактивного отчёта по замеру WT-20 в одну самодостаточную страницу.

Страница НЕ содержит чисел, набранных руками: всё, что в ней показано, читается
из файлов расчёта (`calibration_check.csv`, `unfold_activities.csv`,
`line_activities.csv`, `unfold_spectrum.csv`) и шапок спектров Geant4. Правка
результата — это перезапуск расчёта и пересборка страницы, а не правка текста.

    python analysis/wt20_report.py <каталог расчёта> <каталог шаблонов> [выход.html]
"""
import csv
import io
import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))


def read_csv(path):
    """Читает CSV, сам определяя разделитель.

    В дереве уживаются два: скрипты разложения пишут запятой, скрипты прямой
    задачи — точкой с запятой (в них десятичная запятая недопустима, поэтому
    поля разделены иначе). Ошибка в разделителе не роняет чтение, а тихо
    складывает всю строку в одну ячейку — таблица на странице выходила
    одноколоночной. Поэтому разделитель определяется по первой значащей строке,
    а не задаётся по умолчанию.
    """
    if not os.path.exists(path):
        return []
    with io.open(path, encoding="utf-8") as f:
        raw = [ln for ln in f if ln.strip() and not ln.startswith("#")]
    if not raw:
        return []
    delim = ";" if raw[0].count(";") > raw[0].count(",") else ","
    return [r for r in csv.reader(raw, delimiter=delim) if r]


def read_head(path):
    head = {}
    if not os.path.exists(path):
        return head
    for ln in io.open(path, encoding="utf-8"):
        if not ln.startswith("#"):
            break
        if "=" in ln:
            k, v = ln.lstrip("# ").split("=", 1)
            head[k.strip()] = v.strip()
    return head


def load_spectrum(path, step=1):
    """Спектр разложения. Без прореживания — та же сетка, что у калибровочной
    страницы (2 кэВ у ребиннинга разложения); при step>1 бины суммируются.
    Оператор просил на графике убрать сглаживание, поэтому step по умолчанию 1.
    """
    rows = read_csv(path)
    if not rows:
        return {}, []
    names = rows[0]
    data = [[float(v) for v in r] for r in rows[1:]]
    if step <= 1:
        return names, [[round(v, 4) for v in r] for r in data]
    out = []
    for i in range(0, len(data) - step + 1, step):
        block = data[i:i + step]
        agg = [sum(r[j] for r in block) for j in range(len(names))]
        agg[0] = sum(r[0] for r in block) / len(block)
        out.append([round(v, 4) for v in agg])
    return names, out


def ru(x, nd=1):
    return ("%.*f" % (nd, x)).replace(".", ",")


def data_uri(path):
    """Картинка -> data:URI. Страница обязана быть самодостаточной: внешние
    адреса на ней не грузятся, а file:// не переживёт передачи."""
    import base64
    if not os.path.exists(path):
        return ""
    ext = os.path.splitext(path)[1].lower()
    mime = {".png": "image/png", ".jpg": "image/jpeg",
            ".jpeg": "image/jpeg"}.get(ext, "application/octet-stream")
    with open(path, "rb") as f:
        return "data:%s;base64,%s" % (mime,
                                      base64.b64encode(f.read()).decode())


def load_spectra_pair(d, src_xml=None):
    """Откалиброванные спектры образца и фона для интерактивного графика.

    Читаются напрямую из XML замера через тот же читатель, что использует
    разложение и калибровка (`wt20_unfold`). Поправки берутся из
    `calibration_fitted.csv` того же каталога, что и остальные выходы;
    отсутствие XML или поправок — не отказ, а пустой блок: страница
    собирается без интерактивных спектров.
    """
    if not src_xml or not os.path.exists(src_xml):
        return None
    try:
        sys.path.insert(0, _HERE)
        import numpy as np
        from wt20_unfold import (read_atomspectra_xml, read_correction,
                                 rebin_to_grid)
        import wt20_calibration as C
        import roi_lines as R
    except Exception as e:
        print("не удалось загрузить модули для спектров:", e)
        return None
    spec = read_atomspectra_xml(src_xml)
    if isinstance(spec, (list, tuple)):
        spec = spec[0]
    bg = getattr(spec, "background_embedded", None)
    if bg is None:
        return None
    t_smp, t_bg = float(spec.real_time), float(bg.real_time)
    corr = read_correction(os.path.join(d, "calibration_fitted.csv"))
    if not corr:
        return None
    E_TOP, STEP = 3380.0, 1.0
    edges = np.arange(0.0, E_TOP + STEP, STEP)
    y_smp = rebin_to_grid(np.asarray(spec.counts, float),
                          list(spec.energy_cal), corr.get("sample"), edges)
    y_bg = rebin_to_grid(np.asarray(bg.counts, float),
                         list(bg.energy_cal), corr.get("background"), edges)
    y_bgs = y_bg * (t_smp / t_bg)
    # Якоря и невязки — из calibration_check.csv (уже прочитан выше в main как
    # cal), здесь только формируем список для страницы. Формат тот же, что в
    # wt20_spectra_page: {nuc, E, dE}, разложенный по «sample»/«background».
    anchors = {"sample": [], "background": []}
    pc = os.path.join(d, "calibration_check.csv")
    if os.path.exists(pc):
        for row in csv.reader(io.open(pc, encoding="utf-8")):
            if not row or row[0].startswith(("спектр", "#")):
                continue
            tgt = ("sample" if row[0].startswith("ОБРАЗЕЦ")
                   else "background")
            anchors[tgt].append({"nuc": row[1], "E": float(row[2]),
                                 "dE": float(row[6])})
    xw = next((r for r in R.parse_xml(R.DEFAULT_XML)
               if r["key"] == "XW" and r.get("window")), None)
    e26, _ = C.lib_energy("208tl", 2614.51)
    sums = [{"E": round(e26 + C.lib_energy("208tl", q)[0], 2),
             "lab": "Σ 2615+%d" % round(q)}
            for q in (277.37, 510.77, 583.19)]
    return dict(step=STEP, t_smp=t_smp, t_bg=t_bg,
                xw_zone=(list(xw["window"]) if xw else None),
                sums=sums,
                smp=[round(float(v), 3) for v in y_smp],
                bg=[round(float(v), 3) for v in y_bgs],
                anchors=anchors,
                corr={k: list(map(float, v)) for k, v in corr.items()},
                n_smp=int(y_smp.sum()), n_bg=int(y_bgs.sum()))


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    d, tdir = sys.argv[1], sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(d, "wt20.html")
    # Опциональный XML замера — для интерактивных спектров образца и фона.
    # Позиционный аргумент через -x <путь>, чтобы не ломать существующий вызов.
    src_xml = None
    for i, a in enumerate(sys.argv):
        if a == "-x" and i + 1 < len(sys.argv):
            src_xml = sys.argv[i + 1]

    names, spec = load_spectrum(os.path.join(d, "unfold_spectrum.csv"))
    # Фон вычитается ЗДЕСЬ, до всякого показа: разложение считается по разности
    # «измеренное минус фон», и страница должна показывать ту же величину, что
    # раскладывается. Отдельной компонентой в стопке фон больше не стоит —
    # иначе на рисунке видно одно, а подгоняется другое.
    bg_sum = 0.0
    if names and "фон" in names:
        jb = names.index("фон")
        bg_sum = sum(r[jb] for r in spec)
        for r in spec:
            r[1] -= r[jb]          # измеренное минус фон
            r[2] -= r[jb]          # модель тоже без фона
            del r[jb]
        names = [n for n in names if n != "фон"]
        names[1] = "измерено минус фон"
        names[2] = "сумма компонент"
    meas_sum = sum(r[1] for r in spec) + bg_sum
    acts = read_csv(os.path.join(d, "unfold_activities.csv"))
    lines = read_csv(os.path.join(d, "line_activities.csv"))
    cal = read_csv(os.path.join(d, "calibration_check.csv"))
    head = read_head(os.path.join(tdir, "Tl208.csv"))

    ref = os.path.normpath(os.path.join(_HERE, "..", "reference"))
    draw = os.path.normpath(os.path.join(_HERE, "..", "drawings"))
    # Прямая задача: источник по составу. Файлы не обязательны — страница
    # собирается и без них, просто без соответствующего раздела.
    def opt(name):
        p = os.path.join(d, name)
        return read_csv(p) if os.path.exists(p) else []

    # Шапка прямой задачи: удельная активность и масса лежат в строках с «#»,
    # которые read_csv пропускает как комментарий. Выводятся отдельно —
    # без них раздел показывает выходы линий, но не саму удельную активность.
    src_head = {}
    fp = os.path.join(d, "wt20_source_forward.csv")
    if os.path.exists(fp):
        for ln in io.open(fp, encoding="utf-8"):
            if not ln.startswith("#") or ";" not in ln:
                continue
            k, v = ln.lstrip("# ").rsplit(";", 1)
            src_head[k.strip()] = v.strip()

    spectra = load_spectra_pair(d, src_xml)
    # Разложение методом 2 (матрица отклика × выходы линий) — если есть,
    # кладётся отдельным блоком payload. Формат тот же, что у метода 1:
    # столбец 0 — E, столбец 1 — измерено, 2 — модель, дальше компоненты.
    names_m2, spec_m2 = load_spectrum(os.path.join(d,
                                                   "unfold_matrix_spectrum.csv"))
    acts_m2 = read_csv(os.path.join(d, "unfold_matrix_activities.csv"))
    if spec_m2 and "фон" in names_m2:
        jb = names_m2.index("фон")
        for r in spec_m2:
            r[1] -= r[jb]
            r[2] -= r[jb]
            del r[jb]
        names_m2 = [n for n in names_m2 if n != "фон"]
        names_m2[1] = "измерено минус фон"
        names_m2[2] = "сумма компонент"
    payload = dict(names=names, spec=spec, acts=acts, lines=lines, cal=cal,
                   spectra=spectra,
                   names_m2=names_m2, spec_m2=spec_m2, acts_m2=acts_m2,
                   src_head=src_head,
                   head=head,
                   bg_sum=round(bg_sum, 1),
                   bg_share=round(100.0 * bg_sum / meas_sum, 2) if meas_sum else 0,
                   src_fwd=opt("wt20_source_forward.csv"),
                   peakchk=opt("peak_check.csv"),
                   widthfit=opt("width_fit.csv"),
                   selfabs=opt("wt20_selfabsorption.csv"),
                   scatter=opt("wt20_source_scatter.csv"),
                   surflim=opt("wt20_surface_limit.csv"),
                   winscan=opt("wt20_window_scan.csv"),
                   # Показывается ПРИНЯТЫЙ замер — по красной маркировке марки,
                   # три независимые сканлинии с профилем красноты. Прежде на
                   # странице стояла контрольная картинка отменённого замера по
                   # яркости (одна сканлиния, 6,15 мм), из которой принятый
                   # результат 4,85 мм не следует.
                   img_scan=data_uri(os.path.join(
                       draw, "wt20_tips_scanlines.png")),
                   img_setup=data_uri(os.path.join(draw,
                                                   "nano16pro_wt20_setup.png")))
    tpl = io.open(os.path.join(_HERE, "wt20_report_template.html"),
                  encoding="utf-8").read()
    html = tpl.replace("/*__DATA__*/", json.dumps(payload, ensure_ascii=False))
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    io.open(out, "w", encoding="utf-8").write(html)
    print("записано: %s (%d КБ)" % (out, len(html) // 1024))
    return 0


if __name__ == "__main__":
    sys.exit(main())
