# -*- coding: utf-8 -*-
"""Сборка интерактивного отчёта по замеру WT-20 в одну самодостаточную страницу.

Страница НЕ содержит чисел, набранных руками: всё, что в ней показано, читается
из файлов расчёта (`calibration_check.csv`, `unfold_activities.csv`,
`line_activities.csv`, `unfold_spectrum.csv`) и шапок спектров Geant4. Правка
результата — это перезапуск расчёта и пересборка страницы, а не правка текста.

    python analysis/wt20_report.py <каталог расчёта> <каталог шаблонов>
                                   [-o выход.html] [-x спектр.xml]

Флаг -x — опциональный XML замера для интерактивной панели спектров. Позиционные
и опциональные аргументы разбираются argparse (задача №65: раньше третий
позиционный без строгого парсинга принял путь к XML как out.html и перезаписал
исходный файл замера).
"""
import argparse
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

    Комментарий вида `# текст; с разделителем внутри` `csv.writer` заключает
    в кавычки целиком (иначе разделитель внутри поля сломал бы структуру), и
    СЫРАЯ строка тогда начинается с `"`, а не с `#` — фильтр по первому символу
    такую строку не ловил, она попадала в таблицу отдельной строкой-мусором
    (найдено на живой странице: Таблица 3 «Активности звеньев ряда» схлопнулась
    в один столбец, потому что первой строкой шёл однопольный комментарий).
    Проверка снимает внешнюю кавычку перед тем, как решать, комментарий это.
    """
    if not os.path.exists(path):
        return []
    with io.open(path, encoding="utf-8") as f:
        raw = []
        for ln in f:
            if not ln.strip():
                continue
            probe = ln.lstrip()
            if probe.startswith('"'):
                probe = probe[1:]
            if probe.startswith("#"):
                continue
            raw.append(ln)
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
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("d", help="каталог расчёта (unfold_*.csv, calibration_*)")
    ap.add_argument("tdir", help="каталог шаблонов Geant4 (для чтения шапок)")
    ap.add_argument("-o", "--out",
                    help="выход HTML (по умолчанию <d>/wt20.html)")
    ap.add_argument("-x", "--xml",
                    help="XML замера AtomSpectra для панели спектров")
    args = ap.parse_args()
    d, tdir = args.d, args.tdir
    out = args.out or os.path.join(d, "wt20.html")
    src_xml = args.xml

    # Защита №65: -x должен показывать на .xml, out — на .html.
    # Иначе argparse мог принять xml-путь как позиционный тремя разными
    # способами (через -o=..., через мисклики CLI), и запись HTML в этот путь
    # разрушала бы исходный XML.
    if not out.lower().endswith(".html"):
        raise SystemExit("выходной путь %r не .html — отказ; --out обязателен "
                         "и должен указывать HTML-файл" % out)
    if src_xml and not src_xml.lower().endswith(".xml"):
        raise SystemExit("XML замера %r не .xml — отказ" % src_xml)
    if src_xml and os.path.abspath(src_xml) == os.path.abspath(out):
        raise SystemExit("--out и --xml указывают на один файл — отказ")

    # Рисунки наложения строятся по ПРЯМОЙ задаче: расчёт по составу поверх
    # измеренного (analysis/wt20_forward.py). Формат forward_components.csv
    # совпадает с тем, что ждёт фабрика графика: столбец 0 — энергия,
    # 1 — измерено за вычетом фона, 2 — сумма модели, дальше вклады звеньев.
    names, spec = load_spectrum(os.path.join(d, "forward_components.csv"))
    if names:
        names[1] = "измерено минус фон"
        names[2] = "сумма расчёта"
    bg_sum = 0.0
    meas_sum = sum(r[1] for r in spec) if spec else 0.0
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
    # Рис. 5 — тот же расчёт вторым способом. Столбец «метод2_<возраст>лет»
    # берётся из forward_spectra.csv и подаётся как суммарная кривая; вклады
    # отдельных звеньев вторым способом отдельно не выводятся, поэтому
    # заливка на этом рисунке повторяет разбиение первого способа,
    # отмасштабированное отношением сумм.
    nm_f, sp_f = load_spectrum(os.path.join(d, "forward_spectra.csv"))
    names_m2, spec_m2, acts_m2 = [], [], []
    if nm_f and spec:
        j2 = next((i for i, n in enumerate(nm_f)
                   if n.startswith("метод2_")), None)
        j1 = next((i for i, n in enumerate(nm_f)
                   if n.startswith("метод1_")), None)
        if j2 is not None and j1 is not None and len(sp_f) == len(spec):
            names_m2 = list(names)
            names_m2[2] = "сумма расчёта (метод 2)"
            spec_m2 = []
            for i, row in enumerate(spec):
                tot1 = sp_f[i][j1]
                tot2 = sp_f[i][j2]
                k = (tot2 / tot1) if tot1 > 0 else 0.0
                spec_m2.append([row[0], row[1], round(tot2, 4)]
                               + [round(v * k, 4) for v in row[3:]])

    # Таблицы прямой задачи
    fwd_acts = opt("forward_activities.csv")
    fwd_bands = opt("forward_bands.csv")
    # Согласие способов: отношение сумм по полосам, считается здесь по тем же
    # данным, что и рисунки, чтобы число на странице и кривая не разошлись.
    m2m1 = []
    if nm_f and sp_f:
        j2 = next((i for i, n in enumerate(nm_f)
                   if n.startswith("метод2_")), None)
        j1 = next((i for i, n in enumerate(nm_f)
                   if n.startswith("метод1_")), None)
        if j1 is not None and j2 is not None:
            for lo, hi, lab in ((50, 72, "K-серия вольфрама"),
                                (72, 100, "K-серия дочерних"),
                                (150, 300, "полоса 238,63"),
                                (500, 650, "полоса 583,19"),
                                (850, 1000, "полоса 911,20"),
                                (2500, 2700, "полоса 2614,51"),
                                (3100, 3300, "сумм-пик 2614+583")):
                s1 = sum(r[j1] for r in sp_f if lo <= r[0] < hi)
                s2 = sum(r[j2] for r in sp_f if lo <= r[0] < hi)
                if s1 > 0:
                    m2m1.append(["%d–%d кэВ, %s" % (lo, hi, lab),
                                 "%.4f" % (s2 / s1)])
    # Возраст ряда. Файл `forward_age.csv` держит НЕСКОЛЬКО разнородных
    # блоков в одном файле: таблицу пар (заголовок «E_числитель_кэВ»),
    # чувствительность к ширине окна (заголовок «×ПШПВ» — на странице не
    # показывается, диагностика Б2 уже свёрнута в честный диапазон), разбор
    # состава окон (заголовок «E_кэВ») и итоговую строку («итог_точка_лет; …»).
    # Общей ширины у блоков нет, поэтому раздел определяется по первой ячейке
    # строки — простой автомат состояний. Каждый новый заголовок ОБЯЗАН явно
    # переключить `section`, иначе следующий блок молча утечёт в предыдущий
    # (найдено при добавлении блока ширины окна: без ветки «×ПШПВ» его строки
    # дописывались в конец age_pairs и ломали таблицу на странице).
    age_pairs, age_comp, age_head = [], [], {}
    section = None
    for r in opt("forward_age.csv"):
        if not r:
            continue
        tag = r[0]
        if tag == "E_числитель_кэВ":
            section = "pairs"
        elif tag == "×ПШПВ":
            section = None            # диагностика Б2, на странице не нужна
        elif tag == "E_кэВ":
            section = "comp"
        elif tag.startswith("итог_"):
            for i in range(0, len(r) - 1, 2):
                age_head[r[i]] = r[i + 1]
            section = None
            continue
        if section == "pairs":
            age_pairs.append(r)
        elif section == "comp":
            age_comp.append(r)

    payload = dict(names=names, spec=spec, acts=acts, lines=lines, cal=cal,
                   spectra=spectra,
                   names_m2=names_m2, spec_m2=spec_m2, acts_m2=acts_m2,
                   fwd_acts=fwd_acts, fwd_bands=fwd_bands, m2m1=m2m1,
                   src_head=src_head,
                   head=head,
                   bg_sum=round(bg_sum, 1),
                   bg_share=round(100.0 * bg_sum / meas_sum, 2) if meas_sum else 0,
                   src_fwd=opt("wt20_source_forward.csv"),
                   age_pairs=age_pairs, age_comp=age_comp, age_head=age_head,
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
