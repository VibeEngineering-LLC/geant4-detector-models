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


def load_spectrum(path, step=6):
    """Спектр разложения. Прореживается: странице хватает шага 12 кэВ."""
    rows = read_csv(path)
    if not rows:
        return {}, []
    names = rows[0]
    data = [[float(v) for v in r] for r in rows[1:]]
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


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    d, tdir = sys.argv[1], sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(d, "wt20.html")

    names, spec = load_spectrum(os.path.join(d, "unfold_spectrum.csv"))
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

    payload = dict(names=names, spec=spec, acts=acts, lines=lines, cal=cal,
                   head=head,
                   src_fwd=opt("wt20_source_forward.csv"),
                   selfabs=opt("wt20_selfabsorption.csv"),
                   scatter=opt("wt20_source_scatter.csv"),
                   surflim=opt("wt20_surface_limit.csv"),
                   winscan=opt("wt20_window_scan.csv"),
                   img_scan=data_uri(os.path.join(ref,
                                                  "wt20-pitch-scanline.jpg")),
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
