import sys
import os
import subprocess
import csv
import json
import tempfile
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from adaptive_smooth import presmooth
SMOOTH = float(os.environ.get("PAGE_SMOOTH_TARGET", "0"))   # 0 — без сглаживания; иначе целевая относительная ошибка слоя
WIDTHS = {}

sys.stdout.reconfigure(encoding="utf-8")

def smear(path_in, path_out, emin, emax):
    # свёртка через spec_smear.py; возвращает (y, sig) в отсчётах/ч/кэВ
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "spec_smear.py")
    r = subprocess.run([sys.executable, script, path_in, path_out, "--emin", str(emin), "--emax", str(emax)],
                       capture_output=True, text=True, env={**os.environ, "PYTHONIOENCODING": "utf-8"})
    if r.returncode != 0:
        print("Ошибка в spec_smear.py:", r.stderr)
        sys.exit(2)
    with open(path_out, "r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return [float(x["counts_per_h_per_keV"]) for x in rows], [float(x["sigma_counts_per_h_per_keV"]) for x in rows]

def main():
    if len(sys.argv) < 4:
        print("Usage: make_page_data.py <merged_prefix> <lines_prefix> <out_data_js> [--emin 20] [--emax 10000]")
        sys.exit(1)

    args = sys.argv[1:]
    merged_prefix = args[0]
    lines_prefix = args[1]
    out_data_js = args[2]

    emin = 20
    emax = 10000

    lang = "ru"
    i = 3
    while i < len(args):
        if args[i] == "--emin":
            emin = int(args[i+1])
            i += 2
        elif args[i] == "--emax":
            emax = int(args[i+1])
            i += 2
        elif args[i] == "--lang":
            lang = args[i+1]
            i += 2
        else:
            i += 1

    path_to_spec_smear = os.path.join(os.path.dirname(__file__), "spec_smear.py")

    # Шаг 1: Смазка
    total_file = f"{merged_prefix}_total_all.csv"
    smeared_total = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".csv")
    try:
        result = subprocess.run([
            sys.executable, path_to_spec_smear,
            total_file, smeared_total.name,
            "--emin", str(emin),
            "--emax", str(emax)
        ], capture_output=True, text=True, env={**os.environ, "PYTHONIOENCODING": "utf-8"})
        if result.returncode != 0:
            print("Ошибка в spec_smear.py:", result.stderr)
            sys.exit(2)
    finally:
        smeared_total.close()

    # Шаг 2: Сбор данных
    series = []
    with open(smeared_total.name, "r") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    E0 = float(rows[0]["E_keV"])
    y_vals = [float(row["counts_per_h_per_keV"]) for row in rows]
    y_vals_total = list(y_vals)
    sig_vals = [float(row["sigma_counts_per_h_per_keV"]) for row in rows]

    # Полный спектр
    series.append({
        "id": "total",
        "label": "Full spectrum" if lang == "en" else "Полный спектр",
        "group": "total",
        "y": [float("%.4g" % y) for y in y_vals],
        "sig": [float("%.4g" % s) for s in sig_vals]
    })

    # Компоненты
    components = ["neutron", "proton", "mup", "mum", "em", "ep", "gamma", "k40"]
    labels = {
        "ru": {"neutron": "нейтроны", "proton": "протоны", "mup": "мюоны +", "mum": "мюоны −",
               "em": "электроны", "ep": "позитроны", "gamma": "гамма-кванты", "k40": "K-40 в людях"},
        "en": {"neutron": "neutrons", "proton": "protons", "mup": "muons +", "mum": "muons −",
               "em": "electrons", "ep": "positrons", "gamma": "gamma rays", "k40": "K-40 in people"}
    }[lang]
    tmpd = tempfile.mkdtemp()
    for comp in components:
        comp_file = f"{merged_prefix}_total_{comp}.csv"
        if not os.path.exists(comp_file):
            continue
        dd = np.loadtxt(comp_file, delimiter=",", skiprows=2, ndmin=2)
        pre = os.path.join(tmpd, comp + "_pre.csv")
        WIDTHS[comp] = presmooth(dd[:, 1], dd[:, 2], pre, SMOOTH)[emin:emax]
        y_c, _ = smear(pre, os.path.join(tmpd, comp + "_sm.csv"), emin, emax)
        series.append({"id": comp, "label": labels[comp], "group": "component", "y": [float("%.4g" % v) for v in y_c]})

    # Категории
    cat_file = f"{merged_prefix}_cat_all.csv"
    if os.path.exists(cat_file):
        with open(cat_file, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            cats = list(reader)
        CAT2 = {}
        if os.path.exists(f"{merged_prefix}_cat2_all.csv"):   # дисперсии категорий (cat_sigma.py)
            for r2 in csv.DictReader(open(f"{merged_prefix}_cat2_all.csv", encoding="utf-8")):
                CAT2[(int(r2["pg"]), int(r2["proc"]), int(r2["vol"]))] = np.array([float(r2[k]) for k in r2 if k.startswith("bin_")])

        sums = []
        for cat in cats:
            s = sum(float(cat[f"bin_{i}"]) for i in range(emin, emax))
            sums.append((s, cat))

        sums.sort(key=lambda t: t[0], reverse=True)
        top_cats = sums[:12]

        temp_dir = tempfile.TemporaryDirectory()
        try:
            total_y = [0.0] * (emax - emin)
            for s, cat in top_cats:
                # Создаем временный файл
                tmp_file = os.path.join(temp_dir.name, f"{cat['pg']}_{cat['proc']}_{cat['vol']}.csv")
                vals = [float(cat[f"bin_{i}"]) for i in range(len(cat) - 7)]
                var = CAT2.get((int(cat["pg"]), int(cat["proc"]), int(cat["vol"])))
                wid = presmooth(vals, var, tmp_file, SMOOTH)
                WIDTHS["origin_%s_%s_%s" % (cat["pg"], cat["proc"], cat["vol"])] = wid[emin:emax]
                for i in range(emin, emax):
                    total_y[i - emin] += vals[i]

                # Смазка
                smeared_file = tmp_file.replace(".csv", "_smeared.csv")
                result = subprocess.run([
                    sys.executable, path_to_spec_smear,
                    tmp_file, smeared_file,
                    "--emin", str(emin),
                    "--emax", str(emax)
                ], capture_output=True, text=True, env={**os.environ, "PYTHONIOENCODING": "utf-8"})
                if result.returncode != 0:
                    print("Ошибка в spec_smear.py:", result.stderr)
                    sys.exit(2)

                with open(smeared_file, "r", encoding="utf-8") as f:
                    reader = csv.DictReader(f)
                    rows = list(reader)
                y_vals = [float(row["counts_per_h_per_keV"]) for row in rows]
                label = f"{cat['pgname']} / {cat['procname']} / {cat['volname']}"
                series.append({
                    "id": f"origin_{cat['pg']}_{cat['proc']}_{cat['vol']}",
                    "label": label,
                    "group": "origin",
                    "y": [float("%.4g" % y) for y in y_vals]
                })
        finally:
            temp_dir.cleanup()

        # прочие происхождения = полный спектр минус сумма показанных категорий (не меньше нуля)
        shown = [sr["y"] for sr in series if sr["group"] == "origin"]
        rest = [max(0.0, t - sum(col)) for t, col in zip(y_vals_total, zip(*shown))]
        series.append({"id": "other", "label": ("other origins" if lang == "en" else "прочие происхождения"), "group": "origin", "y": [float("%.4g" % v) for v in rest]})

    # Шаг 3: Округление и подготовка
    E0 = emin + 0.5

    # Шаг 4: Линии
    lines = []
    unlisted = []
    all_candidates = []

    lines_file = f"{lines_prefix}_lines.csv"
    if os.path.exists(lines_file):
        with open(lines_file, "r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        all_candidates = [(float(r["E_keV"]), r["reaction"]) for r in rows]   # для поиска кандидата рядом с непойманным пиком
        rows = [r for r in rows if r["found"] == "1"]
        rows.sort(key=lambda x: float(x["signif"]), reverse=True)
        for row in rows[:80]:
            lines.append({
                "E": float(row["E_keV"]),
                "reaction": row["reaction"],
                "net": float(row["net_cph"]),
                "sig": float(row["sigma_cph"]),
                "signif": float(row["signif"]),
                "origin": row["origin"]
            })

    unlisted_file = f"{lines_prefix}_unlisted.csv"
    if os.path.exists(unlisted_file):
        with open(unlisted_file, "r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        rows.sort(key=lambda x: float(x["signif"]), reverse=True)
        for row in rows[:30]:
            E = float(row["E_keV"])
            near = min(all_candidates, key=lambda c: abs(c[0] - E), default=None) if all_candidates else None
            tol = max(3, E * 0.01)
            unlisted.append({
                "E": E,
                "net": float(row["net_cph"]),
                "sig": float(row["sigma_cph"]),
                "signif": float(row["signif"]),
                "origin": row["origin"],
                "nearReaction": near[1] if near and abs(near[0] - E) <= tol else None,
                "nearDE": round(near[0] - E, 1) if near and abs(near[0] - E) <= tol else None
            })

    # Шаг 5: Запись
    for sr in series:   # ширина сглаживания слоя (кэВ): медиана и максимум по 20–3000 кэВ
        w = WIDTHS.get(sr["id"])
        if w is not None and SMOOTH > 0:
            sr["sm"] = [int(np.median(w[:2980])), int(np.max(w[:2980])), round(float(np.mean(w[:2980] > 1)), 2)]
    data = {
        "E0": E0,
        "dE": 1,
        "unit": "counts per hour per keV" if lang == "en" else "отсчётов в час на кэВ",
        "series": series,
        "lines": lines,
        "unlisted": unlisted,
        "candidates": [{"E": e, "reaction": r} for e, r in all_candidates]   # все кандидаты из базы, вне зависимости от значимости
    }

    with open(out_data_js, "w", encoding="utf-8") as f:
        f.write("window.SPECTRA = ")
        json.dump(data, f, ensure_ascii=False, separators=(',', ': '))
        f.write(";")

    size_kb = os.path.getsize(out_data_js) / 1024
    num_series = len(series)
    print(f"Размер файла: {size_kb:.1f} KB, количество серий: {num_series}")

if __name__ == "__main__":
    main()
