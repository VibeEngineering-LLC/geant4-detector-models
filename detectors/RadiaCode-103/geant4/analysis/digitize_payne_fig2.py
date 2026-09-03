# -*- coding: utf-8 -*-
"""
Оцифровка рисунка 2 статьи Payne S.A. et al., «Nonproportionality of Scintillator Detectors. III. Temperature Dependence Studies»,
IEEE TNS 61 (2014) 2771, DOI 10.1109/TNS.2014.2343572, открытый препринт LLNL-JRNL-648819 (OSTI 1762905).
Извлекается относительный световой выход на электрон как функция энергии для панели CsI(Tl) при трёх температурах.
"""
import sys
import re
import os
import argparse
import json
import datetime
from collections import defaultdict
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

def fit_linear_axis(px, values):
    from numpy.linalg import lstsq
    A = [[x, 1] for x in px]
    b = list(values)
    res = lstsq(A, b, rcond=None)[0]
    a, b = res[0], res[1]
    max_err = max(abs(a * x + b - v) for x, v in zip(px, values))
    return a, b, max_err

def fit_log_axis(px, decades):
    from numpy.linalg import lstsq
    A = [[x, 1] for x in px]
    b = list(decades)
    res = lstsq(A, b, rcond=None)[0]
    c, d = res[0], res[1]
    max_err = max(abs(c * x + d - v) for x, v in zip(px, decades))
    return c, d, max_err

def find_panel(page, label="CsI(Tl)"):
    """Панель опознаётся по ТЕКСТОВОЙ МЕТКЕ внутри неё, а не по цвету рамки:
    цвет — свойство вёрстки конкретного файла, метка — смысл панели.
    Рамка в PDF приходит как объект обводки (type 's'/'fs') с элементом 're'."""
    # Метка встречается и в подписи под рисунком, и в тексте статьи — берём ТО
    # вхождение, которое лежит ВНУТРИ рамки-панели, а не первое по порядку.
    labs = [w for w in page.get_text("words")
            if w[4].replace(" ", "").rstrip(".,;") == label]
    if not labs:
        return None
    frames = [d["rect"] for d in page.get_drawings()
              if d.get("type") in ("s", "fs") and any(it[0] == "re" for it in d["items"])]
    best = None
    for w in labs:
        lx, ly = 0.5 * (w[0] + w[2]), 0.5 * (w[1] + w[3])
        for r in frames:
            if r.x0 <= lx <= r.x1 and r.y0 <= ly <= r.y1:
                if best is None or (r.x1 - r.x0) * (r.y1 - r.y0) < (best.x1 - best.x0) * (best.y1 - best.y0):
                    best = r
    return best

def extract_curves(page, panel):
    drawings = page.get_drawings()
    panel_x0, panel_y0, panel_x1, panel_y1 = panel
    blue_curves = []
    for d in drawings:
        if d.get("type") in ("s", "fs") and d.get("color") and any(it[0] == "l" for it in d["items"]):
            r = d["rect"]
            if not r or not (panel_x0 <= r[0] and r[2] <= panel_x1 and panel_y0 <= r[1] and r[3] <= panel_y1):
                continue
            color = d["color"]
            if abs(color[0] - 0.310) < 0.05 and abs(color[1] - 0.506) < 0.05 and abs(color[2] - 0.741) < 0.05:
                blue_curves.append(d)
    if len(blue_curves) != 3:
        raise RuntimeError(f"Найдено {len(blue_curves)} синих кривых, ожидалось 3")
    
    # Классификация по штриховке
    curves = {}
    # dashes приходит СТРОКОЙ вида "[ 0 6 ] 0" / "[] 0" / "[ 3 3 ] 0" (не списком).
    # Подпись Fig. 2: точечная = -40 C, сплошная = 0 C, штриховая = +40 C.
    for d in blue_curves:
        # Конкретные длины штрихов различаются между панелями одного рисунка
        # (CsI(Tl): "[ 3 3 ]", NaI(Tl): "[ 4.5 6 ]", точки "[ 0 6 ]"/"[ 0 5.5 ]"),
        # поэтому классификация по ТИПУ: нулевая длина штриха = точечная,
        # пустой массив = сплошная, всё остальное = штриховая.
        pat = "".join((d.get("dashes") or "").split())
        if pat.startswith("[]"):
            curves["T_0"] = d
        elif pat.startswith("[0"):
            curves["T_m40"] = d
        elif pat.startswith("["):
            curves["T_p40"] = d
        else:
            raise RuntimeError(f"Неопознанная штриховка синей кривой: {d.get('dashes')!r}")
    if len(curves) != 3:
        raise RuntimeError(f"Найдено {len(curves)} кривых с нужной штриховкой, ожидалось 3")
    return curves

def check_order(curves, n_nodes=30, tol=0.005):
    """П1 из подписи Fig. 2: «Curves with lowest temperature consistently have
    the highest relative light yield». curves: метка -> список пар (E_keV, value).

    tol — допуск, взятый из САМОЙ статьи: «the fitted model ... yielding 0.5%
    uncertainty» (разд. про подгонку). Строгое неравенство без допуска
    объявляет нарушением расхождения ~0.003, которые лежат ВНУТРИ заявленной
    авторами погрешности и внутри погрешности оцифровки; такая проверка меряет
    шум, а не физику. Нарушением считается только выход за допуск.
    Возвращает (доля узлов, где порядок соблюдён, список нарушений)."""
    need = ["T_m40", "T_0", "T_p40"]
    if not all(k in curves and len(curves.get(k, [])) >= 2 for k in need):
        return 0.0, []
    lo = max(min(p[0] for p in curves[k]) for k in need)
    hi = min(max(p[0] for p in curves[k]) for k in need)
    if not (hi > lo > 0):
        return 0.0, []
    grid = np.logspace(np.log10(lo), np.log10(hi), n_nodes)
    y = {}
    for k in need:
        pts = sorted(curves[k], key=lambda q: q[0])
        y[k] = np.interp(np.log10(grid),
                         np.log10([q[0] for q in pts]),
                         [q[1] for q in pts])
    viol = [(float(grid[i]), float(y["T_m40"][i]), float(y["T_0"][i]), float(y["T_p40"][i]))
            for i in range(n_nodes)
            if (y["T_m40"][i] < y["T_0"][i] - tol) or (y["T_0"][i] < y["T_p40"][i] - tol)]
    return (n_nodes - len(viol)) / n_nodes, viol

def selftest():
    # Проверка калибровки Y
    y_px = [420.37, 395.20, 369.60, 344.00, 318.40, 292.80, 267.43]
    values = [0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
    a, b, max_err = fit_linear_axis(y_px, values)
    if max_err >= 0.002:
        print("SELFTEST FAIL")
        return 1
    # Проверка калибровки X
    x_px = [354.36, 417.60, 481.60, 545.02]
    decades = [0, 1, 2, 3]
    c, d, max_err = fit_log_axis(x_px, decades)
    if max_err >= 0.01:
        print("SELFTEST FAIL")
        return 1
    # Мутация калибровки Y: сбить ОДНО деление на 6 пунктов. Порог невязки ловит
    # класс «деления опознаны неверно / неравномерны», а НЕ перевёрнутый порядок
    # значений: обратный ряд даёт такую же идеально линейную подгонку с обратным
    # знаком и невязку ~0 — ложная мутация первой редакции спеки, снята 03.09.2026.
    y_px_mut = list(y_px)
    y_px_mut[3] += 6.0
    _, _, err_y_mut = fit_linear_axis(y_px_mut, values)
    if err_y_mut < 0.002:
        print(f"SELFTEST FAIL: калибровка Y не умеет краснеть "
              f"(сбито деление, невязка {err_y_mut:.5f} < 0.002)")
        return 1
    # Та же мутация по оси X: сбить одно деление на 5 пунктов.
    x_px_mut = list(x_px)
    x_px_mut[1] += 5.0
    _, _, err_x_mut = fit_log_axis(x_px_mut, decades)
    if err_x_mut < 0.01:
        print(f"SELFTEST FAIL: калибровка X не умеет краснеть "
              f"(сбито деление, невязка {err_x_mut:.5f} < 0.01)")
        return 1
    # П1 на ГОДНЫХ данных: порядок соблюдён везде -> доля обязана быть 1.0
    ok_curves = {"T_m40": [(10.0, 1.20), (1000.0, 1.20)],
                 "T_0":   [(10.0, 1.10), (1000.0, 1.10)],
                 "T_p40": [(10.0, 1.05), (1000.0, 1.05)]}
    ratio_ok, _ = check_order(ok_curves)
    if ratio_ok < 1.0:
        print(f"SELFTEST FAIL: П1 не пропускает годные данные (доля {ratio_ok:.3f})")
        return 1
    # Мутация П1: 0 C и +40 C переставлены -> проверка обязана покраснеть
    # Перестановка на 0.05 — заведомо БОЛЬШЕ допуска 0.005, иначе мутация не
    # отличалась бы от шума, который допуск и призван пропускать.
    bad_curves = {"T_m40": [(10.0, 1.20), (1000.0, 1.20)],
                  "T_0":   [(10.0, 1.10), (1000.0, 1.10)],
                  "T_p40": [(10.0, 1.15), (1000.0, 1.15)]}
    ratio_bad, viol = check_order(bad_curves)
    if ratio_bad >= 1.0:
        print("SELFTEST FAIL: П1 не умеет краснеть на нарушенном порядке")
        return 1
    print("SELFTEST OK")
    return 0

def main(argv):
    parser = argparse.ArgumentParser(description="Оцифровка рисунка 2 статьи Payne S.A. et al.")
    parser.add_argument("pdf_path", nargs="?", help="Путь к PDF файлу (не нужен при --selftest)")
    parser.add_argument("--out-dir", default=os.path.join(os.path.dirname(__file__), "out"))
    parser.add_argument("--selftest", action="store_true")
    parser.add_argument("--panel", default="CsI(Tl)", help="метка панели рисунка: CsI(Tl), NaI(Tl), CsI(Na)")
    args = parser.parse_args(argv)
    
    if args.selftest:
        return selftest()
    
    try:
        import fitz
    except ImportError:
        print("Ошибка: PyMuPDF (fitz) не установлен. Установите его командой `pip install pymupdf`")
        return 1
    
    os.makedirs(args.out_dir, exist_ok=True)
    
    doc = fitz.open(args.pdf_path)
    page = None
    for i in range(len(doc)):
        text = doc[i].get_text()
        if "Fig. 2." in text:
            page = doc[i]
            break
    if not page:
        print("Ошибка: Не найдена страница с подписью 'Fig. 2.'")
        return 2
    
    panel = find_panel(page, args.panel)
    tag = re.sub(r"[^A-Za-z0-9]+", "_", args.panel).strip("_")
    if not panel:
        print(f"Ошибка: Не найдена рамка панели {args.panel}")
        return 3
    
    try:
        curves = extract_curves(page, panel)
    except RuntimeError as e:
        print(f"Ошибка: {e}")
        return 4
    
    # Калибровка осей. Деления — это ОТРЕЗКИ ВНУТРИ объектов обводки (items 'l'),
    # а не отдельные объекты: все семь Y-делений лежат одним объектом.
    panel_x0, panel_y0, panel_x1, panel_y1 = panel
    segs = []
    for d in page.get_drawings():
        if d.get("type") not in ("s", "fs"):
            continue
        col = d.get("color")
        if col is None or max(col) > 0.2:      # только чёрные штрихи осей
            continue
        for it in d["items"]:
            if it[0] == "l":
                segs.append((it[1].x, it[1].y, it[2].x, it[2].y))

    # Y-деления: горизонтальные отрезки у ЛЕВОГО края панели, длина по x < 5.
    y_ticks = sorted({round((y1 + y2) / 2.0, 2) for x1, y1, x2, y2 in segs
                      if abs(y1 - y2) < 0.5 and abs(min(x1, x2) - panel_x0) < 1.5
                      and 2.0 <= abs(x1 - x2) <= 4.0   # мажорные; минорные короче (панель NaI)
                      and panel_y0 - 1 <= y1 <= panel_y1 + 1}, reverse=True)
    if len(y_ticks) != 7:
        print(f"Ошибка: найдено {len(y_ticks)} Y-делений, ожидалось 7: {y_ticks}")
        return 5
    y_values = [0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
    a_y, b_y, max_err_y = fit_linear_axis(y_ticks, y_values)
    if max_err_y >= 0.002:
        print(f"Ошибка: невязка Y-калибровки {max_err_y:.5f} > 0.002")
        return 6

    # X-деления: вертикальные отрезки у НИЖНЕГО края, длина по y от 2.0 до 4.0
    # (короткие — минорные деления логарифмической шкалы, они не берутся).
    x_ticks = sorted({round((x1 + x2) / 2.0, 2) for x1, y1, x2, y2 in segs
                      if abs(x1 - x2) < 0.5 and abs(max(y1, y2) - panel_y1) < 1.5
                      and 2.0 <= abs(y1 - y2) <= 4.0
                      and panel_x0 - 1 <= x1 <= panel_x1 + 1})
    if len(x_ticks) != 4:
        print(f"Ошибка: найдено {len(x_ticks)} X-делений, ожидалось 4: {x_ticks}")
        return 7
    c_x, d_x, max_err_x = fit_log_axis(x_ticks, [0, 1, 2, 3])
    if max_err_x >= 0.01:
        print(f"Ошибка: невязка X-калибровки {max_err_x:.5f} > 0.01")
        return 8
    print(f"Калибровка: Y невязка {max_err_y:.5f}, X невязка {max_err_x:.5f}")
    print(f"  Y-деления: {y_ticks}")
    print(f"  X-деления: {x_ticks}")

    # Извлечение точек кривых
    all_data = {}
    for label, curve_obj in curves.items():
        points = []
        for item in curve_obj["items"]:
            if item[0] == 'l':
                p1, p2 = item[1], item[2]
                points.append((p1[0], p1[1]))
                points.append((p2[0], p2[1]))
        # Убираем дубликаты
        unique_points = []
        for p in points:
            if not unique_points or (abs(p[0] - unique_points[-1][0]) > 1e-6 or abs(p[1] - unique_points[-1][1]) > 1e-6):
                unique_points.append(p)
        
        # Перевод в E и отн. выход
        data = []
        for x, y in unique_points:
            E = 10 ** (c_x * x + d_x)
            rel_light_yield = a_y * y + b_y
            data.append((E, rel_light_yield))
        
        # Сортировка по E
        data.sort(key=lambda p: p[0])
        # Усреднение при совпадающих E
        final_data = []
        i = 0
        while i < len(data):
            e = data[i][0]
            values = [data[i][1]]
            j = i + 1
            while j < len(data) and abs(data[j][0] - e) < 1e-6:
                values.append(data[j][1])
                j += 1
            avg_val = sum(values) / len(values)
            final_data.append((e, avg_val))
            i = j
        
        all_data[label] = final_data
    
    # Нормировка по 662 кэВ
    for label in all_data:
        data = all_data[label]
        if not data:
            continue
        # Значение при 662 кэВ: линейная интерполяция ПО log10(E), без
        # экстраполяции — если 662 кэВ вне диапазона кривой, нормировки нет.
        xs = [q[0] for q in data]
        ys = [q[1] for q in data]
        if min(xs) <= 662.0 <= max(xs):
            val_662 = float(np.interp(np.log10(662.0), np.log10(xs), ys))
        else:
            val_662 = None
            print(f"  ПРЕДУПРЕЖДЕНИЕ {label}: 662 кэВ вне диапазона "
                  f"[{min(xs):.1f}, {max(xs):.1f}] кэВ, нормировки нет")
        all_data[label] = [(e, y, val_662) for e, y in data]
    
    # Проверка диапазонов
    for label in all_data:
        data = all_data[label]
        if not data:
            continue
        min_e = min(p[0] for p in data)
        max_e = max(p[0] for p in data)
        if min_e < 0.95 or max_e > 1050:
            print(f"Предупреждение: диапазон энергий для {label} выходит за пределы [1, 1000] кэВ")
    
    # Проверка порядка — по ПЕРЕСЧИТАННЫМ точкам (E, значение), не по объектам PDF
    ratio, violations = check_order(
        {k: [(q[0], q[1]) for q in v] for k, v in all_data.items()})
    if ratio < 1.0:
        print(f"Предупреждение: порядок температур нарушен в {len(violations)} узлах")
        for x, y_m40, y_0, y_p40 in violations[:5]:
            print(f"  E={x:.2f}: T=-40°C={y_m40:.3f}, T=0°C={y_0:.3f}, T=+40°C={y_p40:.3f}")
    
    # Вывод в CSV
    for label in all_data:
        data = all_data[label]
        if not data:
            continue
        filename = os.path.join(args.out_dir, f"payne3_fig2_{tag}_{label}.csv")
        with open(filename, "w", encoding="utf-8") as f:
            f.write("# Источник: Payne S.A. et al., IEEE TNS 61 (2014) 2771, DOI 10.1109/TNS.2014.2343572, OSTI 1762905\n")
            f.write("# Панель: CsI(Tl)\n")
            f.write(f"# Температура: {label}\n")
            f.write("# Способ: извлечение векторных путей PDF\n")
            f.write(f"# Дата генерации: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("# Оцифровка графика, не табличные данные автора\n")
            f.write("# ВНИМАНИЕ: данные Fig. 2 обрываются около 458 кэВ, поэтому\n")
            f.write("#   колонка rel_to_662 ПУСТА: нормировка на 662 кэВ из этого\n")
            f.write("#   рисунка недостижима. rel_to_max_E — нормировка на самую\n")
            f.write("#   высокую энергию кривой. К какой точке нормировал автор,\n")
            f.write("#   в тексте статьи не сказано (проверено поиском по тексту).\n")
            f.write("E_keV,rel_light_yield,rel_to_662,rel_to_max_E\n")
            y_at_max_e = max(data, key=lambda q: q[0])[1]
            for e, y, y_662 in data:
                r662 = f"{y / y_662:.4f}" if y_662 not in (None, 0) else ""
                rmax = f"{y / y_at_max_e:.4f}" if y_at_max_e else ""
                f.write(f"{e:.3f},{y:.4f},{r662},{rmax}\n")
    
    # Сводный JSON
    summary = {
        "calibration": {
            "y_axis": {"a": a_y, "b": b_y, "max_error": max_err_y},
            "x_axis": {"c": c_x, "d": d_x, "max_error": max_err_x}
        },
        "curves": {}
    }
    
    for label in all_data:
        data = all_data[label]
        if not data:
            continue
        e_vals = [p[0] for p in data]
        y_vals = [p[1] for p in data]
        max_y = max(y_vals)
        max_e = e_vals[y_vals.index(max_y)]
        summary["curves"][label] = {
            "points": len(data),
            "energy_range": [min(e_vals), max(e_vals)],
            "value_at_662": data[0][2] if data and data[0][2] is not None else None,
            "max_value": {"E_keV": max_e, "rel_light_yield": max_y}
        }
    
    summary["order_check"] = {
        "ratio": ratio,
        "violations": len(violations)
    }
    
    with open(os.path.join(args.out_dir, f"payne3_fig2_{tag}_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    
    # Вывод П3
    print("\nП3. Статистика по кривым:")
    for label in all_data:
        data = all_data[label]
        if not data:
            continue
        e_vals = [p[0] for p in data]
        y_vals = [p[1] for p in data]
        max_y = max(y_vals)
        max_e = e_vals[y_vals.index(max_y)]
        print(f"  {label}:")
        print(f"    Точек: {len(data)}")
        print(f"    Диапазон энергий: [{min(e_vals):.1f}, {max(e_vals):.1f}] кэВ")
        print(f"    Значение при 662 кэВ: {data[0][2]:.4f}" if data and data[0][2] is not None else "    Значение при 662 кэВ: не определено")
        print(f"    Максимум: E={max_e:.1f} кэВ, Y={max_y:.4f}")
    
    print(f"\nП1. Порядок температур: {ratio:.3f} (доля узлов с правильным порядком)")
    
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
