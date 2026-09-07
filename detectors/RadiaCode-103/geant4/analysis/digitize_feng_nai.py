# -*- coding: utf-8 -*-
import sys
import os
import json
import datetime
import numpy as np
from scipy import ndimage
from scipy.interpolate import interp1d
import argparse
from collections import defaultdict

sys.stdout.reconfigure(encoding="utf-8")

try:
    import fitz
    from digitize_payne_fig2 import fit_linear_axis, fit_log_axis
except ImportError as e:
    print("ОШИБКА: Не найден модуль digitize_payne_fig2", file=sys.stderr)
    sys.exit(2)

# Подпись рисунка приводится ДОСЛОВНО рядом с номером: ссылка одним номером
# непроверяема — в этой же статье под номером 8 идёт фотография установки,
# а не график (инцидент W-070).
FIG_CAPTIONS = {
    "fig11": "Fig. 11. Comparison of the light yield non-linearity of "
             "LaBr3(Ce), LaBr3(Ce,Sr), and NaI(Tl) crystals for 3-400 keV "
             "Compton electrons (page 8)",
    "fig16": "Fig. 16. The light yield non-linearity of NaI(Tl) crystal for "
             "Compton electrons, gamma rays and X-rays (page 9)",
}

def _fail(msg):
    print("SELFTEST FAIL: " + msg, file=sys.stderr)
    sys.exit(1)


def selftest():
    """Самопроверка на синтетике: PDF не читается.

    Каждая проверка ВЫЗЫВАЕТ функцию этого модуля и оценивает её возврат.
    Проверка, повторяющая вычисление своими силами или испытывающая стороннюю
    библиотеку, осталась бы зелёной при любой ошибке в проверяемой функции.
    Половина проверок — мутации: они показывают, что мера умеет краснеть,
    и что на годном входе она НЕ срабатывает.
    """
    # 1. Линейная ось: точный закон обязан дать нулевую невязку.
    px = [0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0]
    vals = [1.2 - 0.1 * i for i in range(9)]
    _, _, err = fit_linear_axis(px, vals)
    if err >= 1e-9:
        _fail("подгонка линейной оси не воспроизводит точный закон")

    # 2. Логарифмическая ось: позиции строятся ПО закону из значений.
    energies = [10.0, 20.0, 50.0, 100.0, 200.0, 400.0]
    decades = [float(np.log10(e)) for e in energies]
    span = decades[-1] - decades[0]
    pos = [(d - decades[0]) / span * 500.0 for d in decades]
    _, _, err = fit_log_axis(pos, decades)
    if err >= 1e-9:
        _fail("подгонка логарифмической оси не воспроизводит точный закон")

    # 3. Мутация калибровки: две метки сопоставлены не своим значениям.
    #    Реверс всего ряда мутацией НЕ является — при равномерных позициях он
    #    остаётся точной прямой (с другим знаком) и невязку не поднимает.
    bad = list(vals)
    bad[2], bad[5] = bad[5], bad[2]
    _, _, err = fit_linear_axis(px, bad)
    if err < 0.004:
        _fail("проверка калибровки не умеет краснеть")

    # 4. Отсев усов: ус 3x40 рядом с маркером 18x18 — принят ровно маркер.
    mask = np.zeros((60, 60), dtype=bool)
    mask[10:50, 10:13] = True          # ус погрешности
    mask[20:38, 30:48] = True          # тело маркера
    centers, _ = extract_markers(mask)
    if len(centers) != 1:
        _fail("отсев усов не работает: принято центров " + str(len(centers)))
    if abs(centers[0][0] - 38.5) > 2 or abs(centers[0][1] - 28.5) > 2:
        _fail("центр маркера смещён: " + str(centers[0]))

    # 5. Мутация: один ус без маркера не должен дать ни одной точки.
    mask = np.zeros((60, 60), dtype=bool)
    mask[10:50, 10:13] = True
    centers, _ = extract_markers(mask)
    if len(centers) != 0:
        _fail("ус принимается за маркер")

    # 6. Расщепление: два слипшихся квадрата 18x18 — два центра.
    comp = np.zeros((40, 60), dtype=bool)
    comp[10:28, 10:28] = True
    comp[10:28, 28:46] = True
    pieces = split_overlapping_components(comp, 18.0)
    if len(pieces) != 2:
        _fail("расщепление не работает: кусков " + str(len(pieces)))
    if abs(pieces[0][0] - pieces[1][0]) < 12:
        _fail("расщепление дало слишком близкие центры")

    # 7. Мутация: одиночный маркер резать нельзя (мера не должна срабатывать
    #    на годном входе — иначе она портит данные, а не защищает их).
    comp = np.zeros((40, 60), dtype=bool)
    comp[10:28, 20:38] = True
    pieces = split_overlapping_components(comp, 18.0)
    if len(pieces) != 1:
        _fail("расщепление режет одиночный маркер: кусков " + str(len(pieces)))

    # 8. Сверка: наборы, разведённые на 5 %, обязаны быть признаны разными.
    d1 = [(10.0, 1.00), (30.0, 1.10), (100.0, 1.05), (300.0, 1.00)]
    d2 = [(e, y * 1.05) for e, y in d1]
    med, mx, _ = compare_datasets(d1, d2)
    if med <= 1.0:
        _fail("сверка не умеет краснеть: медиана " + str(round(float(med), 3)))

    # 9. Мутация в обратную сторону: одинаковые наборы обязаны сойтись.
    med, mx, _ = compare_datasets(d1, list(d1))
    if med > 1e-9 or mx > 1e-9:
        _fail("сверка бракует совпадающие наборы")

    print("SELFTEST OK")
    sys.exit(0)

def get_figure_data(doc, page_idx, bbox, zoom):
    page = doc[page_idx]
    matrix = fitz.Matrix(zoom, zoom)
    pixmap = page.get_pixmap(matrix=matrix, clip=bbox, alpha=False)
    buf = np.frombuffer(pixmap.samples, dtype=np.uint8)
    img = buf.reshape((pixmap.height, pixmap.width, pixmap.n))[:, :, :3].astype(np.int16)
    return img

def find_border(img):
    """Границы поля графика — ВНУТРЕННИЕ, за вычетом всей толщины рамки.

    Рамка нарисована в несколько пикселей толщиной. Если отступить от неё
    только на один пиксель, остаток рамки останется в поле и будет принят за
    метку в КАЖДОМ столбце (наблюдалось: 1271 метка длиной 2 пикселя).
    """
    dark = np.max(img, axis=2) < 100
    rows = np.sum(dark, axis=1) / img.shape[1]
    cols = np.sum(dark, axis=0) / img.shape[0]
    row_indices = np.where(rows > 0.8)[0]
    col_indices = np.where(cols > 0.8)[0]
    if len(row_indices) < 2 or len(col_indices) < 2:
        return None

    def inner(indices):
        s = set(int(i) for i in indices)
        a = int(indices[0])
        while a + 1 in s:
            a += 1
        b = int(indices[-1])
        while b - 1 in s:
            b -= 1
        return a, b

    top, bottom = inner(row_indices)
    left, right = inner(col_indices)
    return (top, bottom, left, right)

def find_ticks(img, axis):
    dark = np.max(img, axis=2) < 100
    if axis == 'x':
        ticks = []
        for col in range(img.shape[1]):
            n = 0
            row = img.shape[0] - 1
            while row > 0 and dark[row, col]:
                n += 1
                row -= 1
            if n > 0:
                ticks.append((col, n))
        return ticks
    elif axis == 'y':
        ticks = []
        for row in range(img.shape[0]):
            n = 0
            col = 0
            while col < img.shape[1] and dark[row, col]:
                n += 1
                col += 1
            if n > 0:
                ticks.append((row, n))
        return ticks

def merge_ticks(ticks, threshold=3):
    if not ticks:
        return []
    sorted_ticks = sorted(ticks, key=lambda x: x[0])
    merged = [sorted_ticks[0]]
    for tick in sorted_ticks[1:]:
        if abs(tick[0] - merged[-1][0]) <= threshold:
            # Склеиваем
            avg_pos = (tick[0] + merged[-1][0]) // 2
            merged[-1] = (avg_pos, max(tick[1], merged[-1][1]))
        else:
            merged.append(tick)
    return merged
def select_large_ticks(ticks):
    """Отбор подписанных (крупных) меток по длине штриха.

    Мелкие метки лог-шкалы вдвое-втрое короче подписанных. Порог — середина
    между самой короткой и самой длинной найденной длиной. Без этого отбора
    в калибровку попадают все метки подряд, и сопоставление их значениям
    становится бессмысленным.
    """
    if not ticks:
        return []
    lengths = [n for _, n in ticks]
    thr = 0.5 * (min(lengths) + max(lengths))
    return [t for t in ticks if t[1] >= thr]


def calibrate_axis(ticks, values):
    positions = np.array([t[0] for t in ticks])
    log_values = np.log10(values)
    a, b, err = fit_log_axis(positions, log_values)
    return a, b, err

def find_reference_lines(mask):
    """Строки, где доля пикселей заданного цвета выше 0.3.

    Порог низкий, потому что пунктирная линия заполняет примерно половину
    ширины поля. Соседние строки одной линии объединяются в её центр.
    """
    frac = mask.sum(axis=1) / float(mask.shape[1])
    rows = np.where(frac > 0.3)[0]
    if rows.size == 0:
        return []
    groups = []
    start = prev = int(rows[0])
    for r in rows[1:]:
        r = int(r)
        if r - prev <= 1:
            prev = r
            continue
        groups.append((start + prev) // 2)
        start = prev = r
    groups.append((start + prev) // 2)
    return groups

def color_mask(img, name):
    """Булевы маски цветов рисунка.

    Разности каналов считаются на знаковом типе: на uint8 выражение b - r
    переполняется и даёт ложные срабатывания.
    """
    r = img[:, :, 0].astype(np.int32)
    g = img[:, :, 1].astype(np.int32)
    b = img[:, :, 2].astype(np.int32)
    if name == "magenta":
        return (r > 150) & (b > 150) & (g < 120)
    if name == "red":
        return (r > 150) & (g < 90) & (b < 90)
    if name == "blue":
        return (b > 110) & (b - r > 60) & (b - g > 60)
    if name == "black":
        return (r < 90) & (g < 90) & (b < 90)
    print("ОТКАЗ: неизвестный цвет " + str(name), file=sys.stderr)
    sys.exit(2)


def split_overlapping_components(component, w_marker):
    """Разрезание слипшихся маркеров по провалам плотности.

    N найденных провалов дают N+1 кусок. Провал засчитывается, только если
    отстоит от краёв компоненты не меньше чем на 0.7 ширины маркера, — иначе
    одиночный маркер с неровным краем был бы разрезан надвое.
    """
    rows, cols = np.where(component)
    if not rows.size:
        return []
    profile = np.bincount(cols, minlength=component.shape[1]).astype(float)
    smoothed = np.convolve(profile, np.ones(3) / 3.0, mode="same")
    lo, hi = int(cols.min()), int(cols.max())
    cuts = []
    for p in range(lo + 1, hi):
        if smoothed[p] <= smoothed[p - 1] and smoothed[p] <= smoothed[p + 1]:
            if p - lo >= 0.7 * w_marker and hi - p >= 0.7 * w_marker:
                if not cuts or p - cuts[-1] >= 0.7 * w_marker:
                    cuts.append(p)
    bounds = [lo] + cuts + [hi + 1]
    out = []
    for i in range(len(bounds) - 1):
        sub = np.zeros_like(component)
        sub[:, bounds[i]:bounds[i + 1]] = component[:, bounds[i]:bounds[i + 1]]
        if sub.sum() == 0:
            continue
        rr, cc = np.where(sub)
        out.append((float(cc.mean()), float(rr.mean())))
    return out


def find_legend_box(img, pad=10):
    """Область легенды — по скоплению ТЕКСТА, а не по её рамке.

    Рамка легенды на этих рисунках тонкая и светло-серая, отдельной связной
    чёрной компоненты не образует. Зато подписи легенды — буквы: их bbox
    заполнен меньше чем на три четверти, тогда как маркер-квадрат заполнен
    почти целиком, а маркер-треугольник имеет другой цвет. Объединённый bbox
    букв, лежащих в правой нижней четверти поля, и есть легенда; он
    расширяется на pad, чтобы захватить образцы маркеров слева от подписей.
    """
    mask = color_mask(img, "black")
    labeled, n = ndimage.label(mask, structure=np.ones((3, 3)))
    h, w = mask.shape
    boxes = []
    for i in range(1, n + 1):
        rows, cols = np.where(labeled == i)
        r0, r1 = int(rows.min()), int(rows.max())
        c0, c1 = int(cols.min()), int(cols.max())
        bh, bw = r1 - r0 + 1, c1 - c0 + 1
        if bh < 0.015 * h or bh > 0.12 * h or bw > 0.40 * w:
            continue
        if len(rows) / float(bh * bw) > 0.75:
            continue          # плотная заливка — это маркер, а не буква
        if 0.5 * (r0 + r1) < 0.40 * h or 0.5 * (c0 + c1) < 0.35 * w:
            continue          # легенда занимает правую нижнюю часть поля
        boxes.append((r0, r1, c0, c1))
    if len(boxes) < 5:
        return None           # текста слишком мало, чтобы это была легенда
    r0 = min(b[0] for b in boxes)
    r1 = max(b[1] for b in boxes)
    c0 = min(b[2] for b in boxes)
    c1 = max(b[3] for b in boxes)
    return (max(0, r0 - pad), min(h, r1 + pad + 1),
            max(0, c0 - 6 * pad), min(w, c1 + pad + 1))


def extract_markers(mask, min_pixels=12):
    """Центры маркеров по булевой маске цвета.

    Ус погрешности того же цвета связан с маркером в одну компоненту, поэтому
    центр берётся по ШИРОКОЙ части: строки, где ширина не меньше половины
    максимальной. Ус — вертикальная линия в 2-4 пикселя, тело маркера шире.
    Возвращает (список центров, число разрезов слипшихся компонент).
    """
    labeled, n = ndimage.label(mask, structure=np.ones((3, 3)))
    raw = []
    for i in range(1, n + 1):
        comp = (labeled == i)
        if comp.sum() < min_pixels:
            continue
        rows, cols = np.where(comp)
        # Штрих оси прилегает к рамке поля; маркер данных — нет.
        if (rows.min() == 0 or cols.min() == 0
                or rows.max() == mask.shape[0] - 1
                or cols.max() == mask.shape[1] - 1):
            continue
        widths = {}
        for rr in np.unique(rows):
            cc = cols[rows == rr]
            widths[int(rr)] = int(cc.max() - cc.min() + 1)
        wmax = max(widths.values())
        if wmax < 5:
            continue
        h = int(rows.max() - rows.min() + 1)
        if h <= 4 and wmax / float(max(h, 1)) > 3.0:
            continue
        if wmax > 0.6 * mask.shape[1]:
            continue
        raw.append((comp, widths, wmax))
    if not raw:
        return [], 0
    w_marker = float(np.median([w for _, _, w in raw]))
    centers = []
    splits = 0
    for comp, widths, wmax in raw:
        if wmax > 1.6 * w_marker:
            pieces = split_overlapping_components(comp, w_marker)
            if len(pieces) > 1:
                splits += len(pieces) - 1
                centers.extend(pieces)
                continue
        wide = [rr for rr, ww in widths.items() if ww >= 0.5 * wmax]
        rows, cols = np.where(comp)
        sel = np.isin(rows, wide)
        centers.append((float(cols[sel].mean()), float(rows[sel].mean())))
    centers.sort()
    return centers, splits

def compare_datasets(data1, data2):
    # Строим сетку
    e_min = max(data1[0][0], data2[0][0])
    e_max = min(data1[-1][0], data2[-1][0])
    if e_min >= e_max:
        return 0, 0, 0
    grid = np.logspace(np.log10(e_min), np.log10(e_max), 25)
    # Интерполируем
    f1 = interp1d([x[0] for x in data1], [x[1] for x in data1], kind='linear', fill_value='extrapolate')
    f2 = interp1d([x[0] for x in data2], [x[1] for x in data2], kind='linear', fill_value='extrapolate')
    y1 = f1(grid)
    y2 = f2(grid)
    diff = np.abs(y1 - y2) / y1 * 100
    return np.median(diff), np.max(diff), grid[np.argmax(diff)]

def match_ticks_to_values(positions, expected, field_size, is_log,
                          must_cover=None, slope_sign=None):
    """Сопоставление найденных меток подписанным значениям — по расстояниям.

    Длина штриха НЕ отличает подписанную метку от неподписанной: на рис. 11
    значения 20 и 400 подписаны, но нарисованы короткими штрихами (3 пикселя
    против 17 у соседей). Поэтому подмножество подбирается перебором: верным
    считается то, при котором позиции ложатся на прямую — в log10(E) для
    логарифмической оси, в самой величине для линейной.

    Условие "все подписанные значения лежат внутри поля" снимает
    неоднозначность сдвига на декаду: наборы (10,50,100,200) и
    (20,100,200,400) дают ОДИНАКОВЫЕ расстояния, но во втором случае
    подписанное значение 10 оказалось бы за левым краем рисунка.

    Возвращает список (невязка, набор значений, a, b), отсортированный по
    невязке; пустой список означает, что ни один набор не подошёл.
    """
    from itertools import combinations
    if not positions or len(positions) > len(expected):
        return []
    fit = fit_log_axis if is_log else fit_linear_axis
    tol_px = 0.02 * field_size
    out = []
    for combo in combinations(expected, len(positions)):
      # Порядок значений вдоль оси заранее не известен, перебираются оба;
      # неверный отсеивается требованием знака наклона (slope_sign).
      for seq in (combo, combo[::-1]):
        vals = [float(np.log10(v)) for v in seq] if is_log else [float(v) for v in seq]
        a, b, err = fit(positions, vals)
        if slope_sign is not None and a * slope_sign <= 0:
            continue
        if a == 0:
            continue
        if must_cover is not None:
            # Диапазон, который ось ОБЯЗАНА покрывать, взят из подписи рисунка.
            # Это и снимает неоднозначность сдвига на декаду: у варианта с
            # удвоенными значениями левый край оси оказывается выше нижней
            # границы данных, заявленной автором.
            v0, v1 = b, a * field_size + b
            lo_ax, hi_ax = (v0, v1) if v0 <= v1 else (v1, v0)
            lo_need = float(np.log10(must_cover[0])) if is_log else float(must_cover[0])
            hi_need = float(np.log10(must_cover[1])) if is_log else float(must_cover[1])
            if lo_ax > lo_need or hi_ax < hi_need:
                continue
        inside = True
        for v in expected:
            t = float(np.log10(v)) if is_log else float(v)
            p = (t - b) / a
            if p < -tol_px or p > field_size + tol_px:
                inside = False
                break
        if inside:
            out.append((err, tuple(seq), a, b))
    out.sort(key=lambda z: z[0])
    return out

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdf")  # не обязателен: самопроверка PDF не читает
    parser.add_argument('--out-dir', default=os.path.join(os.path.dirname(__file__), 'out'))
    parser.add_argument('--zoom', type=float, default=6.0)
    parser.add_argument('--figures', default='fig11,fig16')
    parser.add_argument('--selftest', action='store_true')
    args = parser.parse_args()

    if args.selftest:
        selftest()
    if not args.pdf:
        print("ОТКАЗ: не задан --pdf (обязателен вне режима --selftest)", file=sys.stderr)
        sys.exit(2)

    os.makedirs(args.out_dir, exist_ok=True)

    doc = fitz.open(args.pdf)
    figures = args.figures.split(',')
    results = {}
    meta = {}

    for fig in figures:
        if fig == 'fig11':
            page_idx = 7
            bbox = fitz.Rect(45.354, 86.455, 288.679, 271.219)
            expected_ticks_x = [10, 20, 50, 100, 200, 400]
            cover_x = (3.0, 400.0)     # из подписи: "for 3-400 keV Compton electrons"
            expected_ticks_y = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
            ref_lines = {'magenta': 1.0, 'black_dotted': 0.9}
            marker_name = "blue"
        elif fig == 'fig16':
            page_idx = 8
            bbox = fitz.Rect(306.604, 86.457, 549.932, 271.568)
            expected_ticks_x = [10, 20, 50, 100, 200, 400, 1000]
            cover_x = (10.0, 1000.0)   # подписанные значения оси на рис. 16
            expected_ticks_y = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]
            ref_lines = {'red': 1.0, 'blue_dotted': 0.9, 'black_dotted': 1.1}
            marker_name = "black"
        else:
            print(f"Неизвестный рисунок {fig}", file=sys.stderr)
            sys.exit(2)

        # Рендер
        img = get_figure_data(doc, page_idx, bbox, args.zoom)
        border = find_border(img)
        if not border:
            print("Ошибка: не найдена рамка графика", file=sys.stderr)
            sys.exit(2)
        top, bottom, left, right = border
        field = img[top+1:bottom, left+1:right]

        # Метки осей
        ticks_x = find_ticks(field, 'x')
        ticks_y = find_ticks(field, 'y')

        # Слияние меток
        merged_x = merge_ticks(select_large_ticks(ticks_x))
        merged_y = merge_ticks(select_large_ticks(ticks_y))

        # Сопоставление меток значениям — по взаимным расстояниям.
        mx = match_ticks_to_values([t[0] for t in merged_x], expected_ticks_x,
                                   field.shape[1], is_log=True, must_cover=cover_x,
                                   slope_sign=1)
        my = match_ticks_to_values([t[0] for t in merged_y], expected_ticks_y,
                                   field.shape[0], is_log=False, slope_sign=-1)
        for nm, mt, tol in (("X", mx, 0.005), ("Y", my, 0.004)):
            if not mt:
                print(f"Ошибка: метки оси {nm} не сопоставились ни одному набору значений",
                      file=sys.stderr)
                sys.exit(2)
            if mt[0][0] >= tol:
                print(f"Ошибка: невязка по оси {nm} {mt[0][0]:.5f} >= {tol}", file=sys.stderr)
                sys.exit(2)
            if len(mt) > 1 and mt[1][0] < 3 * max(mt[0][0], 1e-9):
                print(f"Ошибка: сопоставление меток оси {nm} неоднозначно — "
                      f"{list(mt[0][1])} и {list(mt[1][1])} дают сравнимую невязку",
                      file=sys.stderr)
                sys.exit(2)
        err_x, combo_x, a, b = mx[0]
        err_y, combo_y, c, d = my[0]
        print(f"  ось X: метки {list(combo_x)}, невязка {err_x:.5f} декады")
        print(f"  ось Y: метки {list(combo_y)}, невязка {err_y:.5f}")

        # Проверка калибровки Y по опорным линиям — независимо от меток осей
        ref_measured = {}
        for lname, lval in ref_lines.items():
            cname = lname.split("_")[0]
            rows_found = find_reference_lines(color_mask(field, cname))
            if not rows_found:
                print("Ошибка: не найдена опорная линия " + lname, file=sys.stderr)
                sys.exit(2)
            vals = [c * rr + d for rr in rows_found]
            best = min(vals, key=lambda v: abs(v - lval))
            ref_measured[lname] = float(best)
            print(f"  репер {lname}: ожидалось {lval}, измерено {best:.4f}, "
                  f"расхождение {best - lval:+.4f}")
            if abs(best - lval) > 0.006:
                print("Ошибка: опорная линия " + lname + " не подтверждает шкалу",
                      file=sys.stderr)
                sys.exit(2)

        # Выделение маркеров
        legend = find_legend_box(field)
        mmask = color_mask(field, marker_name)
        if legend is not None:
            lr0, lr1, lc0, lc1 = legend
            mmask[lr0:lr1, lc0:lc1] = False
            print(f"  легенда исключена: строки {lr0}-{lr1}, столбцы {lc0}-{lc1}")
        else:
            print("  ВНИМАНИЕ: рамка легенды не найдена, её образцы могут попасть в данные")
        # Горизонтальная линия цвета маркера — это опорная линия рисунка,
        # а не данные: её звенья иначе становятся десятками ложных точек.
        for rr in find_reference_lines(mmask):
            mmask[max(0, rr - 5):rr + 6, :] = False
        markers, n_splits = extract_markers(mmask)
        if not markers:
            print("Ошибка: не найдено маркеров цвета " + marker_name, file=sys.stderr)
            sys.exit(2)
        print(f"  маркеров принято: {len(markers)}, разрезов слипшихся: {n_splits}")
        final_positions = markers

        # Перевод в физические величины
        data = []
        for x, y in final_positions:
            e = 10 ** (a * x + b)
            yield_val = c * y + d
            data.append((e, yield_val))

        data.sort(key=lambda x: x[0])
        results[fig] = data
        meta[fig] = {
            "field_box": [int(top), int(bottom), int(left), int(right)],
            "axis_x": {"a": float(a), "b": float(b), "residual_decades": float(err_x),
                       "ticks_used": [float(v) for v in combo_x]},
            "axis_y": {"c": float(c), "d": float(d), "residual": float(err_y),
                       "ticks_used": [float(v) for v in combo_y]},
            "reference_lines_expected": dict(ref_lines),
            "reference_lines_measured": dict(ref_measured),
            "legend_box": list(legend) if legend else None,
            "marker_color": marker_name,
            "n_points": len(data),
            "n_splits": int(n_splits),
        }
    # Сверка двух рисунков — это и есть проверка самой оцифровки
    exit_code = 0
    if "fig11" in results and "fig16" in results:
        med_diff, max_diff, max_pos = compare_datasets(results["fig11"], results["fig16"])
        print(f"Сверка рисунков: медианное расхождение {med_diff:.2f} %, "
              f"максимальное {max_diff:.2f} % при {max_pos:.1f} кэВ")
        if med_diff > 1.0 or max_diff > 3.0:
            print("РАСХОЖДЕНИЕ СВЕРХ ДОПУСКА: хотя бы одна оцифровка неверна",
                  file=sys.stderr)
            exit_code = 3

    # Вывод CSV
    for fig, data in results.items():
        filename = os.path.join(args.out_dir, f"feng2024_{fig}_NaI_compton.csv")
        with open(filename, 'w', encoding='utf-8') as f:
            print("# Feng et al., arXiv:2312.16658", file=f)
            print("# " + FIG_CAPTIONS[fig], file=f)
            print("# Маркер этого набора: " + meta[fig]["marker_color"], file=f)
            print("# Способ: оцифровка растрового изображения", file=f)
            print("# Множитель рендера: " + str(args.zoom), file=f)
            print("# Дата генерации: " + datetime.date.today().isoformat(), file=f)
            print("# Оцифровка графика, НЕ табличные данные автора", file=f)
            print("E_keV,rel_light_yield", file=f)
            for e, y in data:
                print(f"{e:.2f},{y:.4f}", file=f)

    # Вывод JSON
    summary = {}
    for fig, data in results.items():
        entry = dict(meta[fig])
        entry["energy_range_keV"] = [data[0][0], data[-1][0]]
        top_pt = max(data, key=lambda x: x[1])
        entry["max_yield"] = {"E_keV": top_pt[0], "value": top_pt[1]}
        summary[fig] = entry
    if 'fig11' in results and 'fig16' in results:
        med_diff, max_diff, max_pos = compare_datasets(results['fig11'], results['fig16'])
        summary["comparison"] = {
            "median_difference_percent": med_diff,
            "max_difference_percent": max_diff,
            "max_difference_energy_keV": max_pos
        }

    with open(os.path.join(args.out_dir, "feng2024_NaI_digitize.json"), 'w', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    # Вывод в stdout
    for fig, data in results.items():
        m = meta[fig]
        print("")
        print(f"Рисунок {fig} (маркер: {m['marker_color']})")
        print(f"  Поле графика: {m['field_box']}")
        print(f"  Ось X: метки {m['axis_x']['ticks_used']}, невязка {m['axis_x']['residual_decades']:.6f} декады")
        print(f"  Ось Y: метки {m['axis_y']['ticks_used']}, невязка {m['axis_y']['residual']:.6f}")
        print("  Опорные линии (проверка шкалы, независимая от меток):")
        for line, val in m["reference_lines_expected"].items():
            got = m["reference_lines_measured"].get(line)
            print(f"    {line}: ожидалось {val}, измерено {got:.4f}, расхождение {got - val:+.4f}")
        print(f"  Точек принято: {len(data)}, разрезов слипшихся: {m['n_splits']}")
        print(f"  Диапазон энергий: {data[0][0]:.1f} - {data[-1][0]:.1f} кэВ")
        top_pt = max(data, key=lambda x: x[1])
        print(f"  Максимум кривой: E={top_pt[0]:.1f} кэВ, Y={top_pt[1]:.4f}")

    if 'fig11' in results and 'fig16' in results:
        print("\nСверка двух рисунков:")
        print(f"  Медианное расхождение: {med_diff:.2f}%")
        print(f"  Максимальное расхождение: {max_diff:.2f}% при {max_pos:.1f} кэВ")

    # Требует толкования
    print("\nТРЕБУЕТ ТОЛКОВАНИЯ:")
    if 'fig11' in results:
        max_yield = max(results['fig11'], key=lambda x: x[1])
        print(f"  Максимум кривой fig11: E={max_yield[0]:.1f} кэВ, Y={max_yield[1]:.4f}")
        print("    Сравните с заявлением авторов о максимуме в 15.5% при 14 кэВ")
    if 'fig16' in results:
        max_yield = max(results['fig16'], key=lambda x: x[1])
        print(f"  Максимум кривой fig16: E={max_yield[0]:.1f} кэВ, Y={max_yield[1]:.4f}")
    if 'fig11' in results and 'fig16' in results:
        med_diff, max_diff, max_pos = compare_datasets(results['fig11'], results['fig16'])
        if med_diff > 1.0:
            print(f"  Расхождение выше допуска: медианное {med_diff:.2f}%")

    return exit_code

if __name__ == "__main__":
    sys.exit(main())
