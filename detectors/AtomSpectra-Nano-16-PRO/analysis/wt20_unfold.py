# -*- coding: utf-8 -*-
"""Разложение измеренного спектра по шаблонам нуклидов ряда тория.

Метод шаблонный: измеренный спектр представляется суммой откликов ОТДЕЛЬНЫХ
нуклидов, посчитанных Монте-Карло в той же геометрии, плюс измеренный фон.
Шаблон нормирован на ОДИН распад, поэтому коэффициенты разложения — сразу
активности в беккерелях, без промежуточной «эффективности по линии».

    измеренное(E) = Σ aᵢ · t · Tᵢ(E) + b · фон(E)

где aᵢ — активность нуклида i, Бк; t — время набора, с; Tᵢ — шаблон, отсчёты на
распад в канале E; b — множитель фона (ожидается около единицы).

Побочные пики отдельными компонентами НЕ вводятся: вылет аннигиляционных
квантов, обратное рассеяние, характеристический рентген вольфрама и суммирование
каскада возникают в переносе сами и лежат внутри шаблона своего нуклида. Тем
методика отличается от разложения на аналитические компоненты, где каждый такой
пик приходится перечислять руками.

Приборное разрешение навешивается ЗДЕСЬ: ПШПВ² = f0 + f1·E подбирается вместе с
активностями (перебором по сетке), потому что паспортного разрешения у прибора
нет, а измеренные ширины отдельных линий на мультиплетах ненадёжны.

    python analysis/wt20_unfold.py <спектр.xml> <каталог шаблонов> [каталог вывода]
"""
import csv
import io
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.environ.get("SPECTRAVIBE_ROOT")
if not _ROOT:
    raise SystemExit("не задана переменная окружения SPECTRAVIBE_ROOT")
sys.path.insert(0, os.path.join(_ROOT, "scripts"))

from gamma.io.atomspectra_xml import read_atomspectra_xml        # noqa: E402

E_MAX = 3700.0          # верх энергетической сетки, кэВ
E_STEP = 2.0            # шаг сетки разложения, кэВ
E_FIT = (60.0, 3300.0)  # окно подгонки

# Удельная активность Th-232, Бк/г. ВЫЧИСЛЕНА: A = ln2·N_A/(T½·M) при
# T½ = 4,41797·10¹⁷ с и M = 232,038054 а.е.м. (API МАГАТЭ, ENSDF).
SPEC_ACT_TH232 = 0.6931472 * 6.02214076e23 / (4.41796963644288e17 * 232.038054)

# Порядок и подписи компонент — по месту в ряду.
ORDER = [
    ("Th232", "Th-232", "#7a5c3a"),
    ("Ra228", "Ra-228", "#a08a5c"),
    ("Ac228", "Ac-228", "#d81b8c"),
    ("Th228", "Th-228", "#8a6d3b"),
    ("Ra224", "Ra-224", "#c98b1e"),
    ("Rn220", "Rn-220", "#6b8f3a"),
    ("Po216", "Po-216", "#9bb06a"),
    ("Pb212", "Pb-212", "#b07d2a"),
    ("Bi212", "Bi-212", "#2f6b34"),
    ("Tl208", "Tl-208", "#c8cf7a"),
    ("Po212", "Po-212", "#8fa0a8"),
]

ORDER_MAP = [(k, lab) for k, lab, _ in ORDER]

# Ветвление Bi-212: доля распадов, идущих через Tl-208. Число берётся из
# библиотеки МАГАТЭ при запуске (поле decay_% в файле линий), а не пишется сюда.
LIB = os.path.normpath(os.path.join(_HERE, "..", "reference", "nuclide-lines"))


def branch_to_tl208():
    """Доля распадов Bi-212 по альфа-ветви (на Tl-208), из файла линий Tl-208.

    В библиотеке МАГАТЭ поле `decay_%` строк Tl-208 — это доля РОДИТЕЛЬСКОГО
    распада, приводящая к этому нуклиду. Берётся оттуда, а не пишется числом:
    35,94 % — величина оценённая и может уточняться.
    """
    p = os.path.join(LIB, "212bi_gammas.csv")
    for r in csv.DictReader(io.open(p, encoding="utf-8")):
        d, pc = r.get("decay"), r.get("decay_%")
        if d and d.strip().upper().startswith("A") and pc:
            try:
                return float(pc) / 100.0
            except ValueError:
                pass
    raise SystemExit("в %s нет альфа-ветви Bi-212" % p)


def read_template(path):
    """Спектр Geant4: шапка + E_keV,counts. -> (dict шапки, E[], counts[])."""
    head, e, c = {}, [], []
    for ln in io.open(path, encoding="utf-8"):
        ln = ln.rstrip("\n")
        if ln.startswith("#"):
            if "=" in ln:
                k, v = ln.lstrip("# ").split("=", 1)
                head[k.strip()] = v.strip()
            continue
        if not ln or ln.startswith("E_keV"):
            continue
        a, b = ln.split(",")
        e.append(float(a))
        c.append(float(b))
    return head, np.array(e), np.array(c)


def poly(coefs, x):
    out = np.zeros_like(np.asarray(x, dtype=float))
    for k, c in enumerate(coefs):
        out = out + c * np.asarray(x, dtype=float) ** k
    return out


def read_correction(path):
    """Поправки калибровки из wt20_calibration.py. -> dict имя->коэффициенты."""
    out = {}
    if not os.path.exists(path):
        return out
    for row in csv.reader(io.open(path, encoding="utf-8")):
        if not row or row[0].startswith("#") or row[0] == "спектр":
            continue
        out[row[0]] = [float(v) for v in row[1:] if v not in ("", None)]
    return out


def rebin_to_grid(counts, cal, corr, grid_edges):
    """Отсчёты по каналам -> отсчёты по энергетической сетке.

    Канал считается равномерно заполненным по своему энергетическому интервалу
    [E(ch−½), E(ch+½)]; его отсчёты делятся между ячейками сетки по перекрытию.
    Так сохраняется ПОЛНОЕ число отсчётов — проверяется в конце.
    """
    n = len(counts)
    ch = np.arange(n)
    lo = poly(cal, ch - 0.5)
    hi = poly(cal, ch + 0.5)
    if corr:
        lo = poly(corr, lo)
        hi = poly(corr, hi)
    out = np.zeros(len(grid_edges) - 1)
    g0, g1 = grid_edges[0], grid_edges[-1]
    step = grid_edges[1] - grid_edges[0]
    for i in range(n):
        c = counts[i]
        if c <= 0:
            continue
        a, b = lo[i], hi[i]
        if b <= g0 or a >= g1 or b <= a:
            continue
        a, b = max(a, g0), min(b, g1)
        ka = int((a - g0) / step)
        kb = int((b - g0) / step)
        if ka == kb:
            out[ka] += c
            continue
        w = b - a
        for k in range(ka, min(kb + 1, len(out))):
            left = max(a, g0 + k * step)
            right = min(b, g0 + (k + 1) * step)
            if right > left:
                out[k] += c * (right - left) / w
    return out


def broaden(raw_e, raw_c, grid_centres, f0, f1):
    """Свёртка линейного спектра с гауссовым разрешением ПШПВ² = f0 + f1·E."""
    out = np.zeros(len(grid_centres))
    step = grid_centres[1] - grid_centres[0]
    lo0 = grid_centres[0] - 0.5 * step
    for e, c in zip(raw_e, raw_c):
        if c <= 0:
            continue
        fw = math.sqrt(max(f0 + f1 * e, 1.0))
        s = fw / 2.3548
        k0 = int((e - 4 * s - lo0) / step)
        k1 = int((e + 4 * s - lo0) / step) + 1
        k0 = max(k0, 0)
        k1 = min(k1, len(out))
        if k1 <= k0:
            continue
        x = grid_centres[k0:k1]
        g = np.exp(-0.5 * ((x - e) / s) ** 2)
        ssum = g.sum()
        if ssum <= 0:
            continue
        out[k0:k1] += c * g / ssum
    return out


def nnls_fit(A, y, w):
    """Взвешенный МНК с неотрицательными коэффициентами."""
    from scipy.optimize import nnls
    Aw = A * w[:, None]
    yw = y * w
    x, rnorm = nnls(Aw, yw)
    return x, rnorm


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    src, tdir = sys.argv[1], sys.argv[2]
    outdir = sys.argv[3] if len(sys.argv) > 3 else os.path.dirname(src)
    os.makedirs(outdir, exist_ok=True)

    spec = read_atomspectra_xml(src)
    if isinstance(spec, (list, tuple)):
        spec = spec[0]
    bg = getattr(spec, "background_embedded", None)
    t_smp = float(spec.real_time)

    corr = read_correction(os.path.join(outdir, "calibration_fitted.csv"))
    if not corr:
        print("ВНИМАНИЕ: поправок калибровки нет — работаем по заводской шкале")

    edges = np.arange(0.0, E_MAX + E_STEP, E_STEP)
    centres = 0.5 * (edges[:-1] + edges[1:])

    y = rebin_to_grid(np.asarray(spec.counts, float), list(spec.energy_cal),
                      corr.get("sample"), edges)
    print("образец: %d отсчётов в сетке из %d (в файле %d)"
          % (int(y.sum()), len(centres), int(np.sum(spec.counts))))

    if bg is not None:
        ybg = rebin_to_grid(np.asarray(bg.counts, float),
                            list(bg.energy_cal), corr.get("background"), edges)
        ybg = ybg * (t_smp / float(bg.real_time))     # к времени образца
        print("фон: %.0f отсчётов, приведён к %.0f с (было %.0f с)"
              % (ybg.sum(), t_smp, bg.real_time))
    else:
        ybg = np.zeros_like(y)

    # --- шаблоны --------------------------------------------------------------
    templates, names, colours, stamps = [], [], [], set()
    head0 = {}
    for key, label, colour in ORDER:
        p = os.path.join(tdir, "%s.csv" % key)
        if not os.path.exists(p):
            continue
        head, e, c = read_template(p)
        head0 = head0 or head
        n = float(head["N_primaries"])
        if c.sum() == 0:
            print("  %s: ни одного события — компонента пропущена" % label)
            continue
        stamps.add(head.get("src_sha1", "?"))
        templates.append((key, label, colour, e, c / n, n))
        names.append(label)
        colours.append(colour)
    if len(stamps) > 1:
        raise SystemExit("шаблоны разных ревизий: %s" % ", ".join(stamps))
    print("шаблонов %d, штамп %s" % (len(templates), stamps.pop() if stamps
                                     else "?"))

    m = (centres >= E_FIT[0]) & (centres <= E_FIT[1])
    yy = y[m]
    wgt = 1.0 / np.sqrt(np.maximum(yy, 1.0))          # пуассоновские веса

    # --- сборка компонент подгонки -------------------------------------------
    # Свободная подгонка ПОНУКЛИДНО не годится: шаблоны Th-232, Th-228, Ra-224
    # набраны единицами событий (2, 10 и 100 при 200 тыс. розыгрышей — у этих
    # звеньев почти нет гамма-выхода), и подгонка тянет их множители в сотни
    # килобеккерелей, подгоняя ими шум. Поэтому ряд собирается в ДВЕ подцепочки,
    # внутри которых равновесие обеспечено самими периодами полураспада:
    #
    #   A1 — верхняя часть: Th-232 -> Ra-228 -> Ac-228 (сигнал даёт Ac-228);
    #   A2 — нижняя часть: Th-228 -> Ra-224 -> Rn-220 -> Po-216 -> Pb-212 ->
    #        Bi-212 -> (Tl-208 | Po-212). Все периоды от 0,15 с до 3,6 суток,
    #        равновесие достигается за недели.
    #
    # Ветвление Bi-212 -> Tl-208 берётся из библиотеки МАГАТЭ, а не пишется
    # числом. Отношение A1/A2 остаётся СВОБОДНЫМ: именно оно показывает,
    # нарушено ли равновесие между Ra-228 и Th-228.
    br_tl = branch_to_tl208()
    print("ветвление Bi-212 -> Tl-208: %.2f %% (МАГАТЭ)" % (100.0 * br_tl))
    GROUP = {
        "A1 (Ra-228 -> Ac-228)": {"Th232": 1.0, "Ra228": 1.0, "Ac228": 1.0},
        "A2 (Th-228 -> Tl-208)": {"Th228": 1.0, "Ra224": 1.0, "Rn220": 1.0,
                                  "Po216": 1.0, "Pb212": 1.0, "Bi212": 1.0,
                                  "Tl208": br_tl, "Po212": 1.0 - br_tl},
    }
    tmap = {k: (e, c) for k, _, _, e, c, _ in templates}

    def build(f0, f1):
        cols = []
        for gname, members in GROUP.items():
            acc = np.zeros(len(centres))
            for k, wgt_k in members.items():
                if k in tmap:
                    e, c = tmap[k]
                    acc += wgt_k * broaden(e, c, centres, f0, f1)
            cols.append(acc[m] * t_smp)
        return np.column_stack(cols)

    best = None
    for f1 in np.arange(1.0, 5.01, 0.2):              # ПШПВ² = f0 + f1·E
        for f0 in np.arange(-600.0, 601.0, 100.0):
            if f0 + f1 * E_FIT[0] <= 4.0:
                continue
            A = build(f0, f1)
            # ФОН НЕ ПОДГОНЯЕТСЯ: он измерен тем же прибором и приведён к
            # времени образца, его вклад известен. Свободный множитель фона в
            # пробной подгонке уходил на 4,8 — подгонка затыкала фоном нехватку
            # континуума, то есть лечила симптом.
            resid = yy - ybg[m]
            x, _ = nnls_fit(A, resid, wgt)
            r = (A @ x + ybg[m] - yy) * wgt
            chi2 = float((r ** 2).sum()) / max(1, len(yy) - len(x))
            if best is None or chi2 < best[0]:
                best = (chi2, f0, f1, x, A)

    chi2, f0, f1, x, A = best
    print("\nПШПВ² = %.0f + %.2f·E  (ПШПВ 662 кэВ = %.1f кэВ), хи²/n = %.1f"
          % (f0, f1, math.sqrt(max(f0 + f1 * 661.657, 1.0)), chi2))

    model = A @ x + ybg[m]
    print("\n--- активности по разложению ---")
    rows = []
    for i, gname in enumerate(GROUP):
        print("  %-24s %10.0f Бк" % (gname, x[i]))
        rows.append((gname, x[i]))
    if x[1] > 0:
        print("  отношение A1/A2 = %.3f (равновесие ряда -> 1,000)"
              % (x[0] / x[1]))

    # --- удельная активность и сверка с номиналом этикетки ------------------
    # Масса пачки берётся ИЗ ШАПКИ ШАБЛОНА, а не пересчитывается здесь: считает
    # её геометрия по построенным телам, и второй счёт в другом месте — это
    # ровно тот случай, когда числа расходятся молча.
    mass_g = float(head0.get("wt20_mass_g", "0").split()[0])
    a_meas = 0.5 * (x[0] + x[1])            # ряд в целом, если равновесен
    print("\n--- удельная активность ---")
    if mass_g > 0:
        print("  масса пачки %.1f г (из шапки шаблона)" % mass_g)
        print("  A(ряд) / масса = %.0f Бк/кг сплава"
              % (1000.0 * a_meas / mass_g))
        # Номинал этикетки: 2 % масс. ThO2, доля тория в ThO2 0,878809,
        # удельная активность Th-232 4072 Бк/г (вычислена из T1/2 МАГАТЭ).
        th_g = mass_g * 0.02 * 0.878809
        a_nom = th_g * SPEC_ACT_TH232
        print("  ПО ЭТИКЕТКЕ (2 %% ThO2): тория %.3f г -> %.0f Бк"
              % (th_g, a_nom))
        print("  измерено / номинал = %.3f" % (a_meas / a_nom if a_nom else 0))
        print("  эквивалентное содержание ThO2 = %.2f %% масс."
              % (2.0 * a_meas / a_nom if a_nom else 0))

    with io.open(os.path.join(outdir, "unfold_activities.csv"), "w",
                 encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["компонента", "активность_Бк"])
        for gname, a in rows:
            w.writerow([gname, "%.5g" % a])
        w.writerow(["ПШПВ_662_кэВ", "%.3g" % math.sqrt(f0 + f1 * 661.657)])
        w.writerow(["хи2_на_канал", "%.4g" % chi2])

    # Вклад ОТДЕЛЬНЫХ нуклидов при найденных активностях — для рисунка.
    names = []
    parts = []
    for gi, (gname, members) in enumerate(GROUP.items()):
        for k, wk in members.items():
            if k not in tmap:
                continue
            e, c = tmap[k]
            v = wk * broaden(e, c, centres, f0, f1)[m] * t_smp * x[gi]
            if v.sum() <= 0:
                continue
            names.append(dict(ORDER_MAP).get(k, k))
            parts.append(v)
    names.append("фон")
    parts.append(ybg[m])

    # Пишется через io.open с явным UTF-8: np.savetxt кодирует шапку системной
    # кодировкой, и кириллические имена компонент выходили нечитаемыми.
    tab = np.column_stack([centres[m], yy, model, *parts])
    with io.open(os.path.join(outdir, "unfold_spectrum.csv"), "w",
                 encoding="utf-8", newline="") as f:
        f.write("E_keV,измерено,модель," + ",".join(names) + "\n")
        for row in tab:
            f.write(",".join("%.6g" % v for v in row) + "\n")
    print("\nзаписано: %s" % os.path.join(outdir, "unfold_spectrum.csv"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
