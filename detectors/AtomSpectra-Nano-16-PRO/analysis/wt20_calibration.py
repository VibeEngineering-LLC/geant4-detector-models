# -*- coding: utf-8 -*-
"""Проверка и уточнение энергетической калибровки — ОТДЕЛЬНО для образца и фона.

Зачем отдельно. В файле замера лежат ДВА спектра: сам образец (01.06.2024) и
встроенный фон, набранный раньше и дольше. Калибровка в файле записана одна на
оба, но снимались они в разное время и коэффициенты записаны разные — значит
проверять надо каждый по своим линиям. Вычитание фона при несведённых шкалах
сдвигает пики и рождает ложные структуры на разностном спектре.

Опорные линии берутся не по памяти, а из библиотеки МАГАТЭ, выкачанной в
`reference/nuclide-lines/*.csv` (поле `energy`, `intensity`).

Штатная семилинейная проверка SpectraVibe для ОБРАЗЦА не годится: она построена
на линиях ЕРН (K-40, Bi-214, Pb-214), а в спектре чистого тория их нет — из семи
она нашла пять, и две из них перепутала (351,93 приписано пику 341,5, то есть
Ac-228 338,32). Для образца берутся линии ряда тория, для фона — ЕРН.

    python analysis/wt20_calibration.py <файл.xml> [каталог вывода]
"""
import csv
import io
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.normpath(os.path.join(_HERE, "..", "reference", "nuclide-lines"))

_ROOT = os.environ.get("SPECTRAVIBE_ROOT")
if not _ROOT:
    raise SystemExit("не задана переменная окружения SPECTRAVIBE_ROOT")
sys.path.insert(0, os.path.join(_ROOT, "scripts"))

from gamma.io.atomspectra_xml import read_atomspectra_xml        # noqa: E402
from gamma.peaks.centroid_gost import (                          # noqa: E402
    gost_centroid_graphoanalytic, gost_centroid_weighted_mean,
    gost_select_pedestal_method)
from gamma.peaks.search import estimate_fwhm_at_peak             # noqa: E402

# Собственное разрешение прибора: ПШПВ(662) = 41,60 кэВ по записи Cs-137 (тот
# же закон, что во всех скриптах этого прибора). Нужен только как НАЧАЛЬНОЕ
# приближение для поиска; фактическая ширина меряется по самому пику.
FWHM_662 = 41.60


def fwhm_keV(e):
    return FWHM_662 * math.sqrt(max(e, 1.0) / 661.657)


# --- опорные линии -----------------------------------------------------------
# (нуклид, энергия) — энергия сверяется с библиотекой МАГАТЭ при запуске, чтобы
# набранный от руки список не разошёлся с файлами reference/nuclide-lines.
# Опорной берётся только СИЛЬНАЯ и ОДИНОЧНАЯ линия. Пробный расширенный набор
# показал, за что это правило: Ac-228 338,32 и Bi-212 1620,5 слабы и сидят в
# мультиплетах — центроида уезжала на 15 кэВ, а два ГОСТ-метода расходились
# в 7 ПШПВ. В фоне 583,19 и 609,32 разнесены на 26 кэВ при ПШПВ 40 и в один
# канал не разделяются вовсе: обе «нашлись» в канале 1533.
ANCHORS_SAMPLE = [("212pb", 238.63), ("208tl", 583.19), ("228ac", 911.20),
                  ("208tl", 2614.51)]
ANCHORS_BG = [("214pb", 351.93), ("214bi", 609.32), ("40k", 1460.82),
              ("208tl", 2614.51)]


def lib_energy(nuclide, want, tol=0.5):
    """Энергия линии из библиотеки МАГАТЭ. Отказ, если её там нет."""
    path = os.path.join(LIB, "%s_gammas.csv" % nuclide)
    best = None
    for r in csv.DictReader(io.open(path, encoding="utf-8")):
        try:
            e = float(r["energy"])
            i = float(r["intensity"])
        except (TypeError, ValueError):
            continue
        if abs(e - want) <= tol and (best is None or i > best[1]):
            best = (e, i)
    if best is None:
        raise SystemExit("в %s нет линии %.2f кэВ" % (path, want))
    return best


def poly(coefs, x):
    return sum(c * x ** k for k, c in enumerate(coefs))


def ch_of_energy(coefs, e, n_ch):
    """Обратное преобразование E->канал перебором по монотонному участку."""
    lo, hi = 0.0, float(n_ch - 1)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if poly(coefs, mid) < e:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def measure_line(counts, coefs, e_lib, search_fwhm=1.5):
    """Центроида линии в КАНАЛАХ и её энергия по действующей калибровке.

    ROI берётся УЗКИЙ (±1 ПШПВ). Умолчание ГОСТ-модуля — ±2,5 ПШПВ, и на
    жёсткой линии 2614,5 кэВ такое окно захватывает комптоновскую подложку
    вместе с одиночным вылетом: центроида уезжала на 97 кэВ вниз, то есть
    ДАЛЬШЕ окна поиска, в котором пик искали. Признак ошибки — центроида,
    вышедшая за пределы своего же ROI; теперь это проверяется явно.
    """
    n = len(counts)
    ch0 = ch_of_energy(coefs, e_lib, n)
    # ширина в каналах: ПШПВ по энергии, делённая на dE/dch в этой точке
    dE = poly(coefs, ch0 + 0.5) - poly(coefs, ch0 - 0.5)
    if dE <= 0:
        return None
    fw_ch = fwhm_keV(e_lib) / dE
    lo = max(0, int(ch0 - search_fwhm * fw_ch))
    hi = min(n - 1, int(ch0 + search_fwhm * fw_ch))
    if hi - lo < 4:
        return None
    top = lo + int(np.argmax(counts[lo:hi + 1]))
    fw_meas = estimate_fwhm_at_peak(counts, top, fw_ch)
    if not fw_meas or not math.isfinite(fw_meas) or fw_meas <= 1:
        fw_meas = fw_ch
    # ПШПВ пика не может быть в разы уже приборной: такое значение означает,
    # что «пик» — одиночный выброс, а не линия.
    if fw_meas < 0.4 * fw_ch or fw_meas > 3.0 * fw_ch:
        fw_meas = fw_ch
    ped = gost_select_pedestal_method(counts, top, fw_meas, roi_half_fwhm=1.0)
    net = np.asarray(ped.counts_net, dtype=float)
    if net.size < 5 or net.sum() <= 0:
        return None
    # ЦЕНТРОИДА — ВЗВЕШЕННЫМ СРЕДНИМ (ГОСТ 26874-86 §3.3.2). Графоаналитический
    # метод того же ГОСТа здесь непригоден: он подгоняет параболу к ln N, и на
    # линии 2614,5 кэВ, сидящей на крутом комптоновском спаде, вершина параболы
    # уезжала ЗА ПРЕДЕЛЫ своего же окна — 6084 при окне 6146…6539. Расхождение
    # двух методов проверяется явно и печатается.
    cen = gost_centroid_weighted_mean(net, channel_offset=ped.roi_lo)
    ch = float(cen.n_c)
    if not (ped.roi_lo <= ch <= ped.roi_hi):
        return None
    try:
        alt = gost_centroid_graphoanalytic(net, channel_offset=ped.roi_lo)
        disagree = abs(float(alt.n_c) - ch) / max(fw_meas, 1.0)
    except ValueError:
        # Логарифмическая подгонка требует хотя бы двух пар точек выше
        # полувысоты; на слабой линии их может не быть. Это не отказ замера,
        # а отказ ПЕРЕКРЁСТНОЙ проверки — так и отмечается.
        disagree = float("nan")
    return dict(ch=ch, e_obs=poly(coefs, ch), fwhm_ch=float(fw_meas),
                fwhm_keV=fw_meas * dE, area=float(net.sum()), top=top,
                roi=(ped.roi_lo, ped.roi_hi), disagree=disagree)


def fit_cal(points, degree):
    """Поправка В ПРОСТРАНСТВЕ ЭНЕРГИЙ: E_ист = f(E_по действующей шкале).

    Подгонять заново полином канал->энергия неправильно: заводская шкала
    четвёртой степени уже несёт нелинейность тракта, снятую по многим точкам,
    а у нас опорных линий три-четыре. Дрейф между замерами — это сдвиг и
    масштаб, поэтому поправка берётся линейной (или квадратичной, если точек
    хватает) ПОВЕРХ заводской шкалы. Проверено на фоне: свежая линейная шкала
    в каналах давала невязку 9,8 кэВ, поправка в энергиях — на порядок меньше.

    points: [(E_по действующей шкале, E_библиотечная)].

    Степень ограничена так, чтобы осталась хотя бы одна степень свободы: на
    трёх точках парабола проходит через них ТОЧНО, невязка выходит нулевой и
    «проверка» перестаёт что-либо проверять.
    """
    x = np.array([p[0] for p in points], float)
    y = np.array([p[1] for p in points], float)
    deg = max(1, min(degree, len(points) - 2))
    c = np.polyfit(x, y, deg)[::-1]
    return list(map(float, c)), deg


def report(name, counts, coefs, anchors, out_rows):
    print("\n=== %s ===" % name)
    print("действующая калибровка:", ", ".join("%.6g" % c for c in coefs))
    pts, res = [], []
    for nuc, want in anchors:
        e_lib, inten = lib_energy(nuc, want)
        m = measure_line(counts, coefs, e_lib)
        if not m:
            print("  %-8s %8.2f кэВ — линия не выделена" % (nuc, e_lib))
            continue
        d = m["e_obs"] - e_lib
        frac = d / fwhm_keV(e_lib)
        print("  %-8s %8.2f кэВ (I=%5.2f %%): канал %8.2f, найдено %8.2f, "
              "Δ = %+6.2f кэВ (%+.3f ПШПВ), ПШПВ %5.1f кэВ, площадь %d%s"
              % (nuc, e_lib, inten, m["ch"], m["e_obs"], d, frac,
                 m["fwhm_keV"], int(m["area"]),
                 "" if not (m["disagree"] >= 0.2) else
                 "  ВНИМАНИЕ: методы центроиды расходятся на %.2f ПШПВ"
                 % m["disagree"]))
        pts.append((m["e_obs"], e_lib))
        res.append(d)
        out_rows.append((name, nuc, e_lib, inten, m["ch"], m["e_obs"], d,
                         m["fwhm_keV"], int(m["area"])))
    if len(pts) < 3:
        print("  опорных точек мало — калибровка не уточняется")
        return coefs, res, None
    new, deg = fit_cal(pts, 2)
    res2 = [poly(new, x) - e for x, e in pts]
    print("  до поправки:   макс |Δ| = %.2f кэВ, СКО %.2f"
          % (max(abs(x) for x in res), float(np.std(res))))
    print("  поправка E_ист = %s (степень %d, свободных степеней %d)"
          % (" + ".join("%.8g·E^%d" % (c, k) if k else "%.8g" % c
                        for k, c in enumerate(new)), deg,
             len(pts) - deg - 1))
    print("  после:         макс |Δ| = %.2f кэВ, СКО %.2f"
          % (max(abs(x) for x in res2), float(np.std(res2))))
    return new, res, new


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    src = sys.argv[1]
    outdir = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(src)
    spec = read_atomspectra_xml(src)
    if isinstance(spec, (list, tuple)):
        spec = spec[0]
    bg = getattr(spec, "background_embedded", None)

    rows = []
    smp_counts = np.asarray(spec.counts, float)
    new_smp, _, _ = report("ОБРАЗЕЦ (%s, %.0f с)"
                           % (spec.sample_id or "без имени", spec.real_time),
                           smp_counts, list(spec.energy_cal), ANCHORS_SAMPLE,
                           rows)
    new_bg = None
    if bg is not None:
        bg_counts = np.asarray(bg.counts, float)
        new_bg, _, _ = report("ФОН (встроенный, %.0f с)" % bg.real_time,
                              bg_counts, list(bg.energy_cal), ANCHORS_BG, rows)
        print("\nкалибровка образца и фона в файле %s"
              % ("СОВПАДАЕТ" if list(spec.energy_cal) == list(bg.energy_cal)
                 else "РАЗНАЯ — вычитать без сведения шкал нельзя"))

    os.makedirs(outdir, exist_ok=True)
    p = os.path.join(outdir, "calibration_check.csv")
    with io.open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["спектр", "нуклид", "E_библ_кэВ", "I_%", "канал",
                    "E_найдено_кэВ", "невязка_кэВ", "ПШПВ_кэВ", "площадь"])
        for r in rows:
            w.writerow(["%.4g" % x if isinstance(x, float) else x for x in r])
    print("\nзаписано: %s" % p)

    # Уточнённые калибровки — отдельным файлом: их читает разложение спектра.
    p2 = os.path.join(outdir, "calibration_fitted.csv")
    with io.open(p2, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["# поправка в ЭНЕРГИЯХ поверх заводской шкалы: "
                    "E_ист = a0 + a1*E + a2*E^2, где E — энергия по шкале, "
                    "записанной в самом файле замера"])
        w.writerow(["спектр", "a0", "a1", "a2"])
        for nm, c in (("sample", new_smp), ("background", new_bg)):
            if c:
                w.writerow([nm] + ["%.10g" % x for x in
                                   (list(c) + [0, 0, 0])[:3]])
    print("записано: %s" % p2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
