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
# Опорные линии ОБРАЗЦА берутся из каталога конструктора ROI
# (reference/roi/wizard_lines_iaea.xml) автоматическим отбором, а не набираются
# руками. Критерии отбора:
#   сила     — выход не меньше MIN_YIELD;
#   одиночность — в пределах 1,5 ПШПВ нет чужой линии с выходом больше
#               MIN_NEIGHBOUR доли собственного.
# Отбор снимает и находку S5 внешнего аудита: линия 911,20 держит рядом 964,77
# и 968,97 (суммарный выход сравним с собственным) и в опорные не проходит —
# прежде её пулл впечатывался в поправку шкалы (~6 кэВ в полосе 900-1500).
# Слабый сосед допустИм: у 238,63 сосед 240,99 несёт лишь 9 % собственного
# выхода, и его пулл меньше погрешности центроиды.
#
# Для ФОНА каталог образца не годится (там ЕРН, не ряд тория) — набор прежний,
# проверенный: 351,93 / 609,32 / 1460,82 / 2614,51.
MIN_YIELD = 15.0        # %, порог силы опорной линии
MIN_NEIGHBOUR = 0.20    # доля собственного выхода, с которой сосед мешает

_NUC2LIB = {"Pb212": "212pb", "Bi212": "212bi", "Tl208": "208tl",
            "Ac228": "228ac", "Ra224": "224ra", "Rn220": "220rn",
            "Th228": "228th", "Th232": "232th"}


def anchors_from_catalog():
    """Опорные линии образца — отбором из каталога конструктора ROI.

    Возвращает [(имя_библиотеки | "XW", E, ROI | None)]. Гамма-якоря проходят
    автоотбор из каталога (сила + одиночность); мягкий якорь — суммарный пик
    K-серии вольфрама — строится из МК-шаблонов, см. xray_w_anchor_from_mc.
    """
    sys.path.insert(0, _HERE)
    import roi_lines as R
    cat = R.parse_xml(R.DEFAULT_XML)
    lines = [r for r in cat if r["key"] in _NUC2LIB]
    out = []
    for r in lines:
        if r["yield_pct"] < MIN_YIELD:
            continue
        fw = fwhm_keV(r["E"])
        clean = True
        for q in lines:
            if q is r:
                continue
            if (abs(q["E"] - r["E"]) <= 1.5 * fw
                    and q["yield_pct"] > MIN_NEIGHBOUR * r["yield_pct"]):
                clean = False
                break
        if clean:
            out.append((_NUC2LIB[r["key"]], r["E"], None))
    if len(out) < 3:
        raise SystemExit("отбор опорных линий из каталога дал меньше трёх: %r"
                         % out)
    # Мягкий якорь — ХРИ ВОЛЬФРАМА, суммарный пик K-серии (директива оператора
    # 07.08.2026: «этот пик должен калиброваться по ХРИ вольфрама», «там
    # суммарный пик, нужно складывать пики в области с весами»).
    #
    # Опора берётся ИЗ МОНТЕ-КАРЛО-ШАБЛОНА, а не из каталога и не из
    # библиотеки линий. Причина в постановке: W не радиоактивен, его K-серия —
    # флуоресценция, возбуждаемая в самом источнике, и в схемах распада её
    # нет; при этом расчёт Geant4 её СЧИТАЕТ ПЕРЕНОСОМ, вместе с весами линий,
    # самопоглощением в металле над K-краем 69,5 кэВ и эффективностью
    # регистрации. Центроида профиля шаблона в окне и есть та самая
    # «взвешенная сумма пиков области», к которой привязывается шкала.
    #
    # Прежний якорь XKD (комплекс ХРИ дочерних, опора 82,05 кэВ по
    # библиотечным выходам) СНЯТ: он давал ложное подтверждение. Алгоритм
    # искал пик около 82,05, находил на 82,12 и печатал невязку 0,08 кэВ —
    # «шкала сошлась», тогда как поправка в этой области ошибалась почти на
    # 6 кэВ (проверка 07.08.2026 по разложению мягкого края: показанию
    # прибора 82,8 отвечает Pb Kα1 77,11). Якорь подтверждал сходимость,
    # которой нет, потому что его опора — средневзвешенная по БИБЛИОТЕЧНЫМ
    # выходам, без переноса и эффективности, — не совпадает с центроидой
    # того, что прибор в этом окне действительно видит.
    e_ref, win = xray_anchor_from_mc()
    if e_ref is not None:
        out.insert(0, ("XRAY", e_ref, win))
    return out


# Окно комплекса характеристического рентгена, кэВ. Накрывает ВЕСЬ комплекс
# целиком, а не одну серию: K-серию вольфрама (Kα 58,0-59,3, Kβ 67,0-69,1),
# K-серию дочерних торона — таллия, свинца, висмута (Kα 72,8-77,1, Kβ 85-93) —
# и K-серию актиния (Kα 89,9-93,4).
#
# Резать окно по одной серии НЕЛЬЗЯ (директива оператора 07.08.2026: «это
# нельзя отдельно от прочих ХРИ считать»). При ПШПВ около 12 кэВ в этой
# области серии перекрываются: линия W Kβ 69,1 и линия Tl Kα2 72,8 разнесены
# на 3,7 кэВ, то есть втрое меньше ширины аппаратной линии. Граница между ними
# существует только на бумаге; всякий рез внутри комплекса отсекает часть
# одной серии и оставляет часть другой, и центроида обрезка не отвечает уже
# ничему. Комплекс берётся целиком, а веса внутри него даёт расчёт.
XRAY_WINDOW = (50.0, 100.0)


def xray_anchor_from_mc(tdir=None, model_csv=None):
    """Опора мягкого якоря: центроида комплекса ХРИ по расчёту Geant4.

    Весы для сложения линий области берутся из расчёта, а не назначаются: они
    несут и относительные интенсивности линий внутри каждой серии, и
    самопоглощение в вольфраме над K-краем 69,5 кэВ, и эффективность
    регистрации, и — главное — соотношение между сериями разных элементов.

    Есть два источника весов, и они образуют итерацию.

    1. СУММАРНАЯ МОДЕЛЬ (`model_csv`, колонка «модель» в unfold_spectrum.csv)
       — если разложение уже считалось. Тогда веса серий отвечают найденным
       активностям, то есть тому составу, который в источнике действительно
       есть. Это правильный источник, и он используется, когда доступен.

    2. РАВНЫЕ веса шаблонов — первый проход, когда активностей ещё нет.
       Приближение грубое: центроида комплекса зависит от соотношения серий
       (68,6 кэВ у Pb-212 против 73,8 у поверхностного Pb-212), а оно как раз
       и определяется разложением.

    Отсюда порядок: калибровка на равных весах -> разложение -> калибровка на
    суммарной модели -> разложение. Два круга, дальше опора не движется.

    Возвращает (E_центроиды, окно) либо (None, None), если данных нет.
    """
    lo, hi = XRAY_WINDOW
    # --- источник 1: суммарная модель предыдущего разложения ---------------
    if model_csv is None:
        model_csv = os.environ.get("WT20_MODEL_SPECTRUM", "")
    if model_csv and os.path.exists(model_csv):
        e, v = [], []
        rd = csv.reader(io.open(model_csv, encoding="utf-8"))
        hdr = next(rd, None)
        if hdr and "модель" in hdr:
            j = hdr.index("модель")
            for row in rd:
                if not row:
                    continue
                try:
                    e.append(float(row[0]))
                    v.append(float(row[j]))
                except (ValueError, IndexError):
                    continue
        if e:
            e, v = np.array(e), np.array(v)
            m = (e >= lo) & (e < hi)
            if v[m].sum() > 0:
                return float((e[m] * v[m]).sum() / v[m].sum()), XRAY_WINDOW
    # --- источник 2: шаблоны с равными весами (первый проход) --------------
    if tdir is None:
        tdir = os.environ.get("WT20_TEMPLATES", "")
    if not tdir or not os.path.isdir(tdir):
        return None, None
    acc_e, acc_c = None, None
    for fn in sorted(os.listdir(tdir)):
        if not fn.endswith(".csv") or "_" in fn:
            continue                      # только основные спектры нуклидов
        e, c = [], []
        for ln in io.open(os.path.join(tdir, fn), encoding="utf-8"):
            if ln.startswith("#") or ln.startswith("E_keV") or not ln.strip():
                continue
            a, b = ln.split(",")
            e.append(float(a))
            c.append(float(b))
        if not e:
            continue
        e, c = np.array(e), np.array(c)
        m = (e >= lo) & (e < hi)
        if c[m].sum() <= 0:
            continue
        cn = c[m] / c[m].sum()
        if acc_e is None:
            acc_e, acc_c = e[m], cn
        elif len(acc_e) == len(cn):
            acc_c = acc_c + cn
    if acc_e is None:
        return None, None
    return float((acc_e * acc_c).sum() / acc_c.sum()), XRAY_WINDOW


ANCHORS_BG = [("214pb", 351.93, None), ("214bi", 609.32, None),
              ("40k", 1460.82, None), ("208tl", 2614.51, None)]


def sum_anchor():
    """Якорь по пикам суммирования каскада Tl-208 (оператор, 07.08.2026).

    В спектре образца у 3200 кэВ лежит бугор истинных совпадений: 2614,51 +
    583,19 = 3197,70 и 2614,51 + 510,77 = 3125,28. При ПШПВ ~90 кэВ пара не
    разрешается, поэтому опора — средневзвешенная сумма по выходам вторых
    квантов каскада (85,0 и 22,6 %; вероятность совпадения принята
    одинаковой — полные эффективности на 511 и 583 кэВ близки). Энергии и
    выходы — из библиотеки МАГАТЭ, не числами в коде. Третий сумм-пик
    2614,51 + 277,37 = 2891,88 слаб (6,6 %) и в якорь не входит.
    """
    e26, _ = lib_energy("208tl", 2614.51)
    e58, i58 = lib_energy("208tl", 583.19)
    e51, i51 = lib_energy("208tl", 510.77)
    e_ref = (((e26 + e58) * i58 + (e26 + e51) * i51) / (i58 + i51))
    return ("SUM", e_ref, (3090.0, 3290.0))


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


def measure_line(counts, coefs, e_lib, search_fwhm=1.5, roi_keV=None,
                 broad=False):
    """Центроида линии в КАНАЛАХ и её энергия по действующей калибровке.

    ROI берётся УЗКИЙ (±1 ПШПВ). Умолчание ГОСТ-модуля — ±2,5 ПШПВ, и на
    жёсткой линии 2614,5 кэВ такое окно захватывает комптоновскую подложку
    вместе с одиночным вылетом: центроида уезжала на 97 кэВ вниз, то есть
    ДАЛЬШЕ окна поиска, в котором пик искали. Признак ошибки — центроида,
    вышедшая за пределы своего же ROI; теперь это проверяется явно.

    roi_keV=(lo, hi) заменяет окно поиска окном КАТАЛОГА — для якоря ХРИ
    вольфрама, где стандартное ±1,5 ПШПВ шире, чем расстояние до соседнего
    горба ХРИ Pb/Bi.
    """
    n = len(counts)
    ch0 = ch_of_energy(coefs, e_lib, n)
    # ширина в каналах: ПШПВ по энергии, делённая на dE/dch в этой точке
    dE = poly(coefs, ch0 + 0.5) - poly(coefs, ch0 - 0.5)
    if dE <= 0:
        return None
    fw_ch = fwhm_keV(e_lib) / dE
    if roi_keV is not None:
        lo = int(ch_of_energy(coefs, roi_keV[0], n))
        hi = int(math.ceil(ch_of_energy(coefs, roi_keV[1], n)))
        lo, hi = max(0, lo), min(n - 1, hi)
    else:
        lo = max(0, int(ch0 - search_fwhm * fw_ch))
        hi = min(n - 1, int(ch0 + search_fwhm * fw_ch))
    if broad:
        # Широкая слабая структура (сумм-бугор ~200 отсчётов на 500 каналов):
        # ГОСТ-центроида узкого пика неприменима — argmax на сырых каналах
        # ловит шум. Берётся взвешенная центроида ВСЕГО окна с линейной
        # подложкой по краевым полосам (15 % ширины окна с каждого края).
        seg = counts[lo:hi + 1].astype(float)
        w15 = max(3, (hi - lo) // 7)
        xl, xr = np.arange(lo, lo + w15), np.arange(hi - w15 + 1, hi + 1)
        yl, yr = seg[:w15].mean(), seg[-w15:].mean()
        cl, cr = xl.mean(), xr.mean()
        ch_ax = np.arange(lo, hi + 1)
        base = yl + (yr - yl) * (ch_ax - cl) / max(cr - cl, 1.0)
        # Отрицательное нетто отбрасывается, а не суммируется со знаком.
        # Линейная подложка по краевым полосам — приближение, и на вогнутом
        # участке (в окне W-серии счёт сперва спадает до 53 кэВ, затем растёт
        # к склону горба дочерних) хорда проходит ВЫШЕ середины окна. Со
        # знаком такие каналы тянут центроиду в сторону и обнуляют сумму —
        # структура объявлялась невыделенной там, где она есть. Вес канала
        # ниже подложки физически равен нулю: он не содержит превышения.
        net = np.clip(seg - base, 0.0, None)
        s = float(net.sum())
        if s <= 0:
            return None
        ch = float((ch_ax * net).sum() / s)
        if not (lo + w15 <= ch <= hi - w15):
            return None
        return dict(ch=ch, e_obs=poly(coefs, ch), fwhm_ch=fw_ch,
                    fwhm_keV=fw_ch * dE, area=s, top=int(round(ch)),
                    roi=(lo, hi), disagree=float("nan"))
    if hi - lo < 4:
        return None
    top = lo + int(np.argmax(counts[lo:hi + 1]))
    # Максимум на краю окна — это не пик, а склон соседней структуры:
    # так «пик ХРИ W» на 58-67 кэВ оказался склоном горба ХРИ дочерних, и
    # центроида уезжала на 76 кэВ, за пределы собственного окна.
    if top <= lo + 1 or top >= hi - 1:
        return None
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
    for nuc, want, roi in anchors:
        if nuc in ("XRAY", "SUM"):
            # опора не из библиотеки линий: XW — центроида суммарного пика
            # K-серии вольфрама по МК-шаблонам (флуоресценции нет в схемах
            # распада, но расчёт её считает переносом); SUM —
            # средневзвешенная сумм-пиков каскада Tl-208
            e_lib, inten = want, float("nan")
        else:
            e_lib, inten = lib_energy(nuc, want)
        # XW и SUM — широкие составные структуры, а не одиночные линии:
        # центроида берётся взвешенной по ВСЕМУ окну с линейной подложкой
        # (broad), иначе argmax ловит одну компоненту серии вместо суммы.
        m = measure_line(counts, coefs, e_lib, roi_keV=roi,
                         broad=(nuc in ("SUM", "XRAY")))
        if not m:
            print("  %-8s %8.2f кэВ — линия не выделена" % (nuc, e_lib))
            continue
        d = m["e_obs"] - e_lib
        frac = d / fwhm_keV(e_lib)
        inten_s = ("серия, каталог" if math.isnan(inten)
                   else "I=%5.2f %%" % inten)
        print("  %-8s %8.2f кэВ (%s): канал %8.2f, найдено %8.2f, "
              "Δ = %+6.2f кэВ (%+.3f ПШПВ), ПШПВ %5.1f кэВ, площадь %d%s"
              % (nuc, e_lib, inten_s, m["ch"], m["e_obs"], d, frac,
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

    anchors_sample = anchors_from_catalog()
    anchors_sample.append(sum_anchor())
    print("опорные линии образца из каталога конструктора ROI:")
    for nuc, e, roi in anchors_sample:
        print("  %-6s %8.2f кэВ%s" % (nuc, e,
              "  (ROI %.1f-%.1f из каталога)" % roi if roi else ""))

    rows = []
    smp_counts = np.asarray(spec.counts, float)
    new_smp, _, _ = report("ОБРАЗЕЦ (%s, %.0f с)"
                           % (spec.sample_id or "без имени", spec.real_time),
                           smp_counts, list(spec.energy_cal), anchors_sample,
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
            w.writerow(["%.6g" % x if isinstance(x, float) else x for x in r])
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
