# -*- coding: utf-8 -*-
"""Поиск пиков, опознание линий и проверка калибровки по найденному.

ЗАЧЕМ ЭТОТ МОДУЛЬ ПОЯВИЛСЯ. До него весь разбор спектра шёл по заранее
объявленному списку линий (`VLINES` в kit_recalc): взять окно вокруг
табличной энергии и снять площадь. У такого подхода есть слепая зона, и она
принципиальная — **он не может сообщить о том, чего в списке нет**. Пик в
спектре есть, а мы на него не смотрим; хуже того, мы не знаем, что не смотрим.
Замечание оператора: «ты теряешь часть пиков». Верно, и вот чем это лечится.

Порядок здесь тот же, что у SpectraLine/ЛСРМ и в пакете
gamma-spectrum-analysis (`gamma/peaks/search.py`, `gamma/identification/*`,
`gamma/calibration/anchor_recalibration.py`), потому что порядок в этой задаче
не произволен:

  1. ПОИСК — фильтром Марискотти (вторая производная гауссианы) в шкале
     КАНАЛОВ, до всякого опознания. Значимость — по Карри: отклик фильтра,
     делённый на его же пуассоновскую сигму;
  2. ОПОЗНАНИЕ — по окну ЛСРМ δE(E) вокруг библиотечной линии. Библиотека
     берётся не из справочника, а из спектра испускания ТОГО ЖЕ прогона
     Geant4 (правило проекта: числа из той же базы, что и транспорт);
  3. КАЛИБРОВКА — по невязкам опознанных линий. Правило ЛСРМ: если невязка
     хоть на одном якоре больше 0,3·ПШПВ, запасённая калибровка признаётся
     негодной и E(N) пересчитывается по якорям.

ЧТО ВЗЯТО У НИХ И ПОЧЕМУ ИМЕННО ТАК:

  полосы по ПШПВ   На NaI ПШПВ меняется от ~4 каналов на 30 кэВ до ~90 на
                   2614 — в двадцать раз. Одно ядро фильтра на весь спектр
                   либо замыливает мягкий край, либо рубит жёсткие пики на
                   несколько ложных. Спектр режется на полосы, внутри которых
                   ПШПВ меняется меньше чем в 1,2 раза, каждая свёртывается
                   своим ядром, отклики сшиваются;
  ширина ядра      1,5·ПШПВ на узких пиках (классика Марискотти) и 1,0·ПШПВ
                   при ПШПВ ≥ 15 каналов: широкое ядро сливает соседей,
                   отстоящих на одну ПШПВ, а на NaI это сплошь и рядом;
  подложка         интерполируется ЛОГ-ЛИНЕЙНО между полками слева и справа.
                   Среднее арифметическое на падающем комптоновском склоне
                   завышает подложку на cosh(d/τ)−1 ≈ 3 %, а это и есть вся
                   высота слабого пика — он молча теряется;
  отсев узких      пик с измеренной шириной меньше 0,3 от ожидаемой — не пик,
                   а шумовой выброс в один-два канала.

ЗАЧЕМ СВЕРКА С НИМИ. Если задан `SPECTRAVIBE_ROOT`, тот же спектр прогоняется
через их `mariscotti_search`, и списки сравниваются. Это не украшение: своя
реализация того же метода, дающая другой ответ, — первый признак ошибки, и
дешевле поймать её сверкой, чем спорить с числом в отчёте.

    python detectors/Gamma-1S/analysis/peaksearch.py
    python detectors/Gamma-1S/analysis/peaksearch.py --one Th-232
"""
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "common", "py"))
import paths  # noqa: E402

sys.path.insert(0, HERE)
import becqmoni as bm  # noqa: E402
import deconv as dc  # noqa: E402
import kit_recalc as kr  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))

SIGMA_THR = 3.0        # порог значимости по Карри, в сигмах
MIN_SEP = 1.0          # минимальное расстояние между пиками, в ПШПВ
BAND_RATIO = 1.2       # во сколько раз ПШПВ может меняться внутри полосы
MIN_FWHM_RATIO = 0.3   # отсев шумовых выбросов по измеренной ширине
K_ID = 0.5             # окно опознания: |E_найд − E_лин| ≤ K_ID·ПШПВ(E)
MIN_YIELD = 0.005      # порог выхода линии, чтобы считать её «сильной»
CAL_TOL = 0.3          # доля ПШПВ, выше которой калибровку надо пересчитать


# --- 1. Поиск: фильтр Марискотти -------------------------------------------

def _kernel(fwhm_ch, half):
    """Ядро второй производной гауссианы, нулевой суммы, единичной нормы."""
    sig = fwhm_ch / 2.3548
    x = np.arange(-half, half + 1, dtype=float)
    g = ((x / sig) ** 2 - 1.0) * np.exp(-0.5 * (x / sig) ** 2)
    g -= g.mean()                      # истинная вторая производная даёт нуль
    n = math.sqrt(float((g * g).sum()))
    return g / n if n > 0 else g


def _bands(fw, max_ratio=BAND_RATIO):
    """Разбиение на полосы, внутри которых ПШПВ меняется меньше max_ratio."""
    out, start, lo, hi = [], 0, float(fw[0]), float(fw[0])
    for i in range(1, len(fw)):
        v = float(fw[i])
        nlo, nhi = min(lo, v), max(hi, v)
        if nhi / nlo > max_ratio:
            out.append((start, i, math.sqrt(lo * hi)))
            start, lo, hi = i, v, v
        else:
            lo, hi = nlo, nhi
    out.append((start, len(fw), math.sqrt(lo * hi)))
    return out


def significance(counts, fw):
    """Отклик фильтра в сигмах Карри по всем каналам."""
    n = len(counts)
    sig = np.zeros(n)
    for a, b, f in _bands(fw):
        # Узкое ядро на широких пиках: широкое сливает соседей в одну ПШПВ.
        half = max(int(math.ceil((1.0 if f >= 15.0 else 1.5) * f)), 3)
        k = _kernel(f, half)
        lo, hi = max(0, a - half), min(n, b + half)
        seg = counts[lo:hi]
        resp = -np.convolve(seg, k, mode="same")
        var = np.convolve(np.maximum(seg, 1.0), k ** 2, mode="same")
        s = resp / np.sqrt(np.maximum(var, 1e-12))
        sig[a:b] = s[a - lo:a - lo + (b - a)]
    return sig


def _maxima(sig, mask):
    out, i, n = [], 1, len(sig)
    while i < n - 1:
        if not mask[i] or sig[i] <= sig[i - 1]:
            i += 1
            continue
        j = i
        while j < n - 1 and sig[j + 1] == sig[i]:
            j += 1
        if j < n - 1 and sig[j] > sig[j + 1]:
            out.append(i)
        i = j + 1
    return out


def _thin(cand, sig, fw, factor=MIN_SEP):
    """Прореживание: сильные пики выживают, слабые ближе factor·ПШПВ гибнут."""
    keep = []
    for ch in sorted(cand, key=lambda c: -sig[c]):
        if all(abs(ch - a) >= factor * max(fw[ch], fw[a]) for a in keep):
            keep.append(ch)
    return sorted(keep)


def _baseline(counts, ch, fw):
    """Подложка под пиком — ЛОГ-линейная интерполяция между полками.

    Среднее арифметическое двух полок на падающем экспоненциальном континууме
    завышает подложку; на комптоновском склоне это те самые проценты, из
    которых состоит слабый пик.
    """
    n = len(counts)
    ll = max(0, ch - int(round(3 * fw)))
    lh = max(ll + 1, ch - int(round(2 * fw)))
    rl = min(n - 1, ch + int(round(2 * fw)))
    rh = min(n, ch + int(round(3 * fw)))
    bl = counts[ll:lh].mean() if lh > ll else 0.0
    br = counts[rl:rh].mean() if rh > rl else 0.0
    cl, cr = 0.5 * (ll + lh - 1), 0.5 * (rl + rh - 1)
    if cr <= cl:
        return 0.5 * (bl + br)
    t = (ch - cl) / (cr - cl)
    if bl > 0 and br > 0:
        return float(math.exp(math.log(bl) + (math.log(br) - math.log(bl)) * t))
    return float(bl + (br - bl) * t)


def search(counts, fw, thr=SIGMA_THR):
    """Найденные пики: [(канал, чистая высота, ПШПВ, значимость, площадь)]."""
    counts = np.asarray(counts, dtype=float)
    n = len(counts)
    if n < 50:
        return []
    fw = np.asarray(fw, dtype=float)
    sig = significance(counts, fw)
    em_lo, em_hi = int(math.ceil(2 * fw[0])), int(math.ceil(2 * fw[-1]))
    if em_lo:
        sig[:em_lo] = 0.0
    if em_hi:
        sig[-em_hi:] = 0.0
    cand = _thin(_maxima(sig, sig > thr), sig, fw)
    out = []
    for ch in cand:
        f = float(fw[ch])
        b = _baseline(counts, ch, f)
        h = counts[ch] - b
        if h <= 0:
            continue
        # Отсев шумовых выбросов по ИЗМЕРЕННОЙ ширине на полувысоте.
        if h >= 10.0 and f >= 4.0:
            half = b + 0.5 * h
            rng = max(2, int(round(2.5 * f)))
            r = ch
            while r < min(n - 1, ch + rng) and counts[r] >= half:
                r += 1
            left = ch
            while left > max(0, ch - rng) and counts[left] >= half:
                left -= 1
            if max(1.0, float(r - left)) < MIN_FWHM_RATIO * f:
                continue
        out.append((int(ch), float(h), f, float(sig[ch]),
                    2.507 * (f / 2.3548) * h))
    return out


# --- 2. Опознание по спектру испускания того же прогона ----------------------

def emit_lines(base, min_yield=MIN_YIELD):
    """Линии нуклида/цепочки: [(E, выход на распад)] из *_emit.csv прогона."""
    p = os.path.join(BUILD, base + "_emit.csv")
    if not os.path.exists(p):
        return []
    emit, N = kr.load_hist(p)
    if not N:
        return []
    peaks = []
    for e in sorted(emit):
        c = emit[e]
        if c <= 0:
            continue
        if peaks and e - peaks[-1][1] <= 3.0:
            w = peaks[-1][0] + c
            peaks[-1] = (w, (peaks[-1][1] * peaks[-1][0] + e * c) / w)
        else:
            peaks.append((c, e))
    return sorted((e, c / N) for c, e in peaks if c / N >= min_yield)


def identify(found, lines, k=K_ID, span=None):
    """Сопоставление найденного с библиотечным окном δE = k·ПШПВ(E).

    Возвращает (совпадения, неопознанные пики, слитые, невидимые, вне
    диапазона). Три последних списка — то, ради чего всё затевалось: они
    показывают, чего мы НЕ видим, причём с РАЗНЫМИ причинами. По прежней
    схеме (съём в заранее объявленных линиях) такого вопроса нельзя было
    даже задать.

    span — рабочий диапазон прибора (E_min, E_max) в кэВ.
    """
    used, pairs = set(), []
    for E, y in sorted(lines, key=lambda t: -t[1]):
        win = k * dc.fwhm(E)
        best = None
        for i, f in enumerate(found):
            if i in used or f["E"] is None:
                continue
            d = abs(f["E"] - E)
            if d <= win and (best is None or d < best[0]):
                best = (d, i)
        if best:
            used.add(best[1])
            pairs.append((E, y, found[best[1]], best[0]))
    # Ненайденные линии делятся на два разных случая, и валить их в одну кучу
    # нельзя: «слита с соседом» — предел разрешения прибора, тут алгоритм не
    # виноват и чинится деконволюцией; «не видна» — линия слаба или тонет в
    # континууме, и это вопрос к чувствительности. Признак слияния — рядом,
    # ближе предела Рэлея (1 ПШПВ), стоит НАЙДЕННЫЙ пик, отданный другой линии.
    got = {round(p[0], 3) for p in pairs}
    merged, unseen, outside = [], [], []
    for E, y in lines:
        if round(E, 3) in got:
            continue
        if span and not (span[0] <= E <= span[1]):
            # Третий случай, и его тоже нельзя валить в «не видно»: линия
            # ВНЕ рабочего диапазона прибора. У ряда тория это 2,8 / 12,2 /
            # 16,1 кэВ с выходами до 0,53 — самые сильные линии списка, и
            # без этой оговорки они каждый раз возглавляли бы «потери».
            outside.append((E, y))
            continue
        near = [f for f in found
                if f["E"] is not None and abs(f["E"] - E) <= dc.fwhm(E)]
        (merged if near else unseen).append((E, y))
    extra = [f for i, f in enumerate(found) if i not in used]
    key = lambda t: -t[1]                                        # noqa: E731
    return (pairs, extra, sorted(merged, key=key),
            sorted(unseen, key=key), sorted(outside, key=key))


# --- 3. Калибровка по невязкам опознанного ----------------------------------

def calibration_check(pairs, tol=CAL_TOL):
    """Невязки найденных центроид против библиотечных энергий.

    Правило ЛСРМ: если хоть на одном якоре |ΔE| > tol·ПШПВ(E), запасённая
    калибровка негодна и E(N) надо пересчитать по якорям.
    """
    rows = []
    for E, y, f, d in pairs:
        fw = dc.fwhm(E)
        rows.append(dict(E=E, y=y, ch=f["ch"], E_found=f["E"], dE=f["E"] - E,
                         frac=(f["E"] - E) / fw, sig=f["sig"]))
    worst = max((abs(r["frac"]) for r in rows), default=0.0)
    return rows, worst, worst > tol


def refit(pairs, deg=1):
    """Пересчёт E(N) по якорям: [(канал, E_библ)] -> коэффициенты E = Σ c_k·N^k.

    Степень поднимается только при нужде — правило ЛСРМ и здравый смысл:
    полином высокой степени по десятку якорей уводит края.
    """
    if len(pairs) < 3:
        return None, "якорей меньше трёх — пересчёт небезопасен"
    ch = np.array([p[2]["ch"] for p in pairs], dtype=float)
    E = np.array([p[0] for p in pairs], dtype=float)
    best = None
    for d in range(1, min(deg + 3, 4) + 1):
        if len(ch) < d + 2:
            break
        c = np.polyfit(ch, E, d)
        res = np.abs(np.polyval(c, ch) - E)
        fw = np.array([dc.fwhm(e) for e in E])
        worst = float((res / fw).max())
        if best is None or worst < best[1]:
            best = (list(c[::-1]), worst, d)
        if worst <= CAL_TOL:
            break
    return best, ""


# --- сверка со SpectraVibe ---------------------------------------------------

def cross_check(counts, fw, thr=SIGMA_THR):
    """Тот же спектр их реализацией. -> (список каналов, пояснение)."""
    root = os.environ.get("SPECTRAVIBE_ROOT")
    if not root:
        return None, "SPECTRAVIBE_ROOT не задан — сверка пропущена"
    sp_scripts = os.path.join(root, "scripts")
    if not os.path.isdir(sp_scripts):
        return None, "нет %s — сверка пропущена" % sp_scripts
    if sp_scripts not in sys.path:
        sys.path.insert(0, sp_scripts)
    try:
        from gamma.peaks.search import mariscotti_search
    except ImportError as exc:
        return None, "их модуль не импортируется: %s" % exc
    fwa = np.asarray(fw, dtype=float)
    pk = mariscotti_search(counts, lambda c: float(fwa[min(int(c),
                                                           len(fwa) - 1)]),
                           sigma_threshold=thr, filter_narrow_peaks=True,
                           min_fwhm_ratio=MIN_FWHM_RATIO)
    return [p.channel for p in pk], ""


# --- разбор одной записи -----------------------------------------------------

def analyze(sp, bg, base, subtract_bg=True):
    """Полный разбор записи: поиск, опознание, проверка калибровки."""
    ch = np.arange(len(sp.n), dtype=float)
    en = sp.energy(ch)
    y = sp.n.astype(float)
    if subtract_bg and bg is not None:
        bch = np.arange(len(bg.n), dtype=float)
        y = y - np.interp(en, bg.energy(bch), bg.n.astype(float)) * (
            sp.live / bg.live)
        y = np.maximum(y, 0.0)
    # ПШПВ в каналах: закон в энергии, переведённый через саму калибровку
    fw_ch = np.empty(len(ch))
    for i in range(len(ch)):
        e = max(en[i], 1.0)
        w = dc.fwhm(e)
        lo = np.interp(e - w / 2, en, ch)
        hi = np.interp(e + w / 2, en, ch)
        fw_ch[i] = max(1.0, abs(hi - lo))
    raw = search(y, fw_ch)
    found = [dict(ch=c, E=float(np.interp(c, ch, en)), h=h, fw=f, sig=s, area=a)
             for c, h, f, s, a in raw]
    lines = emit_lines(base)
    # Рабочий диапазон: снизу — где кончается краевая маска фильтра, сверху —
    # верх шкалы. Ниже паспортных 50 кэВ прибор и не заявлен.
    em = int(math.ceil(2 * fw_ch[0]))
    span = (max(50.0, float(en[min(em, len(en) - 1)])), float(en[-1]))
    pairs, extra, merged, unseen, outside = identify(found, lines, span=span)
    rows, worst, need = calibration_check(pairs)
    return dict(found=found, lines=lines, pairs=pairs, extra=extra,
                merged=merged, unseen=unseen, outside=outside, span=span,
                cal_rows=rows, cal_worst=worst, cal_need=need, fw_ch=fw_ch,
                counts=y, energy=en)


def _run():
    only = None
    if "--one" in sys.argv:
        only = sys.argv[sys.argv.index("--one") + 1]
    print("Поиск пиков фильтром Марискотти, опознание по спектру испускания\n"
          "того же прогона, проверка калибровки по невязкам якорей.\n")
    for geom, mask, nuc, aspec, dpct, d0, mass, vol in kr.VOLUME_RECORDS:
        if only and nuc != only:
            continue
        kd = paths.kit_dir(geom)
        files = sorted(str(p) for p in kd.rglob(mask)) if kd else []
        if not files:
            continue
        _lines, ckey = kr.VLINES[nuc]
        base = kr.RUNBASE.get((geom, ckey))
        if not base:
            continue
        sp, bg = bm.read(files[0])
        r = analyze(sp, bg, base)
        print("=" * 78)
        print("%s — %s: найдено пиков %d, линий в библиотеке прогона %d, "
              "опознано %d" % (geom, nuc, len(r["found"]), len(r["lines"]),
                               len(r["pairs"])))
        cc, why = cross_check(r["counts"], r["fw_ch"])
        if cc is None:
            print("   сверка со SpectraVibe: %s" % why)
        else:
            mine = {f["ch"] for f in r["found"]}
            both = sum(1 for c in cc if any(abs(c - m) <= 2 for m in mine))
            print("   сверка со SpectraVibe: у них %d, у нас %d, совпало %d"
                  % (len(cc), len(mine), both))
        print("   калибровка: худшая невязка %.2f ПШПВ — %s"
              % (r["cal_worst"], "ПЕРЕСЧЁТ НУЖЕН" if r["cal_need"] else "годна"))
        for row in sorted(r["cal_rows"], key=lambda x: x["E"]):
            print("      %8.1f кэВ  найдено %8.1f  ΔE %+6.1f  (%+.2f ПШПВ)"
                  % (row["E"], row["E_found"], row["dE"], row["frac"]))
        if r["merged"]:
            print("   СЛИТЫ с соседом — предел разрешения, лечится "
                  "деконволюцией (%d):" % len(r["merged"]))
            for E, yv in r["merged"][:6]:
                print("      %8.1f кэВ, выход %.3f" % (E, yv))
        if r["unseen"]:
            print("   НЕ ВИДНЫ вовсе — вопрос к чувствительности (%d):"
                  % len(r["unseen"]))
            for E, yv in r["unseen"][:6]:
                print("      %8.1f кэВ, выход %.3f" % (E, yv))
        if r["outside"]:
            print("   вне рабочего диапазона %.0f…%.0f кэВ (%d), сильнейшая "
                  "%.1f кэВ с выходом %.3f"
                  % (r["span"][0], r["span"][1], len(r["outside"]),
                     r["outside"][0][0], r["outside"][0][1]))
        if r["extra"]:
            print("   Пики БЕЗ линии в библиотеке прогона (%d):"
                  % len(r["extra"]))
            for f in sorted(r["extra"], key=lambda f: -f["sig"])[:8]:
                print("      %8.1f кэВ, значимость %5.1f, площадь %8.0f"
                      % (f["E"], f["sig"], f["area"]))


if __name__ == "__main__":
    _run()
