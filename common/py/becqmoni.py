"""Чтение спектров BecqMoni XML комплекта поверки Гамма-1С и работа с пиками.

В одном файле лежат два спектра: образец (<EnergySpectrum>) и фон той же
геометрии (<BackgroundEnergySpectrum>). У каждого своя энергетическая
калибровка и своё время, поэтому ROI считается в ЭНЕРГИИ и переводится в
каналы отдельно для каждого спектра — перебивать спектры в общую сетку не
нужно и не надо, это только размазало бы статистику.

Площадь пика — канонический способ [Будыка §7.7; Gilmore §5.5]: трапеция по
левому и правому фоновым окнам, ширина окон и ROI задаются в долях ПШПВ.
ПШПВ НЕ берётся из заголовка (там неясно, полином по каналу или по энергии),
а меряется прямо по пику: по точкам полувысоты.

ПРОВЕРЕНО на Cs-137 комплекта (маринелли, 1025 Бк на дату измерения):
восстановленная отсюда эффективность против записанной ЛСРМ в .efr —
   ROI ±1,00 ПШПВ  ->  1,901e-2  против 1,871e-2   (+1,6 %)
   ROI ±1,25 ПШПВ  ->  1,988e-2                    (+6,3 %)
   ROI ±1,50 ПШПВ  ->  2,024e-2                    (+8,1 %)
То есть чтение спектра, времена и калибровка верны, а ШИРИНА ОКНА — главная
систематика этой процедуры, до 8 %. Гауссиана столько не объясняет (±1,0 ПШПВ
содержит 98,2 % площади, ±1,5 — 99,98 %), значит остаток даёт недовычет
подложки на широком окне. Для сверки с ЛСРМ уместно ±1,0, для МИА — ±1,25,
где окно содержит практически весь пик и согласуется с полной расчётной
эффективностью МК.
"""
import math
import re
import xml.etree.ElementTree as ET

import numpy as np


def _floats(node, tag):
    out = []
    for c in node.iter(tag):
        for v in c.iter():
            if v.text and v is not c:
                try:
                    out.append(float(v.text.strip()))
                except ValueError:
                    pass
    return out


class Spectrum:
    def __init__(self, counts, cal, live, real, label):
        self.n = np.asarray(counts, dtype=float)
        self.cal = list(cal)                 # E = c0 + c1*ch + c2*ch^2 + ...
        self.live = live
        self.real = real
        self.label = label

    def energy(self, ch):
        ch = np.asarray(ch, dtype=float)
        return sum(c * ch ** k for k, c in enumerate(self.cal))

    def channel(self, E):
        """Обратное преобразование численно: калибровка монотонна на диапазоне."""
        ch = np.arange(len(self.n), dtype=float)
        return float(np.interp(E, self.energy(ch), ch))

    def counts_between(self, E0, E1):
        a, b = self.channel(E0), self.channel(E1)
        lo, hi = int(math.floor(a)), int(math.ceil(b))
        lo, hi = max(0, lo), min(len(self.n), hi)
        return float(self.n[lo:hi].sum()), hi - lo


def read(path):
    """-> (образец, фон). Фон может быть None."""
    root = ET.parse(path).getroot()
    out = []
    for tag in ("EnergySpectrum", "BackgroundEnergySpectrum"):
        node = None
        for e in root.iter():
            if e.tag.endswith(tag):
                node = e
                break
        if node is None:
            out.append(None)
            continue
        cal = [float(c.text) for c in node.iter()
               if c.tag.endswith("Coefficient")]
        pulses = None
        for e in node.iter():
            if e.tag.endswith("Spectrum") and e is not node:
                pulses = [float(x.text) for x in e if x.text is not None]
        if pulses is None:
            pulses = []
            for e in node.iter():
                if e.tag.endswith("Pulses"):
                    pulses = [float(x.text) for x in e if x.text is not None]
        # ВАЖНО: MeasurementTime — реальное время, LiveTime — живое. Для
        # скорости счёта нужно ЖИВОЕ: на Cs-137 они расходятся на 2 %.
        live = real = None
        for e in node.iter():
            if e.tag.endswith("MeasurementTime"):
                real = float(e.text)
            elif e.tag.endswith("LiveTime"):
                live = float(e.text)
        if live is None:
            live = real
        out.append(Spectrum(pulses, cal, live, real, tag))
    return out[0], out[1]


def fwhm_at(sp, E0, half_win=None):
    """ПШПВ по точкам полувысоты вокруг пика, кэВ. Грубая оценка окна — 12 % от E."""
    hw = half_win if half_win else max(3 * 0.06 * E0, 20.0)
    ch = np.arange(len(sp.n), dtype=float)
    en = sp.energy(ch)
    m = (en > E0 - hw) & (en < E0 + hw)
    if m.sum() < 5:
        return None
    x, y = en[m], sp.n[m].astype(float)
    # подложка — прямая по краям окна
    k = max(2, len(x) // 8)
    xb = np.r_[x[:k], x[-k:]]
    yb = np.r_[y[:k], y[-k:]]
    a, b = np.polyfit(xb, yb, 1)
    net = y - (a * x + b)
    i0 = int(np.argmax(net))
    top = net[i0]
    if top <= 0:
        return None
    half = 0.5 * top

    def cross(idx, step):
        i = idx
        while 0 <= i + step < len(net) and net[i + step] > half:
            i += step
        if not (0 <= i + step < len(net)):
            return None
        y1, y2 = net[i], net[i + step]
        t = (y1 - half) / (y1 - y2) if y1 != y2 else 0.0
        return x[i] + t * (x[i + step] - x[i])

    lo, hi = cross(i0, -1), cross(i0, +1)
    if lo is None or hi is None:
        return None
    return abs(hi - lo)


def peak_area(sp, E0, fwhm, roi=1.25, side=1.0):
    """Площадь пика с трапецеидальной подложкой.

    roi  — полуширина области пика в долях ПШПВ (1,25 ПШПВ = 2,94 сигма, 99,7 %);
    side — ширина каждого фонового окна в долях ПШПВ, вплотную к ROI.
    Возвращает (площадь, погрешность, подложка).
    """
    w = roi * fwhm
    g, ng = sp.counts_between(E0 - w, E0 + w)
    bl, nl = sp.counts_between(E0 - w - side * fwhm, E0 - w)
    br, nr = sp.counts_between(E0 + w, E0 + w + side * fwhm)
    if nl == 0 or nr == 0:
        return None
    # трапеция: средняя плотность подложки по двум окнам
    dens = 0.5 * (bl / nl + br / nr)
    bg = dens * ng
    var = g + ng * ng * (bl / nl ** 2 + br / nr ** 2) / 4.0
    return g - bg, math.sqrt(max(var, 1.0)), bg


def broaden(hist, fwhm_at_662=49.9, emax=3200.0, bin_keV=1.0):
    """Уширить МОДЕЛЬНЫЙ спектр до разрешения прибора: массив на сетке 1 кэВ.

    Зачем это обязательно. В расчёте линии острые, в измерении — шириной в
    десятки кэВ. Если площадь из модели брать узким окном, а из измерения —
    окном в доли ПШПВ, то БЛЕНДЫ учитываются по-разному, и результат врёт.
    Пример из этой работы: Ac-228 911,2 и 968,97 разнесены на 58 кэВ при
    ПШПВ(911) = 58 кэВ, то есть в NaI это ОДИН пик. Узкое окно по модели брало
    только 911, окно ±1 ПШПВ по измерению — оба, и активность тория по этой
    линии выходила завышенной в полтора раза.

    ПШПВ(E) = fwhm_at_662 * sqrt(E/661,657); форма пика NaI — гауссиана
    [ЛСРМ §8.4.2.1].
    """
    n = int(emax / bin_keV) + 1
    out = np.zeros(n)
    grid = np.arange(n) * bin_keV
    for E0, c in hist.items():
        if c <= 0 or E0 <= 0:
            continue
        sig = fwhm_at_662 * math.sqrt(max(E0, 8.0) / 661.657) / 2.3548
        lo = max(0, int((E0 - 5 * sig) / bin_keV))
        hi = min(n, int((E0 + 5 * sig) / bin_keV) + 1)
        if hi <= lo:
            continue
        g = np.exp(-0.5 * ((grid[lo:hi] - E0) / sig) ** 2)
        s = g.sum()
        if s > 0:
            out[lo:hi] += c * g / s
    return out


def area_broadened(arr, E0, fwhm, roi=1.0, side=1.0, bin_keV=1.0):
    """Площадь пика в уширенном модельном спектре — ТЕМ ЖЕ окном и полками,
    что и в измеренном (peak_area). Возвращает (площадь, подложка)."""
    def win(a, b):
        i0, i1 = int(round(a / bin_keV)), int(round(b / bin_keV))
        i0, i1 = max(0, i0), min(len(arr), i1)
        return arr[i0:i1].sum(), max(1, i1 - i0)

    w = roi * fwhm
    g, ng = win(E0 - w, E0 + w)
    bl, nl = win(E0 - w - side * fwhm, E0 - w)
    br, nr = win(E0 + w, E0 + w + side * fwhm)
    dens = 0.5 * (bl / nl + br / nr)
    return g - dens * ng, dens * ng


def net_rate(sample, bg, E0, fwhm, **kw):
    """Скорость счёта в пике за вычетом ИЗМЕРЕННОГО фона, имп/с."""
    a = peak_area(sample, E0, fwhm, **kw)
    if a is None:
        return None
    rs, ds = a[0] / sample.live, a[1] / sample.live
    if bg is None:
        return rs, ds, 0.0, 0.0
    b = peak_area(bg, E0, fwhm, **kw)
    if b is None:
        return rs, ds, 0.0, 0.0
    rb, db = b[0] / bg.live, b[1] / bg.live
    return rs - rb, math.hypot(ds, db), rb, db
