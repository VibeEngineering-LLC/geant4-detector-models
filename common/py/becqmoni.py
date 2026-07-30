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


def peak_find(sp, E0, half_win=None, roi=1.25):
    """(центроида, ПШПВ, высота) вокруг ожидаемой энергии E0, кэВ.

    ЗАЧЕМ ЦЕНТРОИДА, А НЕ ТОЛЬКО ШИРИНА. Окно площади надо ставить на РЕАЛЬНОЕ
    положение пика, а не на табличную энергию линии: у сцинтилляционного тракта
    калибровка уходит (у одного из приборов в этом хозяйстве — на 4,3 кэВ при
    662). Если окно сдвинуто относительно пика, часть площади срезается, а
    полки фона захватывают склон пика — и то и другое смещает результат в одну
    сторону, вниз. Прежняя fwhm_at находила максимум и ВЫБРАСЫВАЛА его
    положение, оставляя только ширину; площадь при этом бралась по номиналу.

    Центроида считается первым моментом по вычтенной подложке в пределах
    ±roi·ПШПВ вокруг найденного максимума — то есть по той же области, которой
    потом берётся площадь.
    """
    hw = half_win if half_win else max(3 * 0.06 * E0, 20.0)
    ch = np.arange(len(sp.n), dtype=float)
    en = sp.energy(ch)
    m = (en > E0 - hw) & (en < E0 + hw)
    if m.sum() < 5:
        return None
    x, y = en[m], sp.n[m].astype(float)
    k = max(2, len(x) // 8)
    a, b = np.polyfit(np.r_[x[:k], x[-k:]], np.r_[y[:k], y[-k:]], 1)
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
    fw = abs(hi - lo)
    w = roi * fw
    sel = (x > x[i0] - w) & (x < x[i0] + w) & (net > 0)
    if not sel.any():
        return None
    cen = float((x[sel] * net[sel]).sum() / net[sel].sum())
    return cen, fw, float(top)


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


def escape_cuts(E0, w, side, fwhm, escapes, esc_roi=0.5):
    """Полосы исключения из ЛЕВОЙ полки — там, где сидит пик вылета.

    ОДНА реализация на измеренную и на расчётную сторону. Правило, записанное в
    двух местах, расходится: в этом репозитории такое ловилось трижды, и оба
    раза цена была вывод, посчитанный по двум разным определениям.

    Полка слева — `[E−(roi+side)·ПШПВ, E−roi·ПШПВ]`; вылет с энергией `E_выл`
    даёт особенность в `E0 − E_выл`. Справа особенности нет по построению —
    вылет всегда НИЖЕ линии, — поэтому заражение односторонне и симметричным
    усреднением полок не сокращается.

    РЕЗАТЬ ПО ПЕРЕСЕЧЕНИЮ, А НЕ ПО ПОПАДАНИЮ ЦЕНТРА. Первая версия срабатывала
    только когда центр `E0 − E_выл` лежал внутри полки. Если центр чуть снаружи,
    а полоса `±esc_roi·ПШПВ` перекрывает край полки, заражение оставалось, а
    список полос выходил пустым — отказ, замаскированный под норму. Границы
    полосы энергий 54,3…217,4 кэВ как раз и выведены по попаданию ЦЕНТРА, то
    есть краевой случай лежит ровно на её концах (найдено аудитом).

    ПОЛОСЫ ОБЪЕДИНЯЮТСЯ. У CsI два вылета (иод 28,6 и цезий 30,97), их центры
    разнесены на 2,37 кэВ, а полуширина `0,5·ПШПВ` на 122 кэВ — 10,7 кэВ:
    полосы перекрываются почти целиком. Потребитель, вычитающий каждую полосу
    по отдельности, вычел бы перекрытие ДВАЖДЫ и завысил площадь. Возвращаются
    непересекающиеся полосы, обрезанные по полке, — тогда способ применения
    (маской или вычитанием) на результат не влияет.
    """
    lo, hi = E0 - w - side * fwhm, E0 - w
    raw = []
    for e in escapes or ():
        c = E0 - e
        a, b = c - esc_roi * fwhm, c + esc_roi * fwhm
        a, b = max(a, lo), min(b, hi)
        if b > a:
            raw.append((a, b))
    out = []
    for a, b in sorted(raw):
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return tuple(out)


def peak_area(sp, E0, fwhm, roi=1.25, side=1.0, escapes=(), esc_roi=0.5,
              keep_frac=0.25, detail=None):
    """Площадь пика с трапецеидальной подложкой.

    roi  — полуширина области пика в долях ПШПВ (1,25 ПШПВ = 2,94 сигма, 99,7 %);
    side — ширина каждого фонового окна в долях ПШПВ, вплотную к ROI.
    Возвращает (площадь, погрешность, подложка).

    `escapes` — энергии вылета материала кристалла: каналы вокруг `E0 − E_выл`
    исключаются из оценки полки. Границы полосы заражения и порядок величины —
    в docstring `area_broadened`. Правило применяется к ИЗМЕРЕННОЙ стороне тем
    же кодом, что к расчётной, намеренно: полагаться на сокращение смещения в
    отношении нельзя, потому что сокращение держится на том, что модель верно
    воспроизводит вылет, а это и есть проверяемое.

    Откаты здесь ДВА, а не четыре: `both` и `right`. При пустой правой полке
    функция возвращает `None` ещё до вычета — это поведение старше правки и
    менять его тут не место. Режим объявляется через `detail`, а не молчит.
    """
    w = roi * fwhm
    g, ng = sp.counts_between(E0 - w, E0 + w)
    lo_l, hi_l = E0 - w - side * fwhm, E0 - w
    bl, nl = sp.counts_between(lo_l, hi_l)
    br, nr = sp.counts_between(E0 + w, E0 + w + side * fwhm)
    if nl == 0 or nr == 0:
        return None

    # Полосы приходят НЕПЕРЕСЕКАЮЩИМИСЯ и обрезанными по полке (escape_cuts),
    # поэтому вычитание по одной корректно. Первая версия вычитала сырые
    # интервалы и на CsI, где полосы двух вылетов перекрываются почти целиком,
    # снимала перекрытие ДВАЖДЫ: bl и nl занижались, плотность полки падала,
    # площадь завышалась (найдено аудитом).
    cut = escape_cuts(E0, w, side, fwhm, escapes, esc_roi)
    nl_full = nl
    for c0, c1 in cut:
        cg, cn = sp.counts_between(c0, c1)
        bl -= cg
        nl -= cn

    left_ok = nl >= max(1, keep_frac * nl_full)
    if left_ok:
        dens = 0.5 * (bl / nl + br / nr)
        var = g + ng * ng * (bl / nl ** 2 + br / nr ** 2) / 4.0
        mode = "both"
    else:
        # Только правая полка: цена — чувствительность к наклону континуума,
        # но заражённая полка смещает сильнее (5,6 % пика на 81 кэВ).
        # Откатов `left` и `no_shelf` здесь НЕТ: при nr == 0 функция вышла выше
        # с None. Прежний docstring обещал четырёхступенчатую лестницу, которой
        # тут не было, — обещание снято (найдено аудитом).
        dens = br / nr
        var = g + ng * ng * (br / nr ** 2)
        mode = "right"
    bg = dens * ng
    if detail is not None:
        detail.update(mode=mode, n_left=nl, n_left_full=nl_full, n_right=nr,
                      cut=cut)
    return g - bg, math.sqrt(max(var, 1.0)), bg


# Якоря калибровки ФОНА, В ПОРЯДКЕ НАДЁЖНОСТИ, а не интенсивности. Порядок
# взят из методики ранжирования якорей (gamma/identification/anchor_ranks.py
# пакета gamma-spectrum-analysis): на детекторе средней разрешающей силы самая
# интенсивная линия — далеко не самая узнаваемая.
#
#   2614,5 Tl-208 — эталонный маяк, конкурентов рядом нет;
#   1460,8 K-40   — чистая область;
#    351,9 Pb-214 — узнаваема, но зона перегружена;
#    609,3 Bi-214 — сильная, НО в ториевом фоне рядом 583,2 Tl-208, и на NaI
#                   поиск отдаёт один пик на двоих. Проверено: в фоне ториевой
#                   записи этот якорь один даёт невязку 0,49 ПШПВ при том, что
#                   1460,8 и 2614,5 стоят точно. Поэтому он последний и
#                   отбрасывается первым.
BG_ANCHORS = (2614.511, 1460.822, 351.932, 609.320)


def recal_background(bg, anchors=BG_ANCHORS, fwhm_at_662=49.9, tol=0.30,
                     max_shift=40.0):
    """Проверить и при нужде пересчитать калибровку ФОНА по его линиям.

    ЗАЧЕМ. Калибровка фона может отличаться от калибровки пробы, даже если
    прибор тот же: другой сеанс, другая температура, другое усиление. Правило
    ЛСРМ прямое — применять шаги калибровки к фоновому спектру НЕЗАВИСИМО,
    не предполагая согласованности.

    В комплекте это не теория. У фоновых записей калибровка ЛИНЕЙНАЯ из двух
    коэффициентов, тогда как у проб — из четырёх-пяти, и одна на все записи.
    Поиск пиков в самом фоне показал: его собственные природные линии стоят
    ниже своих энергий на 10…15 кэВ, худшая невязка 0,43 ПШПВ при пороге ЛСРМ
    0,30. То есть ошибка усиления около 0,6 % плюс сдвиг.

    ЧЕМ ЭТО ВРЕДИТ. Фон вычитается из пробы по ШКАЛЕ ЭНЕРГИЙ. Если структуры
    фона стоят на 10 кэВ ниже, чем должны, они вычитаются не оттуда, откуда
    надо: под пиком остаётся лишнее, рядом — выеденная яма. На сильных
    источниках комплекта это мелочь, на слабых линиях и в МДА — нет.

    КАК. По каждому якорю ищется пик (bm.peak_find) в окне ±max_shift кэВ,
    центроида переводится обратно в КАНАЛ по текущей калибровке, и E(N)
    подгоняется заново по парам (канал, истинная энергия). Степень 1 при двух
    якорях, 2 при трёх и более — выше нельзя: якорей мало, и полином уведёт
    края. Если худшая невязка меньше tol·ПШПВ, калибровка признаётся годной и
    НЕ трогается (то же правило пропуска, что у ЛСРМ).

    Возвращает (спектр, диагностика). Спектр — тот же объект, если правка не
    нужна или невозможна, иначе новый с исправленной калибровкой.
    """
    if bg is None:
        return bg, {"reason": "фона нет"}
    ch = np.arange(len(bg.n), dtype=float)
    pairs, miss = [], []
    for E in anchors:
        fw = fwhm_at_662 * math.sqrt(E / 661.657)
        f = peak_find(bg, E, half_win=max_shift)
        if not f or abs(f[0] - E) > max_shift:
            miss.append(E)
            continue
        pairs.append((float(bg.channel(f[0])), E, (f[0] - E) / fw))
    if len(pairs) < 2:
        return bg, {"reason": "якорей меньше двух (найдено %d)" % len(pairs),
                    "missing": miss}
    worst = max(abs(p[2]) for p in pairs)
    diag = {"n_anchors": len(pairs), "missing": miss,
            "worst_before": worst, "dropped": [],
            "shifts": [(p[1], p[2]) for p in pairs]}
    if worst <= tol:
        diag["reason"] = "невязка %.2f ПШПВ ≤ порога %.2f — правка не нужна" \
                         % (worst, tol)
        diag["applied"] = False
        return bg, diag

    def fit(ps):
        deg = 1 if len(ps) < 3 else 2
        x = np.array([p[0] for p in ps])
        y = np.array([p[1] for p in ps])
        cf = np.polyfit(x, y, deg)[::-1]      # к виду c0 + c1*N + c2*N²
        sp2 = Spectrum(bg.n, list(cf), bg.live, bg.real, bg.label)
        res = []
        for c, E, _ in ps:
            fw = fwhm_at_662 * math.sqrt(E / 661.657)
            res.append(abs(float(sp2.energy(np.array([c]))[0]) - E) / fw)
        return sp2, res, deg

    # Отбраковка сбитого якоря. Якорь, который после подгонки всё равно уходит
    # дальше порога, не «уточняет» калибровку — он её портит: чаще всего это
    # не сдвиг шкалы, а чужой пик, притянутый поиском (классика — 609,3 рядом
    # с 583,2). Отбрасываем худший и пробуем снова, пока якорей не меньше двух.
    use = list(pairs)
    out, res, deg = fit(use)
    while max(res) > tol and len(use) > 2:
        i = max(range(len(use)), key=lambda k: res[k])
        diag["dropped"].append((use[i][1], res[i]))
        use.pop(i)
        out, res, deg = fit(use)
    diag["applied"] = True
    diag["degree"] = deg
    diag["n_used"] = len(use)
    diag["worst_after"] = max(res)
    diag["reason"] = ("невязка %.2f ПШПВ > порога %.2f — пересчёт по %d "
                      "якорям%s, стало %.2f"
                      % (worst, tol, len(use),
                         "" if not diag["dropped"] else
                         " (отброшены: %s)" % ", ".join(
                             "%.1f" % e for e, _ in diag["dropped"]),
                         max(res)))
    return out, diag


def read_checked(path, **kw):
    """read() + ПРОВЕРКА калибровки фона по его линиям. -> (проба, фон, диаг).

    Отдельная функция, а не флаг у read(): read() остаётся сырым читателем
    формата, а здесь начинается обработка. Все анализирующие скрипты берут
    записи через неё — иначе правило «калибровка фона проверяется независимо»
    держалось бы на памяти автора каждого скрипта.
    """
    sp, bg = read(path)
    bg, diag = recal_background(bg, **kw)
    return sp, bg, diag


def broaden(hist, fwhm_at_662=49.9, emax=3200.0, bin_keV=1.0, fwhm_of=None):
    """Уширить МОДЕЛЬНЫЙ спектр до разрешения прибора: массив на сетке 1 кэВ.

    Зачем это обязательно. В расчёте линии острые, в измерении — шириной в
    десятки кэВ. Если площадь из модели брать узким окном, а из измерения —
    окном в доли ПШПВ, то БЛЕНДЫ учитываются по-разному, и результат врёт.
    Пример из этой работы: Ac-228 911,2 и 968,97 разнесены на 58 кэВ при
    ПШПВ(911) = 58 кэВ, то есть в NaI это ОДИН пик. Узкое окно по модели брало
    только 911, окно ±1 ПШПВ по измерению — оба, и активность тория по этой
    линии выходила завышенной в полтора раза.

    По умолчанию ПШПВ(E) = fwhm_at_662 * sqrt(E/661,657); форма пика NaI —
    гауссиана [ЛСРМ §8.4.2.1]. Закон корня — приближение одной опорной точки:
    на этом приборе он завышает ширину на 583 кэВ и занижает на 2614. Кто
    работает подгонкой, а не окном, обязан передать откалиброванный закон
    fwhm_of(E) — тот же, каким подгоняет измерение (см. deconv.fwhm).
    """
    n = int(emax / bin_keV) + 1
    out = np.zeros(n)
    grid = np.arange(n) * bin_keV
    for E0, c in hist.items():
        if c <= 0 or E0 <= 0:
            continue
        if fwhm_of is not None:
            sig = fwhm_of(max(E0, 8.0)) / 2.3548
        else:
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


def area_broadened(arr, E0, fwhm, roi=1.0, side=1.0, bin_keV=1.0,
                   escapes=(), esc_roi=0.5, keep_frac=0.25, detail=None):
    """Площадь пика в уширенном модельном спектре — ТЕМ ЖЕ окном и полками,
    что и в измеренном (peak_area). Возвращает (площадь, подложка).

    ПОЛКА НЕ ДОЛЖНА СОДЕРЖАТЬ ИЗВЕСТНЫХ ОСОБЕННОСТЕЙ. Правило выведено на
    депозитном окне `[E−30, E−10]`, где в полку попадал пик вылета K-рентгена
    иода (`E−28,6`), но оно относится к ЛЮБОЙ полочной конвенции, включая эту,
    заданную в долях ПШПВ (указание аудитора 30.07.2026). Полка слева здесь —
    `[E−2·ПШПВ, E−1·ПШПВ]`, и вылет лежит внутри неё тогда и только тогда, когда
    `ПШПВ(E) ≤ E_выл ≤ 2·ПШПВ(E)`. При `ПШПВ = k·√E` это замкнутая полоса
    энергий `[(E_выл/2k)², (E_выл/k)²]`: для иода 28,6 кэВ — **54,3…217,4 кэВ**,
    для цезия 30,97 — 63,7…254,9 (CsI). В неё попадают пять линий комплекта:
    59,5; 81,0; 88,0; 122,1; 165,9 — весь мягкий край.

    Порядок величины смещения: вылет/пик на 81 кэВ = 5,64 %; вклад целиком
    садится в левую полку шириной 1 ПШПВ, а вычитается по окну шириной 2 ПШПВ,
    и симметричное усреднение делит вклад надвое — итого около 5,6 % пика.
    Справа особенности нет по построению (вылет всегда НИЖЕ линии), поэтому
    заражение односторонне и усреднением не сокращается.

    `escapes` — энергии вылета материала кристалла, кэВ (NaI: 28,6; CsI: 28,6 и
    30,97). Каналы в пределах `±esc_roi·ПШПВ` от `E0 − E_выл` из оценки полки
    ИСКЛЮЧАЮТСЯ. Плотность считается как `сумма/число каналов`, поэтому
    исключение само по себе смещения не вносит.

    ЛЕСТНИЦА ОТКАТА, объявляемая наружу через `detail`, а не молчаливая:
      `both`      — обе полки, вылет исключён, с каждой стороны осталось не
                    меньше `keep_frac` каналов;
      `right`     — слева осталось слишком мало: берётся только правая полка
                    (ценой чувствительности к наклону континуума);
      `no_shelf`  — не осталось ни одной годной полки, подложка не вычитается.
    Молчаливый откат был бы тем же дефектом, что и заражённая полка: число
    меняется, а причина не видна.
    """
    def win(a, b, cut=()):
        i0, i1 = int(round(a / bin_keV)), int(round(b / bin_keV))
        i0, i1 = max(0, i0), min(len(arr), i1)
        if i1 <= i0:
            return 0.0, 0
        seg = arr[i0:i1]
        n = i1 - i0
        if cut:
            idx = np.arange(i0, i1) * bin_keV
            mask = np.ones(n, dtype=bool)
            for c0, c1 in cut:
                mask &= ~((idx >= c0) & (idx <= c1))
            return float(seg[mask].sum()), int(mask.sum())
        return float(seg.sum()), n

    w = roi * fwhm
    g, ng = win(E0 - w, E0 + w)

    lo_l, hi_l = E0 - w - side * fwhm, E0 - w
    cut = escape_cuts(E0, w, side, fwhm, escapes, esc_roi)
    nl_full = max(1, int(round(hi_l / bin_keV)) - int(round(lo_l / bin_keV)))
    bl, nl = win(lo_l, hi_l, cut)
    br, nr = win(E0 + w, E0 + w + side * fwhm)

    left_ok = nl >= max(1, keep_frac * nl_full)
    if left_ok and nr:
        dens = 0.5 * (bl / nl + br / nr)
        mode = "both"
    elif nr:
        dens = br / nr
        mode = "right"
    elif left_ok:
        dens = bl / nl
        mode = "left"
    else:
        dens = 0.0
        mode = "no_shelf"
    if detail is not None:
        detail.update(mode=mode, n_left=nl, n_left_full=nl_full, n_right=nr,
                      cut=tuple(cut))
    return g - dens * ng, dens * ng


def net_rate(sample, bg, E0, fwhm, sync=True, **kw):
    """Скорость счёта в пике за вычетом ИЗМЕРЕННОГО фона, имп/с.

    САМА СКОРОСТЬ от времён не зависит: R = S/t_s − B/t_b. А вот её
    НЕОПРЕДЕЛЁННОСТЬ зависит, и здесь есть выбор.

    sync=True (по умолчанию) — обе половины пары приводятся к ОБЩЕМУ времени
    T = min(t_пробы, t_фона), и погрешность считается так, будто обе набирались
    ровно T. Смысл правила: **нельзя заявлять точность лучше, чем даёт короткое
    плечо пары**. Указание оператора.

    Когда фон длиннее пробы (обычный случай: у комплекта фон 54 000 с при
    пробах 300…11 000 с), T = t_пробы, и формула СОВПАДАЕТ с обычной — вклад
    фона в дисперсию и так мал. Разница появляется там, где проба длиннее
    фона: в комплекте это три записи — K-40 в «Денте» и Петри и Th-232 в
    Петри, где фон короче на 5…15 %. Там синхронизация даёт погрешность
    больше обычной, то есть более осторожную оценку, и это правильно: лишнее
    время пробы не восполняет недобор статистики фона.

    sync=False — прежняя формула, каждое плечо со своим временем. Оставлена,
    чтобы разницу можно было измерить, а не обсуждать.
    """
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
    if sync:
        T = min(sample.live, bg.live)
        # ВНИМАНИЕ, тонкость. Просто поделить счёты на общее время нельзя:
        # дисперсия СКОРОСТИ от базы времени не зависит вовсе, и такая
        # «синхронизация» не изменила бы ничего. Приведение к общему времени
        # T — это ПРОРЕЖИВАНИЕ более длинного набора: у пуассоновской величины,
        # прореженной с долей k, дисперсия падает в k раз (не в k²). Отсюда
        #     D(R) = σ_s²/(T·t_s) + σ_b²/(T·t_b),
        # что при t_s = t_b совпадает с обычной формулой, а при любом
        # неравенстве даёт большую, то есть осторожную оценку.
        var = a[1] ** 2 / (T * sample.live) + b[1] ** 2 / (T * bg.live)
        return rs - rb, math.sqrt(max(var, 0.0)), rb, db
    return rs - rb, math.hypot(ds, db), rb, db
