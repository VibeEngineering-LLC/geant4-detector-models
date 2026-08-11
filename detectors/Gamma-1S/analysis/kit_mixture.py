"""Состав смесевых источников комплекта — из самих спектров.

В описи комплекта смесь названа «№SRC-04 Am-Ti-Eu-Cs», но поле COMMENT в .spe
хранит ОДНУ строку, поэтому паспортная активность указана лишь для одного
нуклида из четырёх: Am-241 для Маринелли, Cs-137 для Денты и Петри. Полного
состава в файлах нет.

Здесь состав определяется по спектру: ищем значимые пики над измеренным фоном
и сопоставляем их с линиями кандидатов. Кандидаты берутся широким списком,
включающим и то, что могло скрываться за «Ti» в названии.

Метод: сглаженный спектр минус фон (по живому времени), поиск локальных
максимумов с превышением над локальной подложкой более чем 4 сигма, затем
сопоставление найденных энергий с линиями кандидатов по допуску 0,6 ПШПВ.
"""
import glob
import math
import os
import sys

import numpy as np

# Корни путей — из переменных окружения (common/py/paths.py), чтобы в коде не
# было ни одного пути, привязанного к конкретной машине.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import becqmoni as bm  # noqa: E402

KIT = str(paths.ref("Gamma-1S"))

# Состав установлен по .efr: секции «#1 (Ti-44)» с линиями 67,9 / 78,4 / 1157,0
# и «#1 (Eu-152)» стоят в паспорте геометрии Петри. То есть «Ti» в названии
# источника — действительно ТИТАН-44, а не опечатка.
#
# Ti-44 (60 лет) -> Sc-44 (T1/2=4,042 ч, IAEA NUBASE2020 -- см. amticseu-
# remarks.md §11 о более ранней ошибке 3,97 ч по памяти) -> Ca-44.
# Sc-44 испускает позитрон в 94 %
# случаев и следом квант 1157,0 кэВ. Отсюда в спектре:
#   511,0 — аннигиляция (а НЕ Na-22, как можно подумать по силе линии);
#   1668,0 — СУММ-ПИК 511 + 1157: позитрон и квант испускаются в истинном
#            совпадении, и в маринельке этот пик заметен.
CAND = {
    "Am-241": [59.54],
    "Ti-44/Sc-44": [67.87, 78.32, 1157.02],
    "аннигиляция": [511.00],
    "сумм-пик 511+1157": [1668.02],
    "Eu-152": [121.78, 244.70, 344.28, 411.12, 443.97, 778.90, 867.38,
               964.06, 1085.84, 1112.07, 1408.01],
    "Cs-137": [661.657],
    # естественный фон и возможные спутники — если проявятся
    "K-40": [1460.82],
    "Ra-226 ряд": [186.2, 295.22, 351.93, 609.32, 1120.29, 1764.49],
    "Th-232 ряд": [238.63, 338.32, 583.19, 911.20, 968.97, 2614.51],
}

# Линии для самокалибровки: сильные, изолированные и заведомо присутствующие.
ANCHORS = [59.54, 511.00, 661.657, 1157.02, 1408.01]


def centroid(en, net, i, half):
    """Центр тяжести пика по чистым отсчётам в окне ±half каналов.

    Канал максимума систематически смещён вниз на падающем континууме —
    на этом уже обожглись: энергии выходили на 0,5–1,7 % ниже истинных, и
    сопоставление линий превращалось в гадание.
    """
    lo, hi = max(0, i - half), min(len(net), i + half + 1)
    w = np.clip(net[lo:hi], 0, None)
    if w.sum() <= 0:
        return en[i]
    return float((en[lo:hi] * w).sum() / w.sum())


def find_peaks(s, b, emin=40.0, emax=2800.0, nsig=4.0):
    """-> список (энергия, чистые отсчёты, значимость)."""
    ch = np.arange(len(s.n), dtype=float)
    en = s.energy(ch)
    scale = s.live / b.live if b else 0.0
    bg = np.interp(en, b.energy(np.arange(len(b.n), dtype=float)), b.n) * scale \
        if b else np.zeros_like(s.n)
    net = s.n - bg
    out = []
    i = 2
    while i < len(net) - 3:
        E = en[i]
        if E < emin or E > emax:
            i += 1
            continue
        fw = 0.075 * E * math.sqrt(661.657 / max(E, 1.0))     # ПШПВ прибора
        half = max(1, int(round(0.5 * fw / (en[1] - en[0]))))
        lo, hi = max(0, i - 3 * half), min(len(net), i + 3 * half + 1)
        if net[i] != net[lo:hi].max() or net[i] <= 0:
            i += 1
            continue
        # подложка по краям окна
        k = max(1, half)
        base = 0.5 * (net[max(0, i - 3 * half):max(1, i - 2 * half)].mean()
                      + net[min(len(net) - 1, i + 2 * half):hi].mean())
        area = net[max(0, i - half):min(len(net), i + half + 1)].sum() \
            - base * (2 * half + 1)
        gross = s.n[max(0, i - half):min(len(net), i + half + 1)].sum()
        var = gross + (bg[max(0, i - half):min(len(net), i + half + 1)].sum()
                       * max(scale, 1e-9))
        sig = area / math.sqrt(max(var, 1.0))
        if sig >= nsig and area > 0:
            out.append((centroid(en, net - base, i, half), area, sig))
            i += 2 * half
        else:
            i += 1
    return out


def selfcalib(peaks):
    """Поправка к энергетической шкале по опорным линиям.

    Записанная в файлах смесей калибровка не сходится: цезий читается как 657
    вместо 661,7, америций как 53 вместо 59,5. Строим поправку E_ист = a + b*E,
    сопоставляя опорные линии с ближайшими сильными пиками. Начальное
    приближение — по двум крайним опорным, затем уточнение по всем, что нашлись.
    """
    if len(peaks) < 3:
        return 0.0, 1.0, []
    strong = sorted(peaks, key=lambda p: -p[2])[:8]
    # грубо: самый мягкий сильный пик = 59,54; самый близкий к 660 = 661,657
    lo = min(strong, key=lambda p: p[0])
    cs = min(strong, key=lambda p: abs(p[0] - 660.0))
    b = (661.657 - 59.54) / (cs[0] - lo[0])
    a = 59.54 - b * lo[0]
    for _ in range(3):
        pairs = []
        for L in ANCHORS:
            cand = [p for p in peaks if abs(a + b * p[0] - L) <= 0.05 * L]
            if cand:
                pairs.append((max(cand, key=lambda p: p[2])[0], L))
        if len(pairs) < 2:
            break
        x = np.array([p[0] for p in pairs])
        y = np.array([p[1] for p in pairs])
        b, a = np.polyfit(x, y, 1)
    return a, b, pairs


def match(peaks, s):
    """Сопоставление найденных пиков с линиями кандидатов."""
    hits = {}
    unmatched = []
    for E, a, sg in peaks:
        fw = 0.075 * E * math.sqrt(661.657 / max(E, 1.0))
        best, bd = None, 1e9
        for nuc, lines in CAND.items():
            for L in lines:
                if abs(L - E) < bd and abs(L - E) <= 0.6 * fw:
                    bd, best = abs(L - E), (nuc, L)
        if best:
            hits.setdefault(best[0], []).append((best[1], E, a, sg))
        else:
            unmatched.append((E, a, sg))
    return hits, unmatched


if __name__ == "__main__":
    files = sorted(glob.glob(os.path.join(KIT, "*", "*Am-Ti*.xml"))
                   + glob.glob(os.path.join(KIT, "*", "*AmTiCsEu*.xml")))
    if not files:
        raise SystemExit("не найдены файлы смесей")
    for p in files:
        geom = os.path.basename(os.path.dirname(p))
        s, b = bm.read(p)
        pk = find_peaks(s, b)
        a, k, anchors = selfcalib(pk)
        pk = [(a + k * e, ar, sg) for e, ar, sg in pk]
        hits, un = match(pk, s)
        print("\n===== %s : %s" % (geom, os.path.basename(p)))
        print("живое %.0f с, значимых пиков %d; шкала исправлена по %d опорным:"
              " E' = %+.1f %+.4f*E" % (s.live, len(pk), len(anchors), a, k))
        for nuc in sorted(hits, key=lambda k: -max(x[3] for x in hits[k])):
            ls = ", ".join("%.1f (%.0f сигма)" % (x[1], x[3]) for x in hits[nuc])
            print("   %-12s %s" % (nuc, ls))
        if un:
            print("   не опознаны: " +
                  ", ".join("%.1f (%.0f)" % (e, g) for e, _, g in un))
