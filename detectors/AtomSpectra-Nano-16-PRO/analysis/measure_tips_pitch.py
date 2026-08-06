# -*- coding: utf-8 -*-
"""Шаг укладки электродов WT-20 ПО КРАСНЫМ КОНЧИКАМ (маркировка марки).

Зачем именно по кончикам. Голый вольфрам полирован: на нём блик, и порог по
яркости берёт не тело прутка, а светлую полосу — на этом уже сорвался замер по
яркости (получилось 6,15 мм вместо 8,1). Красная маркировка WT-20 матовая и
кроет цилиндр целиком, поэтому в канале «краснота» R − (G+B)/2 сегмент — это
ПОЛНЫЙ силуэт прутка, а не блик. Ширина сегмента = диаметр 3,2 мм, расстояние
между центрами = шаг укладки.

Метод — по скиллу `photo-metrology`: сканлиния, билинейная выборка профиля,
сглаживание, порог (max+min)/2, отбраковка краевых и узких сегментов, масштаб
по эталону НА ТОЙ ЖЕ ЛИНИИ, контрольная картинка, три независимые линии.

Отношение «шаг / ширина» инвариантно к наклону сканлинии, пока прутки
параллельны, поэтому строго перпендикулярную линию проводить не требуется.

    python analysis/measure_tips_pitch.py <фото> x0 y0 x1 y1 [ещё тройки линий]

Без координат берётся зашитая тройка сканлиний для `6971217310.webp`.
"""
import os
import sys

import numpy as np
from PIL import Image

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(_HERE, "..", "drawings",
                                    "wt20_tips_scanlines.png"))
D_ROD = 3.20        # мм, диаметр по этикетке — единственная привязка масштаба

# Три сканлинии поперёк красных кончиков на 6971217310.webp (координаты
# оригинала). Разнесены вдоль прутков, чтобы проверить повторяемость.
# Три сканлинии сдвинуты ВДОЛЬ прутков на ±8 px. Больше сдвигать нельзя:
# прутки в лотке выдвинуты на РАЗНУЮ длину, и линия, отодвинутая на 20 px,
# промахивалась мимо коротких меток — находила 3–4 сегмента вместо девяти.
# Это не шум метода, а свойство кадра, и оно ограничивает базу проверки.
DEFAULT_LINES = [((560, 60), (880, 500)),
                 ((566, 53), (886, 493)),
                 ((554, 67), (874, 507))]


def profile(arr, p0, p1):
    """Билинейная выборка вдоль отрезка. -> (значения, шаг в px)."""
    h, w = arr.shape
    n = int(np.hypot(p1[0] - p0[0], p1[1] - p0[1])) * 2
    xs = np.linspace(p0[0], p1[0], n)
    ys = np.linspace(p0[1], p1[1], n)
    x0 = np.clip(xs.astype(int), 0, w - 2)
    y0 = np.clip(ys.astype(int), 0, h - 2)
    fx, fy = xs - x0, ys - y0
    v = (arr[y0, x0] * (1 - fx) * (1 - fy) + arr[y0, x0 + 1] * fx * (1 - fy)
         + arr[y0 + 1, x0] * (1 - fx) * fy + arr[y0 + 1, x0 + 1] * fx * fy)
    return v, np.hypot(p1[0] - p0[0], p1[1] - p0[1]) / (n - 1)


def segments(mask):
    out, start = [], None
    for i, m in enumerate(mask):
        if m and start is None:
            start = i
        elif not m and start is not None:
            out.append((start, i - 1))
            start = None
    if start is not None:
        out.append((start, len(mask) - 1))
    return out


def refine_edges(prof, s, lo_bg, hi_bg):
    """Границы сегмента по ПОЛУВЫСОТЕ красного плато над НУЛЁМ красноты.

    Ноль здесь физичен: «краснота» R − (G+B)/2 равна нулю у всего нейтрального
    и отрицательна у синего лотка, поэтому граница красной метки — там, где
    сигнал спадает вдвое от плато, а не там, где он переходит в минус.

    Чем плохи два других порога, оба испробованы:
      * (max+min)/2 по всему профилю — профиль несимметричен (метка +190,
        лоток −155), середина лежит на +18, у самого подножия; сегмент
        забирает размытую кромку, и ширина эталона завышается;
      * полувысота относительно СОСЕДНЕГО фона — тот же эффект, потому что
        соседний фон отрицателен.
    В webp красный канал вдобавок подвыбран по цветности, кромка размыта на
    несколько пикселей, поэтому конвенцию границы надо объявлять явно: здесь
    это ПШПВ красной метки.
    """
    peak = prof[s[0]:s[1] + 1].max()
    half = 0.5 * peak
    i = s[0]
    while i > 0 and prof[i] > half:
        i -= 1
    j = s[1]
    while j < len(prof) - 1 and prof[j] > half:
        j += 1
    return i, j


def analyse(red, p0, p1):
    prof, step = profile(red, p0, p1)
    prof = np.convolve(prof, np.ones(5) / 5.0, mode="same")
    # грубая сегментация: всё, что заметно выше нуля, — красная метка
    thr = 0.5 * prof.max()
    raw = segments(prof > thr)
    raw = [s for s in raw if s[0] > 2 and s[1] < len(prof) - 3]
    if len(raw) < 3:
        return None
    # уточнение границ по полувысоте с локальным фоном из соседних провалов
    segs = []
    for k, s in enumerate(raw):
        lo = prof[raw[k - 1][1]:s[0]].min() if k else prof[:s[0]].min()
        hi = (prof[s[1]:raw[k + 1][0]].min() if k + 1 < len(raw)
              else prof[s[1]:].min())
        segs.append(refine_edges(prof, s, lo, hi))
    w = np.array([s[1] - s[0] + 1 for s in segs], float)
    big = np.median(w[w >= 0.5 * w.max()])
    segs = [s for s, ww in zip(segs, w) if ww >= 0.5 * big]
    w = np.array([s[1] - s[0] + 1 for s in segs], float) * step
    c = np.array([0.5 * (s[0] + s[1]) for s in segs], float) * step
    if len(c) < 3:
        return None
    gaps = np.diff(c)
    scale = D_ROD / np.median(w)          # мм/px
    return dict(prof=prof, thr=thr, step=step, segs=segs, widths=w,
                gaps=gaps, scale=scale, pitch=np.median(gaps) * scale,
                p0=p0, p1=p1)


def main():
    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        path = os.path.join(os.path.expanduser("~"), "Мой диск", "Дозиметрия",
                            "Руководства", "Атом Нано 16 про", "WT-20",
                            "6971217310.webp")
    coords = [float(v) for v in sys.argv[2:]]
    lines = ([((coords[i], coords[i + 1]), (coords[i + 2], coords[i + 3]))
              for i in range(0, len(coords) - 3, 4)]
             if len(coords) >= 4 else DEFAULT_LINES)

    im = Image.open(path).convert("RGB")
    a = np.asarray(im, dtype=float)
    # КАНАЛ «КРАСНОТА»: у красной маркировки он велик, у синего лотка и белого
    # фона — около нуля или отрицателен, у голого вольфрама тоже около нуля.
    red = a[:, :, 0] - 0.5 * (a[:, :, 1] + a[:, :, 2])
    print("кадр %d x %d, сканлиний %d" % (im.size[0], im.size[1], len(lines)))

    res, pitches = [], []
    for k, (p0, p1) in enumerate(lines, 1):
        r = analyse(red, p0, p1)
        if not r:
            print("линия %d: сегментов мало — не годится" % k)
            continue
        res.append(r)
        pitches.append(r["pitch"])
        print("линия %d: сегментов %d, ширина медиана %.1f px = %.2f мм, "
              "масштаб %.4f мм/px" % (k, len(r["segs"]), np.median(r["widths"]),
                                      D_ROD, r["scale"]))
        print("   шаги, мм: %s"
              % ", ".join("%.2f" % (g * r["scale"]) for g in r["gaps"]))
        print("   ШАГ (медиана) = %.2f мм, разброс %.2f…%.2f"
              % (r["pitch"], r["gaps"].min() * r["scale"],
                 r["gaps"].max() * r["scale"]))
    if not res:
        raise SystemExit("ни одна сканлиния не дала сегментов")

    print("\nИТОГ по %d линиям: %s -> медиана %.2f мм, разброс между линиями "
          "%.2f мм" % (len(pitches), ", ".join("%.2f" % p for p in pitches),
                       float(np.median(pitches)),
                       float(max(pitches) - min(pitches))))
    # погрешность: порог даёт ±1 px на каждую границу сегмента
    err = 2.0 * res[0]["scale"] * res[0]["step"]
    print("погрешность метода ~±%.2f мм (порог ±1 px на границу)" % err)

    # --- контрольная картинка (по скиллу — обязательна) ---------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 6.4))
    ax1.imshow(np.asarray(im))
    for k, r in enumerate(res, 1):
        (x0, y0), (x1, y1) = r["p0"], r["p1"]
        ax1.plot([x0, x1], [y0, y1], "-", lw=1.2, color="#00b0ff")
        n = int(np.hypot(x1 - x0, y1 - y0)) * 2
        for s in r["segs"]:
            t = 0.5 * (s[0] + s[1]) / (n - 1)
            ax1.plot(x0 + (x1 - x0) * t, y0 + (y1 - y0) * t, "y+", ms=9)
    ax1.set_title("сканлинии и найденные центры красных кончиков",
                  fontsize=10)
    ax1.axis("off")
    for k, r in enumerate(res, 1):
        ax2.plot(np.arange(len(r["prof"])) * r["step"], r["prof"],
                 lw=1.0, label="линия %d" % k)
    ax2.axhline(res[0]["thr"], color="r", ls="--", lw=1.0, label="порог")
    ax2.set_xlabel("px вдоль линии")
    ax2.set_ylabel("краснота  R − (G+B)/2")
    ax2.set_title("профиль красноты; сегмент = красный кончик", fontsize=10)
    ax2.legend(fontsize=8)
    fig.suptitle("WT-20: шаг укладки по красным кончикам, масштаб — диаметр "
                 "%.1f мм на той же линии" % D_ROD, fontsize=11)
    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=130)
    print("контрольная картинка: %s" % OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
