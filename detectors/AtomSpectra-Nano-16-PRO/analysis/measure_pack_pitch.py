# -*- coding: utf-8 -*-
"""Шаг укладки электродов в пенале — замер по фотографии пачки.

ЗАЧЕМ. Шаг укладки входит в модель дважды: он задаёт, какая часть активности
лежит под кристаллом, и он же задаёт самоэкранирование соседним стержнем.
Ни этикетка, ни каталог его не сообщают. Первое допущение «лежат вплотную»
оператор оспорил, и оно оказалось неверным: лоток разводит стержни рёбрами.

КАК. Масштаб задаёт сам стержень — его диаметр известен по этикетке (3,2 мм),
поэтому измеряется БЕЗРАЗМЕРНОЕ отношение «период укладки / ширина стержня», и
никакой линейки в кадре не требуется. Порядок:

  1. направление стержней ищется перебором угла: профиль, усреднённый вдоль
     стержней, при верном угле имеет наибольший разброс;
  2. стержень от пенала отделяется по СИНЕВЕ (B − R): у синего пластика она
     велика, у полированного вольфрама около нуля;
  3. период берётся по автокорреляции профиля, доля металла — по порогу на
     среднем уровне; шаг/диаметр = 1 / доля металла.

ГРАНИЦЫ. Метод даёт ОТНОШЕНИЕ, а не миллиметры: если этикеточные 3,2 мм
неверны, вместе с ними уедет и шаг. Перспектива кадра меняет масштаб поперёк
поля, поэтому участок берётся коротким. Точность ~10 %: порог на «синеве»
двигает долю металла на несколько процентов, а край стержня в кадре размыт.
Проверка линейкой отменяет этот замер.

    python analysis/measure_pack_pitch.py <фото> [x0 x1 y0 y1] [диаметр, мм]

x0..y1 — доли размера кадра, вырезка с укладкой в пенале (по умолчанию весь).
"""
import sys

import numpy as np
from PIL import Image


def profile_across(img_rgb, ang_deg, channel):
    """Профиль поперёк стержней под заданным углом. channel: 'blue'|'lum'."""
    a = np.asarray(img_rgb, dtype=float)
    v = a[:, :, 2] - a[:, :, 0] if channel == "blue" else a.mean(axis=2)
    h, w = v.shape
    yy, xx = np.mgrid[0:h, 0:w]
    t = np.radians(ang_deg)
    s = xx * np.cos(t) + yy * np.sin(t)
    idx = np.rint(s - s.min()).astype(int)
    n = idx.max() + 1
    cnt = np.bincount(idx.ravel(), minlength=n)
    tot = np.bincount(idx.ravel(), weights=v.ravel(), minlength=n)
    ok = cnt > 0.6 * cnt.max()          # только полностью заполненные строки
    if ok.sum() < 40:
        return None
    return tot[ok] / cnt[ok]


def find_angle(img_rgb):
    best = None
    for ang in np.arange(-89, 89, 0.5):
        p = profile_across(img_rgb, ang, "blue")
        if p is None:
            continue
        if best is None or p.std() > best[1]:
            best = (ang, p.std(), p)
    return best


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    path = sys.argv[1]
    box = [float(v) for v in sys.argv[2:6]] if len(sys.argv) >= 6 else \
        [0.0, 1.0, 0.0, 1.0]
    diam = float(sys.argv[6]) if len(sys.argv) > 6 else 3.20

    im = Image.open(path).convert("RGB")
    w, h = im.size
    crop = im.crop((int(box[0] * w), int(box[2] * h),
                    int(box[1] * w), int(box[3] * h)))
    print("кадр %d x %d, вырезка %s" % (w, h, crop.size))

    ang, score, p = find_angle(crop)
    print("направление поперёк стержней %.1f град, контраст %.2f, "
          "длина профиля %d" % (ang, score, len(p)))

    p0 = p - p.mean()
    ac = np.correlate(p0, p0, "full")[len(p0) - 1:]
    ac /= ac[0]
    peaks = [(k, ac[k]) for k in range(6, len(ac) - 1)
             if ac[k] > ac[k - 1] and ac[k] > ac[k + 1] and ac[k] > 0.2]
    if not peaks:
        raise SystemExit("периодичность не найдена — не тот участок кадра")
    print("автокорреляция:", ", ".join("%d px: %.2f" % (k, v)
                                       for k, v in peaks[:6]))
    period = peaks[0][0]
    metal = float((p < p.mean()).mean())
    ratio = 1.0 / metal
    print("период %d px, доля металла %.3f" % (period, metal))
    print("ШАГ / ДИАМЕТР = %.2f" % ratio)
    print("при диаметре %.2f мм: шаг %.2f мм, зазор между стержнями %.2f мм"
          % (diam, diam * ratio, diam * (ratio - 1.0)))
    print("масштаб кадра %.2f px/мм, пачка из 10 стержней %.1f мм"
          % (period / (diam * ratio), diam * (9 * ratio + 1)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
