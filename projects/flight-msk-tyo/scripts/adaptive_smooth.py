"""Адаптивное сглаживание спектра по погрешности: ширина окна растёт, пока относительная ошибка среднего по окну не станет <= target."""
import numpy as np


def adaptive_smooth(y, var, target=0.15, wmax=601):
    """y — скорость счёта на канал, var — дисперсия каждого канала (независимые). Возвращает (сглаженный ряд, ширина окна на канал)."""
    n = len(y)
    cy = np.concatenate([[0.0], np.cumsum(y)]); cv = np.concatenate([[0.0], np.cumsum(var)])
    out = np.array(y, dtype=float); wid = np.ones(n, dtype=int); done = np.zeros(n, dtype=bool)
    idx = np.arange(n)
    for w in range(1, wmax + 1, 2):
        h = w // 2
        lo = np.clip(idx - h, 0, n); hi = np.clip(idx + h + 1, 0, n)
        s = cy[hi] - cy[lo]; v = cv[hi] - cv[lo]
        ok = (~done) & ((s > 0) & (np.sqrt(v) <= target * s) | (w >= wmax))
        out[ok] = s[ok] / np.maximum(hi[ok] - lo[ok], 1); wid[ok] = w; done |= ok
        if done.all():
            break
    return out, wid


def presmooth(vals, var, path, target):
    """Сглаживает слой (если target>0 и есть дисперсия) и пишет в формате total для spec_smear. Возвращает ширины окна на канал."""
    vals = np.asarray(vals, dtype=float)
    wid = np.ones(len(vals), dtype=int)
    if target > 0 and var is not None:
        vals, wid = adaptive_smooth(vals, np.asarray(var, dtype=float), target)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("#\nbin_keV,sumw,sumw2,light_sumw,light_sumw2\n")
        f.writelines("%d,%.10g,0,0,0\n" % (i, v) for i, v in enumerate(vals))
    return wid
