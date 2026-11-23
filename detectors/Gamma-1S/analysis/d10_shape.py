"""D10 (06.11.2026): форма пика «узкая компонента + плечи» при неизменной ПШПВ. Ключи GS_SH_W0/GS_SH_P/GS_SH_R; W0=0 — выключено."""
import math
import os

SH_W0, SH_P, SH_R = (float(os.environ.get(k, v)) for k, v in (("GS_SH_W0", "0"), ("GS_SH_P", "1"), ("GS_SH_R", "2")))
_cache = {}


def sh_f(w):
    """ПШПВ смеси (1−w)·N(1) + w·N(SH_R) в единицах σ узкой компоненты (бисекция, кэш по w)."""
    if w not in _cache:
        h = lambda x: (1 - w) * math.exp(-.5 * x * x) + w / SH_R * math.exp(-.5 * (x / SH_R) ** 2)
        lo, hi = 0.0, 10.0 * SH_R
        for _ in range(60):
            mid = .5 * (lo + hi)
            lo, hi = (mid, hi) if h(mid) > .5 * h(0.0) else (lo, mid)
        _cache[w] = 2 * lo
    return _cache[w]
