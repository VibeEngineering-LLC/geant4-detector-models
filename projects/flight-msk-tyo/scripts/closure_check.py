"""Замыкающая проверка: реплей (пустой мир -> трубка -> прибор) против прямого счёта фотонов на приборе. Использование: closure_check.py <каталог closure> [префикс реплея, по умолчанию replay]. Код 1, если полоса 20-10000 кэВ расходится > 5 % или иная полоса > 3 сигм и > 3 %."""
import sys, numpy as np
d = sys.argv[1]
def load(n):
    a = np.genfromtxt(f"{d}/{n}_total.csv", delimiter=",", skip_header=2)
    return a[:-1, 1], a[:-1, 2]
(r, vr), (q, vq) = load(sys.argv[2] if len(sys.argv) > 2 else "replay"), load("direct")
bad = 0
for lo, hi in [(20, 100), (100, 300), (300, 1000), (1000, 3000), (20, 10000)]:
    a, b = r[lo:hi].sum(), q[lo:hi].sum()
    s = np.sqrt(vr[lo:hi].sum() + vq[lo:hi].sum())
    z = (a - b) / s if s > 0 else 0.0
    ok = abs(a - b) <= 0.05 * b if hi == 10000 else (abs(z) <= 3 or abs(a - b) <= 0.03 * b)
    bad += not ok
    print(f"{lo}-{hi} keV: replay={a:.4f} direct={b:.4f} 1/s diff={100*(a-b)/b:+.2f}% z={z:+.1f} {'ok' if ok else 'FAIL'}")
sys.exit(1 if bad else 0)
