"""Сливает дисперсии категорий (*_cat2.csv этапа II) в один файл: pg,proc,vol,bin_i (скорости счёта^2)."""
import glob, os, sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
d, out = sys.argv[1], sys.argv[2]
T = {os.path.basename(m)[:-9]: float([l for l in open(m, encoding="utf-8") if l.startswith("T_sim_s=")][0].split("=")[1]) for m in glob.glob(os.path.join(d, "*_meta.txt"))}
tt = {}
for b, t in T.items():
    tt[b.rsplit("_", 1)[0]] = tt.get(b.rsplit("_", 1)[0], 0) + t
acc = {}
for b, tf in T.items():
    for l in list(open(os.path.join(d, b + "_cat2.csv"), encoding="utf-8"))[1:]:
        p = l.strip().split(",")
        k = (int(p[0]), int(p[1]), int(p[2]))
        acc[k] = acc.get(k, 0) + np.array(p[3:], dtype=float) * (tf / tt[b.rsplit("_", 1)[0]]) ** 2
f = open(out, "w", encoding="utf-8", newline="\n")
f.write("pg,proc,vol," + ",".join("bin_%d" % i for i in range(len(next(iter(acc.values()))))) + "\n")
for k, v in acc.items():
    f.write("%d,%d,%d," % k + ",".join("%.6g" % x for x in v) + "\n")
print(len(acc), "категорий")
