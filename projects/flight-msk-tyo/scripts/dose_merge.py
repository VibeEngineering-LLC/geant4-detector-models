"""Слияние *_dose.csv фантома: скорости с весом T_f/T_comp, дисперсии с (T_f/T_comp)^2. Выход: results/final/dose_all.csv (comp,layer,class,abs_Gy_s,eq_Sv_s,sig_eq)."""
import csv, glob, os, sys
sys.stdout.reconfigure(encoding="utf-8")
d, out = sys.argv[1], sys.argv[2]
T = {}
for m in glob.glob(os.path.join(d, "*_meta.txt")):
    b = os.path.basename(m)[:-9]
    T[b] = float([l for l in open(m, encoding="utf-8") if l.startswith("T_sim_s=")][0].split("=")[1])
tc = {}
for b, t in T.items():
    tc[b.rsplit("_", 1)[0]] = tc.get(b.rsplit("_", 1)[0], 0) + t
acc = {}
for b, tf in T.items():
    c = b.rsplit("_", 1)[0]; w = tf / tc[c]
    for r in csv.DictReader(open(os.path.join(d, b + "_dose.csv"), encoding="utf-8")):
        k = (c, int(r["layer"]), int(r["class"]))
        a = acc.setdefault(k, [0.0, 0.0, 0.0])
        a[0] += float(r["abs_Gy_per_s"]) * w; a[1] += float(r["eq_Sv_per_s"]) * w; a[2] += (float(r["sigma_eq_Sv_per_s"]) * w) ** 2
with open(out, "w", encoding="utf-8", newline="\n") as f:
    f.write("comp,layer,class,abs_Gy_s,eq_Sv_s,sig_eq\n")
    f.writelines("%s,%d,%d,%.6g,%.6g,%.6g\n" % (c, l, k, a[0], a[1], a[2] ** 0.5) for (c, l, k), a in sorted(acc.items()))
print(len(acc), "строк;", len(T), "файлов;", sorted(tc))
