# #CAL-0: центроид гауссиана+линейная подложка (sa10_th232.fit_peak) по каждому файлу отдельно, шкала — полином файла.
import sys, numpy as np; sys.stdout.reconfigure(encoding="utf-8")
import sa10_th232 as s
out = {}
for tag, path in (("образец", s.SAMPLE_XML), ("фон", s.BKG_XML)):
    d = s.read_atomspectra_xml(path); e = s.channel_to_energy(np.arange(len(d["counts"])), d["coeffs"])
    for E0, hw in ((238.632, 30), (351.93, 25), (583.187, 45), (609.31, 40), (1460.822, 80), (1764.49, 90), (2614.511, 150)):
        f = s.fit_peak(e, d["counts"], E0, hw)
        out[(tag, E0)] = f if f and "mu" in f else None
print("%9s | %-26s | %-26s" % ("E библ", "образец: μ (±) ПШПВ", "фон: μ (±) ПШПВ"))
for E0 in sorted({k[1] for k in out}):
    row = []
    for tag in ("образец", "фон"):
        f = out[(tag, E0)]
        row.append("%8.2f (%.2f) Δ%+6.2f %5.1f" % (f["mu"], f["mu_err"], f["mu"] - E0, f["fwhm"]) if f else "не сошлась")
    print("%9.2f | %-26s | %-26s" % (E0, *row))
