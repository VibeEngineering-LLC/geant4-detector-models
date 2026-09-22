"""Проверка целостности .psp (входы в трубку): размер файла кратен 52, доли входов внутрь, границы y, |u|=1, NaN, сверка числа записей с .meta.
Использование: python psp_check.py <каталог>; код возврата 1, если любая проверка не пройдена."""
import sys, glob, os, numpy as np
sys.stdout.reconfigure(encoding="utf-8")
dt = np.dtype([("evt", "<i4"), ("pdg", "<i4")] + [(k, "<f4") for k in ("E", "phi", "y", "ux", "uy", "uz", "t", "Ep", "wp")] + [(k, "<i2") for k in ("proc", "vol", "ip", "pad")])
bad = 0
for f in sorted(glob.glob(os.path.join(sys.argv[1], "*.psp"))):
    n_bytes = os.path.getsize(f)
    a = np.fromfile(f, dtype=dt)
    meta = dict(l.strip().split("=", 1) for l in open(f + ".meta", encoding="utf-8") if "=" in l)
    inward = float((a["ux"] * np.cos(a["phi"]) + a["uz"] * np.sin(a["phi"]) < 0).mean()) if len(a) else 1.0
    norm = np.sqrt(a["ux"] ** 2 + a["uy"] ** 2 + a["uz"] ** 2)
    ok = (n_bytes % 52 == 0 and len(a) == int(meta["psp_records"]) and inward >= 0.9995 and (len(a) == 0 or (abs(a["y"]).max() <= 150.001 and abs(norm - 1).max() < 1e-5)) and np.isfinite(a["E"]).all() and (a["E"] > 0).all() and (a["t"] >= 0).all())
    bad += not ok
    print("%s records=%d meta=%s inward=%.6f ymax=%.3f %s" % (os.path.basename(f), len(a), meta["psp_records"], inward, abs(a["y"]).max() if len(a) else 0, "OK" if ok else "FAIL"))
sys.exit(1 if bad else 0)
