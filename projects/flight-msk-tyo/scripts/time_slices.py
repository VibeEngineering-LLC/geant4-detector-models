"""Срезы спектрограммы по фазам рейса: спектры 20-300 кэВ и пик 66 кэВ по времени. Использование: time_slices.py <spectrogram.npz> <out_prefix>"""
import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
d = np.load(sys.argv[1]); t = d["t"] - d["t"][0]; C = d["C"]; E = d["E"]; dt = d["dt"]
PH = [("Дананг, в самолёте до руления", 0, 2300), ("руление, взлёт", 2300, 2650), ("набор высоты", 2650, 3500), ("эшелон", 3600, 4500), ("снижение", 4500, 5400), ("после посадки, руление и стоянка", 5600, 6850), ("конец записи (счёт вырос до 4-5 отсч./с)", 6900, 7395)]
pk = (E >= 61) & (E < 72); sb = ((E >= 50) & (E < 58)) | ((E >= 76) & (E < 87)); w = pk.sum() / sb.sum()
fig, ax = plt.subplots(2, 1, figsize=(10, 9)); rows = []
for name, a, b in PH:
    m = (t >= a) & (t < b); T = dt[m].sum(); c = C[m].sum(0); ed = np.arange(24, 200, 2)
    h = np.array([c[(E >= x) & (E < x + 2)].sum() for x in ed]) / T / 2; ax[0].step(ed, h, where="post", label=f"{name} ({T:.0f} с)")
    ex = c[pk].sum() - w * c[sb].sum(); rows.append((name, T, c[(E >= 20) & (E < 3000)].sum() / T, ex / T, np.sqrt(c[pk].sum() + w * w * c[sb].sum()) / T))
ax[0].set_yscale("log"); ax[0].set_xlabel("E, кэВ (шкала прибора)"); ax[0].set_ylabel("отсч./с/кэВ"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
tb = np.arange(0, 7300, 100); tot = [C[(t >= x) & (t < x + 100)][:, (E >= 20) & (E < 3000)].sum() / dt[(t >= x) & (t < x + 100)].sum() for x in tb]
pe = [(C[(t >= x) & (t < x + 100)][:, pk].sum() - w * C[(t >= x) & (t < x + 100)][:, sb].sum()) / dt[(t >= x) & (t < x + 100)].sum() for x in tb]
ax[1].plot(tb, tot, label="счёт 20–3000 кэВ"); ax[1].plot(tb, np.array(pe) * 20, label="пик 66 кэВ, ×20"); ax[1].set_xlabel("время от начала записи, с"); ax[1].legend(); ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig(sys.argv[2] + ".png", dpi=80)
open(sys.argv[2] + ".md", "w", encoding="utf-8", newline="\n").write("| фаза | время, с | счёт 20–3000 кэВ, отсч./с | пик 66 кэВ, отсч./с |\n|---|---|---|---|\n" + "".join(f"| {n} | {T:.0f} | {c:.2f} | {e:.3f} ± {s:.3f} |\n" for n, T, c, e, s in rows))
print(open(sys.argv[2] + ".md", encoding="utf-8").read())
