"""Сырой спектр (несвёрнутый) -> формат lines_report: E_keV, отсч./ч/кэВ, сигма. Использование: raw_to_cph.py <total_all.csv> <out.csv> (бины 20..9999 кэВ)."""
import sys, numpy as np
d = np.genfromtxt(sys.argv[1], delimiter=",", skip_header=2)[20:10000]
o = open(sys.argv[2], "w", encoding="utf-8", newline="\n")
o.write("E_keV,counts_per_h_per_keV,sigma_counts_per_h_per_keV\n")
for b, r, v in zip(d[:, 0], d[:, 1], d[:, 2]):
    o.write(f"{b + 0.5},{r * 3600},{np.sqrt(v) * 3600}\n")
