"""Спектрограмма АтомНано (txt FORMAT 3 + срезы по 5 строк: время мс, широта, долгота, длительность с, строка счётов через табуляцию) -> npz.
Использование: read_spectrogram.py <spectrogram.txt> <out.npz>"""
import sys, numpy as np
L = open(sys.argv[1], encoding="utf-8", errors="replace").read().split("\n")
cal = [float(L[k]) for k in range(11, 15)]
n = int(L[9]); i = 15 + n
t, la, lo, dt, rows = [], [], [], [], []
while i + 4 < len(L) and L[i].strip():
    t.append(int(L[i]) / 1000.0); la.append(float(L[i + 1])); lo.append(float(L[i + 2])); dt.append(float(L[i + 3]))
    rows.append(np.array(L[i + 4].strip("\t").split("\t"), float)[:n]); i += 5
ch = np.arange(n); E = sum(cal[k] * ch ** k for k in range(4))
np.savez_compressed(sys.argv[2], t=np.array(t), lat=np.array(la), lon=np.array(lo), dt=np.array(dt), C=np.array(rows, dtype=np.int32), E=E, cal=cal)
print("срезов:", len(t), "суммарно с:", sum(dt), "отсчётов:", int(np.array(rows).sum()), "калибровка:", cal)
