"""Чтение измеренных спектров АтомНано (txt FORMAT 3 и XML waterfall-viewer): каналы, калибровка, живое время."""
import re, numpy as np
def load(path):
    s = open(path, encoding="utf-8", errors="replace").read()
    if s.startswith("FORMAT"):
        L = s.split("\n"); t = float(L[8]); c = [float(L[11 + k]) for k in range(4)]
        n = int(L[9]); y = np.array([float(x) for x in L[15:15 + n] if x.strip()])
    else:
        t = float(re.search(r"<LiveTime>([\d.]+)", s).group(1))
        c = [float(x) for x in re.findall(r"<Coefficient>([-\d.e+]+)", s)]
        y = np.array([float(x) for x in re.findall(r"<DataPoint>(\d+)", s)])
    ch = np.arange(len(y)); E = sum(c[k] * ch ** k for k in range(len(c)))
    return E, y, t, c
