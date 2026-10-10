# #GS-70 / #CAL-0: проверка шкалы спектра по реперам — гаусс + линия в окне, шкала из файла. argv: xml
# Ba-137m K-рентген: Kα2 31,817 (I 3,64), Kα1 32,194 (6,63), Kβ 36,4 (2,39 сумм.) — % на распад Cs-137 (ENSDF/LNHB)
import sys, xml.etree.ElementTree as ET, numpy as np
from scipy.optimize import curve_fit
sys.stdout.reconfigure(encoding="utf-8")
r = ET.parse(sys.argv[1]).getroot(); es = r.find(".//EnergySpectrum")
y = np.array([float(v.text) for v in es.find("Spectrum").findall("DataPoint")]); c = [float(v.text) for v in r.iter("Coefficient")]
ch = np.arange(len(y)); E = sum(ci * ch ** i for i, ci in enumerate(c))
print("коэффициенты файла:", c, "живое", es.find("LiveTime").text)
def g(x, a, mu, s, b0, b1): return a * np.exp(-0.5 * ((x - mu) / s) ** 2) + b0 + b1 * (x - mu)
for name, E0, hw in (("Ba K (Kα+Kβ)", 32.5, 14), ("Cs-137", 661.657, 75), ("K-40", 1460.822, 120), ("Tl-208", 2614.511, 160)):
    m = np.abs(E - E0) < hw; x, yy = E[m], y[m]
    p, cv = curve_fit(g, x, yy, p0=[yy.max() - yy.min(), E0, hw / 3, yy.min(), 0], sigma=np.sqrt(np.maximum(yy, 1)), maxfev=20000)
    fw = 2.3548 * abs(p[2]); dmu = np.sqrt(cv[1, 1])
    print(f"{name}: центроид {p[1]:.2f} ± {dmu:.2f} кэВ (опорн. {E0}) невязка {p[1]-E0:+.2f} кэВ = {(p[1]-E0)/fw:+.3f} ПШПВ; ПШПВ {fw:.2f} кэВ")
