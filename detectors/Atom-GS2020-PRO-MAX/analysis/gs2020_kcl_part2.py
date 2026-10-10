# Дрейф KCl 88,3 ч (оператор 09.10 «Разобрать дрейф»): 2-я часть = накопление 88,3 ч − 41,5 ч, тот же формат XML → тег kcl через GS2020_KCL_XML
import sys, os, xml.etree.ElementTree as ET, numpy as np
sys.stdout.reconfigure(encoding="utf-8")
R = r"<DOSIM>\Спектры\Atom GS2020 PRO MAX\Референсы"
A, B = (os.path.join(R, "KCl ч 1 л 1085 г (%s).xml" % t) for t in ("41,5 ч, 05-07.10", "88,3 ч, 05-09.10"))
OUT = r"<WORKDIR>\GEANT4\results\kcl_parts\KCl ч 1 л 1085 г часть 2 (88,3 − 41,5 ч).xml"
es = lambda t: t.getroot().find(".//EnergySpectrum")
ta, tb = ET.parse(A), ET.parse(B)
ea, eb = es(ta), es(tb)
ca, cb = (np.array([int(x.text) for x in e.find("Spectrum").findall("DataPoint")]) for e in (ea, eb))
d = cb - ca
neg = np.flatnonzero(d < 0)
print("каналов %d; отрицательных разностей %d%s" % (len(d), len(neg), (" (первые: %s)" % list(zip(neg[:10], d[neg[:10]]))) if len(neg) else ""))
if len(neg): raise SystemExit("ОТКАЗ: 88,3 ч не накопление 41,5 ч — разность отрицательна")
for k in ("LiveTime", "MeasurementTime", "ValidPulseCount", "TotalPulseCount"):
    va, vb = float(ea.find(k).text), float(eb.find(k).text)
    eb.find(k).text = ("%.1f" if k == "LiveTime" else "%d") % (vb - va); print(k, va, vb, eb.find(k).text)
for x, v in zip(eb.find("Spectrum").findall("DataPoint"), d): x.text = str(int(v))
tb.getroot().find(".//StartTime").text = ta.getroot().find(".//EndTime").text
tb.write(OUT, encoding="utf-8", xml_declaration=True)
la, lb = float(ea.find("LiveTime").text), float(eb.find("LiveTime").text)
print("скорость: часть 1 %.3f, часть 2 %.3f имп/с; записано %s" % (ca.sum() / la, d.sum() / lb, OUT))
