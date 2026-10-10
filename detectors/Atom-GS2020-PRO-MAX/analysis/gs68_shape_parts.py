# #GS-68: форма пика (ПШПВ, FW1/5, FW1/10, асимметрия) по сырым спектрам на шкале их файлов; окно 2,5 ПШПВ, боковины 50 кэВ.
import sys, xml.etree.ElementTree as ET, numpy as np
from peak_shape_measure import measure_shape, compare_to_gauss
sys.stdout.reconfigure(encoding="utf-8")
R = "<DOSIM>/Спектры/Atom GS2020 PRO MAX/Референсы/"; K = "<WORKDIR>/GEANT4/results/kcl_parts/"
SP = [("KCl ч1 41,5ч", R + "KCl ч 1 л 1085 г (41,5 ч, 05-07.10).xml", 1460.8), ("KCl ч2 46,8ч", K + "KCl ч 1 л 1085 г часть 2 (88,3 − 41,5 ч).xml", 1460.8),
      ("KCl сумма 88,3ч", R + "KCl ч 1 л 1085 г (88,3 ч, 05-09.10).xml", 1460.8), ("KCl сумма выровн.", K + "KCl ч 1 л 1085 г сумма выровненная (88,3 ч).xml", 1460.8),
      ("KCl 740мл 9,3ч", R + "KCl ч 740 мл 829 г (9,3 ч, 27-28.09).xml", 1460.8),
      ("Th-232 2614", R + "Калибровка Th-232 (без вычета фона).xml", 2614.5), ("Th-232 911", R + "Калибровка Th-232 (без вычета фона).xml", 911.2)]
print("спектр | имп/с | E верш | ПШПВ | FW1/5 изб% | FW1/10 изб% | асим 1/2 | асим 1/10 | лево1/10 | право1/10")
for name, path, pk in SP:
    es = ET.parse(path).getroot().find(".//EnergySpectrum")
    y = np.array([float(x.text) for x in es.find("Spectrum").findall("DataPoint")])
    cal = ET.parse(path).getroot().find(".//EnergyCalibration")  # у 740 мл шкалы в файле нет → шкала ч1 (только отношения ширин)
    c = [float(v.text) for v in cal.find("Coefficients")] if cal is not None else c
    E = sum(ci * (np.arange(len(y)) ** i) for i, ci in enumerate(c))
    lt = float(es.find("LiveTime").text)
    hw = 2.5 * (80.0 if pk < 2000 else 107.0) * (1 if pk > 1000 else 0.75)
    s = measure_shape(E, y, pk, side_keV=50.0, half_win_keV=hw); g = compare_to_gauss(s)
    print(f"{name} | {y.sum()/lt:.0f} | {s['E_вершины']:.1f} | {s['ПШПВ']:.2f} | {g['FW1/5_изб_%']:.1f} | {g['FW1/10_изб_%']:.1f} | "
          f"{s['асимметрия_1/2']:.3f} | {s['асимметрия_1/10']:.3f} | {s['лево_1/10']:.1f} | {s['право_1/10']:.1f}")
