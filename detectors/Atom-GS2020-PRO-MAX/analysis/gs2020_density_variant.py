# -*- coding: utf-8 -*-
"""#GS-64: вариант референсного GDML (mx_oisn10) с другой плотностью, составом и уровнем пробы.
Запуск: python gs2020_density_variant.py <in.gdml> <out.gdml> <rho> <top_mm> "G4_K:0.524447,G4_Cl:0.475553" """
import sys, re, math
sys.stdout.reconfigure(encoding="utf-8")
src, dst, rho, top, spec = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
fr = [(k, float(v)) for k, v in (p.split(":") for p in spec.split(","))]
if abs(sum(v for _, v in fr) - 1.0) > 1e-6:
    raise SystemExit("ОТКАЗ: сумма долей %.6f != 1" % sum(v for _, v in fr))
t = open(src, encoding="utf-8").read()
mat = '<material name="Epoxy_crumb" state="solid"><D value="%.4f" unit="g/cm3"/>%s</material>' % (
    rho, "".join('<fraction n="%.6f" ref="%s"/>' % (v, k) for k, v in fr))
t, n1 = re.subn(r'<material name="Epoxy_crumb" state="solid">.*?</material>', mat, t, flags=re.S)
t, n2 = re.subn(r'<zplane z="42\.5" rmin="0" rmax="70\.5"/></polycone><polycone name="Headspace_solid"(.*?)<zplane z="42\.5"',
                lambda m: '<zplane z="%.2f" rmin="0" rmax="70.5"/></polycone><polycone name="Headspace_solid"%s<zplane z="%.2f"' % (top, m.group(1), top), t, flags=re.S)
if (n1, n2) != (1, 1):
    raise SystemExit("ОТКАЗ: материал найден %d раз, уровень %d раз" % (n1, n2))
V = (math.pi * (70.5**2 - 43**2) * 55.5 + math.pi * 70.5**2 * (top - 3.0)) / 1000.0
open(dst, "w", encoding="utf-8", newline="").write(t)
print("%s: rho=%.4f верх z=%.2f мм V=%.1f см3 масса=%.1f г состав=%s" % (dst, rho, top, V, rho * V, spec))
