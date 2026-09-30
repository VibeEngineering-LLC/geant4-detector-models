# -*- coding: utf-8 -*-
r"""Исправление геометрии уровня пробы (27.09, дефекты W-153): над пробой — воздух, а не Vessel_PP; плотность =
масса / объём пробы при этом уровне. Запуск: python gs2020_fix_level.py <in.gdml> <out.gdml> [масса_г=1052]"""
import math, re, sys
sys.stdout.reconfigure(encoding="utf-8")
src, dst = sys.argv[1], sys.argv[2]
M = float(sys.argv[3]) if len(sys.argv) > 3 else 1052.0
t = open(src, encoding="utf-8").read()
top = float(re.search(r'<polycone name="SourceMatrix_solid".*?<zplane z="([-0-9.]+)" rmin="0" rmax="70.5"/></polycone>', t, re.S).group(1))
V = (math.pi * (70.5**2 - 43**2) * 55.5 + math.pi * 70.5**2 * (top - 3.0)) / 1000.0
rho = M / V
assert t.count('<material name="Epoxy_crumb" state="solid"><D value="0.7987" unit="g/cm3"/>') == 1
t = t.replace('<material name="Epoxy_crumb" state="solid"><D value="0.7987" unit="g/cm3"/>',
              '<material name="Epoxy_crumb" state="solid"><D value="%.4f" unit="g/cm3"/>' % rho)
air = ('<polycone name="Headspace_solid" startphi="0" deltaphi="360" aunit="deg" lunit="mm">'
       '<zplane z="%.1f" rmin="0" rmax="70.5"/><zplane z="52.5" rmin="0" rmax="70.5"/></polycone>' % top)
assert t.count("</solids>") == 1 and t.count('<volume name="LV_VesselWall">') == 1
t = t.replace("</solids>", air + "</solids>")
t = t.replace('<volume name="LV_VesselWall">', '<volume name="LV_Headspace"><materialref ref="G4_AIR"/><solidref ref="Headspace_solid"/></volume><volume name="LV_VesselWall">')
t = t.replace('</physvol></volume><volume name="World">', '</physvol><physvol name="PV_LV_Headspace"><volumeref ref="LV_Headspace"/><position name="pos_Headspace" x="0" y="0" z="0" unit="mm"/></physvol></volume><volume name="World">')
assert t.count("PV_LV_Headspace") == 1, "не вставлен воздух"
open(dst, "w", encoding="utf-8", newline="").write(t)
print("верх пробы z=%.1f мм, объём %.1f см3, плотность %.4f г/см3, воздух z %.1f..52.5 мм -> %s" % (top, V, rho, top, dst))
