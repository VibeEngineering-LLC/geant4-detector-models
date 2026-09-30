# -*- coding: utf-8 -*-
# #GS-40 (оператор 29.09): сосуд воды и KCl — колодец внутр. Ø79, снаружи в полости Ø87 (стенка 4, потолок 2 не замерен); тело
# КОНИЧЕСКОЕ: «донышко внутри 136», у верха Ø141, стенка 2. Объём KCl 740 мл (мерный цилиндр), плотность 829/740. Печать: верх KCl,
# уровень 1 л от внутреннего верха (замер 16,5 мм). Запуск: python gs2020_kcl_well79.py
import math, re, sys
sys.stdout.reconfigure(encoding="utf-8")
SRC = r"C:\g4work\gs2020\kcl\GS2020_marinelli_kcl_740ml.gdml"; DST = r"C:\g4work\gs2020\kcl\GS2020_marinelli_kcl_740ml_well79.gdml"
RB, RT, ZB, ZT, RW, ZW, W = 68.0, 70.5, -52.5, 52.5, 43.5, 3.0, 2.0     # внутр. радиус дна/верха, z дна/верха полости, колодец
ri = lambda z: RB + (RT - RB) * (z - ZB) / (ZT - ZB); ro = lambda z: ri(z) + W if z <= ZT else RT + W
def vol(z):   # объём полости от дна до z, мл: усечённый конус минус колодец
    c = lambda a, b: math.pi * (b - a) / 3 * (ri(a) ** 2 + ri(a) * ri(b) + ri(b) ** 2)
    return (c(ZB, z) - math.pi * RW ** 2 * (min(z, ZW) - ZB)) / 1000.0
level = lambda V: next(z / 1000.0 for z in range(int(ZW * 1000), int(ZT * 1000)) if vol(z / 1000.0) >= V)
top, lit = level(740.0), level(1000.0)
pc = lambda name, zs: '<polycone name="%s" startphi="0" deltaphi="360" aunit="deg" lunit="mm">%s</polycone>' % (name, "".join('<zplane z="%.3f" rmin="%g" rmax="%.3f"/>' % q for q in zs))
wall = pc("VesselWall_solid", [(-54.5, 39.5, ro(-52.5)), (1.0, 39.5, ro(1.0)), (1.0, 0, ro(1.0)), (54.5, 0, RT + W)])
mat = pc("SourceMatrix_solid", [(ZB, RW, ri(ZB)), (ZW, RW, ri(ZW)), (ZW, 0, ri(ZW)), (top, 0, ri(top))]); air = pc("Headspace_solid", [(top, 0, ri(top)), (ZT, 0, RT)])
t = open(SRC, encoding="utf-8").read()
for name, new in (("VesselWall_solid", wall), ("SourceMatrix_solid", mat), ("Headspace_solid", air)):
    t, n = re.subn(r'<polycone name="%s".*?</polycone>' % name, new, t, flags=re.S)
    if n != 1: raise SystemExit("ОТКАЗ: %s заменён %d раз" % (name, n))
open(DST, "w", encoding="utf-8", newline="").write(t)
print("KCl 740 мл: верх z=%.2f (слой над колодцем %.2f мм), плотность %.4f; 1 л: %.1f мм ниже внутр. верха (замер 16,5); полость %.0f мл -> %s" % (top, top - ZW, 829 / 740, ZT - lit, vol(ZT), DST))
