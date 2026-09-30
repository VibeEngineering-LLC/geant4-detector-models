# -*- coding: utf-8 -*-
# #GS-40 v3 (оператор 29.09 днём, чертёж «пластик 3 мм везде»): срез→дно внутри 103, наружная высота 106, срез→верх колодца 46,
# колодец внутри 57; полость Ø136 у дна → Ø141 у верха, колодец внутри Ø79 (прежние замеры), снаружи — аргумент (85 или 87).
# Рамка прежняя: внутренняя поверхность потолка колодца (упор прибора) z = 1,0; крышка 3 мм над срезом (как в прежней модели).
# Запуск: python gs2020_kcl_marinelli_v3.py <Ø колодца снаружи> [write]
import math, re, sys
sys.stdout.reconfigure(encoding="utf-8")
DW = float(sys.argv[1]); W = 3.0; RB, RT, RW = 68.0, float(__import__("os").environ.get("GS_RT", "70.5")), DW / 2   # GS_RT: внутр. радиус у среза
ZWI = 1.0; ZW = ZWI + W; ZBO = ZWI - 57.0; ZB = ZBO + W; ZT = ZB + 103.0    # потолок колодца, дно снаружи/внутри, срез
SRC = r"C:\g4work\gs2020\kcl\GS2020_marinelli_kcl_740ml.gdml"; DST = r"C:\g4work\gs2020\kcl\GS2020_marinelli_kcl_740ml_v3_w%g.gdml" % DW
ri = lambda z: RB + (RT - RB) * (z - ZB) / (ZT - ZB); ro = lambda z: ri(z) + W
def vol(z):   # объём полости от дна до z, мл: усечённый конус минус колодец
    c = lambda a, b: math.pi * (b - a) / 3 * (ri(a) ** 2 + ri(a) * ri(b) + ri(b) ** 2)
    return (c(ZB, z) - math.pi * RW ** 2 * (min(z, ZW) - ZB)) / 1000.0
level = lambda V: next(z / 1000.0 for z in range(int(ZW * 1000), int(ZT * 1000)) if vol(z / 1000.0) >= V)
top, lit = level(740.0), level(1000.0)
pc = lambda name, zs: '<polycone name="%s" startphi="0" deltaphi="360" aunit="deg" lunit="mm">%s</polycone>' % (name, "".join('<zplane z="%.3f" rmin="%g" rmax="%.3f"/>' % q for q in zs))
wall = pc("VesselWall_solid", [(ZBO, 39.5, ro(ZBO)), (ZWI, 39.5, ro(ZWI)), (ZWI, 0, ro(ZWI)), (ZT + W, 0, RT + W)])
print(f"колодец снаружи Ø{DW:g}: стенка колодца {RW - 39.5:.1f}; срез−верх колодца {ZT - ZW:.1f} (замер 46); наружная высота {ZT - ZBO:.1f} (106)")
if len(sys.argv) < 3 or sys.argv[2] != "write": DST = None
mat = pc("SourceMatrix_solid", [(ZB, RW, ri(ZB)), (ZW, RW, ri(ZW)), (ZW, 0, ri(ZW)), (top, 0, ri(top))]); air = pc("Headspace_solid", [(top, 0, ri(top)), (ZT, 0, RT)])
t = open(SRC, encoding="utf-8").read()
for name, new in (("VesselWall_solid", wall), ("SourceMatrix_solid", mat), ("Headspace_solid", air)):
    t, n = re.subn(r'<polycone name="%s".*?</polycone>' % name, new, t, flags=re.S)
    if n != 1: raise SystemExit("ОТКАЗ: %s заменён %d раз" % (name, n))
if DST: open(DST, "w", encoding="utf-8", newline="").write(t)
print("KCl 740 мл: верх z=%.2f (слой над колодцем %.2f мм), плотность %.4f; 1 л: %.1f мм ниже среза (замер 16,5); полость %.0f мл -> %s" % (top, top - ZW, 829 / 740, ZT - lit, vol(ZT), DST))
