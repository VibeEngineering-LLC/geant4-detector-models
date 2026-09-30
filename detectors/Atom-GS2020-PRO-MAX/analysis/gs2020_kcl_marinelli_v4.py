# -*- coding: utf-8 -*-
# #GS-40 v4 (оператор 29.09 вечер «бери небольшой конус - у дна сосуда 85, у дна колодца 83»): колодец снаружи — конус
# Ø85 у дна сосуда (z=ZBO) → Ø83 у дна колодца (потолок, z=ZW); пластик 3 мм везде → изнутри Ø79 → Ø77. Прочее — как v3:
# срез→дно внутри 103, наружная высота 106, срез→верх колодца 46, колодец внутри 57; полость Ø136 → Ø141; крышка 3 мм.
# Потолок колодца 3 мм = 106 − 46 − 57 — вычислен по замерам, точно (оператор 29.09: «вычислен оп замерам. точно»).
# Запуск: python gs2020_kcl_marinelli_v4.py <Ø снаружи у дна сосуда> <Ø снаружи у дна колодца> [write]
import math, re, sys
sys.stdout.reconfigure(encoding="utf-8")
DWB, DWT = float(sys.argv[1]), float(sys.argv[2]); W = 3.0; RB, RT = 68.0, float(__import__("os").environ.get("GS_RT", "70.5"))   # GS_RT: внутр. радиус у среза
ZWI = 1.0; ZW = ZWI + W; ZBO = ZWI - 57.0; ZB = ZBO + W; ZT = ZB + 103.0    # потолок колодца, дно снаружи/внутри, срез
rw = lambda z: DWB / 2 + (DWT - DWB) / 2 * (z - ZBO) / (ZW - ZBO)   # наружный радиус колодца (линейно ZBO → ZW)
SRC = r"C:\g4work\gs2020\kcl\GS2020_marinelli_kcl_740ml.gdml"; DST = r"C:\g4work\gs2020\kcl\GS2020_marinelli_kcl_740ml_v4_w%g_%g.gdml" % (DWB, DWT)
# #GS-44 (29.09): GS_MAT=water — проба вода (G4_WATER) вместо KCl; GS_FILL — объём пробы, мл (по умолчанию 740 — как было)
MAT = __import__("os").environ.get("GS_MAT", "kcl"); FILL = float(__import__("os").environ.get("GS_FILL", "740"))
if MAT not in ("kcl", "water"): raise SystemExit("ОТКАЗ: GS_MAT=%s" % MAT)
if (MAT, FILL) != ("kcl", 740.0): DST = r"C:\g4work\gs2020\kcl\GS2020_marinelli_%s_%gml_v4_w%g_%g.gdml" % (MAT, FILL, DWB, DWT)
ri = lambda z: RB + (RT - RB) * (z - ZB) / (ZT - ZB); ro = lambda z: ri(z) + W
def vol(z):   # объём полости от дна до z, мл: усечённый конус полости минус усечённый конус колодца
    c = lambda f, a, b: math.pi * (b - a) / 3 * (f(a) ** 2 + f(a) * f(b) + f(b) ** 2)
    return (c(ri, ZB, z) - c(rw, ZB, min(z, ZW))) / 1000.0
level = lambda V: next(z / 1000.0 for z in range(int(ZW * 1000), int(ZT * 1000)) if vol(z / 1000.0) >= V)
top, lit = level(FILL), level(1000.0)
pc = lambda name, zs: '<polycone name="%s" startphi="0" deltaphi="360" aunit="deg" lunit="mm">%s</polycone>' % (name, "".join('<zplane z="%.3f" rmin="%g" rmax="%.3f"/>' % q for q in zs))
wall = pc("VesselWall_solid", [(ZBO, rw(ZBO) - W, ro(ZBO)), (ZWI, rw(ZWI) - W, ro(ZWI)), (ZWI, 0, ro(ZWI)), (ZT + W, 0, RT + W)])
print(f"колодец снаружи Ø{2 * rw(ZBO):.2f} у дна сосуда → Ø{2 * rw(ZW):.2f} у дна колодца; изнутри Ø{2 * (rw(ZBO) - W):.2f} → Ø{2 * (rw(ZWI) - W):.2f}; срез−верх колодца {ZT - ZW:.1f} (замер 46); наружная высота {ZT - ZBO:.1f} (106)")
if len(sys.argv) < 4 or sys.argv[3] != "write": DST = None
mat = pc("SourceMatrix_solid", [(ZB, rw(ZB), ri(ZB)), (ZW, rw(ZW), ri(ZW)), (ZW, 0, ri(ZW)), (top, 0, ri(top))]); air = pc("Headspace_solid", [(top, 0, ri(top)), (ZT, 0, RT)])
t = open(SRC, encoding="utf-8").read()
for name, new in (("VesselWall_solid", wall), ("SourceMatrix_solid", mat), ("Headspace_solid", air)):
    t, n = re.subn(r'<polycone name="%s".*?</polycone>' % name, new, t, flags=re.S)
    if n != 1: raise SystemExit("ОТКАЗ: %s заменён %d раз" % (name, n))
if MAT == "water":
    t, n = re.subn(r'(<volume name="LV_SourceMatrix"><materialref ref=")Epoxy_crumb(")', r"\1G4_WATER\2", t)
    if n != 1: raise SystemExit("ОТКАЗ: материал пробы заменён %d раз" % n)
if DST: open(DST, "w", encoding="utf-8", newline="").write(t)
print("%s %g мл: верх z=%.2f (слой над колодцем %.2f мм), плотность %s; 1 л: %.1f мм ниже среза (замер 16,5); полость %.0f мл -> %s" % ("KCl" if MAT == "kcl" else "вода", FILL, top, top - ZW, "%.4f" % (829 / 740) if MAT == "kcl" else "G4_WATER", ZT - lit, vol(ZT), DST))
