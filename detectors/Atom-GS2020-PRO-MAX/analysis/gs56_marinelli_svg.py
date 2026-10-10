# -*- coding: utf-8 -*-
# #GS-56 (оператор 07.10 «маринелли обнови»): схема Маринелли 2 на странице GS2020 под опубликованную модель
# KCl 1 л (генератор gs2020_kcl_marinelli_v4.py 85 83, GS_FILL=1000, GS_RHO=1.085). Координаты — из тех же формул v4,
# не руками. Ось SVG: y = 10 + (ZT − z) мм (срез сосуда — y=10), x = 110 ± r. Правит img/marinelli2_section.svg на месте.
import math, re, sys
sys.stdout.reconfigure(encoding="utf-8")
DWB, DWT, W, RB, RT = 85.0, 83.0, 3.0, 68.0, 70.5
ZWI = 1.0; ZW = ZWI + W; ZBO = ZWI - 57.0; ZB = ZBO + W; ZT = ZB + 103.0
FILL, MASS = 1000.0, 1085
rw = lambda z: DWB / 2 + (DWT - DWB) / 2 * (z - ZBO) / (ZW - ZBO)
ri = lambda z: RB + (RT - RB) * (z - ZB) / (ZT - ZB); ro = lambda z: ri(z) + W
def vol(z):
    c = lambda f, a, b: math.pi * (b - a) / 3 * (f(a) ** 2 + f(a) * f(b) + f(b) ** 2)
    return (c(ri, ZB, z) - c(rw, ZB, min(z, ZW))) / 1000.0
level = lambda V: next(z / 1000.0 for z in range(int(ZW * 1000), int(ZT * 1000)) if vol(z / 1000.0) >= V)
top, lit = level(FILL), level(1000.0)
Y = lambda z: 10 + (ZT - z); P = lambda r, z, s: "%g,%g" % (round(110 + s * r, 2), round(Y(z), 2))
kcl = [P(ri(top), top, -1), P(ri(top), top, 1), P(ri(ZB), ZB, 1), P(rw(ZB), ZB, 1), P(rw(ZW), ZW, 1),
       P(rw(ZW), ZW, -1), P(rw(ZB), ZB, -1), P(ri(ZB), ZB, -1)]
wall = [P(ro(ZT), ZT, -1), P(ro(ZBO), ZBO, -1), P(rw(ZBO) - W, ZBO, -1), P(rw(ZWI) - W, ZWI, -1), P(rw(ZWI) - W, ZWI, 1),
        P(rw(ZBO) - W, ZBO, 1), P(ro(ZBO), ZBO, 1), P(ro(ZT), ZT, 1), P(RT, ZT, 1), P(ri(ZB), ZB, 1), P(rw(ZB), ZB, 1),
        P(rw(ZW), ZW, 1), P(rw(ZW), ZW, -1), P(rw(ZB), ZB, -1), P(ri(ZB), ZB, -1), P(RT, ZT, -1)]
lv_m = 26.5; zm = ZT - (lv_m - 10)
f = r"web\gs2020-th232-page\img\marinelli2_section.svg"
t = open(f, encoding="utf-8").read()
L, Rr = ZT - lit, top - ZW
subs = [
    (r"разрез модели Geant4 v3 \(вода, KCl, цезий\)", "разрез модели Geant4 v4 (KCl 1 л)"),
    (r"разрез модели Geant4 \(вода, KCl, цезий\), размеры в мм", "разрез модели Geant4 (проба KCl 1 л), размеры в мм"),
    (r"<!-- проба KCl 740 мл: верх z=18,37 -->\s*<polygon points=\"[^\"]+\"", '<!-- проба KCl %g мл: верх z=%.2f -->\n  <polygon points="%s"' % (FILL, top, " ".join(kcl))),
    (r"(<!-- стенка сосуда, 3 мм -->\s*<polygon points=\")[^\"]+\"", r'\g<1>%s"' % " ".join(wall)),
    (r"<!-- уровни 1 л: [^>]+-->", "<!-- уровни 1 л: замер 16,5 мм от среза (y 26,5), модель v4 %.1f (y %.2f) -->" % (L, Y(lit))),
    (r'<line x1="[0-9.]+" y1="26.5" x2="[0-9.]+" y2="26.5"', '<line x1="%.2f" y1="26.5" x2="%.2f" y2="26.5"' % (110 - ri(zm), 110 + ri(zm))),
    (r'<line x1="[0-9.]+" y1="24.7" x2="[0-9.]+" y2="24.7"', '<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f"' % (110 - ri(lit), Y(lit), 110 + ri(lit), Y(lit))),
    (r"\(замер\); модель 14,7", "(замер); модель %s" % ("%.1f" % L).replace(".", ",")),
    (r"KCl 740 мл, 829 г;", "KCl %g мл, %d г;" % (FILL, MASS)),
    (r"слой над колодцем 14,4", "слой над колодцем %s" % ("%.1f" % Rr).replace(".", ",")),
    (r"колодец: Ø79 внутри, Ø85 снаружи", "колодец — конус: Ø79→Ø%.0f внутри, Ø85→Ø83 снаружи" % (2 * (rw(ZWI) - W))),
]
for pat, new in subs:
    t, n = re.subn(pat, new, t)
    if n != 1: raise SystemExit("ОТКАЗ: «%s» заменён %d раз" % (pat, n))
open(f, "w", encoding="utf-8", newline="\n").write(t)
print("KCl %g мл: верх z=%.2f, %.1f мм ниже среза, слой над колодцем %.2f; колодец изнутри у потолка Ø%.2f; 1 л модель y=%.2f" % (FILL, top, ZT - top, Rr, 2 * (rw(ZWI) - W), Y(lit)))
