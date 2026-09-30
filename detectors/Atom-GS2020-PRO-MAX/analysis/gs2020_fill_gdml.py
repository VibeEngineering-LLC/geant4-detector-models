# -*- coding: utf-8 -*-
r"""#GS-24: уровень и плотность пробы в сосуде Маринелли под объём и массу (оператор 28.09: KCl 740 мл, 829 г).
Кольцо r 43–70,5 от z −52,5 до 3,0 мм + диск r 0–70,5 до верха; над пробой воздух до 52,5; плотность = масса/объём.
Состав — отдельно, gs2020_matrix_variant.py. Запуск: python gs2020_fill_gdml.py <in_v2.gdml> <out.gdml> <мл> <г>"""
import math, re, sys
sys.stdout.reconfigure(encoding="utf-8")
src, dst, V, M = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
ring = math.pi * (70.5**2 - 43**2) * 55.5 / 1000.0
if not ring < V <= ring + math.pi * 70.5**2 * 49.5 / 1000.0:
    raise SystemExit("ОТКАЗ: объём %.1f мл вне схемы «кольцо %.1f мл + диск»" % (V, ring))
top = 3.0 + (V - ring) * 1000.0 / (math.pi * 70.5**2)
t = open(src, encoding="utf-8").read()
pairs = [('<zplane z="42.5" rmin="0" rmax="70.5"/></polycone>', '<zplane z="%.3f" rmin="0" rmax="70.5"/></polycone>' % top),
         ('lunit="mm"><zplane z="42.5" rmin="0"', 'lunit="mm"><zplane z="%.3f" rmin="0"' % top)]
for old, new in pairs:
    if t.count(old) != 1: raise SystemExit("ОТКАЗ: фрагмент найден %d раз: %s" % (t.count(old), old))
    t = t.replace(old, new)
t, n = re.subn(r'(<material name="Epoxy_crumb" state="solid"><D value=")[\d.]+', r'\g<1>%.4f' % (M / V), t)
if n != 1: raise SystemExit("ОТКАЗ: плотность Epoxy_crumb найдена %d раз" % n)
open(dst, "w", encoding="utf-8", newline="").write(t)
print("кольцо %.1f мл, верх пробы z=%.3f мм, плотность %.4f г/см3 -> %s" % (ring, top, M / V, dst))
