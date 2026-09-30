# -*- coding: utf-8 -*-
"""#GS-20: копия GDML с другим элементным составом Epoxy_crumb (плотность та же).
Запуск: python gs2020_matrix_variant.py <in.gdml> <out.gdml> "G4_C:0.63513,G4_H:0.06345,..." """
import sys, re
sys.stdout.reconfigure(encoding="utf-8")
src, dst, spec = sys.argv[1:4]
fr = [(k, float(v)) for k, v in (p.split(":") for p in spec.split(","))]
if abs(sum(v for _, v in fr) - 1.0) > 1e-6:
    raise SystemExit("ОТКАЗ: сумма долей %.6f != 1" % sum(v for _, v in fr))
t = open(src, encoding="utf-8").read()
pat = r'(<material name="Epoxy_crumb" state="solid"><D value="[\d.]+" unit="g/cm3"/>)(.*?)(</material>)'
body = "".join('<fraction n="%.6f" ref="%s"/>' % (v, k) for k, v in fr)
t2, n = re.subn(pat, lambda m: m.group(1) + body + m.group(3), t, flags=re.S)
if n != 1:
    raise SystemExit("ОТКАЗ: Epoxy_crumb найден %d раз" % n)
open(dst, "w", encoding="utf-8").write(t2)
print("записан", dst, "состав:", spec)
