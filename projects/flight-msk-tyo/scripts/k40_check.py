"""Пересчёт активности K-40 в слое людей из размеров cabin_geom.hh (независимо от C++). Печатает расхождение с meta prod3."""
import re
h = open("src/cabin_geom.hh", encoding="utf-8").read()
g = lambda n: float(re.search(r"\b" + n + r"\s*=\s*([0-9.eE+-]+)", h).group(1))
vol = 2 * g("paxHalfX") * 2 * g("halfLen") * (g("paxTopZ") - g("floorTop")) - 3.141592653589793 * g("tubeR") ** 2 * 2 * g("tubeHalfL")
act = 54.8 * g("paxRho") * vol / 1000.0
m = float([l for l in open("results/stage1/prod3/k40_0.psp.meta", encoding="utf-8", errors="replace").read().replace("\x00", "").splitlines() if l.startswith("k40_activity_Bq")][0].split("=")[1])
print(f"пересчёт={act:.1f} Бк, meta={m:.1f} Бк, расхождение={act-m:+.2f} Бк ({100*(act-m)/m:+.4f} %)")
raise SystemExit(0 if abs(act - m) / m < 1e-3 else 1)
