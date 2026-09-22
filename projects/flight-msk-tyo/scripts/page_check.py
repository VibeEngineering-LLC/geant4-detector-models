"""Проверка data.js: ряд total равен свёрнутому спектру (независимо от make_page_data). Использование: page_check.py <smear_all.csv> <data.js>. Код 1 при расхождении > 4e-4 (округление до 4 значащих)."""
import sys, json, numpy as np
s = open(sys.argv[2], encoding="utf-8").read()
j = json.loads(s[s.index("=") + 1:].rstrip().rstrip(";"))
tot = [x for x in j["series"] if x["id"] == "total"][0]["y"]
d = np.genfromtxt(sys.argv[1], delimiter=",", skip_header=1)
bad = 0
for E in (30.5, 61.5, 511.5, 1460.5, 2223.5, 8900.5):
    i = int(round((E - j["E0"]) / j["dE"])); k = int(np.argmin(np.abs(d[:, 0] - E)))
    rel = abs(tot[i] - d[k, 1]) / d[k, 1]
    bad += rel > 4e-4
    print(f"E={E}: страница={tot[i]:.5g} пересчёт={d[k,1]:.5g} отн.разн.={rel:.1e}")
print("линий на странице:", len(j["lines"]), "серий:", len(j["series"]))
sys.exit(1 if bad else 0)
