# #GS-68: сводка перебора run_gs68_scan.sh по логам. Сортировка по форме ±3 ПШПВ (с крыльями); рядом ±1,5, A/паспорт, χ²/ν.
import sys, re, glob, os
sys.stdout.reconfigure(encoding="utf-8")
D = sys.argv[1] if len(sys.argv) > 1 else "C:/g4work/gs2020/npsm68/scan"; rows = []
for p in glob.glob(os.path.join(D, "*.log")):
    t = open(p, encoding="utf-8", errors="replace").read()
    a = re.search(r"ФОРМА ПИКОВ \(#SHAPE-1[^:]*\): всего ([\d.]+)", t); b = re.search(r"ФОРМА ПИКОВ ±3[^:]*\): всего ([\d.]+); (.*)", t)
    c = re.search(r"ЦЕПОЧКА В РАВНОВЕСИИ.*отношение ([\d.]+); χ²/ν ([\d.]+)", t)
    if not (a and b and c): rows.append((9e9, os.path.basename(p), "НЕТ ИТОГА", "", "", "", "")); continue
    pk = dict(re.findall(r"([\d.]+) → ([\d.]+)", b.group(2)))
    rows.append((float(b.group(1)), os.path.basename(p)[:-4], a.group(1), b.group(1), c.group(1), c.group(2),
                 " ".join(f"{k.split('.')[0]}:{v}" for k, v in pk.items())))
rows.sort()
print("вариант | форма ±1,5 | форма ±3 | A/паспорт | χ²/ν | ±3 по пикам")
for r in rows[: int(os.environ.get("TOP", "40"))]: print(" | ".join(r[1:]))
print(f"всего логов {len(rows)}, без итога {sum(1 for r in rows if r[0] == 9e9)}")
