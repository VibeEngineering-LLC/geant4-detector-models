"""Опись комплекта поверки: что лежит в каждом из 40 XML.

Читает заголовки всех записей и печатает таблицу: геометрия, источник, время,
мёртвое время, масса, объём, комментарий с паспортной активностью. Это база
для пересчёта: дальше по каждой записи считаются площади пиков.

Комментарий ЛСРМ хранит ОДНУ строку, поэтому у смесей там назван лишь один
нуклид из четырёх — состав смесей определяется отдельно, по самим спектрам
(kit_mixture.py).
"""
import glob
import os
import re
import sys

# Корни путей — из переменных окружения (common/py/paths.py), чтобы в коде не
# было ни одного пути, привязанного к конкретной машине.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import becqmoni as bm  # noqa: E402

KIT = str(paths.ref("Gamma-1S"))
GEOM_ORDER = ["Marinelli_1L", "Denta_120mL", "Petri_60mL", "Point_5cm", "Point_25cm"]


def header(path):
    t = open(path, encoding="utf-8", errors="replace").read()
    head = t[:t.find("<Spectrum>")] if "<Spectrum>" in t else t[:4000]
    g = lambda k: (re.search(r"<%s>(.*?)</%s>" % (k, k), head, re.S).group(1).strip()
                   if re.search(r"<%s>(.*?)</%s>" % (k, k), head, re.S) else "")
    note = re.search(r"<Note>(.*?)</Note>", t, re.S)
    return {
        "name": g("Name"),
        "time": g("StartTime")[:10],
        "weight": g("Weight"),
        "volume": g("Volume"),
        "note": (note.group(1).strip().replace("\n", " ") if note else ""),
    }


if __name__ == "__main__":
    rows = []
    for geom in GEOM_ORDER:
        d = os.path.join(KIT, geom)
        if not os.path.isdir(d):
            continue
        for p in sorted(glob.glob(os.path.join(d, "*.xml"))):
            s, b = bm.read(p)
            h = header(p)
            rows.append((geom, h, s, b, os.path.basename(p)))

    print("Комплект поверки Гамма-1С: %d записей\n" % len(rows))
    print("%-13s %-30s %-11s %8s %6s %7s %6s" %
          ("геометрия", "источник", "дата", "живое,с", "мёрт%", "масса,г", "фон,с"))
    for geom, h, s, b, fn in rows:
        dead = 100 * (1 - s.live / s.real) if s.real else 0
        try:
            w = float(h["weight"]) * 1000
        except ValueError:
            w = 0.0
        print("%-13s %-30s %-11s %8.0f %6.2f %7.0f %6.0f" %
              (geom, h["name"][:30], h["time"], s.live, dead, w,
               b.live if b else 0))

    print("\n--- паспортные активности из поля Note ---")
    for geom, h, s, b, fn in rows:
        if h["note"]:
            n = re.sub(r"\(Sum/T[^)]*\)", "", h["note"]).strip()
            print("%-13s %-26s %s" % (geom, h["name"][:26], n[:90]))
