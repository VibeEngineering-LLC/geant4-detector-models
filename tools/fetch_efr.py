"""Скачать и разобрать .efr/.efa комплекта поверки Гамма-1С (репо spectravibe-toolkit).

Файлы в кодировке CP1251. Формат .efr — секции [детектор;геометрия;источник],
затем строки «энергия=эффективность,погрешность%,нуклид,площадь,dплощадь,...».
"""
import base64
import json
import os
import subprocess
import urllib.parse

REPO = "VibeEngineering-LLC/spectravibe-toolkit"
DIR = "detectors/Gamma-1S/efficiency/Gamma-1S_NaI_63x63_USB_SN-01"
HERE = os.path.dirname(os.path.abspath(__file__))
NAMES = [
    "УДС-ГЦ-63х63-USB__SN-01_-_Точечная-25см.efr",
    "УДС-ГЦ-63х63-USB__SN-01_-_Точечная-5см.efr",
    "УДС-ГЦ-63х63-USB__SN-01_-_Маринелли.efr",
    "УДС-ГЦ-63х63-USB__SN-01_-_Маринелли.efa",
    "УДС-ГЦ-63х63-USB__SN-01_-_Дента.efr",
    "УДС-ГЦ-63х63-USB__SN-01_-_Дента.efa",
    "УДС-ГЦ-63х63-USB__SN-01_-_Петри.efr",
    "УДС-ГЦ-63х63-USB__SN-01_-_Петри.efa",
    "УДС-ГЦ-63х63-USB__SN-01_-_Точечная-25см.efa",
]


def fetch(name):
    p = urllib.parse.quote("%s/%s" % (DIR, name))
    out = subprocess.run(
        ["gh", "api", "repos/%s/contents/%s" % (REPO, p), "--jq", ".content"],
        capture_output=True, text=True, check=True)
    raw = base64.b64decode(out.stdout.replace("\n", ""))
    return raw.decode("cp1251")


def parse_efr(text):
    """-> список секций {'source':…, 'geometry':…, 'points':[(E, eff, dpct, nuclide)]}"""
    secs, cur = [], None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            cur = {"header": line[1:-1], "meta": {}, "points": []}
            secs.append(cur)
            continue
        if cur is None:
            continue
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        parts = v.split(",")
        try:
            E = float(k)
        except ValueError:
            cur["meta"][k] = v
            continue
        try:
            cur["points"].append((E, float(parts[0]), float(parts[1]), parts[2]))
        except (ValueError, IndexError):
            pass
    return secs


if __name__ == "__main__":
    os.makedirs(HERE, exist_ok=True)
    summary = {}
    for n in NAMES:
        try:
            t = fetch(n)
        except subprocess.CalledProcessError as e:
            print("!! не скачан:", n, e.stderr[:120])
            continue
        open(os.path.join(HERE, n), "w", encoding="utf-8").write(t)
        if n.endswith(".efr"):
            secs = parse_efr(t)
            pts = [p for s in secs for p in s["points"]]
            summary[n] = pts
            print("%-55s секций %2d, точек %3d" % (n, len(secs), len(pts)))
        else:
            print("%-55s (.efa, %d строк)" % (n, len(t.splitlines())))
    json.dump(summary, open(os.path.join(HERE, "efr_points.json"), "w",
                            encoding="utf-8"), ensure_ascii=False, indent=1)
