"""Скачать ИЗМЕРЕННЫЕ фоновые спектры установки Гамма-1С (репо spectravibe-toolkit).

Берём фон с сосудом Маринелли, залитым водой, и закрытой крышкой — ровно та
геометрия, для которой паспорт задаёт МИА. Файлы .spe ЛСРМ в CP1251.
"""
import base64
import os
import subprocess

REPO = "VibeEngineering-LLC/spectravibe-toolkit"
DIR = "detectors/Gamma-1S/data/averaged_backgrounds"
HERE = os.path.dirname(os.path.abspath(__file__))
NAMES = [
    "bg_2016_marinelli_water_marinelli.spe",
    "bg_2016_marinelli_water_marinelli.provenance.json",
    "bg_2024_marinelli_water_closed_lid_marinelli.spe",
    "bg_2024_marinelli_water_closed_lid_marinelli.provenance.json",
    "MANIFEST.json",
]

for n in NAMES:
    out = subprocess.run(
        ["gh", "api", "repos/%s/contents/%s/%s" % (REPO, DIR, n), "--jq", ".content"],
        capture_output=True, text=True, encoding="utf-8", errors="replace")
    if out.returncode != 0:
        print("!!", n, (out.stderr or "")[:150])
        continue
    raw = base64.b64decode(out.stdout.replace("\n", ""))
    # .spe ЛСРМ — CP1251, но встречаются байты вне таблицы (0x98): заменяем,
    # это только в текстовых подписях, на числа не влияет.
    enc = "utf-8" if n.endswith(".json") else "cp1251"
    txt = raw.decode(enc, errors="replace")
    open(os.path.join(HERE, n), "w", encoding="utf-8").write(txt)
    print("%-58s %7d байт" % (n, len(raw)))
