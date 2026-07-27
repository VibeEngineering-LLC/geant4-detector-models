"""Скачать XML-комплект поверки Маринелли 1 л (BecqMoni): образец + вложенный фон.

Почему XML, а не .spe: в .spe ЛСРМ отсчёты упакованы двоично, и любая
перекодировка их портит. В XML BecqMoni спектр лежит текстом.

Почему через дерево git и SHA, а не через contents/<путь>: в путях кириллица,
и URL-кодирование её ломает. Дерево отдаёт и путь, и SHA блоба — качаем блоб.
"""
import base64
import json
import os
import subprocess

REPO = "VibeEngineering-LLC/spectravibe-toolkit"
WANT = "reference_kits_becqmoni/"          # ВЕСЬ комплект: 5 геометрий
HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "kit")


def api(path, jq=None):
    cmd = ["gh", "api", "repos/%s/%s" % (REPO, path)]
    if jq:
        cmd += ["--jq", jq]
    r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                       errors="replace")
    if r.returncode != 0:
        raise SystemExit("gh api %s: %s" % (path, (r.stderr or "")[:200]))
    return r.stdout


if __name__ == "__main__":
    os.makedirs(HERE, exist_ok=True)
    tree = json.loads(api("git/trees/main?recursive=1"))["tree"]
    hits = [t for t in tree
            if t["type"] == "blob" and WANT in t["path"] and t["path"].endswith(".xml")]
    if not hits:
        raise SystemExit("в дереве нет XML по маске " + WANT)
    for t in hits:
        blob = json.loads(api("git/blobs/" + t["sha"]))
        raw = base64.b64decode(blob["content"])
        # раскладываем по геометриям: .../reference_kits_becqmoni/<геом>/<нукл>/<файл>
        tail = t["path"].split(WANT, 1)[1]
        parts = tail.split("/")
        geom = parts[0] if len(parts) > 1 else "_"
        sub = os.path.join(HERE, geom)
        os.makedirs(sub, exist_ok=True)
        open(os.path.join(sub, parts[-1]), "wb").write(raw)
        print("%-18s %-50s %8d" % (geom, parts[-1], len(raw)))
