# -*- coding: utf-8 -*-
"""Сборка плакатной страницы: src/index.html + out/poster_data.json -> dist/index.html.

Данные встраиваются в тег <script type="application/json"> вместо маркера
/*DATA*/. Страница получается самодостаточной: ни сети, ни внешних файлов.
"""
import os, sys, io, json
sys.stdout.reconfigure(encoding="utf-8")

HERE = os.path.abspath(os.path.dirname(__file__))
SRC = os.path.join(HERE, "src", "index.html")
DATA = os.path.join(HERE, "..", "analysis", "out", "poster_data.json")
CALDATA = os.path.join(HERE, "..", "analysis", "out", "calibration_data.json")
SHDATA = os.path.join(HERE, "..", "analysis", "out", "shield_data.json")
DIST = os.path.join(HERE, "dist")
MARK = "/*DATA*/"
MARK_CAL = "/*CALDATA*/"
MARK_SH = "/*SHIELDDATA*/"


def main():
    html = io.open(SRC, encoding="utf-8").read()
    out = html
    total = 0
    for mark, path, what in ((MARK, DATA, "разложение"),
                             (MARK_CAL, CALDATA, "калибровка"),
                             (MARK_SH, SHDATA, "свинцовый домик")):
        raw = io.open(path, encoding="utf-8").read()
        json.loads(raw)                  # данные обязаны быть валидны до вставки
        if mark not in out:
            raise SystemExit("в шаблоне нет маркера %s" % mark)
        # </script> внутри данных закрыл бы тег раньше времени; в JSON его быть
        # не должно, но проверка дешевле разбирательства с пустой страницей.
        if "</script" in raw.lower():
            raise SystemExit("в данных %s есть </script>" % what)
        out = out.replace(mark, raw)
        total += len(raw)
        print("встроено (%s): %d байт" % (what, len(raw)))
    os.makedirs(DIST, exist_ok=True)
    dst = os.path.join(DIST, "index.html")
    io.open(dst, "w", encoding="utf-8", newline="\n").write(out)
    print("собрано:", os.path.abspath(dst), os.path.getsize(dst), "байт")
    print("данных всего:", total, "байт")


if __name__ == "__main__":
    main()
