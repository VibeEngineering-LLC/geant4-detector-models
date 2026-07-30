# -*- coding: utf-8 -*-
"""Нет исходника, влияющего на прогон, но отсутствующего в ЛЮБОМ отпечатке.

ЗАЧЕМ. Второй, встречный аудит уточнил первую находку задачи 133: расхождение
между `stamp.SRC_LISTS` (питон) и списком в CMakeLists — это ГРОМКИЙ дефект,
он сразу даёт `stale` на реальном прогоне, тихого совпадения там не бывает
(за это отвечает `test_stamp_cmake_match.py`). Настоящая дыра — СОГЛАСНОЕ
УМОЛЧАНИЕ: файл, который реально компилируется в один из бинарников
детектора, но не входит НИ В ОДИН список из `SRC_LISTS`. Тогда обе стороны
(питон и cmake) согласны между собой, вердикт `ok`, а входы фактически
устарели — правку в этом файле штамп не увидит вовсе.

ЧТО ПРОВЕРЯЕТСЯ. Для каждого каталога `geometry/` — что множество
`*.cc`/`*.hh`, лежащих в нём, ЦЕЛИКОМ покрыто объединением записей
`SRC_LISTS`, чей каталог — этот. Проверка «файл есть на диске» простая и
грубая: она не знает, какой файл в какой .exe линкуется, поэтому требует
присутствия в union, а не в каждой отдельной цели. Это нижняя граница, не
доказательство полноты (provenance.cmake сам, макро-файлы .mac, будущие
зависимости уровня выше .cc/.hh проверкой не покрыты — см. задачу 133 п.5
про генератор сборки), но именно класс «новый .cc забыли добавить в список»
она ловит гарантированно и автоматически при каждом запуске.

Запуск:  python common/py/test_stamp_coverage.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths  # noqa: E402
import stamp  # noqa: E402

# набор -> каталог detectors/<...>/geometry, откуда он берёт файлы
LABEL_DETECTOR = {
    "Gamma-1S": "Gamma-1S",
    "Gamma-1S-mucalc": "Gamma-1S",
    "RadiaCode-103": "RadiaCode-103",
    "RadiaCode-103-mucalc": "RadiaCode-103",
    "RadiaCode-103-wallfield": "RadiaCode-103",
}


def main():
    by_detector = {}
    for label, names in stamp.SRC_LISTS.items():
        det = LABEL_DETECTOR.get(label)
        if det is None:
            raise SystemExit(
                "test_stamp_coverage: набор %r из stamp.SRC_LISTS не "
                "приписан ни к одному каталогу geometry/ — впишите его в "
                "LABEL_DETECTOR, иначе проверка покрытия для него не "
                "выполняется." % label)
        by_detector.setdefault(det, set()).update(names)

    fails = []
    for det, covered in sorted(by_detector.items()):
        gdir = paths.geometry(det)
        on_disk = {p.name for p in gdir.glob("*.cc")} | \
                  {p.name for p in gdir.glob("*.hh")}
        missing = sorted(on_disk - covered)
        extra = sorted(covered - on_disk)
        print("%-16s на диске %2d, в отпечатках %2d%s"
              % (det, len(on_disk), len(covered),
                 "" if not missing else "  !! БЕЗ ОТПЕЧАТКА: %s" % missing))
        if extra:
            print("   (в SRC_LISTS есть, файла на диске нет — устаревшая "
                  "запись, не блокер): %s" % extra)
        if missing:
            fails.append((det, missing))

    if fails:
        print("\nПРОВАЛ: файлы вне любого отпечатка провенанса —")
        for det, missing in fails:
            print("  %s: %s" % (det, ", ".join(missing)))
        print("Правка каждого из них НЕ изменит ни один печатаемый src_sha1 "
              "— согласное умолчание из задачи 133, п.2.")
        return 1
    print("\nПокрытие полное: каждый .cc/.hh входит хотя бы в один отпечаток.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
