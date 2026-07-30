# -*- coding: utf-8 -*-
"""Сверка stamp.source_sha1() с provenance.cmake — на КАЖДОМ отпечатке.

ЗАЧЕМ. `stamp.source_sha1()` повторяет алгоритм `provenance.cmake` буква в
букву (см. докстринг там же): `file(SHA1 f)` по каждому файлу, накопление
строк «имя:сумма\\n» в заданном порядке, затем `string(SHA1 ...)` по
накопленному и первые 12 знаков. Две независимые реализации одного правила —
ровно тот класс дефекта, которым этот репозиторий уже наказан пять раз
(method-rules §1): они расходятся молча при первой же правке одной стороны.
Тест запускает НАСТОЯЩИЙ `provenance.cmake` через `cmake -P` на каждом наборе
из `stamp.SRC_LISTS` и сравнивает с тем, что даёт питон, на РЕАЛЬНОМ дереве
исходников — не на синтетических данных.

Требует `cmake` в PATH или переменную окружения G4MODELS_CMAKE с путём к
cmake.exe (Visual Studio Build Tools кладёt его вне обычного PATH — см.
g4setup.ps1). Без cmake тест сообщает об этом и завершается кодом 2
(«нельзя судить»), а не падает и не молчит.

Запуск:  python common/py/test_stamp_cmake_match.py
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths  # noqa: E402
import stamp  # noqa: E402

# набор -> (каталог geometry/, макро-префикс — должен совпадать с CMakeLists)
GEOMETRY_DIRS = {
    "Gamma-1S": ("Gamma-1S", "G1S"),
    "Gamma-1S-mucalc": ("Gamma-1S", "G1SMU"),
    "RadiaCode-103": ("RadiaCode-103", "RC"),
    "RadiaCode-103-mucalc": ("RadiaCode-103", "RCMU"),
    "RadiaCode-103-wallfield": ("RadiaCode-103", "RCWF"),
}


def find_cmake():
    env = os.environ.get("G4MODELS_CMAKE")
    if env and os.path.exists(env):
        return env
    found = shutil.which("cmake")
    if found:
        return found
    return None


def cmake_sha1(cmake_exe, geometry_dir, prefix, names):
    prov = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "cmake", "provenance.cmake")
    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, "out.hh")
        src_list = ";".join(os.path.join(geometry_dir, n) for n in names)
        r = subprocess.run(
            [cmake_exe, "-DSRC_DIR=%s" % geometry_dir, "-DOUT=%s" % out,
             "-DPREFIX=%s" % prefix, "-DSRC_LIST=%s" % src_list,
             "-P", prov],
            capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError("cmake -P упал: %s" % (r.stderr or r.stdout))
        text = open(out, encoding="utf-8").read()
        m = re.search(r'#define %s_SRC_SHA1 "([0-9a-f]+)"' % prefix, text)
        if not m:
            raise RuntimeError("в выходе cmake нет %s_SRC_SHA1: %s"
                               % (prefix, text[:200]))
        return m.group(1)


def main():
    cmake_exe = find_cmake()
    if not cmake_exe:
        print("?? cmake не найден (PATH и G4MODELS_CMAKE пусты) — "
              "судить нельзя, сверка не запущена.")
        return 2

    fails = []
    for label, names in stamp.SRC_LISTS.items():
        detector, prefix = GEOMETRY_DIRS[label]
        gdir = str(paths.geometry(detector))
        py_sha = stamp.source_sha1(gdir, names)
        try:
            cm_sha = cmake_sha1(cmake_exe, gdir, prefix, names)
        except RuntimeError as e:
            fails.append((label, "cmake: %s" % e))
            continue
        ok = py_sha == cm_sha
        print("%-28s python=%s cmake=%s %s"
              % (label, py_sha, cm_sha, "OK" if ok else "!! РАСХОЖДЕНИЕ"))
        if not ok:
            fails.append((label, "python=%s != cmake=%s" % (py_sha, cm_sha)))

    if fails:
        print("\nПРОВАЛ: %d из %d наборов разошлись:" % (len(fails), len(stamp.SRC_LISTS)))
        for label, why in fails:
            print("  %s: %s" % (label, why))
        return 1
    print("\nВСЕ %d наборов совпали." % len(stamp.SRC_LISTS))
    return 0


if __name__ == "__main__":
    sys.exit(main())
