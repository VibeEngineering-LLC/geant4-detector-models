# -*- coding: utf-8 -*-
"""Штамп провенанса и объявление наблюдаемой — ОДНА реализация на репозиторий.

ЗАЧЕМ ДВЕ РАЗНЫЕ ЗАЩИТЫ. Они отвечают на разные вопросы, и ни одна не заменяет
другую:

* **провенанс** — «ИЗ ЧЕГО получено это число»: каким отпечатком исходников
  посчитаны входные спектры, совпадает ли он у всех входов между собой и с
  текущим деревом. Ловит устаревшие входы (три случая за линию Гамма-1С:
  моно-спектры моста, exe в каталоге прогонов, точечные сетки против правок
  торца 29.07);
* **объявление наблюдаемой** — «ЧТО ИМЕННО за число»: пик или полный счёт,
  какое окно, вычтена ли полка и каким окном, размыт спектр или нет. Ловит
  сравнение двух величин, посчитанных по разным правилам (пять случаев за один
  вечер 30.07; см. common/docs/method-rules.md §5).

Штамп провенанса на второй вопрос НЕ отвечает: три расходящиеся таблицы одной
величины могут быть все три свежими.

ПОЧЕМУ ОТПЕЧАТОК, А НЕ mtime. `git checkout`, клон и синхронизация облака
сбрасывают времена в произвольную сторону, а защита, ложно срабатывающая от
`git pull`, будет отключена в первый же день. Отпечаток же приходит из самого
бинарника, считавшего спектр (geometry/provenance.cmake запекает его в exe, а
main.cc печатает в шапку каждого выходного файла) — подделать его расчётом
нельзя.

ГРАММАТИКА. Строки штампа в таблицах results/ — комментарии вида

    #@ ключ = значение

Плоская, машинно проверяемая (tools/check_stamp.py), и не мешает читателям,
пропускающим строки на «#». Ключи с точкой группируются по смыслу:
`src.*` — провенанс, `obs.*` — наблюдаемая.
"""
import hashlib
import os
import subprocess

# Файлы, по которым считается отпечаток. Порядок ЗНАЧИМ: он должен совпадать с
# G1S_PROV_SRC в geometry/CMakeLists.txt, иначе питонная и cmake-сторона дадут
# разные суммы на одном и том же дереве.
SRC_LISTS = {
    "Gamma-1S": ("main.cc", "G1SDetector.cc", "G1SDetector.hh"),
}

NO_STAMP = "БЕЗ-ШТАМПА"


def _sha1_file(path):
    h = hashlib.sha1()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def source_sha1(geometry_dir, names):
    """Отпечаток дерева исходников — тем же способом, что provenance.cmake.

    cmake: `file(SHA1 f)` по каждому файлу, накопление строк «имя:сумма\\n»,
    затем `string(SHA1 ...)` по накопленному и первые 12 знаков. Повторено
    здесь буква в букву; расхождение реализаций сделало бы сторожа генератором
    ложных тревог, а это ровно тот дефект, от которого штамп и заводится.
    """
    acc = "".join("%s:%s\n" % (n, _sha1_file(os.path.join(geometry_dir, n)))
                  for n in names)
    return hashlib.sha1(acc.encode("utf-8")).hexdigest()[:12]


def read_run_stamp(path):
    """{ключ: значение} из шапки расчётного спектра; пустой dict, если нет.

    Пустой результат — законное состояние: спектры, посчитанные до внедрения
    штампа, читаются по-прежнему. Отсутствие штампа НЕ равно совпадению —
    вызывающий обязан различать «сошлось» и «неизвестно», иначе старые входы
    молча пройдут проверку.
    """
    out = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if not ln.startswith("#"):
                break
            body = ln[1:].strip()
            if "=" not in body:
                continue
            k, v = body.split("=", 1)
            k = k.strip()
            if k in ("src_sha1", "git_describe", "build"):
                out[k] = v.strip()
    return out


def check_inputs(paths, geometry_dir=None, names=None):
    """Свести отпечатки входов и (если задано) сравнить с текущим деревом.

    Возвращает (verdict, detail):
      verdict = "ok"        — все входы со штампом, отпечаток один и совпал с деревом;
                "stale"     — отпечаток входов НЕ совпал с деревом;
                "mixed"     — входы посчитаны РАЗНЫМИ отпечатками (худший случай:
                              смесь геометрий внутри одной таблицы);
                "unstamped" — есть входы без штампа, судить нельзя.
    Три разных вердикта вместо булева «годно» — намеренно: «смесь геометрий» и
    «все входы одинаково устарели» лечатся по-разному, а «нет штампа» вообще не
    вывод, а отсутствие вывода.
    """
    seen, missing = {}, []
    for p in paths:
        st = read_run_stamp(p)
        sha = st.get("src_sha1")
        if not sha or sha == NO_STAMP:
            missing.append(os.path.basename(p))
        else:
            seen.setdefault(sha, []).append(os.path.basename(p))
    if len(seen) > 1:
        return "mixed", {"shas": {k: v[:4] for k, v in seen.items()},
                         "unstamped": missing[:4]}
    if missing:
        return "unstamped", {"n": len(missing), "examples": missing[:4],
                             "shas": list(seen)}
    if not seen:
        return "unstamped", {"n": 0, "examples": [], "shas": []}
    got = next(iter(seen))
    if geometry_dir and names:
        cur = source_sha1(geometry_dir, names)
        if cur != got:
            return "stale", {"inputs": got, "tree": cur}
        return "ok", {"sha": got}
    return "ok", {"sha": got, "tree": "не сверялось"}


def git_describe(repo_dir):
    try:
        r = subprocess.run(["git", "-C", repo_dir, "describe", "--always",
                            "--dirty"], capture_output=True, text=True)
        return r.stdout.strip() or NO_STAMP
    except OSError:
        return NO_STAMP


def lines(script, observable, inputs=None, geometry_dir=None, names=None,
          repo_dir=None):
    """Строки `#@ ключ = значение` для шапки таблицы results/.

    script     — имя скрипта-производителя (для «кто это написал»);
    observable — dict объявления наблюдаемой; ОБЯЗАТЕЛЬНЫЕ ключи ниже;
    inputs     — пути расчётных спектров, из которых получена таблица.

    Обязательные ключи наблюдаемой перечислены явно и проверяются: правило,
    которое можно молча не выполнить, не выполняется. Значения — свободный
    текст без запятых (их запрещает csvio, и штамп не должен становиться
    единственным местом, где запятая просачивается в файл).
    """
    need = ("quantity", "area", "window", "shelf", "blurred")
    miss = [k for k in need if k not in observable]
    if miss:
        raise SystemExit(
            "stamp.lines(%s): не объявлены обязательные ключи наблюдаемой: %s.\n"
            "Без них таблица непригодна для сравнения с другой таблицей — "
            "именно это правило нарушалось пять раз за вечер 30.07.2026."
            % (script, ", ".join(miss)))
    out = ["#@ stamp.version = 1", "#@ src.script = %s" % script]
    if repo_dir:
        out.append("#@ src.git = %s" % git_describe(repo_dir))
    if geometry_dir and names:
        out.append("#@ src.tree_sha1 = %s"
                   % source_sha1(geometry_dir, names))
    if inputs:
        verdict, detail = check_inputs(inputs, geometry_dir, names)
        out.append("#@ src.inputs_verdict = %s" % verdict)
        out.append("#@ src.inputs_n = %d" % len(inputs))
        if verdict == "mixed":
            out.append("#@ src.inputs_sha1 = %s"
                       % " | ".join(sorted(detail["shas"])))
        elif verdict == "stale":
            out.append("#@ src.inputs_sha1 = %s" % detail["inputs"])
        elif verdict == "unstamped":
            out.append("#@ src.inputs_sha1 = %s"
                       % (" | ".join(detail["shas"]) or NO_STAMP))
            out.append("#@ src.inputs_unstamped = %d" % detail["n"])
        else:
            out.append("#@ src.inputs_sha1 = %s" % detail["sha"])
    for k in need:
        out.append("#@ obs.%s = %s" % (k, observable[k]))
    for k in sorted(set(observable) - set(need)):
        out.append("#@ obs.%s = %s" % (k, observable[k]))
    bad = [ln for ln in out if "," in ln]
    if bad:
        raise SystemExit(
            "stamp.lines(%s): запятая в штампе — %r.\n"
            "Значения штампа идут в тот же файл, что и данные; замените на "
            "точку с запятой." % (script, bad[0]))
    return out


def read_table_stamp(path):
    """{ключ: значение} из строк `#@` таблицы results/ — для сторожа и отчётов."""
    out = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if not ln.startswith("#"):
                break
            if not ln.startswith("#@"):
                continue
            body = ln[2:].strip()
            if "=" in body:
                k, v = body.split("=", 1)
                out[k.strip()] = v.strip()
    return out
