"""Сплошная ревизия расчётных спектров: чем посчитан и годен ли розыгрыш (R62).

ПОВОД. Разбор R61: `/gps/pos/confine Sample` стоит в большинстве макросов, а
том `Sample` строится ТОЛЬКО в режимах `vessel*`. Запуск такого макроса в
режиме `shield`/`open`/`bare` Geant4 не останавливает: он печатает
предупреждение и МОЛЧА снимает ограничение, после чего источник разыгрывается
по всему заданному телу — включая кристалл. Восемь шаблонов разложения по
нуклидам так и были посчитаны, и заметили это только по «лишнему» горбу на
графике, спустя недели.

Сторож в модели с тех пор есть (main.cc, CheckSourcePlacement — падение на
первой же вершине внутри кристалла, поле `# src_in_crystal` в шапке). Но он
защищает БУДУЩИЕ прогоны. Этот скрипт разбирает УЖЕ ПОСЧИТАННОЕ.

ЧТО ПРОВЕРЯЕТСЯ по каждому CSV:
  1. режим прогона (`# mode`) против тома в `/gps/pos/confine` того макроса,
     которым он посчитан (`# run_args`);
  2. поле `# src_in_crystal` — есть ли и равно ли нулю;
  3. отпечаток исходников (`# src_sha1`) против текущего дерева geometry/ —
     то есть посчитан ли спектр НЫНЕШНЕЙ моделью;
  4. участвует ли файл в опубликованном результате (упоминание имени в
     analysis/, results/, web-th232/, docs/).

Вердикт по файлу — худший из применимых:
  негоден      confine назван, тома в этом режиме нет: розыгрыш ушёл шире пробы
  в кристалле  `src_in_crystal` не ноль
  устарел      посчитан не текущей ревизией исходников
  неизвестно   нет `run_args` или нет `src_in_crystal` — судить не по чему
  годен        всё сошлось

«Неизвестно» — не «годен». Отсутствие поля означает сборку до введения
сторожа, а такие прогоны как раз и подозрительны.

    python detectors/Gamma-1S/analysis/audit_runs.py [--csv отчёт.csv]
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402
import stamp  # noqa: E402

DET = "Gamma-1S"
# Те же файлы, по которым отпечаток считает provenance.cmake при сборке.
SRC_NAMES = ("main.cc", "G1SDetector.cc", "G1SDetector.hh")

# Какие ИМЕНОВАННЫЕ тома существуют в каком режиме. Список короткий намеренно:
# в confine осмысленно попадает только объём пробы, всё остальное — ошибка
# постановки, и её надо увидеть, а не молча разрешить.
#   bare    — только детектор
#   open    — детектор и защита с открытой крышкой
#   shield  — детектор и защита
#   vessel* — то же плюс сосуд с пробой (том `Sample`)
def volume_exists(vol, mode):
    if vol == "Sample":
        return mode.startswith("vessel")
    return None            # неизвестный том — судить не берёмся


def read_header(path):
    """{ключ: значение} из шапки `# ключ = значение`; строки до первых данных."""
    out = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if not ln.startswith("#"):
                break
            body = ln[1:].strip()
            if "=" in body:
                k, v = body.split("=", 1)
                out[k.strip()] = v.strip()
    return out


def parent_of(name):
    """Спутник `X_emit.csv` / `X_chan.csv` -> `X.csv`, иначе None.

    Спутники пишутся тем же прогоном, что и спектр, но своей шапки с режимом и
    аргументами не несут. Судить их отдельно нельзя: без наследования они все
    до одного попадают в «неизвестно» и топят в шуме те файлы, где вопрос
    действительно открыт.
    """
    # Список обязан совпадать с тем, что РЕАЛЬНО пишет main.cc
    # (EndOfRunAction, пять fopen: сам fOut плюс четыре спутника) и с
    # SUFFIXES в web-th232/export_data.py:_grid_main_csvs. Расхождение уже
    # ловилось дважды: _chan.csv (R45) и _shield.csv (R69) — оба раза
    # спутник подходил под шаблон glob основного файла и разбирался как
    # двухколоночный спектр.
    for suf in ("_emit.csv", "_emitx.csv", "_chan.csv", "_shield.csv"):
        if name.endswith(suf):
            return name[:-len(suf)] + ".csv"
    return None


def macro_confine(path):
    """Тома из всех `/gps/pos/confine` макроса, в порядке появления."""
    vols = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            m = re.match(r"\s*/gps/pos/confine\s+(\S+)", ln)
            if m:
                vols.append(m.group(1))
    return vols


def published_names(repo):
    """Текст всех скриптов и результатов одной строкой — для поиска имён файлов.

    Читается один раз: 183 файла против нескольких сотен исходников — это
    десятки тысяч открытий, если искать по каждому имени отдельно.
    """
    buf = []
    roots = [os.path.join(repo, "detectors", DET, d)
             for d in ("analysis", "results", "web-th232", "drivers", "docs")]
    for root in roots:
        for dirpath, _, files in os.walk(root):
            if "__pycache__" in dirpath:
                continue
            for f in files:
                if os.path.splitext(f)[1].lower() not in (
                        ".py", ".md", ".json", ".csv", ".js", ".html", ".mac"):
                    continue
                p = os.path.join(dirpath, f)
                try:
                    with open(p, encoding="utf-8", errors="replace") as fh:
                        buf.append(fh.read())
                except OSError:
                    pass
    return "\n".join(buf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", help="куда сложить построчный отчёт")
    a = ap.parse_args()

    build = str(paths.build(DET))
    geom = str(paths.geometry(DET))
    tree_sha = stamp.source_sha1(geom, SRC_NAMES)
    print("каталог прогонов: %s" % build)
    print("отпечаток текущего дерева geometry/: %s" % tree_sha)

    macros = {}
    for f in os.listdir(build):
        if f.endswith(".mac"):
            macros[f] = macro_confine(os.path.join(build, f))

    pub = published_names(paths.REPO)

    rows, tally = [], {}
    for dirpath, _, files in os.walk(build):
        for f in sorted(files):
            if not f.endswith(".csv"):
                continue
            p = os.path.join(dirpath, f)
            h = read_header(p)
            par = parent_of(f)
            if par and not h.get("run_args"):
                pp = os.path.join(dirpath, par)
                if os.path.exists(pp):
                    ph = read_header(pp)
                    for k in ("mode", "run_args", "src_in_crystal", "src_sha1"):
                        if k in ph:
                            h.setdefault(k, ph[k])
            mode = h.get("mode", "")
            args = h.get("run_args", "")
            mac = args.split()[0] if args else ""
            sic = h.get("src_in_crystal")
            sha = h.get("src_sha1", "")

            vols = macros.get(mac, [])
            confine = vols[0] if vols else ""
            bad_confine = False
            if confine and mode:
                ex = volume_exists(confine, mode)
                bad_confine = (ex is False)

            if bad_confine:
                verdict = "негоден"
            elif sic not in (None, "0"):
                verdict = "в кристалле"
            elif not mac or sic is None:
                verdict = "неизвестно"
            elif sha != tree_sha:
                verdict = "устарел"
            else:
                verdict = "годен"

            rel = os.path.relpath(p, build).replace("\\", "/")
            rows.append((rel, verdict, mode, mac, confine, sic or "-",
                         sha or "-", "да" if f in pub else "нет"))
            tally[verdict] = tally.get(verdict, 0) + 1

    print("\nвсего файлов: %d" % len(rows))
    for k in ("негоден", "в кристалле", "устарел", "неизвестно", "годен"):
        if tally.get(k):
            print("  %-12s %d" % (k, tally[k]))

    # Печатается только то, что требует решения: годные перечислять незачем.
    for want in ("негоден", "в кристалле", "неизвестно", "устарел"):
        sel = [r for r in rows if r[1] == want]
        if not sel:
            continue
        print("\n=== %s (%d) ===" % (want, len(sel)))
        print("%-38s %-10s %-26s %-7s %s"
              % ("файл", "режим", "макрос", "в публ.", "src_sha1"))
        for r in sel[:200]:
            print("%-38s %-10s %-26s %-7s %s" % (r[0][:38], r[2][:10],
                                                 r[3][:26], r[7], r[6][:12]))
        if len(sel) > 200:
            print("... ещё %d" % (len(sel) - 200))

    if a.csv:
        import csv
        with open(a.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["file", "verdict", "mode", "macro", "confine",
                        "src_in_crystal", "src_sha1", "published"])
            w.writerows(rows)
        print("\nотчёт: %s" % a.csv)


if __name__ == "__main__":
    main()
