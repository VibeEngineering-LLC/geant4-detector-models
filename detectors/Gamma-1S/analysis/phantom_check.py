# -*- coding: utf-8 -*-
"""Проверка шаблонов нуклидов на фантомы — энерговыделение там, где его
физически быть не может (R77, директива оператора 09.08.2026).

Повод — R76: сущность K-рентгена показывала структуру выше 500 кэВ, где
рентгена нет вовсе. Дефект был не в физике, а в разметке события: признак
«рентген» ставился, если рентгеновским был ХОТЬ ОДИН вклад в энерговыделение,
и совпадение 75-кэВного рентгена с 2614-кэВной гаммой в окне разрешения
уезжало в рентгеновский слой целиком. Класс дефекта общий, поэтому проверка
здесь сплошная — по всем шаблонам, а не по одному пойманному.

КРИТЕРИЙ — сохранение энергии, не эвристика. В прогоне с nucleusLimits
распадается один нуклид, и энерговыделение в кристалле не может превысить
суммарную энергию квантов, испущенных в ЭТОМ распаде. Верхняя граница
берётся из собственных данных модели: сумма энергий всех линий эмиссионного
спектра (_emit.csv), каждая учтена один раз. Граница заведомо завышена —
испустить все линии сразу один распад не может, разные линии принадлежат
разным ветвям. Тем она и хороша: превышение означает не «редкое совпадение»,
а невозможное событие, то есть дефект.

Второй, более тонкий рубеж — энергия САМОЙ ЖЁСТКОЙ одиночной линии. Выше неё
энерговыделение законно только через каскадное совпадение. Оно реально, но
его доля обязана быть малой и обязана соответствовать схеме распада; здесь
она измеряется и показывается, а не принимается на веру.

Запуск:
    python analysis/phantom_check.py            # каталог сборки из paths
    python analysis/phantom_check.py <каталог>
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths                                              # noqa: E402

# Порог «линия существует»: доля от самой сильной линии эмиссии. Ниже него
# счётчик эмиссии населён единичными отсчётами МК, и брать такие энергии в
# сумму-границу значит поднимать границу шумом.
LINE_REL_MIN = 1e-4
# Порог «энерговыделение есть»: отсчётов в канале депозита. Единичный отсчёт
# на 10^6 розыгрышей — не структура, но и не повод молчать: он показывается
# отдельной колонкой, а не выбрасывается.
DEP_MIN = 1


def read_csv(path):
    """(E_keV, значения по колонкам) из CSV с '#'-шапкой."""
    rows, head = [], None
    with open(path, encoding="utf-8", errors="replace") as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            if head is None:
                head = ln.split(",")
                continue
            p = ln.split(",")
            rows.append([float(x) for x in p])
    if head is None:
        raise ValueError("нет данных: %s" % path)
    cols = {name: [r[i] for r in rows] for i, name in enumerate(head)}
    return head, cols


def emission_bounds(path):
    """(E самой жёсткой линии, сумма энергий всех линий) по _emit.csv."""
    _, c = read_csv(path)
    e, n = c["E_keV"], c["counts"]
    top = max(n) if n else 0.0
    lines = [e[i] for i in range(len(e)) if n[i] > top * LINE_REL_MIN]
    if not lines:
        return 0.0, 0.0
    return max(lines), sum(lines)


def deposit_profile(path, col="counts"):
    _, c = read_csv(path)
    return c["E_keV"], c[col]


def above(e, v, lim):
    """(сумма отсчётов выше lim, максимальная энергия с отсчётом)."""
    s = 0.0
    emax = 0.0
    for i in range(len(e)):
        if v[i] >= DEP_MIN:
            emax = max(emax, e[i])
            if e[i] > lim:
                s += v[i]
    return s, emax


def check(build, key):
    dep = os.path.join(build, "iso_%s.csv" % key)
    emi = os.path.join(build, "iso_%s_emit.csv" % key)
    if not (os.path.exists(dep) and os.path.exists(emi)):
        return None
    e1, cap = emission_bounds(emi)
    e, v = deposit_profile(dep)
    tot = sum(v)
    n_imp, emax = above(e, v, cap)          # выше границы сохранения энергии
    n_sum, _ = above(e, v, e1)              # выше самой жёсткой одиночной линии
    return {
        "nuc": key, "E_line_max": e1, "E_cap": cap, "E_dep_max": emax,
        "total": tot,
        "impossible": n_imp, "impossible_pct": 100.0 * n_imp / max(tot, 1),
        "summing": n_sum, "summing_pct": 100.0 * n_sum / max(tot, 1),
    }


def check_xray(build, keys):
    """Отдельно — сущность K-рентгена: она собирается по всем нуклидам.

    Граница у неё своя и жёсткая: рентген атомной релаксации не бывает выше
    K-края самого тяжёлого атома ряда, поэтому предел — удвоенная (совпадение
    двух рентгеновских квантов в одном событии) максимальная энергия колонки
    x_atomic по всем нуклидам.
    """
    xmax = 0.0
    for k in keys:
        p = os.path.join(build, "iso_%s_emitx.csv" % k)
        if not os.path.exists(p):
            continue
        _, c = read_csv(p)
        e, x = c["E_keV"], c["x_atomic"]
        top = max(x) if x else 0.0
        for i in range(len(e)):
            if x[i] > top * LINE_REL_MIN:
                xmax = max(xmax, e[i])
    tot = imp = 0.0
    edep = 0.0
    for k in keys:
        p = os.path.join(build, "iso_%s_shield.csv" % k)
        if not os.path.exists(p):
            continue
        _, c = read_csv(p)
        if "src_xray" not in c:
            continue
        e, v = c["E_keV"], c["src_xray"]
        tot += sum(v)
        s, m = above(e, v, 2 * xmax)
        imp += s
        edep = max(edep, m)
    return {"nuc": "XRAY", "E_line_max": xmax, "E_cap": 2 * xmax,
            "E_dep_max": edep, "total": tot,
            "impossible": imp, "impossible_pct": 100.0 * imp / max(tot, 1),
            "summing": float("nan"), "summing_pct": float("nan")}


def main():
    build = sys.argv[1] if len(sys.argv) > 1 else str(paths.build("Gamma-1S"))
    keys = ["Th232", "Ac228", "Th228", "Ra224", "Rn220", "Pb212", "Bi212",
            "Tl208"]
    rows = [r for r in (check(build, k) for k in keys) if r]
    if not rows:
        print("не найдено ни одного iso_*.csv в %s" % build)
        return 2
    rows.append(check_xray(build, keys))

    print("каталог: %s" % build)
    print("порог линии эмиссии: %.0e от сильнейшей; порог отсчёта: %d\n"
          % (LINE_REL_MIN, DEP_MIN))
    print("%-7s %9s %9s %9s %12s %12s" %
          ("нуклид", "E линии", "E предел", "E депоз.",
           "невозможно", "суммирование"))
    print("%-7s %9s %9s %9s %12s %12s" %
          ("", "кэВ", "кэВ", "кэВ", "отсч. / %", "отсч. / %"))
    bad = []
    for r in rows:
        smp = ("—" if r["summing"] != r["summing"]
               else "%.0f / %.2f" % (r["summing"], r["summing_pct"]))
        print("%-7s %9.1f %9.1f %9.1f %12s %12s" %
              (r["nuc"], r["E_line_max"], r["E_cap"], r["E_dep_max"],
               "%.0f / %.3f" % (r["impossible"], r["impossible_pct"]), smp))
        if r["impossible"] > 0:
            bad.append(r)

    print("")
    if bad:
        print("ФАНТОМЫ: %d шаблон(ов) дают энерговыделение выше предела "
              "сохранения энергии" % len(bad))
        for r in bad:
            print("   %s: %.0f отсчётов выше %.1f кэВ (предел), максимум "
                  "%.1f кэВ" % (r["nuc"], r["impossible"], r["E_cap"],
                                r["E_dep_max"]))
        return 1
    print("фантомов нет: ни один шаблон не выходит за предел сохранения "
          "энергии")
    return 0


if __name__ == "__main__":
    sys.exit(main())
