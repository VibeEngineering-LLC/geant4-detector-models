# -*- coding: utf-8 -*-
"""Контроль жёсткого края: сколько ППП съедает вылет тормозного излучения.

Два плеча одного прогона (`macros/curve_hard_brems.mac`), одна сборка, один
штамп; различие единственное — `/process/inactivate eBrem`. Отношение плеч
есть мера самого механизма, а не суммы расхождений между постановками.

Зачем. Расхождение с библиотечной кривой растёт с энергией до −15,0 % на
3000 кэВ. Внешний довод (Am6er/BecqMoni, ветка `pie`,
`tools/tccfcalc/README.md`, §5.3): эталонная программа ЛСРМ TCCFCALC не
переносит электрон вовсе и потому не теряет энергию на вылет тормозного;
отключение тормозного в ЕГО расчёте сажает кривую на ЛСРМ. Здесь то же
предположение проверяется на нашем коде и на НАШЕЙ библиотечной кривой,
происхождение которой иное (поле «Calculated from geometry» самого BecqMoni,
а не экспорт ЛСРМ).

Что показал прогон 06.08.2026. Библиотечная кривая легла МЕЖДУ плечами на всех
жёстких узлах: цена тормозного вдвое больше наблюдаемого дефицита (+29,9 %
против 14,0 % на 3000 кэВ). Значит механизм у автора есть, но приближённый, а
версия «тормозного у них нет вовсе» опровергнута: плечо без тормозного
перелетает его кривую.

Узел 661,7 кэВ брался контролем в расчёте на ноль, и ноля там НЕ оказалось:
+1,80 ± 0,63 %. Это не сбой постановки, а свойство строгого съёма — окно
±1,5 кэВ по нерасплывшемуся спектру, поэтому вылет кванта жёстче 1,5 кэВ
выбивает событие из пика, тогда как при приборном разрешении он бы в пике
остался.

    python analysis/brems_control.py <каталог сборки>
"""
import io
import math
import os
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(_HERE, "..", "results"))
REF = os.path.normpath(os.path.join(_HERE, "..", "reference",
                                    "becqmoni-library-curve.csv"))
# Каталог сборки — аргументом или переменной окружения. Абсолютных путей,
# привязанных к машине, в дереве нет (то же правило, что в CMakeLists).
BUILD = None
ARMS = (("on", "spectra_hard_on", "eff_hard_brems_on.csv"),
        ("off", "spectra_hard_off", "eff_hard_brems_off.csv"))


def read_curve(path, cols):
    rows, head = [], None
    for ln in io.open(path, encoding="utf-8"):
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        p = ln.split(",")
        if head is None:
            head = p
            continue
        r = dict(zip(head, p))
        rows.append([float(r[c]) for c in cols])
    return rows


def loglog_interp(xs, ys, x):
    if x <= xs[0] or x >= xs[-1]:
        return None
    for i in range(1, len(xs)):
        if xs[i] >= x:
            t = ((math.log(x) - math.log(xs[i - 1]))
                 / (math.log(xs[i]) - math.log(xs[i - 1])))
            return math.exp(math.log(ys[i - 1])
                            + t * (math.log(ys[i]) - math.log(ys[i - 1])))
    return None


def export(src, out):
    """Площади снимаются ТЕМ ЖЕ export_curve.py, что и опубликованная кривая.

    Своя реализация окна съёма здесь была бы вторым определением наблюдаемой
    и сравнивала бы конвенции, а не физику.
    """
    rc = subprocess.call([sys.executable,
                          os.path.join(_HERE, "export_curve.py"), src, out])
    if rc != 0:
        raise SystemExit("export_curve.py вернул %d на %s" % (rc, src))


def main():
    build = (sys.argv[1] if len(sys.argv) > 1
             else os.environ.get("ASN16_BUILD"))
    if not build:
        raise SystemExit(
            "Не задан каталог сборки, где лежат spectra_hard_on/off.\n"
            "  python analysis/brems_control.py <каталог сборки>\n"
            "либо переменная окружения ASN16_BUILD.")

    curves = {}
    for tag, sub, name in ARMS:
        src = os.path.join(build, sub)
        if not os.path.isdir(src):
            raise SystemExit("нет каталога спектров %s — прогон не сделан" % src)
        out = os.path.join(RES, name)
        export(src, out)
        curves[tag] = {r[0]: (r[1], r[2])
                       for r in read_curve(out, ("E_keV", "eps_peak",
                                                 "d_eps_peak"))}

    lib = read_curve(REF, ("E_keV", "eff_peak", "d_eff_pct"))
    xs = [r[0] for r in lib]
    ys = [r[1] for r in lib]

    es = sorted(curves["on"])
    print("  E, кэВ    с тормозным      без       без/с      библиотека"
          "   с торм./библ.  без торм./библ.")
    for e in es:
        on, d_on = curves["on"][e]
        off, d_off = curves["off"][e]
        gain = 100.0 * (off / on - 1.0)
        # Погрешность отношения: плечи независимы, статистика в каждом своя.
        d_gain = 100.0 * (off / on) * math.hypot(d_on / on, d_off / off)
        lv = loglog_interp(xs, ys, e)
        if lv is None:
            lv = ys[-1] if e >= xs[-1] else ys[0]
        print("%9.1f   %.4e   %.4e   %+6.2f±%.2f %%   %.4e   %+7.2f %%   "
              "%+7.2f %%"
              % (e, on, off, gain, d_gain, lv,
                 100.0 * (on / lv - 1.0), 100.0 * (off / lv - 1.0)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
