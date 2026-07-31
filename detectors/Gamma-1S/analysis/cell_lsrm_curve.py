"""Клетка 2x2 «кривая ЛСРМ x наш алгоритм»: маринелли Th-232 420-7-17.

Дизайн 2x2 (задача 116): расхождение активностей раскладывается на вклад
кривой и вклад алгоритма. Две клетки уже есть: (наша кривая x наш алгоритм)
= 0,796 (kit_recalc) и (ЛСРМ x СпектраЛайн) = 0,59 Base / 0,84 ZBZ+ (сеансы
оператора). Здесь третья: площади и скорости берутся ИЗ НАШЕГО конвейера
(те же, что в kit_recalc_volume.csv), а эффективность — ИЗ ШТАТНОЙ КРИВОЙ
.efa, как считал бы СпектраЛайн без поправки на цепочку (она у оператора
выключена): A_линии = R / (eps_efa(E) * I).

Плотность источника равна плотности аттестации (ОИСН-16, 1,6) — поправка
самопоглощения не нужна, клетка чистая от f-пересчёта.

Побочный контроль — площадь 2614,5 в трёх обработках ОДНОГО спектра:
наша (окно +-1 ПШПВ), сеанс СпектраЛайн 29.07.2026 (54000, оба метода) и
аттестационный .efr (59668). Наша и сеансовая согласны (~54 тыс.),
аттестационная выше обеих на ~10% — это прямо входит в eps_efa узла 2614,5.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
import paths  # noqa: E402
import stamp  # noqa: E402

RESULTS = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results"))
A0 = 3104.0          # Бк: паспорт 1940 Бк/кг x 1,6 кг
TLIVE = 11359.164    # с, живое время записи (свойства спектра)


# Объявление наблюдаемой — что именно за число лежит в таблице. Без него
# таблицу нельзя сравнивать ни с какой другой: за один вечер 30.07.2026
# подмена определения стоила вывода четыре раза (method-rules §5).
OBS = {
    "quantity":
        "скорости счёта нашего конвейера на кривой ЛСРМ — клетка эксперимента 2x2",
    "area":
        "чистая площадь пика за вычетом подложки; съём нашим конвейером",
    "window":
        "окно нашего конвейера в долях ПШПВ",
    "shelf":
        "симметричные полки той же доли ПШПВ",
    "blurred":
        "не применимо — измеренный спектр",
}


def _stamp(inputs=None):
    return stamp.lines("detectors/Gamma-1S/analysis/cell_lsrm_curve.py", OBS,
                       inputs=inputs,
                       geometry_dir=str(paths.geometry("Gamma-1S")),
                       names=stamp.SRC_LISTS["Gamma-1S"],
                       repo_dir=str(paths.REPO))


def efa_lines():
    """{E: (eps, I_доля)} из блока Th-232 маринелльного .efr."""
    p = paths.efficiency_curve("Маринелли", ext="efr")
    out = {}
    inblock = False
    for ln in open(str(p), "rb").read().decode("cp1251").splitlines():
        if ln.startswith("["):
            inblock = "Th232" in ln.replace("-", "")
            continue
        if not inblock or "=" not in ln:
            continue
        key, val = ln.split("=", 1)
        try:
            E = float(key)
        except ValueError:
            continue
        f = val.split(",")
        out[E] = (float(f[0]), float(f[5]) / 100.0)
    return out


def our_rates():
    """{E: (rate_cps, usable)} маринелли Th-232 из kit_recalc_volume.csv."""
    out = {}
    for ln in open(os.path.join(RESULTS, "kit_recalc_volume.csv"),
                   encoding="utf-8"):
        if ln.startswith("#") or ln.startswith("geometry"):
            continue
        f = ln.strip().split(",")
        if f[0] == "Marinelli_1L" and f[1] == "Th-232":
            out[float(f[2])] = (float(f[3]), f[11] == "1")
    return out


if __name__ == "__main__":
    eff = efa_lines()
    rates = our_rates()
    print("Клетка (кривая ЛСРМ x наш алгоритм), маринелли 420-7-17,"
          " паспорт %.0f Бк:\n" % A0)
    print("%9s %8s %11s %7s %9s %8s %s" %
          ("E, кэВ", "имп/с", "eps ЛСРМ", "I", "A, Бк", "A/пасп", ""))
    rows = []
    for E in sorted(rates):
        key = min(eff, key=lambda k: abs(k - E))
        if abs(key - E) > 1.0:
            continue
        eps, I = eff[key]
        r, usable = rates[E]
        A = r / (eps * I)
        note = "" if usable else "  бленд — справочно"
        print("%9.1f %8.3f %11.4e %7.3f %9.0f %8.3f%s"
              % (E, r, eps, I, A, A / A0, note))
        rows.append((E, r, eps, I, A, A / A0, usable))

    out = os.path.join(RESULTS, "cell_lsrm_curve.csv")
    csvio.write(
        out,
        ["E_keV", "rate_cps", "eps_lsrm", "I_frac", "A_Bq", "ratio",
         "usable"],
        [("%.3f" % e, "%.4f" % r, "%.6e" % ep, "%.4f" % i, "%.1f" % a,
          "%.4f" % rt, "%d" % (1 if u else 0))
         for e, r, ep, i, a, rt, u in rows],
        comments=[
            "Клетка 2x2 (кривая ЛСРМ x наш алгоритм): скорости нашего"
            " конвейера, eps из штатного .efa, A=R/(eps*I), без поправки"
            " на цепочку (у оператора выключена).",
            "Паспорт 3104 Бк; плотность source = плотности аттестации,"
            " поправки самопоглощения нет.",
            "Контроль площади 2614,5 (один спектр, три обработки): наша"
            " %.0f, сеанс СпектраЛайн 54000, аттестация .efr 59668."
            % (rates.get(2614.511, (0, 0))[0] * TLIVE),
        ],
        stamp=_stamp())
    print("\nтаблица: %s" % out)
    print("контроль площади 2614,5: наша %.0f, СпектраЛайн 54000,"
          " аттестация .efr 59668 (+%.1f%% к нашей)"
          % (rates.get(2614.511, (0, 0))[0] * TLIVE,
             100 * (59668 / (rates[2614.511][0] * TLIVE) - 1)))
