"""Предсказание просадки узлов Th-228 моделью (путь 2 задачи 140).

ЗАЧЕМ. Узлы аттестованной кривой, снятые с Th-228, лежат на 4,5 % ниже линии,
построенной по одиночным излучателям (`attested_cascade_check.py`). Знак
совпадает с ожидаемым для неучтённого каскадного суммирования, но 2,1 сигмы —
это направление, а не доказательство. Совпадение ВЕЛИЧИНЫ с независимым
предсказанием было бы подтверждением куда более сильным, чем совпадение знака.

ЧТО СЧИТАЕТСЯ. Поправка на каскадное суммирование `C = eps_моно / eps_распад`
для точечной геометрии 5 см:

  eps_моно   — эффективность на чистой моноэнергии (складывать не с чем);
  eps_распад — та же линия в прогоне полного распада Th-228, где часть
               отсчётов уходит из пика в суммарный.

`C > 1` означает, что в распаде пик недобирает. Если аттестация НЕ вводила
поправку, её узел занижен ровно на `1 − 1/C`, и это число сравнивается с
наблюдённой просадкой.

ОБЕ СТОРОНЫ СНИМАЮТСЯ ОДНИМ ПРАВИЛОМ — `peakwin.area`, окно и полка в каналах.
Разные конвенции на числителе и знаменателе дали бы `C`, в которую войдёт
разница правил, а не физика; этот класс дефекта в линии Гамма-1С ловился пять
раз за один вечер.

ЧЕГО СРАВНЕНИЕ НЕ ДОКАЗЫВАЕТ. `C` модели относится к цепочке Th-228 в том
изотопном равновесии, в каком её разыгрывает Geant4. Источник комплекта мог
быть в ином состоянии равновесия, и тогда доли каскадных переходов отличаются.
Совпадение величин — сильный довод, расхождение — не опровержение.
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
import paths  # noqa: E402
import peakwin  # noqa: E402
import stamp  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))
OUT = str(paths.results("Gamma-1S"))

DECAY = "p5_Th228"
# Моно-прогоны того же положения источника и того же 4pi — из набора моста
# (задача 115). Ключ — энергия линии, значение — имя файла.
MONO = {583.187: "bridge4pi_E0583.2.csv",
        727.330: "bridge4pi_E0727.3.csv",
        860.557: "bridge4pi_E0860.6.csv",
        2614.511: "bridge4pi_E2614.5.csv"}
# Наблюдённая просадка узла относительно опоры одиночных излучателей,
# results/attested_cascade_check.csv. Вне 662…1116 кэВ опоры нет — там None.
OBSERVED = {583.187: None, 727.330: -4.31, 860.557: -4.67, 2614.511: None}

OBS = {
    "quantity": "C — поправка на каскадное суммирование для точечной 5 см;"
                " отношение эффективности на моноэнергии к эффективности"
                " той же линии в прогоне полного распада",
    "area": "чистая площадь пика за вычетом полки; ОДНО правило на обе"
            " стороны (common/py/peakwin.area)",
    "window": "+-6 кэВ в каналах; полка [E-25; E-10] — конвенция peakwin",
    "shelf": "односторонняя слева; вычитается на обеих сторонах одинаково",
    "blurred": "нет — депозит-спектры как есть",
    "geometry": "точечный источник z = 91 мм; полный 4pi; моно и распад в"
                " одинаковых условиях",
}


def read_hist(path):
    hist, N = {}, None
    for ln in open(path, encoding="utf-8"):
        if ln.startswith("#"):
            if "N_primaries" in ln:
                N = int(ln.split("=")[1])
            continue
        if ln and ln[0].isdigit():
            e, c = ln.split(",")
            hist[float(e)] = float(c)
    return hist, N


def emitted(E, tol=2.0):
    """Сколько квантов этой энергии испущено в прогоне распада.

    Из `*_emit.csv` того же расчёта, а не из справочника выходов: справочник —
    другой источник, и разница выходов вошла бы в C как физика.
    """
    p = os.path.join(BUILD, "%s_emit.csv" % DECAY)
    if not os.path.exists(p):
        return None
    tot = 0.0
    for ln in open(p, encoding="utf-8"):
        if ln.startswith("#") or not ln.strip():
            continue
        if ln[0].isdigit():
            e, c = ln.split(",")
            if abs(float(e) - E) <= tol:
                tot += float(c)
    return tot


def main():
    dpath = os.path.join(BUILD, "%s.csv" % DECAY)
    if not os.path.exists(dpath):
        raise SystemExit("Нет прогона распада %s в %s." % (DECAY, BUILD))
    dhist, _ = read_hist(dpath)

    print("Предсказание просадки: поправка на суммирование C для точечной"
          " 5 см.\n")
    print("%10s %12s %12s %8s %12s %12s %8s"
          % ("E, кэВ", "площ.распад", "испущено", "C", "занижение",
             "наблюдено", "разница"))
    rows = []
    for E in sorted(MONO):
        mp = os.path.join(BUILD, MONO[E])
        if not os.path.exists(mp):
            print("%10.3f  нет моно-прогона %s — пропущена" % (E, MONO[E]))
            continue
        mhist, mN = read_hist(mp)
        nem = emitted(E)
        if not nem:
            print("%10.3f  в *_emit.csv нет этой линии — пропущена" % E)
            continue
        a_mono = peakwin.area(mhist, E)
        a_dec = peakwin.area(dhist, E)
        eps_mono = a_mono / mN
        eps_dec = a_dec / nem
        C = eps_mono / eps_dec
        under = 100 * (1.0 / C - 1.0)          # занижение узла; %
        # Погрешность C определяется РАСПАДНОЙ стороной: там площадь пика на
        # два порядка меньше, чем в моно. Без этой строки предсказание из
        # сотни отсчётов выглядело бы наравне с наблюдением из миллиона.
        rel = math.sqrt(1.0 / max(a_dec, 1.0) + 1.0 / max(a_mono, 1.0))
        d_under = 100 * rel / C
        obs = OBSERVED.get(E)
        if obs is None:
            diff = ""
        else:
            diff = "%+.2f" % (under - obs)
        print("%10.3f %12.0f %12.0f %8.4f %8.2f+-%.2f %11s %8s"
              % (E, a_dec, nem, C, under, d_under,
                 "—" if obs is None else "%+.2f %%" % obs, diff))
        rows.append((E, a_dec, nem, C, under, d_under, obs))

    print("\nВЫВОД.")
    pair = [(E, u, du, o) for E, _, _, _, u, du, o in rows if o is not None]
    if not pair:
        print("  Нет линий; где есть и предсказание; и наблюдение.")
    else:
        for E, u, du, o in pair:
            n = abs(u - o) / du if du > 0 else float("inf")
            print("  %.1f кэВ: модель %+.2f +- %.2f %%; наблюдается %+.2f %%;"
                  " расхождение %+.2f п.п. = %.1f сигмы"
                  % (E, u, du, o, u - o, n))
        worst = max(abs(u - o) / du for _, u, du, o in pair)
        prec = max(du for _, _, du, _ in pair)
        need = max(abs(o) for _, _, _, o in pair)
        if prec > need:
            # Решает не расхождение, а ТОЧНОСТЬ предсказания: пока она хуже
            # самой проверяемой величины, ни совпадение, ни расхождение не
            # значат ничего. Объявлять «величины расходятся» в таком
            # положении — выдавать шум за результат.
            print("  ТЕСТ НЕ РАЗРЕШЁН. Точность предсказания (%.1f %%) хуже"
                  " самой проверяемой просадки (%.1f %%):" % (prec, need))
            print("  распадный прогон даёт в пике 144…260 отсчётов. Ни"
                  " совпадение; ни расхождение здесь не довод.")
            print("  Нужен прогон распада примерно на порядок длиннее — тот"
                  " же; что требуется задаче 122.")
        elif worst < 2.0:
            print("  Величины СОГЛАСУЮТСЯ (%.1f сигмы). Это уже не совпадение"
                  " знака: независимый расчёт" % worst)
            print("  воспроизводит наблюдённую просадку по величине —"
                  " механизм каскадного суммирования в аттестации"
                  " подтверждён.")
        else:
            print("  Величины РАСХОДЯТСЯ значимо (%.1f сигмы). Знак совпадает;"
                  " но механизм объясняет не всё:" % worst)
            print("  либо суммирование в аттестации учитывалось частично;"
                  " либо изотопное равновесие источника")
            print("  отличается от разыгрываемого моделью; либо просадка"
                  " имеет и другую причину.")
    print("  ОГОВОРКА: C модели относится к цепочке Th-228 в равновесии"
          " Geant4; состояние источника комплекта")
    print("  не установлено. Совпадение — сильный довод; расхождение — не"
          " опровержение.")

    csvio.write(
        os.path.join(OUT, "cascade_predict_p5.csv"),
        ["E_keV", "peak_decay", "emitted", "C_summing", "underest_pct",
         "d_underest_pct", "observed_pct"],
        [("%.3f" % E, "%.0f" % a, "%.0f" % n, "%.4f" % C, "%+.2f" % u,
          "%.2f" % du, "" if o is None else "%+.2f" % o)
         for E, a, n, C, u, du, o in rows],
        comments=[
            "C = eps(моно) / eps(распад) для точечной 5 см: во сколько раз"
            " пик недобирает в полном распаде из-за суммирования.",
            "underest_pct = 1/C - 1 — на столько занижен узел аттестации;"
            " ЕСЛИ поправка при её построении не вводилась.",
            "observed_pct — просадка узла относительно опоры одиночных"
            " излучателей (attested_cascade_check.csv); вне 662…1116 опоры нет.",
            "Обе стороны C сняты ОДНИМ правилом peakwin.area — иначе в C"
            " вошла бы разница конвенций; а не физика.",
            "C модели относится к цепочке в равновесии Geant4; состояние"
            " источника комплекта не установлено.",
        ],
        stamp=stamp.lines(
            "detectors/Gamma-1S/analysis/cascade_predict_p5.py", OBS,
            inputs=[dpath] + [os.path.join(BUILD, MONO[E]) for E in sorted(MONO)
                              if os.path.exists(os.path.join(BUILD, MONO[E]))],
            geometry_dir=str(paths.geometry("Gamma-1S")),
            names=stamp.SRC_LISTS["Gamma-1S"], repo_dir=str(paths.REPO)))
    print("\nтаблица: %s" % os.path.join(OUT, "cascade_predict_p5.csv"))


if __name__ == "__main__":
    main()
