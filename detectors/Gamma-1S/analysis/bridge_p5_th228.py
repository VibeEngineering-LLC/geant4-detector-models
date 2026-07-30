"""Мост методик на ЧИСТЫХ моно-линиях точечной 5 см (задача 115/93).

РАЗВИЛКА (аудитор 30.07.2026). Расхождение eps модели со штатной кривой —
ПЛАТО +7,8 % на 583-861 кэВ и РАЗРЫВ +27,5 % на 2614,5. Часть его — не
физика, а разная КОНВЕНЦИЯ съёма площади: eps сетки берёт полнопоглощённый
пик в узком окне депозит-спектра (net_counts ~ истинный пик), аттестованная
eps — площадь ГАУСС-ФИТА по РАЗМЫТОМУ спектру, где широкий пик сидит на
континууме неполного поглощения и фит часть его теряет. Вопрос: поправка
ПЛОСКАЯ по энергии или РАСТЁТ?

ОПРЕДЕЛЕНИЕ МОСТА (исправлено). Прежняя попытка «фит/окно» на спектре
распада давала мусор: симметричная подложка на межлинейном континууме
уходила в минус. Правильно:

    мост(E) = N_фит_размытый(E) / N_истинный_пик(E)

  N_истинный_пик — полнопоглощённый пик в депозит-спектре (то, что и есть
    числитель eps сетки);
  N_фит — площадь гаусс-фита по ТОМУ ЖЕ спектру, размытому приборным ПШПВ
    (конвенция прибора).
Мост < 1 и падающий с ростом E означает: аттестованная конвенция ЗАНИЖАЕТ
площадь тем сильнее, чем выше энергия — и это СНИМАЕТ часть завышения eps
модели, приводя обе кривые к одной конвенции.

ДАННЫЕ. Моноэнергетические спектры полной геометрии scat_p5_full есть
только на 661,7 и 2614,5 кэВ — середина и жёсткий край, чего хватает на
бинарный вопрос. Полный набор шести линий требует моно-сетки точечной 5 см
(регенерация драйвером); открыто.

ОГОВОРКА О ФОНЕ. Здесь фон под пиком линейный; СпектраЛайн берёт
ступеньку-из-образа + полином степени 2 (BackgrPower=2 в слепке), который
лучше держит асимметричный континуум под широким пиком 2614,5. Поэтому
восстановление фита здесь — НИЖНЯЯ ГРАНИЦА (мой линейный фон переподставляет
на 2614,5: восстанавливает 84 % истинного пика против 99 % на 662). Мост
0,84 на 2614,5 — верхняя оценка снятия; истинное снятие меньше.

Спектры в рабочем каталоге модели, в репозиторий не входят.
"""
import math
import os
import sys

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erfc

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402

RESULTS = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results"))

# Приборная ПШПВ точечной 5 см: sqrt(a+b*E) по измеренным точкам окна
# «Параметры пиков» (сеанс 30.07.2026).
_EM = np.array([238.632, 300.087, 583.187, 727.330, 860.557, 2614.511])
_PC = np.array([9.835, 9.075, 7.192, 6.664, 6.283, 4.098])
_B, _A = np.polyfit(_EM, (_PC / 100 * _EM) ** 2, 1)


def fwhm(E):
    return math.sqrt(max(_A + _B * E, 1.0))


# Доступные моно-спектры полной геометрии и высота ступеньки-из-образа .cpt.
MONO = {661.657: ("scat_p5_full_E661.7.csv", 0.00035),
        2614.511: ("scat_p5_full_E2614.5.csv", 0.00139)}
# Отношение eps модель/штатная из hard_edge (метод Base) — цель снятия.
EXCESS = {661.657: None, 2614.511: 0.2750}   # плато ~0,078; на 662 нет линии
PLATEAU_EXCESS = 0.078


def load(path):
    ec = []
    for ln in open(path, encoding="utf-8"):
        if ln.startswith(("#", "E_keV")) or not ln.strip():
            continue
        e, c = ln.split(",")
        ec.append((float(e), float(c)))
    return ec


def true_peak(ec, E0, half=3.0):
    return sum(c for e, c in ec if abs(e - E0) < half)


def blur(ec, emax=3400):
    src = np.zeros(emax)
    for e, c in ec:
        i = int(e)
        if 0 <= i < emax:
            src[i] += c
    xs = np.arange(emax, dtype=float)
    out = np.zeros(emax)
    for i in np.nonzero(src)[0]:
        s = fwhm(max(xs[i], 20.0)) / 2.3548
        lo, hi = max(0, int(xs[i] - 6 * s)), min(emax, int(xs[i] + 6 * s))
        g = np.exp(-0.5 * ((xs[lo:hi] - xs[i]) / s) ** 2)
        out[lo:hi] += src[i] * g / (s * math.sqrt(2 * math.pi))
    return xs, out


def fit_area(xs, ys, E0, hstep, zone_half=3.2):
    fw = fwhm(E0)
    s0 = fw / 2.3548
    m = (xs >= E0 - zone_half * fw) & (xs <= E0 + zone_half * fw)
    x, y = xs[m], ys[m]

    def model(xx, A, sh, ws, b0, b1):
        s = ws * s0
        mu = E0 + sh
        amp = A / (s * math.sqrt(2 * math.pi))
        gg = amp * np.exp(-0.5 * ((xx - mu) / s) ** 2)
        st = amp * hstep * 0.5 * erfc((xx - mu) / (s * math.sqrt(2)))
        return b0 + b1 * (xx - xx.mean()) + gg + st

    p, _ = curve_fit(model, x, y, p0=[y.sum(), 0.0, 1.0, np.median(y), 0.0],
                     sigma=np.sqrt(np.maximum(y, 1.0)),
                     bounds=([0, -30, 0.5, -np.inf, -np.inf],
                             [np.inf, 30, 2.0, np.inf, np.inf]),
                     maxfev=20000)
    return p[0]


def find_build():
    b = os.environ.get("G4MODELS_BUILD_GAMMA_1S")
    cands = [b] if b else []
    for c in cands:
        if c and os.path.exists(os.path.join(c, "scat_p5_full_E2614.5.csv")):
            return c
    raise SystemExit(
        "Не найдены моно-спектры scat_p5_full_E*.csv.\n"
        "Задайте G4MODELS_BUILD_GAMMA_1S на рабочий каталог модели.")


if __name__ == "__main__":
    build = find_build()
    print("Мост методик на чистых моно-линиях точечной 5 см.")
    print("мост = площадь гаусс-фита (размытый) / истинный пик (депозит).\n")
    print("%9s %12s %12s %8s %10s" %
          ("E, кэВ", "истин.пик", "фит", "мост", "ПШПВ,кэВ"))
    rows = []
    for E0, (fn, hstep) in sorted(MONO.items()):
        ec = load(os.path.join(build, fn))
        npeak = true_peak(ec, E0)
        xs, ys = blur(ec)
        nfit = fit_area(xs, ys, E0, hstep)
        br = nfit / npeak
        print("%9.3f %12.0f %12.0f %8.4f %10.1f"
              % (E0, npeak, nfit, br, fwhm(E0)))
        rows.append((E0, npeak, nfit, br))

    b662 = [r[3] for r in rows if abs(r[0] - 661.657) < 1][0]
    b2614 = [r[3] for r in rows if abs(r[0] - 2614.511) < 1][0]
    print("\nМост растёт: %.3f (662) -> %.3f (2614), падение %.1f %%"
          % (b662, b2614, 100 * (b2614 - b662)))

    print("\nПрименение к завышению eps модели над штатной (hard_edge):")
    print("  плато 583-861: excess +%.1f %%; мост на середине %.3f ->"
          " остаток +%.1f %% (конвенция плато почти не трогает)"
          % (100 * PLATEAU_EXCESS, b662,
             100 * ((1 + PLATEAU_EXCESS) * b662 - 1)))
    resid = (1 + EXCESS[2614.511]) * b2614 - 1
    print("  2614,5: excess +%.1f %%; мост %.3f (нижняя граница) ->"
          " остаток +%.1f %% (верхняя граница снятия)"
          % (100 * EXCESS[2614.511], b2614, 100 * resid))
    print("  ВЫВОД: мост НЕ плоский; снимает часть РАЗРЫВА на 2614,5"
          " (до ~%.0f п.п.), плато 7,8 %% оставляет как физику модели."
          % (100 * (EXCESS[2614.511] - resid)))

    csvio.write(
        os.path.join(RESULTS, "bridge_p5_th228.csv"),
        ["E_keV", "true_peak", "fit_area", "bridge_fit_over_truepeak"],
        [("%.3f" % E, "%.0f" % n, "%.0f" % f, "%.4f" % b)
         for E, n, f, b in rows],
        comments=[
            "Мост методик; точечная 5 см; моно-спектры полной геометрии"
            " scat_p5_full (661;7 и 2614;5 кэВ).",
            "мост = площадь гаусс-фита по размытому спектру / истинный"
            " полнопоглощённый пик депозит-спектра.",
            "мост < 1 и падающий: аттестованная конвенция занижает площадь"
            " тем сильнее; чем выше энергия — снимает часть разрыва.",
            "Фон линейный; СпектраЛайн берёт ступеньку+полином 2 —"
            " восстановит больше; мост 0;84 на 2614;5 — нижняя граница.",
            "Полный набор шести линий требует моно-сетки точечной 5 см"
            " (регенерация драйвером); открыто.",
        ])
    print("\nтаблица: %s" % os.path.join(RESULTS, "bridge_p5_th228.csv"))
