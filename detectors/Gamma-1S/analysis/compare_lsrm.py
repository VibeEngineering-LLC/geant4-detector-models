"""Сверка расчётной кривой ППП Гамма-1С с измеренной кривой ЛСРМ (.efr).

Расчёт: сетка grid/rho1.60_E*.csv (Маринелли 1 л, ОИСН-16 1,6 г/см³).
Эксперимент: 15 точек из УДС-ГЦ-63х63-USB__SN-01_-_Маринелли.efa
(тот же файл в ref/, разбор — fetch_efr.parse_efr).

Выход: таблица точка-к-точке, средневзвешенное отношение МК/эксп (аналог
K_NORM), его разброс, и сверка формы кривой после нормировки.
"""
import glob
import math
import os
import re
import sys

# Корни путей — из переменных окружения (common/py/paths.py), чтобы в коде не
# было ни одного пути, привязанного к конкретной машине.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
import paths  # noqa: E402
import stamp  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import peakwin  # noqa: E402

# Объявление наблюдаемой — что именно за число лежит в сводке. Строки окна
# и полки собираются из констант peakwin: прежде здесь стояло рукописное
# «полка E-30…E-10», а правило мигрировало на [E-25; E-10] — объявление лгало
# (внутренний аудит 31.07.2026).
OBS_SUMMARY = dict(
    {
        "quantity": "сводка сверки расчётной кривой Маринелли с аттестованной"
                    " .efr: нормировка; разброс формы; хи2 на степень свободы",
        "area": "чистая площадь пика за вычетом полки континуума"
                " (правило common/py/peakwin)",
    },
    **peakwin.declare())

# parse_efr живёт в инструментах репозитория, а не среди данных
sys.path.insert(0, str(paths.tools()))
from fetch_efr import parse_efr  # noqa: E402

BUILD = str(paths.build("Gamma-1S"))

if not os.path.isdir(BUILD):
    raise SystemExit(
        "Нет каталога расчётных спектров %s.\n"
        "Они не коммитятся (сотни файлов), а воспроизводятся драйверами:\n"
        "    python detectors/Gamma-1S/drivers/run_grid.py\n"
        "    python detectors/Gamma-1S/drivers/run_all_grids.py\n"
        "Либо укажите G4MODELS_BUILD_GAMMA_1S на готовый каталог."
        % BUILD)
REF = str(paths.ref("Gamma-1S"))

SUMMARY = os.path.join(str(paths.results("Gamma-1S")), "compare_lsrm_summary.csv")


def marinelli_k():
    """Отношение МК/эксперимент по маринелли — ИЗ ФАЙЛА, не из литерала.

    Это число служит масштабом сразу нескольким разборам. Пока оно было
    переписано от руки в трёх скриптах и в отчёте, после каждого пересчёта
    сетки часть копий отставала, а на глаз такое не ловится: и 1,171, и 1,165
    выглядят одинаково правдоподобно.
    """
    try:
        with open(SUMMARY, encoding="utf-8") as fh:
            fh.readline()
            return "%.3f" % float(fh.readline().split(",")[0])
    except (OSError, ValueError, IndexError):
        return "неизвестно (запустите compare_lsrm.py)"


# Окно ППП и полка — единственная реализация common/py/peakwin (полка
# [E−25; E−10], счёт в каналах). Прежняя собственная копия правила держала
# полку E−30, захватывающую пик вылета иода E−28,6 (аудит 31.07.2026).
# Континуум под пиком образуют события с почти полным энерговыделением
# (многократное рассеяние с малой потерей). Справа от E0 при моноэнергетическом
# источнике без наложений отсчётов нет, поэтому полка только левая.


def read_run(path):
    N, E0 = None, None
    hist = {}
    for line in open(path, encoding="utf-8"):
        if line.startswith("#"):
            if "N_primaries" in line:
                N = int(line.split("=")[1])
            elif "E_prim_keV" in line:
                E0 = float(line.split("=")[1])
            continue
        if line and line[0].isdigit():
            e, c = line.split(",")
            hist[float(e)] = int(c)
    det = {}
    net = peakwin.area(hist, E0, detail=det)
    bg = (det["side"] / det["n_side"] * det["n_peak"]
          if det["n_side"] else 0.0)
    # D(bg) = (n/nside)^2 * side = (n/nside)*bg
    var = (det["gross"] + (det["n_peak"] / det["n_side"]) * bg
           if det["n_side"] else det["gross"])
    return E0, net, math.sqrt(max(var, 1.0)), N


def mc_curve(rho_tag="rho1.60"):
    out = {}
    for f in sorted(glob.glob(os.path.join(BUILD, "grid", rho_tag + "_E*.csv"))):
        E0, net, dnet, N = read_run(f)
        if net > 0:
            out[round(E0, 3)] = (net / N, dnet / N)
    return out


def lsrm_points():
    """Точки .efa Маринелли: E -> (eff, dpct, nuclide)."""
    p = paths.efficiency_curve("Маринелли", "efa")
    if p is None:
        raise SystemExit("не найдена измеренная кривая .efa для Маринелли "
                         "в %s" % paths.ref("Gamma-1S"))
    txt = paths.read_text(p)
    pts = []
    for line in txt.splitlines():
        m = re.match(r"^(\d+\.?\d*)=([0-9.E+-]+),([0-9.]+),([\w-]+),", line)
        if m:
            pts.append((float(m.group(1)), float(m.group(2)),
                        float(m.group(3)), m.group(4)))
    return pts


if __name__ == "__main__":
    mc = mc_curve()
    pts = lsrm_points()
    print("%-10s %-8s %12s %12s %8s %8s" %
          ("E, кэВ", "нуклид", "эксп", "МК", "МК/эксп", "+-"))
    logs, ws = [], []
    rows = []
    for E, eff, dpct, nuc in pts:
        key = min(mc, key=lambda k: abs(k - E))
        if abs(key - E) > 1.0:
            print("%-10.1f %-8s %12.3e   нет расчётной точки" % (E, nuc, eff))
            continue
        m, dm = mc[key]
        r = m / eff
        dr = r * math.sqrt((dm / m) ** 2 + (dpct / 100) ** 2)
        rows.append((E, nuc, eff, m, r, dr))
        logs.append(math.log(r))
        ws.append(1.0 / (dr / r) ** 2)
        print("%-10.1f %-8s %12.3e %12.3e %8.3f %8.3f" % (E, nuc, eff, m, r, dr))

    lw = sum(l * w for l, w in zip(logs, ws)) / sum(ws)
    k = math.exp(lw)
    dk = 1.0 / math.sqrt(sum(ws))
    # разброс формы после нормировки
    dev = [math.exp(l - lw) - 1 for l in logs]
    rms = math.sqrt(sum(d * d for d in dev) / len(dev))
    chi2 = sum(w * (l - lw) ** 2 for l, w in zip(logs, ws)) / (len(logs) - 1)
    print("\nсредневзвешенное МК/эксп  = %.3f +- %.3f" % (k, k * dk))
    print("эквивалент K_NORM         = %.3f   (RadiaCode: 0.833)" % (1 / k))
    print("разброс формы (RMS)       = %.1f %%" % (100 * rms))
    print("chi2/dof формы            = %.2f" % chi2)
    print("макс. отклонение формы    = %+.1f %%  и  %+.1f %%"
          % (100 * max(dev), 100 * min(dev)))

    # Сводка кладётся файлом, чтобы это число ЖИЛО В ОДНОМ МЕСТЕ. Раньше оно
    # было переписано от руки в compare_point.py и в отчёте и отставало после
    # каждого пересчёта сетки.
    # Запись — общей реализацией csvio, а не ручным fh.write: ручная запись
    # обходит и сторожа формата, и объявление наблюдаемой, из-за чего эта
    # таблица оставалась несравнимой ни с какой другой.
    csvio.write(
        SUMMARY,
        ["k_mc_over_exp", "d_k", "n_points", "rms_shape", "chi2_dof"],
        [("%.4f" % k, "%.4f" % (k * dk), "%d" % len(logs), "%.4f" % rms,
          "%.3f" % chi2)],
        stamp=stamp.lines(
            "detectors/Gamma-1S/analysis/compare_lsrm.py", OBS_SUMMARY,
            geometry_dir=str(paths.geometry("Gamma-1S")),
            names=stamp.SRC_LISTS["Gamma-1S"], repo_dir=str(paths.REPO)))
