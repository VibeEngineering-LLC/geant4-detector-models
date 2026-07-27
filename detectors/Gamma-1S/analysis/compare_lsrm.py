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
import paths  # noqa: E402

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

# Окно ППП: расчёт без уширения, пик острый; края учитывают утечку в
# соседний канал.
WIN = 6.0        # +- кэВ вокруг E0
BG0, BG1 = 30.0, 10.0   # левая полка континуума: E0-30 .. E0-10

# Континуум под пиком образуют события с почти полным энерговыделением
# (многократное рассеяние с малой потерей). Справа от E0 при моноэнергетическом
# источнике без наложений отсчётов нет, поэтому полка только левая.
SUBTRACT_BG = True


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
    gross = sum(c for e, c in hist.items() if abs(e - E0) <= WIN)
    side = sum(c for e, c in hist.items() if E0 - BG0 <= e <= E0 - BG1)
    nside = BG0 - BG1                      # ширина полки в каналах по 1 кэВ
    n = 2 * WIN + 1
    bg = side / nside * n if SUBTRACT_BG else 0.0
    # D(bg) = (n/nside)^2 * side = (n/nside)*bg; вывод — в export_curves.py
    var = gross + (n / nside) * bg
    return E0, gross - bg, math.sqrt(max(var, 1.0)), N


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
