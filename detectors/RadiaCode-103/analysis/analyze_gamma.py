# -*- coding: utf-8 -*-
"""Кривые эффективности по сетке прогонов.

Из спектров энерговыделения в кристалле извлекает:
  eps_p(E)  — фотопиковая эффективность на испущенный фотон (абсолютная),
  eps_t(E)  — полная эффективность счёта выше порога,
  f_abs     — множитель самопоглощения относительно предельного случая (воздух).

Активность изотопа в пробе:  A [Бк] = S_net / (eps_p * p_gamma * t),
удельная:                    a [Бк/кг] = A / m,  m = 0.20015 л * rho.
"""
import os
import sys

import numpy as np

# Модули прибора лежат в двух каталогах: разбор в analysis/, запуск прогонов
# в drivers/. Импорт через каталог-сосед иначе не находится: python кладёт
# в sys.path только каталог запускаемого файла.
for _d in ("analysis", "drivers"):
    _p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import rcspec

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = rcspec.RESULTS
PEAK_HALF = 4.0      # кэВ, полуширина окна полного поглощения
THRESH = 20.0        # кэВ, нижний порог счёта прибора
VESSEL_L = 0.20015   # л, промеренный объём сосуда


def analyse_dir(d):
    rows = []
    for fn in sorted(os.listdir(d)):
        if not fn.startswith("E") or not fn.endswith(".csv"):
            continue
        meta, hist = rcspec.read_spec(os.path.join(d, fn))
        n = float(meta["N_primaries"])
        e0 = float(meta["E_prim_keV"])
        lo, hi = int(e0 - PEAK_HALF), int(e0 + PEAK_HALF) + 1
        peak = hist[max(0, lo):hi].sum()
        tot = hist[int(THRESH):].sum()
        rows.append(dict(E=e0, N=n, peak=peak, tot=tot,
                         eps_p=peak / n, d_eps_p=np.sqrt(max(peak, 1)) / n,
                         eps_t=tot / n, d_eps_t=np.sqrt(max(tot, 1)) / n,
                         matrix=meta.get("matrix", "?"),
                         rho=float(meta.get("density_gcm3", 0)),
                         vol=float(meta.get("sample_cm3", 0))))
    rows.sort(key=lambda r: r["E"])
    return rows


def fit_loglog(E, eps, deg=5):
    """Классическая аппроксимация ln(eps) полиномом от ln(E)."""
    ok = eps > 0
    x, y = np.log(E[ok]), np.log(eps[ok])
    p = np.polyfit(x, y, deg)
    resid = np.exp(np.polyval(p, x)) / eps[ok] - 1.0
    return p, float(np.sqrt(np.mean(resid ** 2)) * 100)


def main():
    base = rcspec.rdir("gamma")
    configs = sorted(os.listdir(base)) if os.path.isdir(base) else []
    if not configs:
        raise SystemExit("нет результатов в " + base)

    data = {}
    for cfg in configs:
        d = os.path.join(base, cfg)
        if not os.path.isdir(d):
            continue
        rows = analyse_dir(d)
        if rows:
            data[cfg] = rows
            print("%-22s точек %2d  проба %.2f см³  rho %.3f"
                  % (cfg, len(rows), rows[0]["vol"], rows[0]["rho"]))

    # Проверка ДО записи, и это принципиально. Раньше она стояла ниже, и при
    # пустом каталоге расчётов скрипт успевал записать в results/ файл из
    # одного заголовка — то есть затирал закоммиченную таблицу
    # эффективности. Именно так из репозитория и пропали 150 строк для
    # маринельки 200 мл: результат перезаписан отсутствием данных.
    if not data:
        raise SystemExit(
            "ни в одной конфигурации нет разобранных прогонов: каталог\n"
            "расчётов пуст. Готовая таблица в results/ НЕ тронута.\n"
            "Посчитайте сетку (drivers/run_grid.py) или укажите\n"
            "G4MODELS_BUILD_RADIACODE_103 на готовый каталог.")

    # опорная конфигурация — воздух: предельный случай без самопоглощения
    ref = next((k for k in data if "_air_" in k), None)

    out = rcspec.rdir("efficiency.csv")
    with open(out, "w", encoding="utf-8") as f:
        f.write("config,matrix,rho_gcm3,E_keV,eps_p,d_eps_p,eps_t,d_eps_t,f_abs\n")
        for cfg, rows in data.items():
            rref = {r["E"]: r["eps_p"] for r in data[ref]} if ref else {}
            for r in rows:
                fa = r["eps_p"] / rref[r["E"]] if rref.get(r["E"]) else float("nan")
                f.write("%s,%s,%.4f,%.1f,%.6e,%.2e,%.6e,%.2e,%.4f\n"
                        % (cfg, r["matrix"], r["rho"], r["E"], r["eps_p"],
                           r["d_eps_p"], r["eps_t"], r["d_eps_t"], fa))
    print("\nтаблица:", out)

    # Сводка строится только по энергиям, посчитанным ВО ВСЕХ конфигурациях,
    # иначе при незавершённой сетке сравниваются разные точки.
    if not data:
        raise SystemExit(
            "нет ни одного разобранного прогона: каталог расчётов пуст.\n"
            "Посчитайте сетку (drivers/run_grid.py) или укажите\n"
            "G4MODELS_BUILD_RADIACODE_103 на готовый каталог.")
    common = set.intersection(*(set(r["E"] for r in rows) for rows in data.values()))
    show = [e for e in (661.7, 351.9, 1460.8, 2614.5) if e in common]
    if not show:
        show = sorted(common)[-2:]
    print("\nобщих энергий: %d, сводка при %s кэВ"
          % (len(common), ", ".join("%.1f" % e for e in show)))
    hdr = "%-22s %8s" % ("конфигурация", "rms фита")
    for e in show:
        hdr += " %12s %12s" % ("eps_p(%.0f)" % e, "eps_t(%.0f)" % e)
    print(hdr)
    for cfg, rows in data.items():
        E = np.array([r["E"] for r in rows])
        ep = np.array([r["eps_p"] for r in rows])
        et = np.array([r["eps_t"] for r in rows])
        _, rms = fit_loglog(E, ep)
        line = "%-22s %7.1f%%" % (cfg, rms)
        for e in show:
            i = int(np.argmin(np.abs(E - e)))
            line += " %12.4e %12.4e" % (ep[i], et[i])
        print(line)


if __name__ == "__main__":
    main()
