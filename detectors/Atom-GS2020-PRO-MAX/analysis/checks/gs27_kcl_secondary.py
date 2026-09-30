# -*- coding: utf-8 -*-
# #GS-27 (оператор 28.09: «так же нужно разбирать вторичные пики»): остатки KCl (изм.−фон − модель) у ожидаемых вторичных особенностей 1460,8 кэВ и слепой поиск окон с |pull|>4.
import sys, json, numpy as np
sys.stdout.reconfigure(encoding="utf-8")
E0, M = 1460.822, 510.999
FW = lambda E: 0.7880 * E ** 0.6335                       # закон ПШПВ страницы (fwhm_cal)
FEAT = {"обратное рассеяние": E0 / (1 + 2 * E0 / M), "двойной вылет": E0 - 2 * M, "аннигиляция 511": M, "одиночный вылет": E0 - M,
        "край Комптона": E0 - E0 / (1 + 2 * E0 / M), "вылет I Kα": E0 - 28.6, "1460,8": E0, "наложение 2×1460": 2 * E0}
def win(d, E):
    e = np.asarray(d["e"]); s = (e > E - 1.5 * FW(E)) & (e < E + 1.5 * FW(E))
    n, m, v = (np.asarray(d[k])[s].sum() for k in ("net", "model", "var"))
    return n, m, (n - m) / np.sqrt(v)
for tag in ("fit_kcl_bgw", "fit_kcl_m2_bgw"):
    d = json.load(open((sys.argv[1] if len(sys.argv) > 1 else r"C:\g4work\gs2020\kcl") + "\\%s.json" % tag, encoding="utf-8"))   # argv[1]: папка (v3 — kcl_v3w85)
    print("==", tag, "окно ±1,5 ПШПВ: изм.−фон, модель, отношение, pull")
    for k, E in FEAT.items():
        n, m, p = win(d, E); print("  %-20s %7.1f кэВ  %10.0f %10.0f  %6.3f  %+6.1f" % (k, E, n, m, n / m if m else float("nan"), p))
    hits = [(E, *win(d, E)) for E in np.arange(160, 2990, 5.0)]
    bad = [h for h in hits if abs(h[3]) > 4]
    print("  слепой поиск |pull|>4 (шаг 5 кэВ, окно ±1,5 ПШПВ): %d окон; центры: %s" % (len(bad), ", ".join("%.0f(%+.1f)" % (h[0], h[3]) for h in bad[:40])))
