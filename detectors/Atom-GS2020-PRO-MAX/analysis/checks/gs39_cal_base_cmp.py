# -*- coding: utf-8 -*-
# #GS-39: сравнение основ шкалы по форме (cal_shape_base*.json): наклон кэВ/кан у 911 и 2614 (у СпектраЛайн 0,478/0,460–0,479),
# χ²/ν групп X-K, 238, 2614 и общий χ²/ν из лога. Запуск без аргументов; файлы из out_v5 и логи из out_v5_oisn10.
import sys, json, re, os, numpy as np
sys.stdout.reconfigure(encoding="utf-8")
O, L = r"C:\g4work\gs2020\run_marinelli\out_v5", r"C:\g4work\gs2020\run_marinelli\out_v5_oisn10"
V = [("квартика (опубл.)", O + r"\base_shape_2026-09-29\cal_shape.json", L + r"\base_shape_2026-09-29\cal_shape.log")] + \
    [("основа %d" % d, O + r"\cal_shape_base%d.json" % d, L + r"\cal_shape_base%d.log" % d) for d in (1, 2, 3)]
def axis(r, n=8192):
    ch = np.arange(n, dtype=float); e = sum(c * ch ** i for i, c in enumerate(r["file_coeffs"]))
    return e + np.interp(e, r["knots"], r["theta"])
for name, fj, fl in V:
    if not os.path.exists(fj): print(name, "— нет файла"); continue
    D = json.load(open(fj, encoding="utf-8")); log = open(fl, encoding="utf-8").read() if os.path.exists(fl) else ""
    tot = dict(re.findall(r"== (\w+) ==.*?χ²/ν = ([\d.]+)", log, re.S))
    for t in ("sample", "bg", "bgw", "kcl"):
        r = D[t]; e = axis(r); sl = lambda E: np.gradient(e)[int(np.interp(E, e, np.arange(len(e))))]
        g = {k: v.get("chi2_nu", float("nan")) for k, v in r["groups"].items()}
        print("%-18s %-6s наклон 911 %.4f 2614 %.4f | χ²/ν общий %s X-K %.2f 238 %.2f 2614 %.2f | узлы %s" % (name, t, sl(911.2), sl(2614.5),
              tot.get(t, "?"), g.get("X-K", np.nan), g.get("238", np.nan), g.get("2614", np.nan), ",".join("%g" % k for k in r["knots"])))
