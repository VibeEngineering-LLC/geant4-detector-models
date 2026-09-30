# -*- coding: utf-8 -*-
"""#GS-28: сдвиг δ по форме групп на ШКАЛЕ ФАЙЛА (без поправки по реперам) против центроида одиночного гаусса репера."""
import sys, os, json, numpy as np
sys.stdout.reconfigure(encoding="utf-8"); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "cal"))
import gs2020_calib_mplet as M
EXT = float(os.environ.get("GS28_EXT", "15")); own = json.load(open(M.C.OUT_JSON, encoding="utf-8")); ds = np.arange(-15, 15.05, 0.05)
for tag in sys.argv[1:] or ["bg", "bgw", "kcl", "sample"]:
    d = M.C.read(tag); y = d["counts"]; x = M.C.s.channel_to_energy(np.arange(len(y)), own[tag]["file_coeffs"])
    ref = {round(r["E_lib"]): r["mu_file"] - r["E_lib"] for r in own[tag]["refs"]}
    for g in M.GROUPS_DEF:
        gl = [l for l in (M.get_lines_in_window(*c, g["lo"] - EXT, g["hi"] + EXT) for c in g["comp"]) if l]
        m = (x >= g["lo"] - EXT) & (x <= g["hi"] + EXT); xw, yw = x[m], y[m]
        c2 = np.array([(M.solve_for_delta(xw, yw, gl, g["cont"], v)[1] or np.inf) for v in ds]); i = int(np.argmin(c2))
        nu = len(yw) - 6; ok = ds[c2 <= c2[i] * (1 + 1 / nu)]; sd = max((ok[-1] - ok[0]) / 2, 0.05)
        r = next((f"{v:+.2f}" for k, v in ref.items() if abs(k - int(g["name"])) < 3), "—") if g["name"].isdigit() else "—"
        print(f"{tag} {g['name']}: δ_файл(форма) {ds[i]:+.2f} ± {sd:.2f}{' КРАЙ' if abs(ds[i]) > 14.9 else ''}; одиночный гаусс μ−E {r}; χ²/ν {c2[i]:.2f}")
