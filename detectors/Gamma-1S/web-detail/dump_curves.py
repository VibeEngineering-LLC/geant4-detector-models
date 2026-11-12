"""Кривые подгонки fit_kit2024 (Th232 маринелли) для страницы детализации: измерение, фон, модель каждого набора.
Usage: SPECTRAVIBE_ROOT=... python dump_curves.py out.json k24g1t dtF dtP ... (перехват core.unfold, метод не меняется)"""
import json, os, sys
sys.path.insert(0, os.environ.get("FIT_SCRIPTS_DIR", "."))  # folder with fit_kit2024.py
import fit_kit2024 as fk
cap, orig = {}, fk.core.unfold
fk.core.unfold = lambda *a, **k: cap.setdefault("r", orig(*a, **k))
row = next(r for r in fk.SPECTRA if r[1].startswith("Th232_420-7-17_Маринелли"))
fk.write_fwhm_points(); out = {"model": {}, "A": {}}
for p in sys.argv[2:]:
    os.environ["TPL_PREFIX"] = p; cap.clear(); res = fk.fit_one(row, fk.BLUR_2024, verbose=False); r = cap["r"]
    m = (r["cols"] * r["coef"][:, None]).sum(axis=0) + r["bg_scaled"]
    out["model"][p] = [round(float(x), 2) for x in m]; out["A"][p] = {g: v["A_over_passport"] for g, v in res["groups"].items()}
out["e"] = [round(float(x), 2) for x in r["e"]]; out["counts"] = [float(x) for x in r["counts"]]
out["bg"] = [round(float(x), 2) for x in r["bg_scaled"]]; out["live_s"] = res["live_s"]
idx = [i for i, s in enumerate(r["sel"]) if s]; out["sel"] = [idx[0], idx[-1]]
json.dump(out, open(sys.argv[1], "w", encoding="utf-8"), ensure_ascii=False)
print({p: out["A"][p] for p in out["A"]}, len(out["e"]), out["sel"])
