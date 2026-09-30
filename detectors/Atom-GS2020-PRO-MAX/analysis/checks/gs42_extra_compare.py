# -*- coding: utf-8 -*-
r"""#GS-42 приёмка: таблица «было/стало» по подгонкам GS2020 (KCl М1 основной/хвост, KCl М2, торий М1/М2/М2 full, фон воды):
A, отношение к ожидаемой/паспорту, χ²/ν, строка ФОРМА ПИКОВ (из лога), строки ФИЗИКА (#GS-42) и extra_components.
Запуск: python gs42_extra_compare.py KCL_БЫЛО TH_БЫЛО KCL_СТАЛО TH_СТАЛО [--physics]"""
import sys, os, json
sys.stdout.reconfigure(encoding="utf-8")
FITS = [("KCl М1", 0, "fit_kcl_bgw.json", "fit_kcl_main.log", None),
        ("KCl М1 хвост донора", 0, "fit_kcl_bgw_tail.json", "fit_kcl_tail.log", None),
        ("KCl М2", 0, "fit_kcl_m2_bgw.json", "fit_kcl_m2.log", None),
        ("Th М1 цепочка", 1, "fit_m1_bgw_tail0_fwscale_calshape_calsum.json", "m1_bgw.log", "chain"),
        ("Th М2 цепочка", 1, "fit_m2_bgw_tail0.json", "m2_bgw.log", "chain"),
        ("Th М2 full цепочка", 1, "fit_m2_full_noThresh_bgw_tail0.json", "m2f_bgw.log", "chain")]


def one(dirs, fit):
    lab, d, js, lg, sub = fit
    j = json.load(open(os.path.join(dirs[d], js), encoding="utf-8"))
    c = j[sub] if sub else j
    ref = j["passport_Bq"] if sub else j["A_expected_Bq"]
    txt = open(os.path.join(dirs[d], lg), encoding="utf-8").read().splitlines()
    shape = [ln.split("):", 1)[1].strip() for ln in txt if ln.startswith("ФОРМА ПИКОВ")]
    phys = [ln for ln in txt if ln.startswith("ФИЗИКА (#GS-42)")]
    return c["A_Bq"], c["A_Bq"] / ref, c["chi2"] / c["ndof"], (shape[0] if shape else "—"), phys, j.get("extra_components")


def physics(rows):
    for lab, ph, ex in rows:
        print("\n## %s: строк ФИЗИКА (#GS-42) %d" % (lab, len(ph)))
        for ln in ph:
            print("  " + ln.replace("C:\\g4work\\gs2020\\", ""))
        for k, v in sorted((ex or {}).items()):
            if isinstance(v, dict) and "br" in v:
                print("  %s BR %.4f: " % (k, v["br"]) + "; ".join("%s окно %.1f, <150 %.1f, 150–400 %.1f, >400 %.1f" % (
                    p, v[p]["window"], v[p]["lt150"], v[p]["150_400"], v[p]["gt400"]) for p in ("beta", "ib") if p in v))


def main():
    a, b = sys.argv[1:3], sys.argv[3:5]
    print("| подгонка | A было, Бк | A стало, Бк | ΔA/A | отн. было | отн. стало | χ²/ν было | χ²/ν стало | ФОРМА ПИКОВ было | ФОРМА ПИКОВ стало |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    worst, same, rows = 0.0, True, []
    for fit in FITS:
        A0, r0, x0, s0, _, _ = one(a, fit)
        A1, r1, x1, s1, ph, ex = one(b, fit)
        worst = max(worst, abs(A1 - A0) / A0, abs(x1 - x0) / x0)
        same = same and (s0 == s1 or s0 == "—")
        rows.append((fit[0], ph, ex))
        print("| %s | %.3f | %.3f | %+.2e | %.4f | %.4f | %.4f | %.4f | %s | %s |" % (fit[0], A0, A1, (A1 - A0) / A0, r0, r1, x0, x1, s0, s1))
    print("наибольшее относительное расхождение A и χ²/ν: %.3e; строки ФОРМА ПИКОВ совпали (где были в «было»): %s" % (worst, same))
    if "--physics" in sys.argv:
        physics(rows)


if __name__ == "__main__":
    main()
