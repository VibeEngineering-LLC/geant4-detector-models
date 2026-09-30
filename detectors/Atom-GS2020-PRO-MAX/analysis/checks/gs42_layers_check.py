# -*- coding: utf-8 -*-
r"""GS-42 п.2 приёмка: поканальные слои BETA/IB (extra_components фитов + стек страниц).
(а) регрессия 9 JSON против pre_layers_2026-09-29; (б) Σ(*_col[sel]) == band_counts window;
(в) доля BETA/IB в окне на слой, независимый пересчёт из extra_components (без общего кода с export_*).
Запуск: PYTHONIOENCODING=utf-8 python scripts\checks\gs42_layers_check.py"""
import os
import sys
import json
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

KCL_DIR = r"C:\g4work\gs2020\kcl_v3w85"
TH_DIR = r"C:\g4work\gs2020\run_marinelli\out_v5_oisn10"
BACKUP = "pre_layers_2026-09-29"

KCL_FILES = ["fit_kcl_bgw.json", "fit_kcl_bgw_tail.json", "fit_kcl_m2_bgw.json"]
TH_FILES = ["fit_m1_bgw_tail0_fwscale_calshape_calsum.json", "fit_m1_tail0_fwscale_calshape_calsum.json",
            "fit_m2_bgw_tail0.json", "fit_m2_tail0.json",
            "fit_m2_full_noThresh_bgw_tail0.json", "fit_m2_full_noThresh_tail0.json"]


def load(p):
    return json.load(open(p, encoding="utf-8"))


def _is_num_list(v):
    return isinstance(v, list) and (len(v) == 0 or isinstance(v[0], (int, float)))


def compare_regression(cur, backup, path=""):
    """(a) Ключи, БЫВШИЕ в backup, сверяются с cur побитово (числовые массивы — векторно, rel <= 1e-12);
    новые ключи cur (ib_col/beta_col и т.п.) не трогаем. Возвращает (n_сравненных_значений, [расхождения])."""
    n, bad = 0, []
    if isinstance(backup, dict):
        for k, v in backup.items():
            if k not in cur:
                bad.append(path + "." + k + ": ключ пропал"); continue
            nn, bb = compare_regression(cur[k], v, path + "." + k)
            n += nn; bad += bb
        return n, bad
    if _is_num_list(backup):
        b = np.asarray(backup, dtype=float)
        ok_shape = isinstance(cur, list) and len(cur) == len(b)
        a = np.asarray(cur, dtype=float) if ok_shape else None
        if not ok_shape:
            return len(b), [path + ": длина изменилась (%d -> %s)" % (len(b), len(cur) if isinstance(cur, list) else type(cur))]
        rel = np.abs(a - b) / np.maximum(np.abs(b), 1e-300)
        bad_mask = (a != b) & (rel > 1e-12)
        if bad_mask.any():
            bad.append(path + ": %d/%d элементов разошлись, макс отн. |Δ|=%.3e" % (int(bad_mask.sum()), len(b), float(rel[bad_mask].max())))
        return len(b), bad
    if isinstance(backup, list):   # список не-чисел (напр. bool) — поэлементно, списки в JSON здесь не вложены глубоко
        if not isinstance(cur, list) or len(cur) != len(backup):
            return 1, [path + ": длина списка изменилась"]
        n, bad = 0, []
        for i, (a, b) in enumerate(zip(cur, backup)):
            nn, bb = compare_regression(a, b, path + "[%d]" % i)
            n += nn; bad += bb
        return n, bad
    if isinstance(backup, float):
        a = float(cur)
        if a != backup and (backup == 0 or abs(a - backup) / abs(backup) > 1e-12):
            return 1, [path + ": %r != %r" % (cur, backup)]
        return 1, []
    if cur != backup:
        return 1, [path + ": %r != %r" % (cur, backup)]
    return 1, []


def sel_of(j):
    """Окно подгонки: берём готовое поле "sel" (Th-232 м1/м2), а если его нет (KCl — плоский JSON
    без "sel") — тот же (e>=lo)&(e<=hi), что строит mix_unfold_core.unfold."""
    if "sel" in j:
        return np.asarray(j["sel"], dtype=bool)
    e = np.asarray(j["e"], dtype=float)
    return (e >= j["lo"]) & (e <= j["hi"])


def check_new_arrays(j, model, extra, label):
    """(б) Для каждого *_col в extra_components: длина == len(model), Σ по sel == band_counts[...]["window"]
    (отн. <= 1e-9). Возвращает (n_проверенных_массивов, [расхождения])."""
    n_ok, problems = 0, []
    sel = sel_of(j)
    for k, d in (extra or {}).items():
        if not isinstance(d, dict) or "br" not in d:
            continue   # служебные ключи extra_components (amplitude, bands_keV, var_excess_max, chain_only)
        for p in ("beta", "ib"):
            if (p + "_col") not in d:
                continue
            col = np.asarray(d[p + "_col"], dtype=float)
            if len(col) != len(model):
                problems.append("%s/%s/%s_col: длина %d != model %d" % (label, k, p, len(col), len(model)))
                continue
            win_saved = float(d[p]["window"])
            win_calc = float(col[sel].sum())
            rel = abs(win_calc - win_saved) / abs(win_saved) if win_saved != 0 else abs(win_calc)
            n_ok += 1
            if rel > 1e-9:
                problems.append("%s/%s/%s_col: Σ[sel]=%.6f, band_counts.window=%.6f, отн.расхожд=%.3e"
                                 % (label, k, p, win_calc, win_saved, rel))
    return n_ok, problems


def layer_shares(j, label):
    """(в) Доля BETA/IB в окне — прямым суммированием band_counts["window"] по всем нуклидам extra_components
    (независимый пересчёт из сырых JSON, без кода export-скриптов)."""
    chain = j.get("chain")
    model = chain["model"] if chain else j["model"]
    sel = sel_of(j)
    extra = j.get("extra_components")   # extra_components — ключ ВЕРХНЕГО уровня (сосед "chain"), не внутри chain
    model_win = float(np.asarray(model, dtype=float)[sel].sum())
    beta_win = ib_win = 0.0
    nucs = []
    for k, d in (extra or {}).items():
        if not isinstance(d, dict) or "br" not in d:
            continue
        nucs.append(k)
        beta_win += float(d["beta"]["window"]) if "beta" in d else 0.0
        ib_win += float(d["ib"]["window"]) if "ib" in d else 0.0
    pct = (lambda x: 100.0 * x / model_win) if model_win else (lambda x: float("nan"))
    print("  %s: модель в окне %.3f; нуклидов с extra %d (%s); BETA %.3f (%.4f%%); IB %.3f (%.4f%%)"
          % (label, model_win, len(nucs), ",".join(nucs), beta_win, pct(beta_win), ib_win, pct(ib_win)))
    return model_win, beta_win, ib_win


ALL_FILES = [(KCL_DIR, n) for n in KCL_FILES] + [(TH_DIR, n) for n in TH_FILES]


def section_a():
    print("=== (а) РЕГРЕССИЯ (9 JSON против %s) ===" % BACKUP)
    total_n, total_bad = 0, []
    for base, name in ALL_FILES:
        cur = load(os.path.join(base, name))
        bak = load(os.path.join(base, BACKUP, name))
        n, bad = compare_regression(cur, bak, name)
        total_n += n; total_bad += bad
        print("  %s: сравнено значений %d, расхождений %d" % (name, n, len(bad)))
    print("ИТОГО (а): сравнено %d, расхождений %d" % (total_n, len(total_bad)))
    for b in total_bad[:20]:
        print("    ! " + b)
    return total_n, total_bad


def section_c():
    print("\n=== (в) ДОЛЯ BETA/IB В ОКНЕ (независимый пересчёт из extra_components каждого JSON) ===")
    for base, name in ALL_FILES:
        layer_shares(load(os.path.join(base, name)), name)


def section_b():
    print("\n=== (б) ПОКАНАЛЬНЫЕ МАССИВЫ *_col: Σ[sel] == band_counts.window, длина == len(model) ===")
    total_n, total_bad = 0, []
    for base, name in ALL_FILES:
        j = load(os.path.join(base, name))
        chain = j.get("chain")
        model = chain["model"] if chain else j["model"]
        extra = j.get("extra_components")   # extra_components — ключ ВЕРХНЕГО уровня, не внутри chain
        n_ok, bad = check_new_arrays(j, model, extra, name)
        total_n += n_ok; total_bad += bad
        print("  %s: проверено массивов %d, расхождений %d" % (name, n_ok, len(bad)))
    print("ИТОГО (б): проверено %d, расхождений %d" % (total_n, len(total_bad)))
    for b in total_bad[:20]:
        print("    ! " + b)
    return total_n, total_bad


def main():
    n_a, bad_a = section_a()
    n_b, bad_b = section_b()
    section_c()
    ok = not bad_a and not bad_b
    print("\nРЕЗУЛЬТАТ: %s — (а) сравнено %d/расхождений %d; (б) проверено %d/расхождений %d"
          % ("ЗЕЛЕНО" if ok else "КРАСНО", n_a, len(bad_a), n_b, len(bad_b)))
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
