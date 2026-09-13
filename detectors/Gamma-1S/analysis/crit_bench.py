# -*- coding: utf-8 -*-
"""crit_bench.py — сравнительный стенд критериев подгонки шаблонов (#CRIT-1).

Спека: GEANT4/scripts/_spec_crit_bench.md. Ядро (mix_unfold_core) сюда не импортируется и не
копируется: стенд получает готовые r["cols"], counts, bg_scaled, k, n_events, sel.
Мутационные условия (спека п.6): проверка 2 обязана краснеть, если A1 реализовать с обрезкой
max(p-b, 0); проверка 3 — если в cash потерять слагаемое -p*ln(m).
"""
import sys
import os
import json
import argparse
import numpy as np
from scipy.optimize import nnls, minimize, linprog

sys.stdout.reconfigure(encoding="utf-8")

BANDS = [(25, 37), (37, 50), (50, 70), (70, 100), (100, 300), (300, 1500)]
CRIT_ORDER = ["A0", "A1", "A2", "B", "C", "D", "E1", "E2"]


def _nnls_w(C, y, w):
    """Взвешенный NNLS: a >= 0; sd по активным столбцам (a > 0), nan у прижатых к нулю."""
    a, _ = nnls(C * w[:, None], y * w)
    act = a > 0
    sd = np.full(len(a), np.nan)
    if act.any():
        Cw = C[:, act] * w[:, None]
        cov = np.linalg.inv(Cw.T @ Cw)
        sd[act] = np.sqrt(np.diag(cov))
    return a, sd


def _prep(cols, counts, bg_scaled, sel):
    """C = cols[:, sel].T (n x K), p = counts[sel], b = bg_scaled[sel]."""
    sel = np.asarray(sel, dtype=bool)
    C = np.asarray(cols, dtype=float)[:, sel].T
    p = np.asarray(counts, dtype=float)[sel]
    b = np.asarray(bg_scaled, dtype=float)[sel]
    return C, p, b


def fit_A0(cols, counts, bg_scaled, k, n_events, sel):
    """A0 — нынешний: нетто с обрезкой нулём, веса 1/sqrt(max(y, 1))."""
    C, p, b = _prep(cols, counts, bg_scaled, sel)
    y = np.maximum(p - b, 0)
    w = 1.0 / np.sqrt(np.maximum(y, 1.0))
    a, sd = _nnls_w(C, y, w)
    return a, sd, {}


def fit_A1(cols, counts, bg_scaled, k, n_events, sel):
    """A1 — нетто БЕЗ обрезки; var = p + k*b, т.к. var(k*фон) = k^2*фон = k*b."""
    C, p, b = _prep(cols, counts, bg_scaled, sel)
    y = p - b
    var = np.maximum(p + k * b, 1.0)
    w = 1.0 / np.sqrt(var)
    a, sd = _nnls_w(C, y, w)
    return a, sd, {"var": var}


def fit_B(cols, counts, bg_scaled, k, n_events, sel):
    """B — относительные веса 1/max(y, 1), нетто с обрезкой нулём."""
    C, p, b = _prep(cols, counts, bg_scaled, sel)
    y = np.maximum(p - b, 0)
    w = 1.0 / np.maximum(y, 1.0)
    a, sd = _nnls_w(C, y, w)
    return a, sd, {}


def _iter_template_var(C, y, var, C_tmpl, n_events, a0, n_iter=6):
    """Итерации: w = 1/sqrt(var + model_var), model_var = (C_tmpl/n_events) @ a_tmpl**2.
    C_tmpl — первые K столбцов C (шаблоны); столбцы подложки D дисперсии не имеют."""
    n_ev = np.asarray(n_events, dtype=float)[None, :]
    K = C_tmpl.shape[1]
    a = np.asarray(a0, dtype=float)
    sd = None
    for _ in range(n_iter):
        model_var = (C_tmpl / n_ev) @ (a[:K] ** 2)
        w = 1.0 / np.sqrt(var + model_var)
        a, sd = _nnls_w(C, y, w)
    return a, sd


def fit_A2(cols, counts, bg_scaled, k, n_events, sel):
    """A2 — A1 + дисперсия шаблонов var_k = cols_k/n_events_k, 6 итераций, старт с A1."""
    C, p, b = _prep(cols, counts, bg_scaled, sel)
    a1, _, ex1 = fit_A1(cols, counts, bg_scaled, k, n_events, sel)
    a, sd = _iter_template_var(C, p - b, ex1["var"], C, n_events, a1)
    return a, sd, {"iterations": 6}


def cash_stat(a, C, p, b):
    """Статистика Кэша: 2*sum(m - p*ln m), m = C@a + b (фон в модели, не вычитается)."""
    m = C @ a + b
    return 2.0 * np.sum(m - p * np.log(np.maximum(m, 1e-300)))


def cash_grad(a, C, p, b):
    """Аналитический градиент cash по a: 2*C.T @ (1 - p/m)."""
    m = np.maximum(C @ a + b, 1e-300)
    return 2.0 * (C.T @ (1.0 - p / m))


def cash_shifted(a, C, p, b):
    """Сдвинутая поканальная форма Кэша: 2*sum[(m - p) - p*ln(m/p)] (насыщенная модель вычтена
    ВНУТРИ суммы; p = 0 даёт слагаемое 0). Тот же минимум, что у cash_stat, без потери точности."""
    m = np.maximum(C @ a + b, 1e-300)
    ps = np.where(p > 0, p, 1.0)
    return 2.0 * np.sum((m - p) - np.where(p > 0, p * np.log(m / ps), 0.0))


def _cash_hess_analytic(a, act, C, p, b):
    """Аналитический гессиан cash по активным параметрам: 2*C.T @ diag(p/m^2) @ C."""
    m = np.maximum(C @ a + b, 1e-300)
    Ca = C[:, act]
    return 2.0 * (Ca.T @ (Ca * (p / m ** 2)[:, None]))


def _hess_cd(gradf, a, act):
    """Численный гессиан по активным параметрам: центральные разности градиента gradf(a),
    шаг 1e-3*max(a_k, 1)."""
    idx = np.flatnonzero(act)
    H = np.zeros((len(idx), len(idx)))
    for j, jj in enumerate(idx):
        h = 1e-3 * max(a[jj], 1.0)
        ap, am = a.copy(), a.copy()
        ap[jj] += h
        am[jj] -= h
        H[:, j] = (gradf(ap)[idx] - gradf(am)[idx]) / (2.0 * h)
    return 0.5 * (H + H.T)


def _cash_hessian(a, act, C, p, b):
    """Численный гессиан cash (спека C)."""
    return _hess_cd(lambda x: cash_grad(x, C, p, b), a, act)


def fit_C(cols, counts, bg_scaled, k, n_events, sel):
    """C — пуассоновская правдоподобность (Cash). Обход координатора 12.09 (часть спеки):
    параметр x = a/sd_A1, минимизируется cash_shifted с аналитическим градиентом; метод L-BFGS-B,
    границы a >= 0, старт A1, options — по спеке. sd — численный гессиан (спека), рядом аналитический."""
    C, p, b = _prep(cols, counts, bg_scaled, sel)
    a1, sd1, _ = fit_A1(cols, counts, bg_scaled, k, n_events, sel)
    res, a = _minimize_cash(a1, _c_scale(a1, sd1), C, p, b)
    K = C.shape[1]
    act = a > 0
    sd, sd_an = np.full(K, np.nan), np.full(K, np.nan)
    if act.any():
        sd[act] = np.sqrt(np.diag(2.0 * np.linalg.inv(_cash_hessian(a, act, C, p, b))))
        sd_an[act] = np.sqrt(np.diag(2.0 * np.linalg.inv(_cash_hess_analytic(a, act, C, p, b))))
    rel = float(np.max(np.abs(sd[act] - sd_an[act]) / sd_an[act])) if act.any() else float("nan")
    extra = {"cash": float(cash_stat(a, C, p, b)), "cash_shifted": float(res.fun), "nit": int(res.nit),
             "success": bool(res.success), "message": str(res.message),
             "sd_analytic": sd_an, "sd_hess_reldiff": rel}
    return a, sd, extra


def _c_scale(a1, sd1):
    """Масштаб параметров C: sd_A1; у прижатых к нулю (sd nan) — max(a, 1)."""
    return np.where(np.isfinite(sd1) & (sd1 > 0), sd1, np.maximum(a1, 1.0))


def _minimize_cash(a_start, scale, C, p, b):
    """L-BFGS-B в масштабе x = a/scale по сдвинутой форме Кэша; отказ — SystemExit (спека)."""
    f = lambda x: cash_shifted(x * scale, C, p, b)
    g = lambda x: scale * cash_grad(x * scale, C, p, b)
    res = minimize(f, np.asarray(a_start, dtype=float) / scale, jac=g, method="L-BFGS-B",
                   bounds=[(0, None)] * len(scale), options={"maxiter": 2000, "ftol": 1e-12})
    if not res.success:
        raise SystemExit("ОТКАЗ оптимизатора C (L-BFGS-B): %s" % res.message)
    return res, np.asarray(res.x, dtype=float) * scale


def make_fit_D(e):
    """D — A2 + свободная подложка (константа и e/1000, обе >= 0 через NNLS). Энергии каналов
    e приходят замыканием, чтобы сохранить единую сигнатуру fit_X(...)."""
    e = np.asarray(e, dtype=float)

    def fit_D(cols, counts, bg_scaled, k, n_events, sel):
        C, p, b = _prep(cols, counts, bg_scaled, sel)
        K = C.shape[1]
        es = e[np.asarray(sel, dtype=bool)]
        C2 = np.hstack([C, np.ones((len(es), 1)), (es / 1000.0)[:, None]])
        a1, _, ex1 = fit_A1(cols, counts, bg_scaled, k, n_events, sel)
        a0 = np.concatenate([a1, [0.0, 0.0]])
        a_all, sd_all = _iter_template_var(C2, p - b, ex1["var"], C, n_events, a0)
        tot = float(np.sum(C2 @ a_all))
        bg_part = float(np.sum(C2[:, K:] @ a_all[K:]))
        extra = {"iterations": 6, "bg_const": float(a_all[K]), "bg_slope": float(a_all[K + 1]),
                 "bg_share_pct": 100.0 * bg_part / tot if tot > 0 else float("nan")}
        return a_all[:K], sd_all[:K], extra
    return fit_D


def _lp_tv(C, y):
    """E1: LP  min Σt  при  t ≥ ±(C a − y),  Σ_i (C a)_i = Σ y  (нормировка),  a ≥ 0, t ≥ 0.
    Решатель highs. Возвращает (a, Σt, res)."""
    n, K = C.shape
    I = np.eye(n)
    A_ub = np.vstack([np.hstack([C, -I]), np.hstack([-C, -I])])
    b_ub = np.concatenate([y, -y])
    A_eq, b_eq = np.concatenate([C.sum(0), np.zeros(n)])[None, :], np.array([y.sum()])
    res = linprog(np.concatenate([np.zeros(K), np.ones(n)]), A_ub=A_ub, b_ub=b_ub,
                  A_eq=A_eq, b_eq=b_eq, bounds=[(0, None)] * (K + n), method="highs")
    if not res.success:
        raise SystemExit("ОТКАЗ LP E1 (highs): %s" % res.message)
    return np.asarray(res.x[:K], dtype=float), float(np.sum(res.x[K:])), res


def c_acceptance(cols, counts, bg_scaled, k, n_events, sel):
    """Приёмка обхода C числами: минимум с трёх стартов (0,95/1,05/1,20 a_A1) и сверка sd
    численного гессиана с аналитическим 2*C.T diag(p/m^2) C. Гейты — у вызывающего."""
    C, p, b = _prep(cols, counts, bg_scaled, sel)
    a1, sd1, _ = fit_A1(cols, counts, bg_scaled, k, n_events, sel)
    scale = _c_scale(a1, sd1)
    fs, As = [], []
    for fac in (0.95, 1.05, 1.20):
        res, a = _minimize_cash(fac * a1, scale, C, p, b)
        fs.append(float(res.fun))
        As.append(a)
    fs, As = np.array(fs), np.array(As)
    f_rel = float((fs.max() - fs.min()) / max(np.abs(fs).max(), 1e-300))
    a_rel = float(np.max((As.max(0) - As.min(0)) / np.maximum(np.abs(As).max(0), 1e-300)))
    a, sd, ex = fit_C(cols, counts, bg_scaled, k, n_events, sel)
    return {"f_starts": fs, "f_rel_spread": f_rel, "a_starts": As, "a_rel_spread": a_rel,
            "sd_numeric": sd, "sd_analytic": ex["sd_analytic"], "sd_rel_diff": ex["sd_hess_reldiff"]}


def make_fit_E1(n_boot=25, seed=20260912):
    """E1 — полная вариация по нормированным формам, точное LP с ограничением Σm = S.
    sd — параметрический бутстрап n_boot реплик p* = Poisson(m + b) с тем же b (в бою 25, в closure 10)."""
    def fit_E1(cols, counts, bg_scaled, k, n_events, sel):
        import time
        C, p, b = _prep(cols, counts, bg_scaled, sel)
        y = np.maximum(p - b, 0.0)
        t0 = time.perf_counter()
        a, sum_t, res = _lp_tv(C, y)
        dt = time.perf_counter() - t0
        rng, m = np.random.default_rng(seed), C @ a
        boot = np.array([_lp_tv(C, np.maximum(rng.poisson(m + b) - b, 0.0))[0] for _ in range(n_boot)])
        sd = boot.std(axis=0, ddof=1)
        sd[~(a > 0)] = np.nan
        extra = {"tv": 0.5 * sum_t / float(y.sum()), "sd_method": "bootstrap%d" % n_boot,
                 "lp_time_s": dt, "lp_nit": int(res.nit), "model_sum_over_S": float(m.sum() / y.sum())}
        return a, sd, extra
    return fit_E1


fit_E1 = make_fit_E1(25)


def metrics_all(a, cols, counts, bg_scaled, k, sel, e, names, passport, live_s, n_params=None):
    """Единые метрики для любого критерия: A/паспорт, chi2_ref/nu (дисперсия A1), shape, полосы.
    nu = n − число подгоняемых параметров ДАННОГО критерия (координатор 12.09): по умолчанию K,
    у D — K + 2 (константа и наклон подложки)."""
    C, p, b = _prep(cols, counts, bg_scaled, sel)
    a = np.asarray(a, dtype=float)
    n, K = C.shape
    nu = n - (K if n_params is None else int(n_params))
    pred, net = C @ a, p - b
    y = np.maximum(net, 0.0)
    chi2_ref = float(np.sum((pred - net) ** 2 / np.maximum(p + k * b, 1.0)))
    sp, sy = float(np.sum(pred)), float(np.sum(y))
    shape = 0.5 * float(np.sum(np.abs(pred / sp - y / sy))) if sp > 0 and sy > 0 else float("nan")
    es = np.asarray(e, dtype=float)[np.asarray(sel, dtype=bool)]
    bands = {}
    for lo, hi in BANDS:
        m = (es >= lo) & ((es < hi) | ((hi == BANDS[-1][1]) & (es <= hi)))
        den = float(np.sum(y[m]))
        bands["%d-%d" % (lo, hi)] = float(np.sum(pred[m])) / den if den > 0 else float("nan")
    aop = {nm: float(a[i] / live_s / passport[nm]) for i, nm in enumerate(names)}
    return {"A_over_passport": aop, "chi2_ref": chi2_ref, "nu": int(nu), "chi2_ref_nu": chi2_ref / nu,
            "shape": shape, "bands": bands}


def _chi2_shape(a, C, y):
    """E2: Σ (C a − y)² / max(y, 1)."""
    return float(np.sum((C @ a - y) ** 2 / np.maximum(y, 1.0)))


def _chi2_shape_grad(a, C, y):
    return 2.0 * (C.T @ ((C @ a - y) / np.maximum(y, 1.0)))


def fit_E2(cols, counts, bg_scaled, k, n_events, sel):
    """E2 — χ² по форме с жёсткой нормировкой Σm = S: SLSQP, старт A1, a ≥ 0; sd — численный
    гессиан целевой по активным параметрам (без учёта ограничения — по спеке), cov = 2·inv(H).
    Параметризация x = a/sd_A1 и целевая, делённая на f(a_A1): SLSQP чувствителен к масштабу."""
    C, p, b = _prep(cols, counts, bg_scaled, sel)
    y, csum = np.maximum(p - b, 0.0), C.sum(0)
    S = float(y.sum())
    a1, sd1, _ = fit_A1(cols, counts, bg_scaled, k, n_events, sel)
    scale, f0 = _c_scale(a1, sd1), max(_chi2_shape(a1, C, y), 1e-300)
    return _e2_solve(a1, scale, f0, C, y, csum, S)


def _synth_rep(rng, cols, n_ev, mu, bg_raw, k, template_noise):
    """Один повтор closure: фон и проба независимо; шум шаблонов — по желанию."""
    b_synth = k * rng.poisson(bg_raw)
    counts_synth = rng.poisson(mu).astype(float)
    cols_fit = cols
    if template_noise:
        cols_fit = rng.poisson(cols * n_ev[:, None]) / n_ev[:, None]
    return b_synth, counts_synth, cols_fit


def _e2_solve(a1, scale, f0, C, y, csum, S):
    K = len(a1)
    cons = [{"type": "eq", "fun": lambda x: (csum @ (x * scale) - S) / S, "jac": lambda x: csum * scale / S}]
    res = minimize(lambda x: _chi2_shape(x * scale, C, y) / f0, a1 / scale,
                   jac=lambda x: scale * _chi2_shape_grad(x * scale, C, y) / f0, method="SLSQP",
                   constraints=cons, bounds=[(0, None)] * K, options={"maxiter": 2000, "ftol": 1e-14})
    if not res.success:
        raise SystemExit("ОТКАЗ оптимизатора E2 (SLSQP): %s" % res.message)
    a = np.asarray(res.x, dtype=float) * scale
    act = a > 0
    sd = np.full(K, np.nan)
    if act.any():
        H = _hess_cd(lambda x: _chi2_shape_grad(x, C, y), a, act)
        sd[act] = np.sqrt(np.diag(2.0 * np.linalg.inv(H)))
    extra = {"chi2_shape": _chi2_shape(a, C, y), "nit": int(res.nit), "success": bool(res.success),
             "message": str(res.message), "model_sum_over_S": float(csum @ a / S)}
    return a, sd, extra


def closure(cols, counts, bg_scaled, k, n_events, sel, a_true, n_rep=100, seed=20260912,
            template_noise=False, fits=None):
    """Замкнутый тест на синтетике. fits — {имя: fit_X}; None — пять критериев без D (D нужен e,
    тогда fits передать явно). Генерация — по точным cols и истинному ожидаемому фону."""
    if fits is None:
        fits = {"A0": fit_A0, "A1": fit_A1, "A2": fit_A2, "B": fit_B, "C": fit_C}
    rng = np.random.default_rng(seed)
    cols, n_ev = np.asarray(cols, dtype=float), np.asarray(n_events, dtype=float)
    a_true = np.asarray(a_true, dtype=float)
    bg_raw = np.asarray(bg_scaled, dtype=float) / k
    mu = cols.T @ a_true + k * bg_raw
    est = {nm: np.zeros((n_rep, len(a_true))) for nm in fits}
    err = {nm: np.zeros((n_rep, len(a_true))) for nm in fits}
    for i in range(n_rep):
        b_synth, counts_synth, cols_fit = _synth_rep(rng, cols, n_ev, mu, bg_raw, k, template_noise)
        for nm, f in fits.items():
            a_hat, sd_hat, _ = f(cols_fit, counts_synth, b_synth, k, n_events, sel)
            est[nm][i], err[nm][i] = a_hat, sd_hat
    return _closure_stats(est, err, a_true, n_rep)


def _closure_stats(est, err, a_true, n_rep):
    """bias_pct, spread_pct (std, ddof=1), coverage_68, pull_rms по повторам с конечным sd."""
    out = {}
    for nm in est:
        rows = []
        for j in range(len(a_true)):
            ah, sh = est[nm][:, j], err[nm][:, j]
            fin = np.isfinite(sh) & (sh > 0)
            dev = ah - a_true[j]
            pull = dev[fin] / sh[fin]
            rows.append({"bias_pct": float(100.0 * (ah.mean() / a_true[j] - 1.0)),
                         "spread_pct": float(100.0 * ah.std(ddof=1) / a_true[j]),
                         "coverage_68": float(np.mean(np.abs(dev) <= np.where(fin, sh, 0.0))),
                         "sd_calib_68": float(np.mean(np.abs(ah - ah.mean()) <= np.where(fin, sh, 0.0))),
                         "pull_rms": float(np.sqrt(np.mean(pull ** 2))) if fin.any() else float("nan"),
                         "n_sd_finite": int(fin.sum()), "n_rep": int(n_rep)})
        out[nm] = rows
    return out


def _flags_closure(tag, cl, n_rep, names):
    flags = []
    for nm, rows in cl.items():
        for j, row in enumerate(rows):
            if not (0.55 <= row["coverage_68"] <= 0.80):
                flags.append("closure %s: %s/%s coverage_68 = %.2f вне [0,55; 0,80]"
                             % (tag, nm, names[j], row["coverage_68"]))
            lim = 2.0 * row["spread_pct"] / np.sqrt(n_rep)
            if abs(row["bias_pct"]) > lim:
                flags.append("closure %s: %s/%s |bias| = %.3f %% > 2*spread/sqrt(n_rep) = %.3f %%"
                             % (tag, nm, names[j], row["bias_pct"], lim))
    return flags


def interpretation_flags(real, cl0, cl1, n_rep, names):
    """Блок «ТРЕБУЕТ ТОЛКОВАНИЯ» — правила спеки, печатается всегда (пустой — «пусто»)."""
    flags = _flags_closure("без шума шаблонов", cl0, n_rep, names)
    flags += _flags_closure("с шумом шаблонов", cl1, n_rep, names)
    exc = real["C"]["extra"]
    if not exc.get("success", False):
        flags.append("C: оптимизатор не сошёлся: %s" % exc.get("message"))
    acc = real["C"].get("acceptance")
    if acc is not None and (acc["f_rel_spread"] > 1e-6 or acc["sd_rel_diff"] > 1e-4):
        flags.append("C приёмка обхода: разброс минимума по стартам %.2e (порог 1e-6), sd числ/аналит %.2e (порог 1e-4)"
                     % (acc["f_rel_spread"], acc["sd_rel_diff"]))
    exd = real["D"]["extra"]
    if exd["bg_share_pct"] > 20:
        flags.append("D: доля подложки в модели %.1f %% > 20 %%" % exd["bg_share_pct"])
    for nm in names:
        x = real["A2"]["metrics"]["A_over_passport"][nm]
        y = real["C"]["metrics"]["A_over_passport"][nm]
        rel = 100.0 * abs(x - y) / max(abs(x), 1e-300)
        if rel > 5.0:
            flags.append("A2 vs C по %s: A/паспорт %.4f против %.4f (расхождение %.1f %% от A2)"
                         % (nm, x, y, rel))
        for en in ("E1", "E2"):
            z = real[en]["metrics"]["A_over_passport"][nm]
            rel = 100.0 * abs(x - z) / max(abs(x), 1e-300)
            if rel > 5.0:
                flags.append("%s vs A2 по %s: A/паспорт %.4f против %.4f (расхождение %.1f %% от A2)"
                             % (en, nm, z, x, rel))
    tv, sh1 = real["E1"]["extra"]["tv"], real["A1"]["metrics"]["shape"]
    if tv > sh1:
        flags.append("E1: TV в минимуме %.6f больше shape у A1 %.6f" % (tv, sh1))
    return flags


def _print_real(real, names):
    print("(1) РЕАЛЬНЫЙ СПЕКТР — A/паспорт (sd, %) по нуклидам; chi2_ref/nu; shape; полосы модель/измерение")
    hdr = "%-4s" % "крит" + "".join("%18s" % nm for nm in names) + "%12s %5s %7s" % ("chi2_ref/nu", "nu", "shape")
    print(hdr + "".join("%9s" % ("%d-%d" % bd) for bd in BANDS))
    for nm in CRIT_ORDER:
        rr, m = real[nm], real[nm]["metrics"]
        cells = ""
        for j, nuc in enumerate(names):
            a, sd = rr["a"][j], rr["sd"][j]
            sdp = 100.0 * sd / a if a > 0 and np.isfinite(sd) else float("nan")
            cells += "%10.4f (%5.2f)" % (m["A_over_passport"][nuc], sdp)
        line = "%-4s%s%12.3f %5d %7.4f" % (nm, cells, m["chi2_ref_nu"], m["nu"], m["shape"])
        print(line + "".join("%9.4f" % m["bands"]["%d-%d" % bd] for bd in BANDS))
    exc, exd = real["C"]["extra"], real["D"]["extra"]
    print("  C: cash = %.3f (сдвинутая форма %.6f), nit = %d, success = %s, sd числ/аналит расх. %.2e"
          % (exc["cash"], exc["cash_shifted"], exc["nit"], exc["success"], exc["sd_hess_reldiff"]))
    acc = real["C"]["acceptance"]
    print("  C приёмка: f(3 старта 0,95/1,05/1,20) = %s, разброс %.2e отн.; a разброс %.2e отн."
          % (np.array2string(acc["f_starts"], precision=9), acc["f_rel_spread"], acc["a_rel_spread"]))
    print("  D: bg_const = %.4g, bg_slope = %.4g, bg_share_pct = %.2f"
          % (exd["bg_const"], exd["bg_slope"], exd["bg_share_pct"]))
    e1, e2 = real["E1"]["extra"], real["E2"]["extra"]
    print("  E1: TV в минимуме = %.6f; shape(E1) по единой метрике = %.6f; shape(A1) = %.6f; Σm/S = %.9f;"
          " LP: %.3f с, nit = %d; sd = %s"
          % (e1["tv"], real["E1"]["metrics"]["shape"], real["A1"]["metrics"]["shape"],
             e1["model_sum_over_S"], e1["lp_time_s"], e1["lp_nit"], e1["sd_method"]))
    print("  E2: chi2_shape = %.3f, nit = %d, success = %s, Σm/S = %.9f"
          % (e2["chi2_shape"], e2["nit"], e2["success"], e2["model_sum_over_S"]))


def _print_closure(title, cl, names):
    print(title)
    print("%-4s %-12s %10s %10s %12s %9s %6s"
          % ("крит", "нуклид", "bias_%", "spread_%", "coverage_68", "pull_rms", "n_sd"))
    for nm in CRIT_ORDER:
        if nm not in cl:
            continue
        for j, row in enumerate(cl[nm]):
            print("%-4s %-12s %10.3f %10.3f %12.2f %9.3f %6d"
                  % (nm, names[j], row["bias_pct"], row["spread_pct"], row["coverage_68"],
                     row["pull_rms"], row["n_sd_finite"]))


def _jsonable(o):
    """numpy → встроенные типы; nan/inf → None (строгий JSON)."""
    if isinstance(o, dict):
        return {str(kk): _jsonable(v) for kk, v in o.items()}
    if isinstance(o, (list, tuple, np.ndarray)):
        return [_jsonable(v) for v in (o.tolist() if isinstance(o, np.ndarray) else o)]
    if isinstance(o, (np.floating, float)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, (np.integer, np.bool_)):
        return o.item()
    return o


def _unpack(r):
    """Вход стенда из словаря r ядра (спека, раздел «Вход стенда»)."""
    cols, sel = np.asarray(r["cols"], dtype=float), np.asarray(r["sel"], dtype=bool)
    e, names = np.asarray(r["e"], dtype=float), list(r["names"])
    counts = np.asarray(r["spec"].counts, dtype=float)
    bg_scaled = np.asarray(r["bg_scaled"], dtype=float)
    return cols, sel, e, names, counts, bg_scaled, list(r["n_events"]), float(r["live_s"])


def run(r, k, passport, out_json, n_rep=100):
    """Отчёт: (1) реальный спектр по шести критериям, (2)/(3) closure без/с шумом шаблонов,
    блок «ТРЕБУЕТ ТОЛКОВАНИЯ»; всё то же — в out_json."""
    cols, sel, e, names, counts, bg_scaled, n_events, live_s = _unpack(r)
    fits = {"A0": fit_A0, "A1": fit_A1, "A2": fit_A2, "B": fit_B, "C": fit_C, "D": make_fit_D(e),
            "E1": fit_E1, "E2": fit_E2}
    fits_cl = dict(fits, E1=make_fit_E1(10))    # в closure бутстрап E1 снижен до 10 реплик
    print("окно: %d каналов из %d, K = %d, live_s = %.1f, k = %.6f, n_rep = %d, seed = 20260912"
          % (sel.sum(), len(sel), len(names), live_s, k, n_rep))
    real = {}
    for nm in CRIT_ORDER:
        a, sd, extra = fits[nm](cols, counts, bg_scaled, k, n_events, sel)
        n_par = len(names) + 2 if nm == "D" else len(names)
        real[nm] = {"a": a, "sd": sd, "extra": extra,
                    "metrics": metrics_all(a, cols, counts, bg_scaled, k, sel, e, names, passport, live_s, n_par)}
    real["C"]["acceptance"] = c_acceptance(cols, counts, bg_scaled, k, n_events, sel)
    _print_real(real, names)
    a_true = real["A1"]["a"]
    print("closure: a_true = амплитуды A1 на реальном спектре = %s распадов"
          % np.array2string(a_true, precision=1, max_line_width=200))
    print("closure: E1 с бутстрапом 10 реплик (в бою 25)")
    return _run_closures(real, a_true, (cols, counts, bg_scaled, k, n_events, sel), names, fits_cl,
                         out_json, n_rep)


def _run_closures(real, a_true, args, names, fits, out_json, n_rep):
    import time
    t0 = time.time()
    cl0 = closure(*args, a_true, n_rep=n_rep, template_noise=False, fits=fits)
    _print_closure("(2) CLOSURE без шума шаблонов, %.0f с" % (time.time() - t0), cl0, names)
    t0 = time.time()
    cl1 = closure(*args, a_true, n_rep=n_rep, template_noise=True, fits=fits)
    _print_closure("(3) CLOSURE с шумом шаблонов, %.0f с" % (time.time() - t0), cl1, names)
    flags = interpretation_flags(real, cl0, cl1, n_rep, names)
    print("ТРЕБУЕТ ТОЛКОВАНИЯ:")
    for fl in flags or ["пусто"]:
        print("  - " + fl)
    out = {"names": names, "n_rep": n_rep, "seed": 20260912, "a_true": a_true,
           "real": real, "closure_no_template_noise": cl0, "closure_template_noise": cl1,
           "requires_interpretation": flags}
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(_jsonable(out), fh, ensure_ascii=False, indent=1)
    print("JSON: %s" % out_json)
    return out


def _synth(bg_raw_level, rng):
    """Синтетика самопроверки: K=2, N=300; гауссианы (площади 0,05 и 0,02, центры 100 и 200,
    sigma 8) + общий континуум 1e-4*exp(-ch/150); n_events [1e6, 1e6]; a_true [2e5, 1e5]; k=0,5."""
    N, k = 300, 0.5
    ch = np.arange(N, dtype=float)
    cont = 1e-4 * np.exp(-ch / 150.0)
    g = lambda mu, area: area * np.exp(-0.5 * ((ch - mu) / 8.0) ** 2) / (8.0 * np.sqrt(2 * np.pi))
    cols = np.array([g(100.0, 0.05) + cont, g(200.0, 0.02) + cont])
    a_true = np.array([2e5, 1e5])
    bg_raw = bg_raw_level + 0.0 * ch
    sel = (ch >= 20) & (ch <= 280)
    counts = rng.poisson(cols.T @ a_true + k * bg_raw).astype(float)
    return dict(cols=cols, counts=counts, bg_scaled=k * bg_raw, k=k, n_events=[1e6, 1e6],
                sel=sel, a_true=a_true, e=ch)


def _fail(msg):
    raise SystemExit("SELFTEST FAIL: " + msg)


def _neyman_bound(S):
    """D1: верхняя оценка сдвига Неймана у весов по НАБЛЮДЁННЫМ отсчётам (A1/A2): ~ -1 отсчёт на
    канал окна; в амплитуде нуклида k — n_win / sum_sel(cols_k), если весь сдвиг лёг на него.
    Измерено 12.09 (40 зёрен, bg_raw=3): A1 -1,2 %/-3,8 %, C -0,09 %/-0,04 % — п.1 спеки для A1
    на этой синтетике невыполним без допуска; C гейтится строго."""
    return S["sel"].sum() / S["cols"][:, S["sel"]].sum(axis=1)


def _check0_cash_forms(S, args):
    """п.0 (часть п.3, до первого вызова оптимизатора): обе формы Кэша модуля при a_A1 совпадают
    с независимыми формулами теста — иначе мутация в объективе C падает как отказ L-BFGS-B, а не как
    SELFTEST FAIL."""
    C, p, b = _prep(args[0], args[1], args[2], args[5])
    a1, _, _ = fit_A1(*args)
    m = C @ a1 + b
    ref_stat = 2.0 * float(np.sum(m - p * np.log(m)))
    ref_sat = 2.0 * float(np.sum(p[p > 0] - p[p > 0] * np.log(p[p > 0])))
    got_stat, got_shift = float(cash_stat(a1, C, p, b)), float(cash_shifted(a1, C, p, b))
    if not abs(got_stat - ref_stat) <= 1e-9 * abs(ref_stat):
        _fail("п.0 cash_stat при a_A1: ожидалось %.6f, получено %.6f" % (ref_stat, got_stat))
    if not abs(got_shift - (ref_stat - ref_sat)) <= 1e-7 * abs(ref_stat - ref_sat):
        _fail("п.0 cash_shifted при a_A1: ожидалось %.6f, получено %.6f" % (ref_stat - ref_sat, got_shift))


def _check1_return(S, args):
    """1. Возврат закладки: C — |a - a_true| <= 3*sd; A1, A2 — 3*sd + допуск Неймана (D1)."""
    res, ney = {}, _neyman_bound(S)
    for nm, f in (("A1", fit_A1), ("A2", fit_A2), ("C", fit_C)):
        a, sd, ex = f(*args)
        res[nm] = (a, sd, ex)
        for j in range(2):
            d, tol = abs(a[j] - S["a_true"][j]), 3 * sd[j] + (0.0 if nm == "C" else ney[j])
            if not (np.isfinite(sd[j]) and d <= tol):
                _fail("п.1 %s нуклид %d: ожидалось |a - a_true| <= 3*sd%s = %.1f, получено a = %.1f "
                      "(a_true = %.1f), |разность| = %.1f"
                      % (nm, j, "" if nm == "C" else " + n_win/sum(cols_k)", tol, a[j], S["a_true"][j], d))
    return res


def _check2_strong_bg(rng):
    """2. A0 != A1 на сильном фоне (bg_raw = 200): разница > 0,5*sd хотя бы у одного нуклида,
    и A1 при этом остаётся в 3*sd от закладки (обрезка нулём в A1 смещает амплитуду)."""
    S2 = _synth(200.0, rng)
    args2 = (S2["cols"], S2["counts"], S2["bg_scaled"], S2["k"], S2["n_events"], S2["sel"])
    a0, _, _ = fit_A0(*args2)
    a1, sd1, _ = fit_A1(*args2)
    diff = np.abs(a0 - a1) / sd1
    if not np.any(diff > 0.5):
        _fail("п.2 ожидалось |A0 - A1| > 0,5*sd хотя бы у одного нуклида, получено %s (в единицах sd)"
              % np.array2string(diff, precision=3))
    # D2: обрезка нулём в A1 меняет амплитуду лишь на 0,06/0,45 % при sd 1,2/3,9 % (измерено 12.09,
    # 40 зёрен) — амплитудный признак её не видит. Детерминированный признак: A1 обязан
    # удовлетворять нормальным уравнениям для НЕобрезанного нетто p-b с весами 1/(p+k*b)
    # по активным столбцам (a>0); у мутанта с обрезкой невязка -0,014/-0,023, у A1 ~1e-15.
    C, p, b = _prep(S2["cols"], S2["counts"], S2["bg_scaled"], S2["sel"])
    r_ = (p - b) - C @ a1
    w2 = 1.0 / np.maximum(p + S2["k"] * b, 1.0)
    kkt = (C.T @ (w2 * r_)) / (np.abs(C.T) @ (w2 * np.abs(r_)))
    if not np.all(np.abs(kkt[a1 > 0]) <= 1e-8):
        _fail("п.2 A1 при bg_raw = 200 не удовлетворяет нормальным уравнениям для необрезанного "
              "нетто: ожидалось |невязка/масштаб| <= 1e-8, получено %s" % np.array2string(kkt, precision=3))


def _check3_cash_min(S, args, res):
    """3. C — минимум правдоподобия. Опорная формула — СВОЯ у теста (не cash модуля), иначе
    мутация в cash не видна; плюс сверка cash модуля с формулой теста при a_A1."""
    C, p, b = _prep(args[0], args[1], args[2], args[5])
    ref = lambda a: 2.0 * float(np.sum((C @ a + b) - p * np.log(C @ a + b)))
    a1s, aCs = res["A1"][0], res["C"][0]
    mod = float(cash_stat(a1s, C, p, b))
    if not abs(mod - ref(a1s)) <= 1e-9 * abs(ref(a1s)):
        _fail("п.3 cash модуля при a_A1 не совпал с формулой теста: ожидалось %.6f, получено %.6f"
              % (ref(a1s), mod))
    if not ref(aCs) <= ref(a1s):
        _fail("п.3 ожидалось cash(a_C) <= cash(a_A1), получено %.6f > %.6f" % (ref(aCs), ref(a1s)))
    # 3б. Приёмка обхода C (координатор 12.09): три старта и гессиан — числами.
    acc = c_acceptance(*args)
    if not acc["f_rel_spread"] <= 1e-6:
        _fail("п.3б минимум C с трёх стартов расходится: ожидалось <= 1e-6 отн., получено %.3e (f = %s)"
              % (acc["f_rel_spread"], np.array2string(acc["f_starts"], precision=9)))
    if not acc["sd_rel_diff"] <= 1e-4:
        _fail("п.3б sd численного гессиана против аналитического: ожидалось <= 1e-4 отн., получено %.3e"
              % acc["sd_rel_diff"])
    return acc


def _check4_metrics(S, args, res):
    """4. Единая метрика: chi2_ref бит в бит при двух вызовах; shape в [0, 1]."""
    passport, nm = {"n0": 1.0, "n1": 1.0}, ["n0", "n1"]
    m1 = metrics_all(res["A1"][0], args[0], args[1], args[2], S["k"], S["sel"], S["e"], nm, passport, 1.0)
    m2 = metrics_all(res["A1"][0], args[0], args[1], args[2], S["k"], S["sel"], S["e"], nm, passport, 1.0)
    if m1["chi2_ref"] != m2["chi2_ref"]:
        _fail("п.4 chi2_ref не воспроизводится: %r != %r" % (m1["chi2_ref"], m2["chi2_ref"]))
    if not (0.0 <= m1["shape"] <= 1.0):
        _fail("п.4 shape вне [0, 1]: получено %r" % m1["shape"])


def _check5_closure(S, args):
    """5. Closure на синтетике (n_rep = 30), все шесть критериев. Гейты: C — по спеке строго
    (|bias| < 3*spread/sqrt(30), coverage_68 в [0,45; 0,90]); A1 — bias с допуском Неймана (D1),
    а вместо coverage вокруг истины — калибровка sd вокруг среднего (D4: coverage_68 A1 измерена
    0,43/0,27 из-за смещения Неймана, критерий спеки для A1 невыполним; число печатается)."""
    fits = {"A0": fit_A0, "A1": fit_A1, "A2": fit_A2, "B": fit_B, "C": fit_C, "D": make_fit_D(S["e"]),
            "E1": make_fit_E1(10), "E2": fit_E2}
    cl = closure(*args, S["a_true"], n_rep=30, template_noise=False, fits=fits)
    ney_pct = 100.0 * _neyman_bound(S) / S["a_true"]
    for nm, key, allow in (("C", "coverage_68", 0.0 * ney_pct), ("A1", "sd_calib_68", ney_pct)):
        for j, row in enumerate(cl[nm]):
            lim = 3.0 * row["spread_pct"] / np.sqrt(30) + allow[j]
            if not abs(row["bias_pct"]) < lim:
                _fail("п.5 %s нуклид %d: ожидалось |bias_pct| < 3*spread/sqrt(30)%s = %.3f, получено %.3f"
                      % (nm, j, "" if nm == "C" else " + допуск Неймана", lim, row["bias_pct"]))
            if not (0.45 <= row[key] <= 0.90):
                _fail("п.5 %s нуклид %d: ожидалось %s в [0,45; 0,90], получено %.3f" % (nm, j, key, row[key]))
    return cl


def _check7_shape(S, args, res, cl):
    """7. E1, E2: shape по единой метрике не больше shape(A1); a в 3*spread (closure) от a_true;
    детерминированно: Σm = S (нормировка), у E1 TV ≡ shape; у E2 — решение KKT точной задачи."""
    C, p, b = _prep(args[0], args[1], args[2], args[5])
    y, S_ = np.maximum(p - b, 0.0), float(np.maximum(p - b, 0.0).sum())
    ma = lambda a: metrics_all(a, args[0], args[1], args[2], S["k"], S["sel"], S["e"], ["n0", "n1"], {"n0": 1.0, "n1": 1.0}, 1.0)
    sh1, out = ma(res["A1"][0])["shape"], {}
    for nm, f in (("E1", fit_E1), ("E2", fit_E2)):
        a, sd, ex = f(*args)
        sh = ma(a)["shape"]
        out[nm] = (a, sd, ex, sh)
        if not sh <= sh1 + 1e-9:
            _fail("п.7 %s: ожидалось shape <= shape(A1) = %.6f, получено %.6f" % (nm, sh1, sh))
        for j in range(2):
            spread = cl[nm][j]["spread_pct"] / 100.0 * S["a_true"][j]
            if not abs(a[j] - S["a_true"][j]) <= 3 * spread:
                _fail("п.7 %s нуклид %d: ожидалось |a - a_true| <= 3*spread = %.1f, получено %.1f"
                      % (nm, j, 3 * spread, abs(a[j] - S["a_true"][j])))
    _check7_exact(C, y, S_, out)
    return out


def _check7_exact(C, y, S_, out):
    """п.7, детерминированная часть. Мутация «убрано ограничение нормировки в E1» даёт LAD-решение
    с Σm ≠ S и TV ≠ shape — амплитудный признак её не обязан видеть, эти два — видят."""
    a1, _, ex1, sh1 = out["E1"]
    if not abs(ex1["model_sum_over_S"] - 1.0) <= 1e-9:
        _fail("п.7 E1 нормировка: ожидалось Σm/S = 1 ± 1e-9, получено %.12f" % ex1["model_sum_over_S"])
    if not abs(ex1["tv"] - sh1) <= 1e-9:
        _fail("п.7 E1: ожидалось TV == shape (единая метрика) ± 1e-9, получено TV = %.9f, shape = %.9f" % (ex1["tv"], sh1))
    a2, _, ex2, _ = out["E2"]
    if not abs(ex2["model_sum_over_S"] - 1.0) <= 1e-6:
        _fail("п.7 E2 нормировка: ожидалось Σm/S = 1 ± 1e-6, получено %.9f" % ex2["model_sum_over_S"])
    if np.all(a2 > 0):   # точное решение ККТ задачи E2 при всех активных параметрах
        W, c = 1.0 / np.maximum(y, 1.0), C.sum(0)
        M = np.block([[2.0 * C.T @ (C * W[:, None]), c[:, None]], [c[None, :], np.zeros((1, 1))]])
        a_ex = np.linalg.solve(M, np.concatenate([2.0 * C.T @ (W * y), [S_]]))[:len(c)]
        rel = float(np.max(np.abs(a2 - a_ex) / np.abs(a_ex)))
        if not rel <= 1e-6:
            _fail("п.7 E2 против точного ККТ-решения: ожидалось <= 1e-6 отн., получено %.3e (a = %s, точное %s)"
                  % (rel, np.array2string(a2, precision=2), np.array2string(a_ex, precision=2)))


def selftest():
    """Самопроверка без файлов (спека, раздел «Самопроверка»). Провал — SystemExit('SELFTEST FAIL: …')."""
    rng = np.random.default_rng(1)
    S = _synth(3.0, rng)
    args = (S["cols"], S["counts"], S["bg_scaled"], S["k"], S["n_events"], S["sel"])
    _check0_cash_forms(S, args)
    res = _check1_return(S, args)
    _check2_strong_bg(rng)
    acc = _check3_cash_min(S, args, res)
    _check4_metrics(S, args, res)
    cl = _check5_closure(S, args)
    e7 = _check7_shape(S, args, res, cl)
    for nm in ("A1", "A2", "C"):
        print("  п.1 %s: a = %s, sd = %s" % (nm, np.array2string(res[nm][0], precision=1),
                                            np.array2string(res[nm][1], precision=1)))
    print("  п.3б C: f(3 старта) = %s, разброс %.2e отн.; a разброс %.2e отн.; sd числ/аналит расх. %.2e"
          % (np.array2string(acc["f_starts"], precision=9), acc["f_rel_spread"], acc["a_rel_spread"], acc["sd_rel_diff"]))
    for nm in ("A1", "C"):
        print("  п.5 %s: bias_pct = %s, spread_pct = %s, coverage_68 = %s, sd_calib_68 = %s"
              % (nm, [round(r_["bias_pct"], 3) for r_ in cl[nm]], [round(r_["spread_pct"], 3) for r_ in cl[nm]],
                 [r_["coverage_68"] for r_ in cl[nm]], [r_["sd_calib_68"] for r_ in cl[nm]]))
    for nm in ("E1", "E2"):
        a, sd, ex, sh = e7[nm]
        print("  п.7 %s: a = %s, sd = %s, shape = %.6f, Σm/S = %.9f, %s; closure bias_pct = %s, coverage_68 = %s"
              % (nm, np.array2string(a, precision=1), np.array2string(sd, precision=1), sh, ex["model_sum_over_S"],
                 "TV = %.6f, LP %.3f с" % (ex["tv"], ex["lp_time_s"]) if nm == "E1" else "chi2_shape = %.3f" % ex["chi2_shape"],
                 [round(r_["bias_pct"], 3) for r_ in cl[nm]], [r_["coverage_68"] for r_ in cl[nm]]))


def main():
    ap = argparse.ArgumentParser(description="crit_bench — стенд сравнения критериев подгонки")
    ap.add_argument("--selftest", action="store_true", help="самопроверка на синтетике, без файлов")
    ns = ap.parse_args()
    if not ns.selftest:
        ap.print_help()
        return 2
    selftest()
    print("SELFTEST OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
