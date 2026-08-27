# -*- coding: utf-8 -*-
"""Экспорт данных для плакатной веб-страницы «Из чего состоит фон RadiaCode-103».

Считает ОБА разложения на действующей модели (реальная комната, кирпич+бетон
вместе) и выгружает в JSON всё, что рисует страница: измеренный спектр,
суммарную модель, вклад каждого нуклида, амплитуды с ошибками, метрики,
полосы и матрицу корреляций.

ПРИВАТНОСТЬ: габариты комнаты и положение прибора в JSON НЕ ПОПАДАЮТ —
страница публичная, а это личные данные оператора (см. .gitignore репозитория).
Из шапок шаблонов берутся только n_events / t_run_s / flux, размеров там нет.
"""
import os, sys, json
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
from scipy.optimize import nnls

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault("G4MODELS_TEMPLATE_FMT", "rc103_field_room_%s.csv")
import fit_two_criteria as ftc

TEMPLATE_FMT = "rc103_field_room_%s.csv"
TL208_BRANCH = 0.3594
NUCS = ftc.NUCS
# Оптимум связанной модели — из прогона fit_physical_chains.py (см. out/).
F_RN, R_TH, F_TN = 0.000, 2.000, 0.000

RU = {"K40": "K-40", "Ra226": "Ra-226", "Pb214": "Pb-214", "Bi214": "Bi-214",
      "Pb212": "Pb-212", "Ac228": "Ac-228", "Bi212": "Bi-212", "Tl208": "Tl-208",
      "mu": "мюоны"}


def main():
    smp = ftc.read_rcxml.read(os.path.join(ftc.MEAS_DIR, ftc.MEAS_NAME))[0]
    cnt = smp.counts[:-1].astype(float)
    ch = np.arange(len(cnt))
    e_meas = np.asarray(sum(c * ch ** i for i, c in enumerate(ftc.CAL_ROOM)))
    live = float(smp.live)

    cps, var = {}, {}
    for nuc in NUCS:
        meta, arr, cnt_mc = ftc.read_template(
            os.path.join(ftc.TEMPLATE_DIR, TEMPLATE_FMT % nuc))
        cps[nuc] = ftc.rcspec.fold(arr, "103")
        var[nuc] = ftc.template_variance(cnt_mc, float(meta.get("t_run_s", 0.0)))
    meta_mu, arr_mu, cnt_mu = ftc.read_template(ftc.MUON_CSV)
    cps["mu"] = ftc.rcspec.fold(arr_mu, "103")
    var["mu"] = ftc.template_variance(cnt_mu, float(meta_mu.get("n_events", 0.0)))

    cols = NUCS + ["mu"]
    A = np.zeros((len(e_meas), len(cols)))
    V = np.zeros_like(A)
    for k, name in enumerate(cols):
        A[:, k] = ftc.fl.rebin_model_to_meas(
            np.arange(len(cps[name])) + 0.5, cps[name], e_meas)
        V[:, k] = ftc.fl.rebin_model_to_meas(
            np.arange(len(var[name])) + 0.5, var[name], e_meas)

    sel = (e_meas >= ftc.E_LO) & (e_meas < ftc.E_HI)
    A, V, y, e = A[sel], V[sel], cnt[sel], e_meas[sel]
    A_counts, VAR_counts = A * live, V * live * live

    # --- понуклидное разложение, оба критерия ---
    out_free = {}
    for crit, w in (("A", 1.0 / np.sqrt(np.maximum(y, 1.0))),
                    ("B", 1.0 / np.maximum(y, 1.0))):
        amp, sd, pred, chi2ndf, shape = ftc.fit(
            A_counts, y, w, cols, f"экспорт-{crit}", "", var_counts=VAR_counts)
        out_free[crit] = {
            "amp": {n: float(a) for n, a in zip(cols, amp)},
            "sd": {n: (float(s) if np.isfinite(s) else None) for n, s in zip(cols, sd)},
            "chi2ndf": float(chi2ndf), "shape": float(shape),
        }

    # --- связанная модель на оптимуме ---
    col_ra = A_counts[:, cols.index("Ra226")] + (1 - F_RN) * (
        A_counts[:, cols.index("Pb214")] + A_counts[:, cols.index("Bi214")])
    col_th = A_counts[:, cols.index("Ac228")] + R_TH * (1 - F_TN) * (
        A_counts[:, cols.index("Pb212")] + A_counts[:, cols.index("Bi212")]
        + TL208_BRANCH * A_counts[:, cols.index("Tl208")])
    Ac = np.column_stack([A_counts[:, cols.index("K40")], col_ra, col_th,
                          A_counts[:, cols.index("mu")]])
    vc_ra = VAR_counts[:, cols.index("Ra226")] + (1 - F_RN) ** 2 * (
        VAR_counts[:, cols.index("Pb214")] + VAR_counts[:, cols.index("Bi214")])
    vc_th = VAR_counts[:, cols.index("Ac228")] + (R_TH * (1 - F_TN)) ** 2 * (
        VAR_counts[:, cols.index("Pb212")] + VAR_counts[:, cols.index("Bi212")]
        + TL208_BRANCH ** 2 * VAR_counts[:, cols.index("Tl208")])
    Vc = np.column_stack([VAR_counts[:, cols.index("K40")], vc_ra, vc_th,
                          VAR_counts[:, cols.index("mu")]])
    wA = 1.0 / np.sqrt(np.maximum(y, 1.0))
    amp_c, sd_c, pred_c, chi2_c, shape_c = ftc.fit(
        Ac, y, wA, ["K40", "Ra226_chain", "Th232_chain", "mu"],
        "экспорт-связанный", "", var_counts=Vc)

    a_k, a_ra, a_th, a_mu = (float(v) for v in amp_c)
    act = {"K40": a_k, "Ra226": a_ra, "Pb214": a_ra * (1 - F_RN),
           "Bi214": a_ra * (1 - F_RN), "Ac228": a_th,
           "Pb212": a_th * R_TH * (1 - F_TN), "Bi212": a_th * R_TH * (1 - F_TN),
           "Tl208": TL208_BRANCH * a_th * R_TH * (1 - F_TN), "mu": a_mu}

    # Вклад каждого нуклида = его активность x его шаблон. Сумма тождественно
    # равна pred_c: столбцы связанной модели — те же шаблоны с теми же весами.
    parts = {n: act[n] * A_counts[:, cols.index(n)] / live for n in cols}
    resid = float(np.max(np.abs(sum(parts.values()) * live - pred_c)))

    step = 2                       # прореживание для веса страницы
    idx = np.arange(0, len(e), step)
    ser = lambda v: [round(float(x), 6) for x in np.asarray(v)[idx]]

    bands = []
    for lo, hi in ftc.BANDS:
        m = (e >= lo) & (e < hi)
        bands.append({"lo": lo, "hi": hi,
                      "meas": round(float(y[m].sum() / live), 4),
                      "model": round(float(pred_c[m].sum() / live), 4)})

    # Вырожденность считаем ТЕМ ЖЕ способом, что degeneracy_report в
    # fit_two_criteria.py:164-166: нормировка столбцов на евклидову норму, мера
    # близости norm.T @ norm (косинус угла, без центрирования), cond — отношение
    # крайних сингулярных чисел. Иначе страница показывает не то число, что
    # печатают отчёты: np.corrcoef центрирует, а деление на сумму столбца даёт
    # другую матрицу — на этом уже разошлись 265 на странице против 246 в отчёте.
    norm_cols = A_counts / np.maximum(np.linalg.norm(A_counts, axis=0), 1e-300)
    C = norm_cols.T @ norm_cols
    sv = np.linalg.svd(norm_cols, compute_uv=False)
    cond_val = float(sv[0] / sv[-1])

    data = {
        "meas": {"file": "фон помещения, 7 суток", "live_s": round(live, 1),
                 "live_h": round(live / 3600.0, 2), "total_counts": int(cnt.sum()),
                 "cps": round(float(y.sum() / live), 3),
                 "e_lo": ftc.E_LO, "e_hi": ftc.E_HI, "n_ch": int(sel.sum())},
        "energy": ser(e),
        "measured": ser(y / live),
        "model": ser(pred_c / live),
        "parts": {n: ser(parts[n]) for n in cols},
        "chain": {
            "amp": {"K40": a_k, "Ra226_chain": a_ra, "Th232_chain": a_th, "mu": a_mu},
            "sd": {k: float(v) for k, v in zip(
                ["K40", "Ra226_chain", "Th232_chain", "mu"], sd_c)},
            "act": act, "f_rn": F_RN, "r_th": R_TH, "f_tn": F_TN,
            "chi2ndf": float(chi2_c), "shape": float(shape_c),
            "balance": round(float(pred_c.sum() / y.sum()), 4),
            "closure_counts": resid,
        },
        "free": out_free,
        "bands": bands,
        "cond": cond_val,
        "corr": [[round(float(v), 3) for v in row] for row in C],
        "corr_names": [RU[n] for n in cols],
        "muon_pdg": ftc.MUON_PDG_PER_S,
        "ru": RU,
    }

    out = os.path.join(HERE, "out", "poster_data.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, separators=(",", ":"))
    print("записано:", out, os.path.getsize(out), "байт")
    print("замыкание суммы вкладов на модель, макс |разность| =", resid, "отсчётов")
    print("баланс модель/измерение =", data["chain"]["balance"])


if __name__ == "__main__":
    main()
