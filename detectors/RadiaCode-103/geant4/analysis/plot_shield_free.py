# -*- coding: utf-8 -*-
"""
plot_shield_free.py

Свободная подгонка шаблонов (K40, Ra226_chain, Th232_chain, mu) на спектр
в домике БЕЗ ОПОРЫ на расчётную активность открытого фона. Повторной
реализации чтения/свёртки/подгонки нет — через predict_shield (ps) и
fit_two_criteria (ftc), как в plot_shield_predict.py.
"""
import os, sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
import fit_two_criteria as ftc
import predict_shield as ps

matplotlib.rcParams["font.size"] = 9
matplotlib.rcParams["figure.dpi"] = 130
matplotlib.rcParams["axes.grid"] = True
matplotlib.rcParams["grid.alpha"] = 0.3

OUT_DIR = os.path.join(HERE, "..", "verify")
os.makedirs(OUT_DIR, exist_ok=True)
OUT_PNG = os.path.join(OUT_DIR, "RC103_shield_free_fit.png")

C_MEAS = "#111111"
C_MODEL = "#d62728"
C_MUON = "#8c8c8c"
C_K40 = "#59a14f"
C_RA = "#4e79a7"
C_TH = "#e15759"


def build():
    cnt_s, e_s, live_s, _ = ps.read_meas(ps.MEAS_SHIELD, None)
    cps_s, var_s, hm_s, posture_s = ps.load_cps(ps.FMT_SHIELD, ps.MUON_SHIELD_CSV)
    if not hm_s:
        sys.exit("нет мюонного шаблона с домиком")

    names_s, A_s, V_s = ps.columns_on_grid(cps_s, var_s, e_s, hm_s)
    sel_s = (e_s >= ftc.E_LO) & (e_s < ftc.E_HI)
    cnt_s = cnt_s[sel_s]
    e_s = e_s[sel_s]
    A_s = A_s[sel_s, :]
    V_s = V_s[sel_s, :]

    y_s = cnt_s
    w_s = 1.0 / np.sqrt(np.maximum(y_s, 1.0))

    title = "Свободная подгонка в домике (без расчётной активности)"
    note = "критерий A (chi2/ndf), веса пуассоновские; активности открытого фона НЕ используются"
    amp, sd, pred, chi2ndf, shape = ftc.fit(
        A_s * live_s, y_s, w_s, names_s, title, note, var_counts=V_s * live_s * live_s)

    i_mu = names_s.index("mu")
    i_k = names_s.index("K40")
    i_ra = names_s.index("Ra226_chain")
    i_th = names_s.index("Th232_chain")

    pred_mu = A_s[:, i_mu] * amp[i_mu]
    pred_k = A_s[:, i_k] * amp[i_k]
    pred_ra = A_s[:, i_ra] * amp[i_ra]
    pred_th = A_s[:, i_th] * amp[i_th]
    pred_tot = pred_mu + pred_k + pred_ra + pred_th
    meas_rate = cnt_s / live_s

    # pred_* уже в имп/с: A_s — шаблоны cps, amp получена подгонкой
    # (A_s*live_s) к counts, поэтому A_s@amp = counts/live_s = rate.
    # Повторное деление на live_s было бы двойным (W-класс единиц).
    return {
        "e_s": e_s, "meas": meas_rate,
        "pred_mu": pred_mu, "pred_k": pred_k,
        "pred_ra": pred_ra, "pred_th": pred_th,
        "pred_tot": pred_tot,
        "amp": dict(zip(names_s, amp)), "sd": dict(zip(names_s, sd)),
        "chi2ndf": chi2ndf, "shape": shape, "posture": posture_s,
    }


def draw(d):
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(11, 8), gridspec_kw={"height_ratios": [3, 1]}, sharex=True)

    ax_top.fill_between(d["e_s"], 0, d["pred_mu"], color=C_MUON, linewidth=0, label="подгонка: мюоны")
    base = d["pred_mu"]
    ax_top.fill_between(d["e_s"], base, base + d["pred_k"], color=C_K40, linewidth=0, label="подгонка: K-40")
    base = base + d["pred_k"]
    ax_top.fill_between(d["e_s"], base, base + d["pred_ra"], color=C_RA, linewidth=0, label="подгонка: Ra226_chain")
    base = base + d["pred_ra"]
    ax_top.fill_between(d["e_s"], base, base + d["pred_th"], color=C_TH, linewidth=0, label="подгонка: Th232_chain")

    ax_top.plot(d["e_s"], d["pred_tot"], color=C_MODEL, lw=1.3, label="сумма подгонки")
    ax_top.step(d["e_s"], d["meas"], where="mid", color=C_MEAS, lw=1.0, label="измерено, домик")

    ax_top.set_yscale("log")
    y_min = max(1e-7, np.min(d["meas"][d["meas"] > 0]))
    y_max = 1.6 * np.max(d["meas"])
    ax_top.set_ylim(y_min, y_max)
    ax_top.set_xlim(ftc.E_LO, ftc.E_HI)
    ax_top.set_ylabel("скорость счёта, имп/с на канал")
    posture = d["posture"]
    p_str = (f"опора {posture[0]} мм, экран {'вверх' if posture[1] == 1.0 else 'вниз'}"
             if posture is not None else "посадка не указана в шаблонах")
    ax_top.set_title(
        "RadiaCode-103: свободная подгонка в домике, БЕЗ расчётной активности "
        f"открытого фона ({p_str})\nchi2/ndf = {d['chi2ndf']:.2f}, невязка формы = {d['shape']:.4f}")
    ax_top.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=3, frameon=False)

    ratio = np.where(d["meas"] > 0, d["pred_tot"] / d["meas"], np.nan)
    ax_bot.step(d["e_s"], ratio, where="mid", color=C_MODEL, lw=1.0)
    ax_bot.axhline(1.0, color="#999999", lw=0.8, linestyle="--")
    ax_bot.set_ylim(0, 2.5)
    ax_bot.set_ylabel("модель/измерение")
    ax_bot.set_xlabel("энергия, кэВ")

    fig.tight_layout()
    fig.savefig(OUT_PNG)
    print(f"Сохранено: {OUT_PNG}")


if __name__ == "__main__":
    dd = build()
    print(f"\nПостановка (P-005): {dd['posture']}")
    print(f"chi2/ndf = {dd['chi2ndf']:.3f}, невязка формы = {dd['shape']:.4f}")
    print("Компонента     | амплитуда   | sigma")
    for n in ["K40", "Ra226_chain", "Th232_chain", "mu"]:
        print(f"{n:<15}| {dd['amp'][n]:11.4f} | {dd['sd'][n]:.4f}")
    draw(dd)
