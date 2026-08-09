# -*- coding: utf-8 -*-
"""Быстрый проход МЕТОДА 2 для Ra-226 (Пилот 1, задача #175/#182).

НЕ полноценный export_data.py-эквивалент -- метод 1 (МК-шаблоны по
нуклидам) здесь не строится вовсе (iso_*.csv звеньев Ra-226 ещё не
прогнаны), и main() в export_data.py требует их безусловно (там всё
одной линейной функцией, run_method1/run_method2 -- вложенные замыкания,
не переиспользуемые отдельно -- см. TODO ниже). Этот скрипт переиспользует
ТОЛЬКО настоящие модульные функции (make_full_response, make_eps_peak_
interp, fit_fwhm_calibration, fit_amplitudes, broaden_and_rebin) и пишет
свою тонкую склейку для измеренного .spe (не BecqMoni XML) и одноамплитудной
подгонки библиотеки+сумм-пиков.

TODO (для полноценного метода 1 + чистой архитектуры): вынести run_method1/
run_method2 из main() в модульные функции, как остальной движок -- тогда
подобный скрипт не нужен, echo export_data.py --method2-only.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault("G4MODELS_SOURCE_CONFIG",
                      os.path.join(HERE, "configs", "ra226.yaml"))
os.environ.setdefault("PYTHONIOENCODING", "utf-8")

sys.path.insert(0, HERE)
import export_data as ed  # noqa: E402

SPECTRAVIBE_ROOT = (r"C:\Users\Дмитрий\Мой диск\Дозиметрия\ИИ\1 Скилы"
                    r"\0_Work\gamma-spectrum-analysis")
sys.path.insert(0, os.path.join(SPECTRAVIBE_ROOT, "scripts"))
from gamma.io.lsrm_spe import read_lsrm_spe  # noqa: E402


def read_pair():
    cfg = ed._CFG
    sample_rel = cfg["source"]["measured_sample_spe_rel"]
    sample_path = os.path.join(str(ed.KIT), *sample_rel.split("/"))
    bg_env = cfg["source"]["measured_background_spe_env"]
    bg_path = os.environ.get(bg_env)
    if not bg_path:
        raise SystemExit(
            "Задайте %s -- путь к фоновому .spe (личные данные, "
            "в репозитории их нет)." % bg_env)
    if not os.path.isfile(sample_path):
        raise SystemExit("Нет образца: " + sample_path)
    if not os.path.isfile(bg_path):
        raise SystemExit("Нет фона: " + bg_path)

    s = read_lsrm_spe(sample_path)
    b = read_lsrm_spe(bg_path)
    ch_s = np.arange(s.n_channels, dtype=float)
    ch_b = np.arange(b.n_channels, dtype=float)
    return ({"counts": np.asarray(s.counts, dtype=np.int64),
             "e_of_ch": np.asarray(s.channel_to_energy(ch_s), dtype=float),
             "live_s": float(s.live_time), "real_s": float(s.real_time)},
            {"counts": np.asarray(b.counts, dtype=np.int64),
             "e_of_ch": np.asarray(b.channel_to_energy(ch_b), dtype=float),
             "live_s": float(b.live_time), "real_s": float(b.real_time)})


def main():
    meas, bg = read_pair()

    bg_on_meas = np.interp(meas["e_of_ch"], bg["e_of_ch"],
                           bg["counts"].astype(float), left=0.0, right=0.0)
    bg_scale_time = meas["live_s"] / bg["live_s"]
    bg_scaled = bg_on_meas * bg_scale_time

    e = meas["e_of_ch"]
    ch_edges = np.concatenate((
        [e[0] - 0.5 * (e[1] - e[0])],
        0.5 * (e[:-1] + e[1:]),
        [e[-1] + 0.5 * (e[-1] - e[-2])],
    ))

    fwhm_cal = ed.fit_fwhm_calibration(meas["counts"], e)
    ed.FWHM_LAW.update({"kind": "power", "k": fwhm_cal["k"], "p": fwhm_cal["p"]})
    print("ПШПВ калибровка: k=%.4f p=%.4f (СКО %.1f%%)"
          % (fwhm_cal["k"], fwhm_cal["p"], fwhm_cal["rms_dev_pct"]))

    eps_peak = ed.make_eps_peak_interp(os.path.join(ed.BUILD, "grid"))
    resp = ed.make_full_response(os.path.join(ed.BUILD, "grid"), ch_edges,
                                 True, eps_peak)

    T = meas["live_s"]
    sel = (e >= ed.E_FIT_LO) & (e <= ed.E_FIT_HI)
    y_sel = meas["counts"][sel].astype(float)
    bgm = bg_scaled[sel]

    # Та же схема, что run_method2 в export_data.py: одна СУММАРНАЯ форма
    # shape_total (веса всех линий + сумм-пиков внутри одного распада
    # ветви, депопуляция F_B вычтена из одиночных линий) -- одна
    # амплитуда на весь спектр, не по столбцу на линию (см. SKILL.md
    # geant4-spectrum-pipeline, раздел TCS/F_B).
    shape_total = np.zeros_like(e)

    # ИСПРАВЛЕНО 09.08.2026 (аудит Б2, коммит df5d178, та же находка, что и
    # в export_data.py.run_method2 -- третья независимая копия формулы,
    # тоже не была синхронизирована): эффективность ПАРТНЁРА каскада в
    # депопуляции -- ПОЛНАЯ (eps_total=shape.sum()), не пиковая.
    depl = {}
    for E1, E2, nuc_key, I1_pct, I2_pct, note, fb_pct in ed.SUM_PEAKS:
        shp1, _, _ = resp(E1)
        shp2, _, _ = resp(E2)
        eps1s_tot = float(shp1.sum())
        eps2s_tot = float(shp2.sum())
        fb_frac_s = fb_pct / 100.0
        k1 = (nuc_key, round(E1, 3))
        k2 = (nuc_key, round(E2, 3))
        depl[k1] = depl.get(k1, 0.0) + (
            (I1_pct / 100.0) * (I2_pct / 100.0) * eps2s_tot / fb_frac_s)
        depl[k2] = depl.get(k2, 0.0) + (
            (I2_pct / 100.0) * (I1_pct / 100.0) * eps1s_tot / fb_frac_s)

    n_lines_used = 0
    for E, I_pct, nuc_key, note in ed.GAMMA_LIBRARY:
        shape, chans, eps = resp(E)
        w = I_pct / 100.0
        w_depl = depl.get((nuc_key, round(E, 3)), 0.0)
        if w_depl > 0:
            w = max(0.0, w - w_depl)
        shape_total += w * shape
        n_lines_used += 1

    n_sum_used = 0
    for E1, E2, nuc_key, I1_pct, I2_pct, note, fb_pct in ed.SUM_PEAKS:
        Esum = E1 + E2
        if Esum > ed.E_FIT_HI:
            continue
        _, _, eps1 = resp(E1)
        _, _, eps2 = resp(E2)
        shape_sum, _, eps_sum_node = resp(Esum)
        w = ((I1_pct / 100.0) * (I2_pct / 100.0) * eps1 * eps2
             / max(eps_sum_node, 1e-30) / (fb_pct / 100.0))
        shape_total += w * shape_sum
        n_sum_used += 1

    coef, dcoef, chi2b, ndofb, model2 = ed.fit_amplitudes(
        y_sel, [shape_total[sel] * T, bgm], ed.SYS_FLOOR)
    A_Bq, dA_Bq = float(coef[0]), float(dcoef[0])
    bg_amp = float(coef[1])

    p = ed._CFG["passport"]
    decay_f = ed.decay_factor_years(p["half_life_years"], p["days_pass_to_meas"])
    A_pass = p["bq_per_kg"] * (p["mass_g"] / 1000.0) * decay_f
    dA_pass = A_pass * p["unc_pct"] / 100.0

    print()
    print("живое время образца: %.1f с, фона: %.1f с (масштаб %.4f)"
          % (meas["live_s"], bg["live_s"], bg_scale_time))
    print("диапазон подгонки: %.0f..%.0f кэВ, %d линий + %d сумм-пиков"
          % (ed.E_FIT_LO, ed.E_FIT_HI, n_lines_used, n_sum_used))
    print("chi2/ndof = %.2f/%d = %.2f  амплитуда фона (ожидается ~1) = %.3f"
          % (chi2b, ndofb, chi2b / ndofb, bg_amp))
    print()
    print("паспорт (распад-корр. на дату измерения): A = %.0f +- %.0f Бк"
          % (A_pass, dA_pass))
    print("метод 2 (быстрый проход)               : A = %.0f +- %.0f Бк"
          % (A_Bq, dA_Bq))
    print("отношение метод2/паспорт: %.3f" % (A_Bq / A_pass))


if __name__ == "__main__":
    main()
