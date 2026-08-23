# -*- coding: utf-8 -*-
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import numpy as np
from scipy.optimize import linprog
import fit_nuclides as fn
import fit_lines as fl
import read_rcxml
import io
import contextlib

def assemble():
    # Таблица линий архивирована 21.08 (results/_attic_table_method_20260821/),
    # fn.BG_PREFIX/WF_PREFIX по умолчанию уже указывают на ионный режим -
    # параметризация по префиксам больше не нужна, единственный активный путь.
    names, cols = fn.merge_by_chain(*fn.load_templates())
    mu, pdg = fn.load_muons()
    if mu is not None:
        names.append("mu")
        cols.append(mu)
    
    smp = read_rcxml.read(fn.MEASURED)[0]
    cnt = smp.counts[:-1].astype(float)
    ch = np.arange(len(cnt))
    e_meas = np.asarray(sum(c * ch ** i for i, c in enumerate(fn.CAL_ROOM)))[:len(cnt)]
    cps_meas = cnt / smp.live
    
    A = np.zeros((len(e_meas), len(cols)))
    for k, c in enumerate(cols):
        A[:, k] = fl.rebin_model_to_meas(np.arange(len(c)) + 0.5, c, e_meas)
    
    return (names, A, e_meas, cps_meas, smp.live)

def fraction_covered(model, meas):
    return float(np.sum(np.minimum(model, meas)) / np.sum(meas))

def form_residual_pct(model, meas):
    # Невязка формы по Морозову: сумма модулей поканальных отклонений к сумме
    # измерения. НЕ зависит от live-time (в отличие от хи2/ndf) — мера ФОРМЫ.
    return float(100.0 * np.sum(np.abs(model - meas)) / np.sum(meas))

def chi2_of(model, meas, live, e_meas):
    # Маска считалась, но не применялась - хи2 суммировался по ВСЕМУ спектру,
    # включая пороговый шум каналов 0-9 (10^5-10^6 отсчётов, README:139-144).
    # Без обрезки эти каналы забивают статистику целиком.
    mask = (e_meas >= 20) & (e_meas < 2830)
    ndf = mask.sum() - 1
    m, y = model[mask], meas[mask]
    var = np.maximum(y * live, 1) / live ** 2
    chi2 = np.sum((m - y) ** 2 / var)
    return (chi2, ndf)

def fit_by_coverage(names, A, e_meas, cps_meas):
    mask = (e_meas >= 20) & (e_meas < 2830)
    M = A[mask]
    y = cps_meas[mask]
    
    n = M.shape[1]
    m = M.shape[0]
    
    # Линейное программирование для минимизации суммы отклонений L1
    A_ub = np.block([[M, -np.eye(m)], [-M, -np.eye(m)]])
    b_ub = np.concatenate([y, -y])
    c = np.concatenate([np.zeros(n), np.ones(m)])
    bounds = [(0, None)] * (n + m)
    
    try:
        res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
    except Exception as e:
        print("Ошибка при решении линейной задачи: %s" % str(e))
        return (None, False)
    
    if not res.success:
        print("Ошибка в решении LP: %s" % res.message)
        return (None, False)
    
    return (res.x[:n], True)

def main():
    print("=== ПОЛНЫЙ РАСПАД (ион), поиск ошибок разложения ===")
    
    names, A, e_meas, cps_meas, live = assemble()
    
    # Выполняем существующий хи2-фит
    with contextlib.redirect_stdout(io.StringIO()):
        amps_dict, a_mu = fn.main()
    
    amp_chi2 = np.array([amps_dict.get(n, 0.0) for n in names])
    model_chi2 = A @ amp_chi2
    
    # Выполняем L1-фит
    amp_l1, ok = fit_by_coverage(names, A, e_meas, cps_meas)
    
    if not ok:
        print("L1-фит не удался для этой конфигурации")
        return 0
    
    model_l1 = A @ amp_l1
    
    # Выводим амплитуды
    print("%-10s %14s %14s" % ("звено", "хи2-фит (линии)", "L1-фит (весь спектр)"))
    for i, name in enumerate(names):
        print("%-10s %14.3f %14.3f" % (name, amp_chi2[i], amp_l1[i]))
    
    # Пересчет метрик
    mask = (e_meas >= 20) & (e_meas < 2830)
    fc_chi2 = fraction_covered(model_chi2[mask], cps_meas[mask])
    fc_l1 = fraction_covered(model_l1[mask], cps_meas[mask])
    
    chi2_val, ndf = chi2_of(model_chi2, cps_meas, live, e_meas)
    chi2_l1, _ = chi2_of(model_l1, cps_meas, live, e_meas)
    
    fr_chi2 = form_residual_pct(model_chi2[mask], cps_meas[mask])
    fr_l1 = form_residual_pct(model_l1[mask], cps_meas[mask])

    print("\nМетрика                         хи2-фит      L1-фит")
    print("%% заполнения формой (весь спектр)  %.1f%%        %.1f%%" % (fc_chi2 * 100, fc_l1 * 100))
    print("невязка формы (весь спектр)        %.1f%%         %.1f%%" % (fr_chi2, fr_l1))
    print("хи2/ndf (по всему спектру)         %.2f         %.2f" % (chi2_val / ndf, chi2_l1 / ndf))
    
    # Невязка по окнам диагностических линий (#FIT-1, 22.08): та же оконная
    # конвенция nsig=2.5, что у fit_lines.line_net_area. Модуль сгенерирован
    # по спеке _spec_coverage_windows.md (IRON MODE ступень 2).
    import rcspec
    import coverage_windows
    print("")
    coverage_windows.print_window_residuals(
        e_meas, cps_meas, [("chi2fit", model_chi2), ("L1fit", model_l1)],
        lambda e0: rcspec.fwhm(e0, fl.MODEL))

    # Полоса-за-полосой сравнение
    bands = ((20, 100), (100, 300), (300, 700), (700, 1500), (1500, 2000), (2000, 2400), (2400, 2830))
    print("\nполоса,кэВ   изм     хи2-модель L1-модель  хи2 м/и   L1 м/и")
    
    for b_low, b_high in bands:
        band_mask = (e_meas >= b_low) & (e_meas < b_high)
        if not np.any(band_mask):
            continue
        meas_band = cps_meas[band_mask]
        model_chi2_band = model_chi2[band_mask]
        model_l1_band = model_l1[band_mask]
        
        ratio_chi2 = np.sum(model_chi2_band) / np.sum(meas_band) if np.sum(meas_band) > 0 else 0
        ratio_l1 = np.sum(model_l1_band) / np.sum(meas_band) if np.sum(meas_band) > 0 else 0
        
        print("%6d-%-6d %7.1f %10.2f %10.2f %8.3f %8.3f" % (
            b_low, b_high,
            np.sum(meas_band),
            np.sum(model_chi2_band),
            np.sum(model_l1_band),
            ratio_chi2,
            ratio_l1
        ))
    
    # Поиск расхождений по каналам. ВАЖНО: diff/rel в УРЕЗАННОМ маской
    # пространстве, индекс из argsort - позиция ВНУТРИ mask, не в полном
    # массиве. Сгенерированная версия применяла такой индекс к полноразмерным
    # маскам полос и к model_chi2/model_l1 напрямую - рассогласование,
    # топ давал не те энергии. Переводим индекс в энергию через e_meas[mask].
    print("\n=== ПОКАНАЛЬНЫЙ ВКЛАД: где расходятся хи2-фит и L1-фит ===")
    e_m = e_meas[mask]
    mc_m, ml_m, y_m = model_chi2[mask], model_l1[mask], cps_meas[mask]
    rel = (ml_m - mc_m) / np.maximum(y_m, 1e-12)
    abs_rel = np.abs(rel)
    sorted_indices = np.argsort(abs_rel)[::-1]

    print("полоса,кэВ   хи2-модель L1-модель  относит. расхождение")
    shown = 0
    for idx in sorted_indices:
        if abs_rel[idx] <= 0:
            continue
        e0 = e_m[idx]
        band = next((b for b in bands if b[0] <= e0 < b[1]), None)
        if band is None:
            continue
        print("%6d-%-6d %10.4f %10.4f %15.3f"
              % (band[0], band[1], mc_m[idx], ml_m[idx], rel[idx]))
        shown += 1
        if shown >= 10:
            break
    
    # Изменение амплитуд по нуклидам
    print("\nИзменение амплитуд по нуклидам (в процентах):")
    amp_changes = []
    for i in range(len(names)):
        if amp_chi2[i] > 0:
            change = 100 * (amp_l1[i] - amp_chi2[i]) / amp_chi2[i]
            amp_changes.append((names[i], change))
    
    amp_changes.sort(key=lambda x: abs(x[1]), reverse=True)
    
    for name, change in amp_changes[:10]:
        print("%-10s %8.1f%%" % (name, change))
    
    print("\nПримечание: этот анализ сравнивает два КРИТЕРИЯ ФИТА на ОДНИХ И ТЕХ ЖЕ шаблонах.")
    print("Он не говорит, какой из них 'правильнее' — большие различия между двумя фитами")
    print("свидетельствуют о том, что амплитуды плохо определены (дегенерация базиса).")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
