# -*- coding: utf-8 -*-
import math, os, sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import io
import contextlib
HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
import fit_two_criteria as ftc

TEMPLATE_FMT = "rc103_field_room_%s.csv"   # действующая модель: реальная комната, кирпич+бетон вместе
TL208_BRANCH = 0.3594
NUCS = ["K40", "Ra226", "Pb214", "Bi214", "Pb212", "Ac228", "Bi212", "Tl208"]

# Режим задаётся окружением, чтобы не плодить копии скрипта (умолчания
# воспроизводят прежнее поведение):
#   G4M_NO_MUON=1     — убрать мюонный шаблон из базиса;
#   G4M_FIX_RTH=1.0   — зафиксировать r_Th вместо сканирования;
#   G4M_FIX_FTN=0.0   — зафиксировать f_Tn;
#   G4M_FRN_LO/HI/N   — границы и число узлов сетки по f_Rn.
NO_MUON = os.environ.get("G4M_NO_MUON", "") == "1"
_FIX_RTH = os.environ.get("G4M_FIX_RTH")
_FIX_FTN = os.environ.get("G4M_FIX_FTN")

GRID_F_RN = np.linspace(float(os.environ.get("G4M_FRN_LO", 0.0)),
                        float(os.environ.get("G4M_FRN_HI", 0.50)),
                        int(os.environ.get("G4M_FRN_N", 26)))
# Сдвиг равновесия Th-228/Ra-228. Верхняя граница поднята с 2.00 до 4.00
# (2026-08-27): на прежней сетке оптимум ложился РОВНО на край, то есть предел
# ставила сетка, а не данные, и отличить настоящий оптимум от упора в границу
# было невозможно. Шаг сохранён прежний, 0.05.
# Границы сетки вынесены в окружение (2026-08-28): при свободной подгонке без
# мюонов оптимум лёг РОВНО на верхний край 4.00, и отличить настоящий оптимум
# от сползания было нельзя, не двигая границу.
GRID_R_TH = (np.array([float(_FIX_RTH)]) if _FIX_RTH
             else np.linspace(float(os.environ.get("G4M_RTH_LO", 0.50)),
                              float(os.environ.get("G4M_RTH_HI", 4.00)),
                              int(os.environ.get("G4M_RTH_N", 71))))
GRID_F_TN = (np.array([float(_FIX_FTN)]) if _FIX_FTN
             else np.linspace(0.0, 0.20, 11))   # эманирование торона; физически мало

def load_all():
    cps = {}
    var = {}
    for nuc in NUCS:
        path = os.path.join(ftc.TEMPLATE_DIR, TEMPLATE_FMT % nuc)
        try:
            meta, arr, cnt_mc = ftc.read_template(path)
        except FileNotFoundError:
            print(f"Файл шаблона не найден: {path}")
            sys.exit(1)
        cps[nuc] = ftc.rcspec.fold(arr, "103")
        var[nuc] = ftc.template_variance(cnt_mc, float(meta.get("t_run_s", 0.0)))
    if not NO_MUON:
        meta_mu, arr_mu, cnt_mu = ftc.read_template(ftc.MUON_CSV)
        cps["mu"] = ftc.rcspec.fold(arr_mu, "103")
        var["mu"] = ftc.template_variance(cnt_mu, float(meta_mu.get("n_events", 0.0)))
    return cps, var

def build_columns(cps, var, f_rn, r_th, f_tn):
    # Строим столбцы для модели
    A_Ra = cps["Ra226"]
    A_Pb214 = cps["Pb214"]
    A_Bi214 = cps["Bi214"]
    A_Th = cps["Ac228"]
    A_Pb212 = cps["Pb212"]
    A_Bi212 = cps["Bi212"]
    A_Tl208 = cps["Tl208"]

    # Дисперсии
    var_Ra = var["Ra226"]
    var_Pb214 = var["Pb214"]
    var_Bi214 = var["Bi214"]
    var_Th = var["Ac228"]
    var_Pb212 = var["Pb212"]
    var_Bi212 = var["Bi212"]
    var_Tl208 = var["Tl208"]

    # Столбцы
    col_Ra_chain = A_Ra + (1 - f_rn) * (A_Pb214 + A_Bi214)
    var_Ra_chain = var_Ra + (1 - f_rn)**2 * (var_Pb214 + var_Bi214)

    col_Th_chain = A_Th + r_th * (1 - f_tn) * (A_Pb212 + A_Bi212 + TL208_BRANCH * A_Tl208)
    var_Th_chain = var_Th + r_th**2 * (1 - f_tn)**2 * (var_Pb212 + var_Bi212 + TL208_BRANCH**2 * var_Tl208)

    names = ["K40", "Ra226_chain", "Th232_chain"]
    cols = [cps["K40"], col_Ra_chain, col_Th_chain]
    vars_ = [var["K40"], var_Ra_chain, var_Th_chain]
    if not NO_MUON:
        names.append("mu")
        cols.append(cps["mu"])
        vars_.append(var["mu"])

    return names, cols, vars_

def bands_by_component(cps, var, f_rn, r_th, f_tn, amp, e_meas, cnt, live):
    """Вклад КАЖДОГО столбца по полосам (диагностика #MU-1, 2026-08-28).

    bands_report печатает только суммарную модель, поэтому по нему нельзя
    сказать, КТО именно перебирает в полосе. Включается G4M_BANDS_PARTS=1.
    """
    names, cols, _ = build_columns(cps, var, f_rn, r_th, f_tn)
    A = np.zeros((len(e_meas), len(cols)))
    for k, c in enumerate(cols):
        A[:, k] = ftc.fl.rebin_model_to_meas(np.arange(len(c)) + 0.5, c, e_meas) * live
    print()
    print("=== ВКЛАД КОМПОНЕНТ ПО ПОЛОСАМ, с^-1 ===")
    print("полоса,кэВ  | измерено |" + "|".join("%11s" % n for n in names)
          + "|      сумма |   м/и")
    # Полосы можно переопределить: G4M_BANDS="1500:2830:50" — окна по 50 кэВ.
    # Нужно, чтобы увидеть СТРУКТУРУ расхождения внутри широкой полосы:
    # ровный сдвиг означает нехватку масштаба компоненты, локальный выброс —
    # дефект линии или функции отклика.
    spec = os.environ.get("G4M_BANDS")
    if spec:
        _lo, _hi, _w = (float(x) for x in spec.split(":"))
        bands = [(_lo + i * _w, _lo + (i + 1) * _w)
                 for i in range(int(round((_hi - _lo) / _w)))]
    else:
        bands = list(ftc.BANDS)
    for lo, hi in bands:
        s = (e_meas >= lo) & (e_meas < hi)
        m = cnt[s].sum() / live
        parts = [amp[k] * A[s, k].sum() / live for k in range(len(names))]
        tot = sum(parts)
        print("%4d-%4d   | %8.5f |" % (int(lo), int(hi), m)
              + "|".join("%11.5f" % v for v in parts)
              + "| %10.5f | %6.3f" % (tot, tot / m if m else float("nan")))


def main():
    # Чтение измерения
    try:
        smp = ftc.read_rcxml.read(os.path.join(ftc.MEAS_DIR, ftc.MEAS_NAME))[0]
    except FileNotFoundError:
        print(f"Файл измерения не найден: {os.path.join(ftc.MEAS_DIR, ftc.MEAS_NAME)}")
        sys.exit(1)
    cnt = smp.counts[:-1].astype(float)
    ch = np.arange(len(cnt))
    e_meas = np.asarray(sum(c * ch ** i for i, c in enumerate(ftc.CAL_ROOM)))
    live = float(smp.live)

    cps, var = load_all()

    def evaluate(f_rn, r_th, f_tn, criterion, quiet=True, title=""):
        names, cols, vars_ = build_columns(cps, var, f_rn, r_th, f_tn)
        A = np.zeros((len(e_meas), len(cols)))
        V = np.zeros_like(A)
        for k, (c, v) in enumerate(zip(cols, vars_)):
            A[:, k] = ftc.fl.rebin_model_to_meas(np.arange(len(c)) + 0.5, c, e_meas)
            V[:, k] = ftc.fl.rebin_model_to_meas(np.arange(len(v)) + 0.5, v, e_meas)
        sel = (e_meas >= ftc.E_LO) & (e_meas < ftc.E_HI)
        A, V, y, e_sel = A[sel], V[sel], cnt[sel], e_meas[sel]
        A_counts, VAR_counts = A * live, V * live * live
        w = 1.0/np.sqrt(np.maximum(y,1.0)) if criterion=="A" else 1.0/np.maximum(y,1.0)
        note = ("веса пуассоновские" if criterion == "A" else "веса относительные") + \
               ", цепочки связаны физикой (см. докстроку)"
        if quiet:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                amp, sd, pred, chi2ndf, shape = ftc.fit(A_counts, y, w, names, "", "", var_counts=VAR_counts)
        else:
            amp, sd, pred, chi2ndf, shape = ftc.fit(A_counts, y, w, names, title, note, var_counts=VAR_counts)
        return amp, sd, pred, chi2ndf, shape, names

    best_A = None
    best_B = None
    min_chi2ndf = float('inf')
    min_shape = float('inf')

    total = len(GRID_F_RN) * len(GRID_R_TH) * len(GRID_F_TN)
    count = 0

    for f_rn in GRID_F_RN:
        for r_th in GRID_R_TH:
            for f_tn in GRID_F_TN:
                count += 1
                if count % 200 == 0:
                    print(f"Прогресс: {count}/{total} узлов")
                amp_A, sd_A, pred_A, chi2ndf_A, shape_A, names_A = evaluate(f_rn, r_th, f_tn, "A")
                if chi2ndf_A < min_chi2ndf:
                    min_chi2ndf = chi2ndf_A
                    best_A = (f_rn, r_th, f_tn, amp_A, sd_A, pred_A, chi2ndf_A, shape_A)
                amp_B, sd_B, pred_B, chi2ndf_B, shape_B, names_B = evaluate(f_rn, r_th, f_tn, "B")
                if shape_B < min_shape:
                    min_shape = shape_B
                    best_B = (f_rn, r_th, f_tn, amp_B, sd_B, pred_B, chi2ndf_B, shape_B)

    print("\n=== ЛУЧШИЙ УЗЕЛ ПО КРИТЕРИЮ A ===")
    f_rn_A, r_th_A, f_tn_A, amp_A, sd_A, pred_A, chi2ndf_A, shape_A = best_A
    print(f"f_Rn = {f_rn_A:.3f}, r_Th = {r_th_A:.3f}, f_Tn = {f_tn_A:.3f}")
    # Повторный вызов БЕЗ подавления stdout — печатает таблицу амплитуд.
    evaluate(f_rn_A, r_th_A, f_tn_A, "A", quiet=False,
             title="Лучший узел по критерию A (цепочки связаны физикой)")

    print("\n=== ВОССТАНОВЛЕННЫЕ АКТИВНОСТИ ===")
    if NO_MUON:
        A_K, A_Ra, A_Th = amp_A
        A_mu = float("nan")
    else:
        A_K, A_Ra, A_Th, A_mu = amp_A
    A_Pb214 = A_Bi214 = A_Ra * (1 - f_rn_A)
    A_Pb212 = A_Bi212 = A_Th * r_th_A * (1 - f_tn_A)
    A_Tl208 = TL208_BRANCH * A_Bi212
    print("Нуклид | Бк/кг")
    print(f"K40     | {A_K:.3f}")
    print(f"Ra226   | {A_Ra:.3f}")
    print(f"Pb214   | {A_Pb214:.3f}")
    print(f"Bi214   | {A_Bi214:.3f}")
    print(f"Pb212   | {A_Pb212:.3f}")
    print(f"Ac228   | {A_Th:.3f}")
    print(f"Bi212   | {A_Bi212:.3f}")
    print(f"Tl208   | {A_Tl208:.3f}")

    _sel = (e_meas >= ftc.E_LO) & (e_meas < ftc.E_HI)
    ftc.bands_report(pred_A, cnt[_sel], e_meas[_sel], live)
    if os.environ.get("G4M_BANDS_PARTS") == "1":
        bands_by_component(cps, var, f_rn_A, r_th_A, f_tn_A, amp_A,
                           e_meas[_sel], cnt[_sel], live)

    print("\n=== ПРОВЕРКА ФИЗИЧНОСТИ ===")
    if A_Pb214 <= A_Ra:
        print("✓ A(Pb214) <= A(Ra226)")
    else:
        print("✗ ДЕФЕКТ МОДЕЛИ: A(Pb214) > A(Ra226)")
    ratio = A_Tl208 / A_Bi212
    if abs(ratio - TL208_BRANCH) < 0.001:
        print("✓ A(Tl208)/A(Bi212) = 0.3594")
    else:
        print(f"✗ A(Tl208)/A(Bi212) = {ratio:.4f}, ожидалось 0.3594")

    if 0.03 <= f_rn_A <= 0.24:
        print("✓ f_Rn в пределах литературы")
    elif f_rn_A < 0.03:
        print("✗ f_Rn ниже литературного диапазона")
    else:
        print("✗ f_Rn выше литературного диапазона")

    if abs(f_tn_A) < 0.01:
        print("✓ f_Tn ≈ 0 (торон не успевает уйти)")
    else:
        print(f"⚠ f_Tn = {f_tn_A:.3f}")

    print("\n=== ГРАНИЦЫ ПАРАМЕТРОВ (профиль) ===")
    f_rn_opt, r_th_opt, f_tn_opt, amp_B, sd_B, pred_B, chi2ndf_B, shape_B = best_B
    min_profile = shape_B

    print("f_Rn:")
    prof_frn = []
    for fr in GRID_F_RN:
        amp_prof, _, _, _, shape_prof, _ = evaluate(fr, r_th_opt, f_tn_opt, "B")
        if shape_prof <= min_profile * 1.01:
            prof_frn.append(fr)
    print(f"  диапазон: [{prof_frn[0]:.3f}, {prof_frn[-1]:.3f}]")

    print("r_Th:")
    prof_rth = []
    for rt in GRID_R_TH:
        amp_prof, _, _, _, shape_prof, _ = evaluate(f_rn_opt, rt, f_tn_opt, "B")
        if shape_prof <= min_profile * 1.01:
            prof_rth.append(rt)
    print(f"  диапазон: [{prof_rth[0]:.3f}, {prof_rth[-1]:.3f}]")

    print("f_Tn:")
    prof_fnt = []
    for ft in GRID_F_TN:
        amp_prof, _, _, _, shape_prof, _ = evaluate(f_rn_opt, r_th_opt, ft, "B")
        if shape_prof <= min_profile * 1.01:
            prof_fnt.append(ft)
    print(f"  диапазон: [{prof_fnt[0]:.3f}, {prof_fnt[-1]:.3f}]")

    print("\n=== ТРЕБУЕТ ТОЛКОВАНИЯ ===")
    if f_rn_A in [GRID_F_RN[0], GRID_F_RN[-1]]:
        print("⚠ f_Rn упёрся в край сетки")
    if r_th_A in [GRID_R_TH[0], GRID_R_TH[-1]]:
        print("⚠ r_Th упёрся в край сетки")
    if f_tn_A in [GRID_F_TN[0], GRID_F_TN[-1]]:
        print("⚠ f_Tn упёрся в край сетки")
    if chi2ndf_A > 10:
        print("⚠ chi2/ndf > 10")
    if any(a < 1e-10 for a in amp_A):
        print("⚠ одна из амплитуд обнулена NNLS")

if __name__ == "__main__":
    main()
