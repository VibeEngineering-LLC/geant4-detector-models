import math, os, sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
import fit_two_criteria as ftc

MATERIALS = ["brick", "concrete"]
TEMPLATE_DIR = os.path.join(HERE, "..", "run_field", "output")
TEMPLATE_FMT = "rc103_field_room_%s_%s.csv"   # % (material, nuclide)
TL208_BRANCH = 0.3594   # доля Bi-212 -> Tl-208 (остальное альфа-канал Po-212)
CHAINS = {
    "K40": (["K40"], {}),
    "Ra226": (["Pb214", "Bi214"], {}),
    "Th232": (["Ac228", "Pb212", "Bi212", "Tl208"], {"Tl208": TL208_BRANCH}),
}

def load_chain_material(chain_name, material):
    """Загружает цепочку нуклидов для заданного материала и возвращает свернутые значения."""
    cps_sum = 0.0
    var_sum = 0.0
    found_any = False
    nuclides, weights = CHAINS[chain_name]
    for nuclide in nuclides:
        path = os.path.join(TEMPLATE_DIR, TEMPLATE_FMT % (material, nuclide))
        if not os.path.exists(path):
            print(f"Предупреждение: файл не найден {path}")
            continue
        found_any = True
        meta, arr, cnt_mc = ftc.read_template(path)
        weight = weights.get(nuclide, 1.0)
        cps_sum += weight * ftc.rcspec.fold(arr, "103")
        var_sum += (weight**2) * ftc.template_variance(cnt_mc, float(meta.get("t_run_s", 0.0)))
    if not found_any:
        return None, None
    return cps_sum, var_sum

def prepare():
    """Собирает шаблоны, подгоняет оба критерия, возвращает всё для печати/рисования."""
    smp = ftc.read_rcxml.read(os.path.join(ftc.MEAS_DIR, ftc.MEAS_NAME))[0]
    cnt = smp.counts[:-1].astype(float)
    ch = np.arange(len(cnt))
    e_meas = np.asarray(sum(c*ch**i for i,c in enumerate(ftc.CAL_ROOM)))
    live = float(smp.live)

    names = []
    cols_raw = []
    vars_raw = []

    for chain_name in CHAINS:
        for material in MATERIALS:
            cps, var = load_chain_material(chain_name, material)
            if cps is not None:
                name = f"{chain_name}_{material}"
                names.append(name)
                cols_raw.append(cps)
                vars_raw.append(var)

    if hasattr(ftc, "MUON_CSV") and os.path.exists(ftc.MUON_CSV):
        meta_mu, raw_mu, cnt_mu = ftc.read_template(ftc.MUON_CSV)
        cps_mu = ftc.rcspec.fold(raw_mu, "103")
        var_mu = ftc.template_variance(cnt_mu, float(meta_mu.get("n_events", 0.0)))
        names.append("mu")
        cols_raw.append(cps_mu)
        vars_raw.append(var_mu)

    A = np.zeros((len(e_meas), len(cols_raw)))
    VAR = np.zeros_like(A)
    for k, (col, var) in enumerate(zip(cols_raw, vars_raw)):
        A[:, k] = ftc.fl.rebin_model_to_meas(np.arange(len(col))+0.5, col, e_meas)
        VAR[:, k] = ftc.fl.rebin_model_to_meas(np.arange(len(var))+0.5, var, e_meas)

    sel = (e_meas >= ftc.E_LO) & (e_meas < ftc.E_HI)
    A = A[sel]
    VAR = VAR[sel]
    cnt = cnt[sel]
    e_meas = e_meas[sel]

    A_counts = A * live
    VAR_counts = VAR * live * live

    cond = ftc.degeneracy_report(A, names)

    # ftc.fit() возвращает (amp, sd, pred, chi2ndf, shape) — позиционно, не словарь.
    weights_a = 1.0/np.sqrt(np.maximum(cnt,1.0))
    amp_a, sd_a, pred_a, chi2ndf_a, shape_a = ftc.fit(
        A_counts, cnt, weights_a, names,
        "A — критерий chi2/ndf (цепочки, раздельные материалы)",
        "веса пуассоновские sigma_i=sqrt(N_i) плюс дисперсия шаблонов (итеративно)",
        var_counts=VAR_counts)
    ftc.bands_report(pred_a, cnt, e_meas, live)
    d_a = dict(zip(names, amp_a))
    ds_a = dict(zip(names, sd_a))

    weights_b = 1.0/np.maximum(cnt,1.0)
    amp_b, sd_b, pred_b, chi2ndf_b, shape_b = ftc.fit(
        A_counts, cnt, weights_b, names,
        "B — критерий невязки формы (цепочки, раздельные материалы)",
        "веса относительные sigma_i=N_i плюс дисперсия шаблонов (итеративно)",
        var_counts=VAR_counts)
    ftc.bands_report(pred_b, cnt, e_meas, live)
    d_b = dict(zip(names, amp_b))
    ds_b = dict(zip(names, sd_b))

    return dict(names=names, e_meas=e_meas, cnt=cnt, live=live, A_counts=A_counts,
                cond=cond, amp_a=amp_a, sd_a=sd_a, pred_a=pred_a, chi2ndf_a=chi2ndf_a,
                shape_a=shape_a, d_a=d_a, ds_a=ds_a, amp_b=amp_b, sd_b=sd_b,
                pred_b=pred_b, chi2ndf_b=chi2ndf_b, shape_b=shape_b, d_b=d_b, ds_b=ds_b)


def main():
    r = prepare()
    (names, e_meas, cnt, live, cond) = (r["names"], r["e_meas"], r["cnt"], r["live"], r["cond"])
    (pred_a, chi2ndf_a, d_a, ds_a) = (r["pred_a"], r["chi2ndf_a"], r["d_a"], r["ds_a"])
    (pred_b, chi2ndf_b, d_b, ds_b) = (r["pred_b"], r["chi2ndf_b"], r["d_b"], r["ds_b"])

    print("\nВЗАИМНАЯ АКТИВНОСТЬ (кирпич/бетон)")
    print("проверка: NNLS сам увидел ли, что кирпич активнее бетона — ожидание оператора: кирпич АКТИВНЕЕ")
    for chain_name in CHAINS:
        for label, d, ds in (("A", d_a, ds_a), ("B", d_b, ds_b)):
            nb, nc = f"{chain_name}_brick", f"{chain_name}_concrete"
            if d.get(nb, 0) > 0 and d.get(nc, 0) > 0:
                r = d[nb] / d[nc]
                er = r * math.sqrt((ds[nb]/d[nb])**2 + (ds[nc]/d[nc])**2)
                print(f"  {chain_name} ({label}): {r:.3f} +- {er:.3f}")

    print("\nРАВНОВЕСИЕ ЦЕПОЧЕК (Th/Ra) ПО МАТЕРИАЛАМ")
    print("в кирпиче Th/Ra > 1 (кирпич); типичный диапазон по независимым региональным исследованиям 0.5-1.7 (источники: Nigeria concrete ~1.49, Nigeria brick ~0.48, Egypt brick ~1.72, Nepal brick ~1.57 — большой региональный разброс, не единое число)")
    for material in MATERIALS:
        for label, d, ds in (("A", d_a, ds_a), ("B", d_b, ds_b)):
            nt, nr = f"Th232_{material}", f"Ra226_{material}"
            if d.get(nt, 0) > 0 and d.get(nr, 0) > 0:
                r = d[nt] / d[nr]
                er = r * math.sqrt((ds[nt]/d[nt])**2 + (ds[nr]/d[nr])**2)
                print(f"  {material} ({label}): {r:.3f} +- {er:.3f}")

    print("\nТРЕБУЕТ ТОЛКОВАНИЯ")
    trig = False
    for n in names:
        if d_a.get(n, 0) == 0 or d_b.get(n, 0) == 0:
            print(f"  {n} обнулён в одном из разложений (A={d_a.get(n,0):.2f}, B={d_b.get(n,0):.2f})")
            trig = True
        elif n in d_a and n in d_b and d_a[n] > 0:
            ratio = d_b[n] / d_a[n]
            if ratio > 2.0 or ratio < 0.5:
                print(f"  {n} расходится между критериями более чем вдвое (A={d_a[n]:.2f}, B={d_b[n]:.2f})")
                trig = True
    if max(chi2ndf_a, chi2ndf_b) > 10:
        print(f"  max(chi2/ndf) = {max(chi2ndf_a, chi2ndf_b):.1f} > 10")
        trig = True
    if not trig:
        print("  (пусто — ни один триггер не сработал)")

    print("\nСравнение вырожденности с прежними 9-параметрическими разложениями:")
    print(f"  cond (7 параметров, эта модель)        = {cond:.1f}")
    print(f"  cond (комната, 9 параметров, прошлый прогон этой сессии) = 246.4")
    print(f"  cond (сфера, 9 параметров, прошлый прогон этой сессии)   = 271.2")
    print(f"  меньше ли вырожденность у этой модели: {cond < 246.4}")

if __name__ == "__main__":
    main()
