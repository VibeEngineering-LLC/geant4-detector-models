import os, sys, pickle, json, argparse
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))              # .../Gamma-1S/analysis
WEB = os.path.join(os.path.dirname(HERE), "web-th232")
os.environ["G4MODELS_SOURCE_CONFIG"] = os.path.join(WEB, "configs", "amticseu.yaml")
if not os.environ.get("SPECTRAVIBE_ROOT"):
    raise SystemExit("ОТКАЗ: задайте SPECTRAVIBE_ROOT (нужен модулю export_ra226_data)")
sys.path.insert(0, HERE)
sys.path.insert(0, WEB)
import crit_bench as cb
import export_ra226_data as erd

def _max_rel(x, ref, name):
    """Максимум относительного расхождения по ненулевым элементам эталона; иначе громкий отказ."""
    x = np.asarray(x, dtype=float)
    ref = np.asarray(ref, dtype=float)
    if x.shape != ref.shape:
        raise SystemExit(f"ОТКАЗ гейта {name}: форма {x.shape} против эталона {ref.shape}")
    nz = ref != 0
    if not nz.any():
        raise SystemExit(f"ОТКАЗ гейта {name}: эталон целиком нулевой")
    if np.any(x[~nz] != 0):
        raise SystemExit(f"ОТКАЗ гейта {name}: ненулевые значения там, где эталон нулевой")
    return float(np.max(np.abs(x[nz] - ref[nz]) / np.abs(ref[nz])))

def load_dump(path):
    if not os.path.exists(path):
        raise SystemExit(f"ОТКАЗ: нет файла дампа {path}")
    with open(path, "rb") as f:
        return pickle.load(f)

def build_cols(d, nodes, pair):
    def resp(E):
        if abs(E - 511.0) < 0.01:
            return (pair[0], {}, pair[1])
        closest = min(nodes.keys(), key=lambda x: abs(x - E))
        if abs(closest - E) > 0.01:
            raise SystemExit(f"ОТКАЗ: не найден узел для энергии {E}")
        return (nodes[closest][0], {}, nodes[closest][1])

    def n_of(E):
        if abs(E - 511.0) < 0.01:
            return d["pair_n"]
        closest = min(d["n_map"].keys(), key=lambda x: abs(x - E))
        if abs(closest - E) > 0.01:
            raise SystemExit(f"ОТКАЗ: не найдено число событий для энергии {E}")
        return d["n_map"][closest]

    acc = {}
    _, w, _, _ = erd.run_method2(d["lib_m2"], d["sum_peaks"], resp, d["e"], d["ch_edges"], d["keys"], var_acc=acc, n_of=n_of)
    cols = np.array([w[key] for key in d["keys"]])
    V = np.array([acc[key] for key in d["keys"]])
    return cols, V

def gate_build(d):
    cols, V = build_cols(d, d["nodes"], d["pair"])
    cols2_ref = d["cols2_ref"]
    V2_ref = d["V2_ref"]

    m_cols = _max_rel(cols, d["cols2_ref"], "cols")
    m_V = _max_rel(V, d["V2_ref"], "V")

    if m_cols > 1e-12 or m_V > 1e-12:
        raise SystemExit(f"ОТКАЗ: расхождение столбцов или дисперсий больше 1e-12: cols={m_cols:.3e}, V={m_V:.3e}")

    print(f"Максимальное расхождение столбцов: {m_cols:.3e}")
    print(f"Максимальное расхождение дисперсий: {m_V:.3e}")

    return cols, V

def fit_A2V(cols, counts, bg_scaled, k, V, sel):
    K = cols.shape[0]
    a1, _, ex1 = cb.fit_A1(cols, counts, bg_scaled, k, np.ones(K), sel)
    C, p, b = cb._prep(cols, counts, bg_scaled, sel)
    a, sd = cb._iter_template_var(C, p - b, ex1["var"], np.asarray(V)[:, sel].T, np.ones(K), a1)
    return a, sd, {"iterations": 6}

def noisy_nodes(d, rng):
    nodes_new = {}
    for E, (s, eps) in d["nodes"].items():
        n = d["n_map"][E]
        s_noisy = rng.poisson(s * n) / n
        nodes_new[E] = (s_noisy, eps)
    pair_new = (rng.poisson(d["pair"][0] * d["pair_n"]) / d["pair_n"], d["pair"][1])
    return nodes_new, pair_new

def fit_real(d, cols, V):
    real = {}
    for name, fit_func in [("A1", cb.fit_A1), ("A2V", fit_A2V), ("E1", cb.fit_E1)]:
        if name == "A2V":
            a, sd, extra = fit_func(cols, d["counts"], d["bg_scaled"], d["k"], V, d["sel"])
        else:
            a, sd, extra = fit_func(cols, d["counts"], d["bg_scaled"], d["k"], np.ones(len(d["keys"])), d["sel"])
        metrics = cb.metrics_all(a, cols, d["counts"], d["bg_scaled"], d["k"], d["sel"], d["e"], d["keys"], d["passport"], d["live_s"])
        real[name] = {"a": a, "sd": sd, "metrics": metrics}
        print(f"\n{name}:")
        for i, key in enumerate(d["keys"]):
            ap = d["passport"][key]
            a_over_passport = metrics["A_over_passport"][key]
            print(f"  {key}: A/паспорт = {a_over_passport:.3f}, 100*sd/a = {100*sd[i]/a[i] if (a[i] > 0 and np.isfinite(sd[i])) else 'н/д'}%")
        print(f"  chi2_ref_nu: {metrics['chi2_ref_nu']:.3f}, shape: {metrics['shape']:.3f}")
    return real

def closure_m2(d, cols, V, a_true, n_rep, seed, node_noise):
    rng = np.random.default_rng(seed)
    bg_raw = d["bg_scaled"] / d["k"]
    mu = cols.T @ a_true + d["k"] * bg_raw
    est = {name: np.zeros((n_rep, len(d["keys"]))) for name in ["A1", "A2V", "E1"]}
    err = {name: np.zeros((n_rep, len(d["keys"]))) for name in ["A1", "A2V", "E1"]}

    for i in range(n_rep):
        b_synth = d["k"] * rng.poisson(bg_raw)
        counts_synth = rng.poisson(mu).astype(float)

        if node_noise:
            cols_fit, V_fit = build_cols(d, *noisy_nodes(d, rng))
        else:
            cols_fit, V_fit = cols, V

        for name, fit_func in [("A1", cb.fit_A1), ("A2V", fit_A2V), ("E1", cb.make_fit_E1(10))]:
            if name == "A2V":
                a, sd, _ = fit_func(cols_fit, counts_synth, b_synth, d["k"], V_fit, d["sel"])
            else:
                a, sd, _ = fit_func(cols_fit, counts_synth, b_synth, d["k"], np.ones(len(d["keys"])), d["sel"])
            est[name][i] = a
            err[name][i] = sd

    return cb._closure_stats(est, err, a_true, n_rep)

def flags(d, real, cl0, cl1, n_rep, names, m1):
    out = cb._flags_closure("без шума узлов", cl0, n_rep, names) + cb._flags_closure("с шумом узлов", cl1, n_rep, names)
    for key in names:
        x = real["A2V"]["metrics"]["A_over_passport"][key]
        z = real["E1"]["metrics"]["A_over_passport"][key]
        rel = 100.0 * abs(z - x) / abs(x)
        if rel > 5.0:
            out.append(f"{key}: E1 {z:.4f} против A2V {x:.4f} — расхождение {rel:.1f} %")
        if m1 is not None:
            rel1 = 100.0 * abs(m1[key] - x) / abs(x)
            if rel1 > 5.0:
                out.append(f"{key}: метод 1 {m1[key]:.4f} против метода 2 A2V {x:.4f} — расхождение {rel1:.1f} %")
    out.append("шум весов депопуляции и пиковых эффективностей не разыгрывается")
    out.append("корреляция соседних каналов после свёртки не учтена (как в методе 1)")
    return out

def selftest():
    N = 400; x = np.arange(N); n = 1e5
    centers = [50, 110, 170, 230, 290, 350]; w = [0.5, 0.3, 0.2]
    S = []
    for i in range(6):
        s = np.exp(-0.5*((x - centers[i])/4.0)**2)
        s *= 0.05 / s.sum()
        S.append(s)
    def assemble(Sl):
        cols = np.array([ sum(w[j]*Sl[j] for j in range(3)), sum(w[j]*Sl[3+j] for j in range(3)) ])
        V = np.array([ sum(w[j]**2*Sl[j]/n for j in range(3)), sum(w[j]**2*Sl[3+j]/n for j in range(3)) ])
        return cols, V
    cols, V = assemble(S)
    a_true = np.array([2e6, 1e6]); k = 0.5; bg_raw = np.full(N, 5.0); bg_scaled = k*bg_raw
    sel = np.zeros(N, dtype=bool); sel[10:391] = True
    rng = np.random.default_rng(1)
    counts = rng.poisson(cols.T @ a_true + k*bg_raw).astype(float)
    a, sd, _ = fit_A2V(cols, counts, bg_scaled, k, V, sel)
    for i in range(len(a_true)):
        if abs(a[i] - a_true[i]) > 3*sd[i]:
            raise SystemExit(f"ОТКАЗ: отклонение A2V от истинного значения превышает 3σ: {a[i]} vs {a_true[i]}")
    for j in range(6):
        nuc = 0 if j < 3 else 1
        ch = centers[j]
        r = V[nuc, ch] / (cols[nuc, ch] / n)
        if abs(r - 1) <= 0.10:
            raise SystemExit(f"ОТКАЗ: дисперсия метода 2 не отличается от cols/n для узла {j}, r={r}")
    # Исправлено координатором 13.09: генерация передавала в _closure_stats нулевые sd — покрытие
    # было нулём по построению при любом коде; и требовала «A1 хуже у каждого» вместо «хотя бы у одного».
    est = {"A1": np.zeros((40, 2)), "A2V": np.zeros((40, 2))}
    err = {"A1": np.zeros((40, 2)), "A2V": np.zeros((40, 2))}
    for i in range(40):
        b_synth = k*rng.poisson(bg_raw)
        counts_synth = rng.poisson(cols.T @ a_true + k*bg_raw).astype(float)
        Sn = [rng.poisson(s*n)/n for s in S]
        cf, Vf = assemble(Sn)
        est["A1"][i], err["A1"][i], _ = cb.fit_A1(cf, counts_synth, b_synth, k, np.ones(2), sel)
        est["A2V"][i], err["A2V"][i], _ = fit_A2V(cf, counts_synth, b_synth, k, Vf, sel)
    st = cb._closure_stats(est, err, a_true, 40)
    cov = {c: [st[c][i]["coverage_68"] for i in range(2)] for c in ("A1", "A2V")}
    print(f"самопроверка 3: coverage_68 A1 = {cov['A1']}, A2V = {cov['A2V']}")
    for i in range(2):
        if not (0.50 <= cov["A2V"][i] <= 0.85):
            raise SystemExit(f"ОТКАЗ: coverage_68 A2V для нуклида {i} = {cov['A2V'][i]:.3f} вне [0.50; 0.85]")
    if not any(cov["A1"][i] < cov["A2V"][i] for i in range(2)):
        raise SystemExit(f"ОТКАЗ: A1 не хуже A2V ни у одного нуклида: A1 {cov['A1']}, A2V {cov['A2V']}")
    print("SELFTEST OK")

def _jsonable(o):
    """Рекурсивно: ndarray -> list, numpy-числа -> float/int, нечисловые float -> None."""
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (float, np.floating)):
        return float(o) if np.isfinite(o) else None
    return o

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dump", type=str)
    parser.add_argument("--json", type=str)
    parser.add_argument("--m1-json", type=str)
    parser.add_argument("--n-rep", type=int, default=100)
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()

    if args.selftest:
        selftest()
        return

    if not args.dump:
        raise SystemExit("ОТКАЗ: укажите --dump")

    d = load_dump(args.dump)
    cols, V = gate_build(d)
    real = fit_real(d, cols, V)
    a_true = real["A1"]["a"]
    cl0 = closure_m2(d, cols, V, a_true, args.n_rep, 20260913, False)
    cl1 = closure_m2(d, cols, V, a_true, args.n_rep, 20260913, True)

    print("\nClosure без шума узлов:")
    for name in ["A1", "A2V", "E1"]:
        print(f"  {name}:")
        for i, key in enumerate(d["keys"]):
            stats = cl0[name][i]
            print(f"    {key}: bias_pct={stats['bias_pct']:.2f}, spread_pct={stats['spread_pct']:.2f}, coverage_68={stats['coverage_68']:.2f}, pull_rms={stats['pull_rms']:.2f}")

    print("\nClosure с шумом узлов:")
    for name in ["A1", "A2V", "E1"]:
        print(f"  {name}:")
        for i, key in enumerate(d["keys"]):
            stats = cl1[name][i]
            print(f"    {key}: bias_pct={stats['bias_pct']:.2f}, spread_pct={stats['spread_pct']:.2f}, coverage_68={stats['coverage_68']:.2f}, pull_rms={stats['pull_rms']:.2f}")

    m1 = None
    if args.m1_json:
        with open(args.m1_json, "r", encoding="utf-8") as f:
            data = json.load(f)
            m1 = {key: data["method1"]["groups"][key]["A_over_passport"] for key in d["keys"]}

    flgs = flags(d, real, cl0, cl1, args.n_rep, d["keys"], m1)

    # Исправлено координатором 13.09: блок толкования обязан печататься в вывод (спека, #SA-4), а не только
    # попадать в JSON; в JSON — сами строки флагов.
    print("\nТРЕБУЕТ ТОЛКОВАНИЯ:")
    for f_line in (flgs or ["пусто"]):
        print(f"  {f_line}")

    result = {
        "real": real,
        "closure_no_node_noise": cl0,
        "closure_node_noise": cl1,
        "requires_interpretation": flgs
    }

    if args.json:
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(_jsonable(result), f, ensure_ascii=False, indent=1)

if __name__ == "__main__":
    main()
