# -*- coding: utf-8 -*-
r"""#GS-42 (29.09.2026): готовые компоненты модели GS2020 сверх шаблонов М1 и линий М2 (#PHYS-1); Geant4 не запускается.
β/e⁻: <папка образца>\beta_merged\beta_<нуклид>.csv — распад иона, фотоны распада убиты при рождении (γ и рентген распада
уже в библиотеке М2), электроны β/конверсии/Оже перенесены; на распад = отсчёты / n_events_processed.
IB (внутреннее тормозное, KUB): <папка>\ib_merged\ib_<нуклид>.csv — фотоны из primary_table; на распад = отсчёты·Y/table_drawn,
Y = Σ dNdk·шаг таблицы. У α-излучателей ALPHA IB нет — явный список, не ошибка.
Отсутствие ожидаемого файла/ключа шапки, чужая геометрия — ОТКАЗ (SystemExit), не тихий ноль. GS_EXTRA=0 — всё выключено.
М1 берёт только IB (β уже в распаде иона), М2 — β + IB. Потребители: fit_gs2020_kcl*.py, fit_gs2020_th232_m*.py."""
import os, sys, glob, contextlib
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from merge_templates_gs2020 import read_chunk   # §33: тот же разбор шаблона, что у сумматора отрезков

ENABLED = os.environ.get("GS_EXTRA", "1") != "0"
IB_OFF = os.environ.get("GS_IB", "0") != "1"   # 01.10: внутреннее тормозное (IB) по умолчанию ВЫКЛЮЧЕНО везде; GS_IB=1 — вернуть (KUB-оценка)
ALPHA = ("Th232", "Th228", "Ra224", "Rn220")    # α-излучатели цепочки Th-232: IB β-распада у них нет
MATCH = ("em_cut_mm", "em_deex", "em_option", "vessel", "src_mode", "sample_matrix", "sample_rho_g_cm3", "npsm_enabled")
BANDS = (("lt150", -1e30, 150.0), ("150_400", 150.0, 400.0), ("gt400", 400.0, 1e30))   # кэВ, по всей шкале спектра


def _fail(msg):
    raise SystemExit("ОТКАЗ (#GS-42): " + msg)


def _int(h, key, path):
    if key not in h:
        _fail("в шапке %s нет ключа %s" % (path, key))
    return int(h[key])


def _read(path, what):
    """Шапка, GDML из строки-комментария и гистограмма count_edep {центр бина 1 кэВ: отсчёты > 0} (как g1s.read_template
    при npsm_enabled=0 — единственный режим шаблонов GS2020)."""
    if not os.path.exists(path):
        _fail("нет %s: %s" % (what, path))
    c = read_chunk(path)
    h = dict(c["header"])
    gdml = [x.split("GDML", 1)[1].strip() for x in c["comments"] if "GDML" in x]
    if h.get("npsm_enabled") != "0":
        _fail("%s: npsm_enabled=%s, ожидается 0 (М1 читает count_edep)" % (path, h.get("npsm_enabled")))
    return h, (gdml[0] if gdml else None), c


def _hist(c, path):
    """Гистограмма count_edep разобранного шаблона: {центр бина: отсчёты > 0}; сетка 1 кэВ с центрами x,5 обязательна
    (сложение с шаблоном М1 по ключам бинов)."""
    hist = {}
    for b, n in zip(c["bins"], c["edep"]):
        k = float(b)
        if n < 0 or abs((k - 0.5) - round(k - 0.5)) > 1e-6:
            _fail("%s: бин %s (счёт %d) — отрицательный счёт или вне сетки 1 кэВ с центрами x,5" % (path, b, n))
        if n > 0:
            hist[k] = float(n)
    return hist


def table_yield(path):
    """Y = Σ dNdk·шаг таблицы primary_table (CSV k_keV,dNdk_per_decay_per_keV; строки '#' — комментарии).
    Шаг обязан быть постоянным и положительным — иначе «Σ·шаг» не определён, ОТКАЗ."""
    if not os.path.exists(path):
        _fail("нет таблицы IB primary_table: %s" % path)
    k, d = [], []
    with open(path, encoding="utf-8") as f:
        for ln in f:
            ln = ln.strip()
            if ln and not ln.startswith("#") and not ln.startswith("k_keV"):
                a, b = ln.split(",")[:2]
                k.append(float(a)); d.append(float(b))
    st = np.diff(k)
    if len(k) < 2 or np.any(st <= 0) or float(np.ptp(st)) > 1e-6 * float(st[0]) or min(d) < 0:
        _fail("таблица %s: строк %d, шаг непостоянен/не возрастает или плотность < 0" % (path, len(k)))
    return float(np.sum(d) * st[0]), float(st[0]), len(k)


def _ref(folder, nuc):
    """Шапка шаблона М1 того же образца для сверки геометрии: mix_<нуклид>_npsmoff.csv, а у нуклида без своего шаблона
    М1 (Ra-228) — первый mix_*_npsmoff.csv с GDML в шапке (шаблон цепочки её не несёт)."""
    own = os.path.join(folder, "mix_%s_npsmoff.csv" % nuc)
    for p in ([own] if os.path.exists(own) else sorted(glob.glob(os.path.join(folder, "mix_*_npsmoff.csv")))):
        h, g, _ = _read(p, "шаблона М1")
        if g:
            return p, h, g
    _fail("в %s нет шаблона М1 mix_*_npsmoff.csv с GDML для сверки геометрии" % folder)


def _check_geo(path, h, gdml, ref):
    """Геометрия и физ-лист β/IB-шаблона обязаны совпасть с шаблоном М1 того же образца (GDML + ключи MATCH)."""
    rp, rh, rg = ref
    bad = ["GDML %s ≠ %s" % (gdml, rg)] if gdml != rg else []
    bad += ["%s %s ≠ %s" % (k, h.get(k), rh.get(k)) for k in MATCH if h.get(k) != rh.get(k)]
    if bad:
        _fail("%s не совпадает с шаблоном М1 %s: %s" % (path, rp, "; ".join(bad)))


def _load_ib(nuc, folder, ref):
    """IB нуклида: n_eff = table_drawn / Y — «распадов» на отсчёт шаблона; table_drawn обязан равняться n_events_processed
    (один разыгранный фотон на событие), иначе нормировка неоднозначна — ОТКАЗ."""
    p = os.path.join(folder, "ib_merged", "ib_%s.csv" % nuc)
    h, g, c = _read(p, "IB-шаблона (β⁻-излучатель %s)" % nuc)
    _check_geo(p, h, g, ref)
    drawn, n = _int(h, "table_drawn", p), _int(h, "n_events_processed", p)
    if drawn != n or drawn <= 0 or "primary_table" not in h:
        _fail("%s: table_drawn=%d, n_events_processed=%d, primary_table=%s" % (p, drawn, n, h.get("primary_table")))
    Y, step, rows = table_yield(h["primary_table"])
    if not Y > 0:
        _fail("%s: интеграл таблицы IB Y=%g не положителен" % (h["primary_table"], Y))
    return {"hist": _hist(c, p), "n_eff": drawn / Y, "file": p, "Y": Y, "drawn": drawn, "table": h["primary_table"],
            "table_step_keV": step, "table_rows": rows}


def _load_beta(nuc, folder, ref):
    """β/e⁻-компонента: n_eff = n_events_processed (распадов иона); шапка обязана нести beta_only=1."""
    p = os.path.join(folder, "beta_merged", "beta_%s.csv" % nuc)
    h, g, c = _read(p, "β-шаблона %s" % nuc)
    _check_geo(p, h, g, ref)
    n = _int(h, "n_events_processed", p)
    if h.get("beta_only") != "1" or n <= 0:
        _fail("%s: beta_only=%s, n_events_processed=%d — не β-шаблон" % (p, h.get("beta_only"), n))
    return {"hist": _hist(c, p), "n_eff": float(n), "file": p, "n": n}


def sums(comp):
    """Σ по спектру на распад (событий с энерговыделением в кристалле): (β, IB); нет части — 0."""
    return tuple(sum(pt["hist"].values()) / pt["n_eff"] if pt else 0.0 for pt in (comp["beta"], comp["ib"]))


def load(nuc, folder, beta=True):
    """β + IB нуклида на распад на сетке шаблонов (бины 1 кэВ по энерговыделению в кристалле). beta=False — только IB
    (метод 1). Возвращает {"nuc", "beta": часть|None, "ib": часть|None, "per_decay": {бин: β/N + IB·Y/drawn}};
    часть — {"hist": {бин: отсчёты}, "n_eff": делитель на распад, "file", …}. Печатает строку ФИЗИКА (#GS-42)."""
    ref = _ref(folder, nuc)
    out = {"nuc": nuc, "beta": _load_beta(nuc, folder, ref) if beta else None,
           "ib": None if nuc in ALPHA else _load_ib(nuc, folder, ref)}
    if IB_OFF and out["ib"]:   # 01.10 (оператор «да, везде»): IB убран из расчёта до числового подтверждения — гистограмма обнуляется, структура прежняя
        out["ib"] = dict(out["ib"], hist={k: 0.0 for k in out["ib"]["hist"]})
    tb = ("%s, распадов %d" % (out["beta"]["file"], out["beta"]["n"]) if beta
          else "не добавляется — уже в шаблоне М1 (распад иона)")
    ti = ("нет — α-излучатель" if out["ib"] is None
          else "%s, Y=%.6e, drawn=%d" % (out["ib"]["file"], out["ib"]["Y"], out["ib"]["drawn"]))
    per = {}
    for pt in (out["beta"], out["ib"]):
        for k, v in (pt["hist"].items() if pt else ()):
            per[k] = per.get(k, 0.0) + v / pt["n_eff"]
    out["per_decay"] = per
    print("ФИЗИКА (#GS-42): %s β: %s; IB: %s; Σ на распад: β %.4e, IB %.4e" % ((nuc, tb, ti) + sums(out)))
    return out


def fold(comp, broaden, ch_edges, fwhm):
    """Та же свёртка разрешения и перевод в каналы, что у шаблонов М1 (mix_unfold_core.unfold → g1s.broaden(hist, n,
    ch_edges, blur·ПШПВ)): передать g1s.broaden и ту же функцию ширины. Возвращает {"beta"|"ib": (столбец на распад,
    дисперсия на распад²)}; дисперсия — столбец/n_eff, та же форма, что у A2 (cols/n_events) и V2 метода 2; нет части — нули."""
    res = {}
    for key in ("beta", "ib"):
        pt = comp[key]
        col = np.zeros(len(ch_edges) - 1) if pt is None else broaden(pt["hist"], 1.0, ch_edges, fwhm) / pt["n_eff"]
        if pt is None:
            var = col.copy()
        elif os.environ.get("GS_TVAR", "exact") == "raw":
            var = col / pt["n_eff"]   # прежнее: верно только без размытия
        else:   # #GS-63: точная дисперсия размытой части Σ N_E·f_i(E)²/n_eff² (mix_unfold_core.template_var)
            import mix_unfold_core
            var = mix_unfold_core.template_var(pt["hist"], pt["n_eff"], ch_edges, fwhm, broaden)
        res[key] = (col, var)
    return res


def node_var_of(build_out, ch_edges, fwhm):
    """#GS-63: var_of для run_method2 — точная дисперсия размытой формы узла сетки на событие² (тот же отбор бинов ≥ 1 кэВ
    и та же свёртка g1s.broaden, что grid_response); GS_TVAR=raw — None (прежняя shape/n)."""
    if os.environ.get("GS_TVAR", "exact") == "raw":
        return None
    import mix_unfold_core as muc
    files = {float(os.path.basename(p)[len("grid_mar_E"):-len(".csv")]): p for p in glob.glob(os.path.join(build_out, "grid_mar_E*.csv"))}
    cache = {}

    def var_of(E):
        k = min(files, key=lambda x: abs(x - E))
        if abs(k - E) > 0.01:
            _fail("нет прогона сетки для %.3f кэВ (ближайший %.3f)" % (E, k))
        if k not in cache:
            hist, n, _ = muc.g1s.read_template(files[k])
            cache[k] = muc.template_var({a: v for a, v in hist.items() if a >= 1.0}, n, ch_edges, fwhm)
        return cache[k]
    return var_of


def band_counts(col, e, sel):
    """Отсчёты столбца (уже × амплитуда) в окне подгонки sel и в полосах BANDS по всей шкале."""
    e = np.asarray(e, dtype=float)
    d = {"window": float(col[np.asarray(sel, dtype=bool)].sum())}
    d.update({nm: float(col[(e >= lo) & (e < hi)].sum()) for nm, lo, hi in BANDS})
    return d


def add_ib(hist, n, items):
    """Шаблон М1 + IB по формуле донора add_ib_template.py: count = ion + BR·Y·(N_ion/N_ib)·ib, нормировка n = N_ion
    не меняется (N_ib/Y = n_eff). items = [(comp, BR), …]; β не добавляется (в М1 он уже в распаде иона)."""
    out = dict(hist)
    for comp, br in items:
        for k, v in (comp["ib"]["hist"].items() if comp["ib"] else ()):
            out[k] = out.get(k, 0.0) + br * n / comp["ib"]["n_eff"] * v
    return out


@contextlib.contextmanager
def m1_ib(g1s, plan):
    """Метод 1: на время контекста g1s.read_template (им одним mix_unfold_core.unfold читает шаблоны) отдаёт add_ib(…).
    plan = {путь шаблона М1: [(comp, BR), …]}; пустой plan — без подмены (прежнее поведение). Файлы mix_* не меняются.
    Оговорка донора: дисперсия шаблона в A2 — col/N_ion и для IB-части (истинная col·Y/N_ib) — оценивается в подгонке."""
    orig, used = g1s.read_template, set()
    todo = {os.path.normcase(os.path.abspath(p)): v for p, v in plan.items()}

    def patched(path):
        hist, n, npsm = orig(path)
        key = os.path.normcase(os.path.abspath(path))
        used.update([key] if key in todo else [])
        return (add_ib(hist, n, todo[key]) if key in todo else hist), n, npsm
    g1s.read_template = patched if todo else orig
    try:
        yield
    finally:
        g1s.read_template = orig
    if set(todo) - used:
        _fail("шаблоны М1 из плана IB не прочитаны подгонкой: %s" % sorted(set(todo) - used))


def m1_var_excess(items, n_m1, coef, var, sel):
    """Наибольшая по окну доля, на которую приближение донора (дисперсия IB-части шаблона М1 = col/N_ion) меняет дисперсию
    канала критерия A2 против точной BR²·col/n_eff: max_sel |coef²·Σ(BR·col/N_ion − BR²·col/n_eff)| / var.
    items = [(BR, столбец IB на распад нуклида, n_eff)]."""
    if not items:
        return 0.0
    ex = sum(br * c / n_m1 - br * br * c / ne for br, c, ne in items) * coef ** 2
    s = np.asarray(sel, dtype=bool)
    return float(np.max(np.abs(ex[s]) / np.asarray(var)[s]))
