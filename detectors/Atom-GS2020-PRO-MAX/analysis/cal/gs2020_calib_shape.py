# -*- coding: utf-8 -*-
r"""#GS-37: калибровка шкалы GS2020 подгонкой ФОРМЫ спектра по откликам Geant4 (сетка метода 2 grid_mar_E*.csv):
одна гладкая поправка d(E) на спектр, совместно по всем группам линий и мультиплетам (ЛСРМ «совместная калибровка»,
Salathe/Bandstra full-spectrum gain). Свёртка — донор g1s.broaden (Гамма-1С), один раз на линию, кэш .npz.
Запуск: python gs2020_calib_shape.py [теги]; синтетика: GS_SHAPE_SYNTH="θ0,θ1,θ2" python gs2020_calib_shape.py sample"""
import os, sys, csv, json, math, glob
import numpy as np
from scipy.optimize import lsq_linear, least_squares
sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import gs2020_calib as C
import fit_gs2020_th232_m1 as m1          # настраивает sys.path донора, ПШПВ, хвост ядра (env gs2020_fit_env.sh)
import mix_unfold_g1s as g1s
if m1.TAIL is not None:
    g1s.TAIL_T = m1.TAIL
LINES_DIR = r"C:\g4work\gs2020\lines"
GRID_DIR = os.environ.get("GS_SHAPE_GRID", r"C:\g4work\gs2020\run_marinelli\out_v5_oisn10")
OUT_JSON = C.SHAPE_JSON   # его же читает C.energy_axis при GS_CAL_SHAPE=1
CACHE = os.path.join(os.path.dirname(C.OUT_JSON), "cal_shape_cache")
DEG = int(os.environ.get("GS_SHAPE_DEG", "2"))
# узлы, кэВ шкалы файла: по узлу на группу, лишние по спектру снимает AUTOKNOTS; GS_SHAPE_KNOTS="" — полином DEG
KNOTS = [float(x) for x in os.environ.get("GS_SHAPE_KNOTS", "80,243,320,600,750,960,1130,1480,1730,2614").split(",") if x.strip()]
NP = len(KNOTS) if KNOTS else DEG + 1   # число параметров шкалы θ
ALL_KNOTS = list(KNOTS)                  # полный набор; KNOTS/NP сужаются по спектру (AUTOKNOTS) в main
SYNTH_DROP = [x.strip() for x in os.environ.get("GS_SHAPE_SYNTH_DROP", "").split(",") if x.strip()]   # группы без линий в синтетике
FINE_STEP = 0.25

# --- ПШПВ ---
def get_fwhm_func():
    fw_csv = os.path.join(m1.OUT, "fwhm_points_gs2020.csv")
    m1.write_fwhm_csv(fw_csv)
    return g1s.make_fwhm(fw_csv)

FWHM = None # Инициализируется в main

# --- Линии ---
def load(nuc, kind):
    path = os.path.join(LINES_DIR, f"{nuc}_{kind}.csv")
    lines = []
    if not os.path.exists(path):
        return lines
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                e = float(row["energy"])
                i = float(row["intensity"])
                lines.append((e, i))
            except (ValueError, KeyError):
                continue
    return lines

# --- Группы ---
GROUPS_DEF = [
    {"name": "X-K", "lo": 55, "hi": 110, "comps": [("Bi K","214pb","x"), ("Pb K","208tl","x"), ("Po K","214bi","x"), ("Th K","228ac","x")], "cdeg": 3},
    # #DBG-1: окна групп не перекрываются (каналы 270–275 учитывались в 238 и 352 дважды)
    {"name": "238", "lo": 205, "hi": 268, "comps": [("Pb-212","212pb","g"), ("Ra-224","224ra","g"), ("Pb-214","214pb","g")], "cdeg": 2},
    {"name": "352", "lo": 272, "hi": 380, "comps": [("Pb-214","214pb","g"), ("Ac-228","228ac","g"), ("Pb-212","212pb","g"), ("Tl-208","208tl","g")], "cdeg": 2},
    {"name": "609", "lo": 550, "hi": 650, "comps": [("Tl-208","208tl","g"), ("Bi-214","214bi","g")], "cdeg": 2},
    {"name": "727", "lo": 690, "hi": 810, "comps": [("Bi-212","212bi","g"), ("Ac-228","228ac","g"), ("Bi-214","214bi","g")], "cdeg": 2},
    {"name": "911", "lo": 870, "hi": 1000, "comps": [("Ac-228","228ac","g"), ("Bi-214","214bi","g")], "cdeg": 2},
    {"name": "1120", "lo": 1060, "hi": 1200, "comps": [("Bi-214","214bi","g"), ("Ac-228","228ac","g")], "cdeg": 2},
    {"name": "1461", "lo": 1300, "hi": 1640, "comps": [("K-40","40k","g"), ("Ac-228","228ac","g"), ("Bi-214","214bi","g")], "cdeg": 2},
    # #GS-37: Ac-228 1588/1630 заходит хвостом на левый край окна; без него узел 1730 у bg/bgw уходил на +16 кэВ
    {"name": "1764", "lo": 1650, "hi": 1880, "comps": [("Bi-214","214bi","g"), ("Ac-228","228ac","g")], "cdeg": 2},
    {"name": "2204", "lo": 2100, "hi": 2320, "comps": [("Bi-214","214bi","g")], "cdeg": 2},
    {"name": "2614", "lo": 2400, "hi": 2850, "comps": [("Tl-208","208tl","g"), ("Bi-214","214bi","g")], "cdeg": 2},
]

# --- Сетка и кэш отклика ---
GRID = {}
def init_grid():
    global GRID
    for path in glob.glob(os.path.join(GRID_DIR, "grid_mar_E*.csv")):
        name = os.path.basename(path)
        try:
            e0_str = name[len("grid_mar_E"):-4]
            e0 = float(e0_str)
            GRID[e0] = path
        except ValueError:
            continue

def line_cum(E):
    if not GRID: init_grid()
    keys = list(GRID.keys())
    if not keys: return None
    # Ближайший узел
    e0 = min(keys, key=lambda k: abs(k - E))
    if abs(E - e0) > 15:
        print(f"  нет узла для {E:.1f} кэВ")
        return None
    
    path = GRID[e0]
    cache_path = os.path.join(CACHE, f"L_{E:.3f}.npz")
    
    if os.path.exists(cache_path):
        data = np.load(cache_path)
        return data["cum"]

    hist, n, _ = g1s.read_template(path)
    # Оставить бины k >= 1.0
    valid_keys = [k for k in hist.keys() if k >= 1.0]
    if not valid_keys: return None
    
    dE = E - e0
    hist_shifted = {k + dE: v for k, v in hist.items() if k >= 1.0}
    
    fine_edges = np.arange(0, g1s.FINE_E_MAX + FINE_STEP, FINE_STEP)
    dens = g1s.broaden(hist_shifted, n, fine_edges, FWHM)
    cum = np.concatenate([[0.0], np.cumsum(dens)])
    
    os.makedirs(CACHE, exist_ok=True)
    np.savez(cache_path, cum=cum)
    return cum

# --- Подготовка групп для спектра ---
def prepare_groups(coeffs, n_channels):
    """Возвращает список активных групп с отфильтрованными линиями и кэшированными откликами."""
    active_groups = []
    
    # Границы каналов в энергии файла
    ch_edges_file = C.s.channel_to_energy(np.arange(n_channels + 1) - 0.5, coeffs)
    
    for g_def in GROUPS_DEF:
        lo, hi = g_def["lo"], g_def["hi"]
        comps_data = []
        all_comps_valid = True
        
        for comp_name, nuc, kind in g_def["comps"]:
            lines = load(nuc, kind)
            # Фильтр по энергии и интенсивности
            # lo - 2*FWHM(lo) <= E <= hi + 2*FWHM(hi)
            e_lo_lim = lo - 2 * FWHM(lo)
            e_hi_lim = hi + 2 * FWHM(hi)
            
            filtered_lines = [(e, i) for e, i in lines if e_lo_lim <= e <= e_hi_lim]
            
            if not filtered_lines:
                print(f"  {g_def['name']}: {comp_name} — нет линий")
                all_comps_valid = False
                continue
                
            max_i = max(i for _, i in filtered_lines)
            threshold = 0.01 * max_i
            final_lines = [(e, i) for e, i in filtered_lines if i >= threshold]
            
            if not final_lines:
                 print(f"  {g_def['name']}: {comp_name} — нет линий после порога")
                 all_comps_valid = False
                 continue

            # Суммарный отклик компоненты
            cum_total = None
            for e_line, i_line in final_lines:
                c = line_cum(e_line)
                if c is not None:
                    if cum_total is None:
                        cum_total = np.zeros_like(c)
                    cum_total += (i_line / 100.0) * c
            
            if cum_total is not None:
                comps_data.append({"name": comp_name, "cum": cum_total, "lines": final_lines})
            else:
                # Если хотя бы одна линия не найдена, считаем компоненту отсутствующей? 
                # По ТЗ: "Компонента без линий — пропускается". Здесь линии есть, но узла нет.
                # Будем считать, что если cum_total None, то компонента не вносит вклад.
                pass

        # #DBG-1: компонента без линий пропускается, группа остаётся; каналы группы — ПО ШКАЛЕ ФАЙЛА, один раз
        # (от θ не зависят: иначе длина вектора невязок меняется и least_squares не работает)
        cf = 0.5 * (ch_edges_file[:-1] + ch_edges_file[1:])
        idx = np.where((cf >= lo) & (cf <= hi))[0]
        if comps_data and len(idx) > 0:
            active_groups.append({
                "name": g_def["name"],
                "lo": lo,
                "hi": hi,
                "cdeg": g_def["cdeg"],
                "comps": comps_data,
                "idx": idx
            })
            
    return active_groups

def solve_lin(M, y, sy, n_comps):
    """#DBG-1: точный BVLS с нормировкой колонок (trf не доходил до минимума при совпадающих линиях 238,6/241/242)."""
    Mw = M / sy[:, None]
    nrm = np.linalg.norm(Mw, axis=0); nrm[nrm <= 0] = 1.0
    lb = np.concatenate([np.zeros(n_comps), -np.inf * np.ones(M.shape[1] - n_comps)])
    q = lsq_linear(Mw / nrm, y / sy, bounds=(lb, np.inf), method="bvls", max_iter=2000).x
    return q / nrm

# --- Модель и невязки ---
def basis(E_file):
    """Базис поправки d(E) = basis(E)·θ: полином по u (GS_SHAPE_DEG) либо кусочно-линейная по узлам GS_SHAPE_KNOTS
    (#GS-37: реальная шкала «волнистая» — квадратичная оставляла остатки групп до края скана ±6 кэВ)."""
    E = np.atleast_1d(np.asarray(E_file, float))
    if KNOTS:
        return np.stack([np.interp(E, KNOTS, np.eye(NP)[j]) for j in range(NP)], axis=1)   # «шляпки», за краями — константа
    u = (E - 1000) / 1000.0
    return np.stack([u ** j for j in range(NP)], axis=1)

def get_dE(E_file, theta):
    return basis(E_file) @ np.asarray(theta, float)

def build_residuals(theta, coeffs, y, active_groups):
    n = len(y)
    ch = np.arange(n)
    
    # Границы каналов в энергии файла
    edges_file = C.s.channel_to_energy(np.arange(n + 1) - 0.5, coeffs)
    
    # Поправка d(E) на границах
    d_edges = get_dE(edges_file, theta)
    edges_true = edges_file + d_edges
    
    # Центры каналов в истинной энергии
    centers_true = 0.5 * (edges_true[:-1] + edges_true[1:])
    
    all_r = []
    total_params = 0
    
    for g in active_groups:
        lo, hi = g["lo"], g["hi"]
        idx = g["idx"]   # фиксированные каналы группы (prepare_groups)
        
        if len(idx) == 0:
            continue
            
        n_ch_g = len(idx)
        y_g = y[idx]
        
        # Матрица модели M для группы
        # Колонки: компоненты, континуум
        comps = g["comps"]
        n_comps = len(comps)
        cdeg = g["cdeg"]
        n_cont = cdeg + 1
        
        M_g = np.zeros((n_ch_g, n_comps + n_cont))
        
        # Компоненты
        for k, comp in enumerate(comps):
            cum_c = comp["cum"]
            fine_edges = np.arange(0, g1s.FINE_E_MAX + FINE_STEP, FINE_STEP)
            
            # Границы каналов группы в истинной энергии
            e_lo_g = edges_true[idx[0]]
            e_hi_g = edges_true[idx[-1] + 1]
            
            # Интерполяция кумуляты на границы каналов
            # interp возвращает значения в точках x. Нам нужны разности.
            # Границы: idx[0], idx[1], ..., idx[-1]+1
            chan_edges_true = edges_true[idx[0] : idx[-1] + 2]
            
            cum_vals = np.interp(chan_edges_true, fine_edges, cum_c)
            col = np.diff(cum_vals)
            M_g[:, k] = col
            
        # Континуум
        centers_g = centers_true[idx]
        t = (centers_g - 0.5 * (lo + hi)) / (hi - lo)
        for j in range(n_cont):
            M_g[:, n_comps + j] = t ** j
            
        # Линейное решение для амплитуд и континуума
        sy = np.sqrt(np.maximum(y_g, 1.0))
        
        # Bounds: comps >= 0, cont unbounded (lb=-inf)
        lb = np.concatenate([np.zeros(n_comps), -np.inf * np.ones(n_cont)])
        ub = np.inf * np.ones(n_comps + n_cont)
        
        try:
            p_opt = solve_lin(M_g, y_g, sy, n_comps)
            r_g = (y_g - M_g @ p_opt) / sy
            all_r.append(r_g)
            total_params += len(p_opt)
        except Exception:
            # Если линейное решение не сходится, возвращаем большие невязки
            all_r.append(y_g / sy)
            total_params += n_comps + n_cont

    if not all_r:
        return np.array([0.0]), 0
        
    r_total = np.concatenate(all_r)
    return r_total, total_params

def group_dchi2(theta, coeffs, y, g):
    """#GS-37: значимость группы целиком — Δχ² модели «линии + континуум» против «только континуум». По компонентам
    по отдельности нельзя: рентген Bi/Pb/Po/Th K почти коллинеарен, каждая компонента «незначима» при сильной группе."""
    r_full, _ = build_residuals(theta, coeffs, y, [g])
    g0 = dict(g); g0["comps"] = [dict(c, cum=np.zeros_like(c["cum"])) for c in g["comps"]]
    r_cont, _ = build_residuals(theta, coeffs, y, [g0])
    return float(np.sum(r_cont ** 2) - np.sum(r_full ** 2))

AUTOKNOTS = os.environ.get("GS_SHAPE_AUTOKNOTS", "1") == "1"   # узел только у значимой группы спектра
# 100, не 25 (28.09): при 25 у KCl оставался узел 600 (группа 609, Δχ² 34) и уплывал на −20 кэВ — вырождение Tl-208 583
# / Bi-214 609 (разность 26 кэВ) при свободных амплитудах; при 100 синтетика и реальные — без выбросов
DCHI2_MIN = float(os.environ.get("GS_SHAPE_DCHI2", "100"))

# #GS-39 (оператор 29.09 «у тория пшпв 2614 теперь меньше чем нужно»): основа шкалы — НЕ полином файла 4-й степени, а его
# приближение степени BASE_DEG по 60–2700 кэВ. Узлы закрепляют только ПОЛОЖЕНИЯ; наклон между узлами и за последним берётся
# из основы, а у квартики файла образца наклон у 2614 — 0,5495 кэВ/кан против 0,46–0,48 у СпектраЛайн: ПШПВ в кэВ
# переводилась в каналы на 13 % уже измеренной. BASE_DEG=0 — прежнее поведение (полином файла как есть).
# 2, не 1 и не 3 (29.09, scripts\checks\gs39_cal_base_cmp.py): общий χ²/ν образца 1,57 (квартика) → 1,97 (1) / 1,23 (2) /
# 1,85 (3), X-K 6,3 / 28,3 / 7,1 / 26,3; у bg, bgw, kcl степень 2 тоже лучшая; наклон у 2614 — 0,486 кэВ/кан
BASE_DEG = int(os.environ.get("GS_SHAPE_BASE_DEG", "2"))
def base_coeffs(coeffs, n_channels):
    if BASE_DEG <= 0 or BASE_DEG >= len(coeffs) - 1: return list(coeffs)
    ch = np.arange(n_channels, dtype=float); e = C.s.channel_to_energy(ch, coeffs); m = (e > 60) & (e < 2700)
    b = np.polyfit(ch[m], e[m], BASE_DEG)[::-1].tolist()
    print("  основа шкалы: степень %d по 60–2700 кэВ, макс. отличие от полинома файла %.1f кэВ" % (BASE_DEG, np.max(np.abs(C.s.channel_to_energy(ch[m], b) - e[m]))))
    return b

def knots_for(theta, coeffs, y, groups):
    """Узел остаётся, если ближайшая к нему (по центру окна) группа значима: Δχ² ≥ DCHI2_MIN."""
    sig = {g["name"]: group_dchi2(theta, coeffs, y, g) for g in groups}
    keep = []
    for k in KNOTS:
        g = min(groups, key=lambda q: abs(0.5 * (q["lo"] + q["hi"]) - k))
        if sig[g["name"]] >= DCHI2_MIN: keep.append(k)
    return keep, sig

# масштаб кривизны узловой шкалы, кэВ. 10, не 3 (28.09): после AUTOKNOTS свободных узлов нет, а 3 кэВ сглаживали
# настоящий излом (синтетика с изломом 1480→1730→2614: −1,06 кэВ на 1765 при 3, −0,003 при 10)
REG_KEV = float(os.environ.get("GS_SHAPE_REG", "10.0"))
def objective(theta, coeffs, y, active_groups):
    r, _ = build_residuals(theta, coeffs, y, active_groups)
    if KNOTS and NP >= 3 and REG_KEV > 0:
        # #GS-37: штраф на вторую разность узлов — узел без значимой группы рядом (Bi-214 1764 у тория, 960/1130 у
        # KCl) «уплывал» на 25 кэВ; хорошо определённые узлы (σ 0,1–0,5 кэВ) штраф с масштабом 3 кэВ не сдвигает
        r = np.concatenate([r, np.diff(np.asarray(theta, float), 2) / REG_KEV])
    return r

# --- Синтетика ---
def generate_synthetic(coeffs, n_channels, theta_true, rng):
    # Подготовка групп (структура линий та же)
    # Для синтетики нам нужно знать структуру групп, но не зависеть от данных.
    # Используем prepare_groups, но он зависит от FWHM и линий. Это ок.
    active_groups = prepare_groups(coeffs, n_channels)
    
    y_syn = np.zeros(n_channels)
    
    ch_edges_file = C.s.channel_to_energy(np.arange(n_channels + 1) - 0.5, coeffs)
    d_edges_true = get_dE(ch_edges_file, theta_true)
    edges_true = ch_edges_file + d_edges_true
    centers_true = 0.5 * (edges_true[:-1] + edges_true[1:])
    
    for g in active_groups:
        lo, hi = g["lo"], g["hi"]
        idx = g["idx"]

        if len(idx) == 0: continue

        n_ch_g = len(idx)
        comps = [] if g["name"] in SYNTH_DROP else g["comps"]   # проверка AUTOKNOTS: группа без линий
        cdeg = g["cdeg"]

        # Моделируем отклик при theta_true
        # Компоненты
        mu_comps = np.zeros(n_ch_g)
        for comp in comps:
            cum_c = comp["cum"]
            fine_edges = np.arange(0, g1s.FINE_E_MAX + FINE_STEP, FINE_STEP)
            chan_edges_true = edges_true[idx[0] : idx[-1] + 2]
            cum_vals = np.interp(chan_edges_true, fine_edges, cum_c)
            col = np.diff(cum_vals)
            mu_comps += col
            
        # Амплитуды: 1e6 / max(1e-12, sum(col)) * (0.2 + 0.8*k/(K-1))
        K = len(comps)
        for k_idx, comp in enumerate(comps):
            cum_c = comp["cum"]
            fine_edges = np.arange(0, g1s.FINE_E_MAX + FINE_STEP, FINE_STEP)
            chan_edges_true = edges_true[idx[0] : idx[-1] + 2]
            cum_vals = np.interp(chan_edges_true, fine_edges, cum_c)
            col = np.diff(cum_vals)
            
            s_col = np.sum(col)
            amp = 1e6 / max(1e-12, s_col)
            if K > 1:
                factor = 0.2 + 0.8 * k_idx / (K - 1)
            else:
                factor = 1.0
            mu_comps += amp * factor * col
            
        # Континуум: 200 * (1 - 0.3*t)
        centers_g = centers_true[idx]
        t = (centers_g - 0.5 * (lo + hi)) / (hi - lo)
        mu_cont = 200.0 * (1.0 - 0.3 * t)
        
        y_syn[idx] += mu_comps + mu_cont
        
    # Шум Пуассона
    y_poisson = rng.poisson(y_syn)
    return y_poisson

# --- Анализ группы после подгонки ---
def analyze_group(theta, coeffs, y, g, active_groups):
    """Локальный остаток и диагностика."""
    n = len(y)
    ch_edges_file = C.s.channel_to_energy(np.arange(n + 1) - 0.5, coeffs)
    d_edges = get_dE(ch_edges_file, theta)
    edges_true = ch_edges_file + d_edges
    centers_true = 0.5 * (edges_true[:-1] + edges_true[1:])
    
    lo, hi = g["lo"], g["hi"]
    idx = g["idx"]

    if len(idx) == 0: return None

    y_g = y[idx]
    comps = g["comps"]
    cdeg = g["cdeg"]
    n_comps = len(comps)
    n_cont = cdeg + 1
    
    # Построение M для текущего theta
    M_g = np.zeros((len(idx), n_comps + n_cont))
    
    for k, comp in enumerate(comps):
        cum_c = comp["cum"]
        fine_edges = np.arange(0, g1s.FINE_E_MAX + FINE_STEP, FINE_STEP)
        chan_edges_true = edges_true[idx[0] : idx[-1] + 2]
        cum_vals = np.interp(chan_edges_true, fine_edges, cum_c)
        M_g[:, k] = np.diff(cum_vals)
        
    centers_g = centers_true[idx]
    t = (centers_g - 0.5 * (lo + hi)) / (hi - lo)
    for j in range(n_cont):
        M_g[:, n_comps + j] = t ** j
        
    sy = np.sqrt(np.maximum(y_g, 1.0))
    
    # Базовое решение
    lb = np.concatenate([np.zeros(n_comps), -np.inf * np.ones(n_cont)])
    ub = np.inf * np.ones(n_comps + n_cont)
    
    try:
        p_base = solve_lin(M_g, y_g, sy, n_comps)
        r_base = (y_g - M_g @ p_base) / sy
        chi2_min = np.sum(r_base**2)
    except:
        return None
        
    # Скан сдвига s
    # Добавляем сдвиг к d(E) только для этой группы? 
    # ТЗ: "скан добавочного сдвига s от −6 до +6 шаг 0,05 (к d этой группы)"
    # Это означает, что мы меняем границы каналов для этой группы локально.
    
    s_vals = np.arange(-12.0, 12.05, 0.05)   # было ±6: остатки упирались в край
    chi2_s = []
    
    for s in s_vals:
        # Локальные границы со сдвигом
        edges_true_shifted = edges_true[idx[0] : idx[-1] + 2] + s
        centers_g_shifted = 0.5 * (edges_true_shifted[:-1] + edges_true_shifted[1:])
        
        M_g_s = np.zeros((len(idx), n_comps + n_cont))
        
        for k, comp in enumerate(comps):
            cum_c = comp["cum"]
            fine_edges = np.arange(0, g1s.FINE_E_MAX + FINE_STEP, FINE_STEP)
            cum_vals = np.interp(edges_true_shifted, fine_edges, cum_c)
            M_g_s[:, k] = np.diff(cum_vals)
            
        t_s = (centers_g_shifted - 0.5 * (lo + hi)) / (hi - lo)
        for j in range(n_cont):
            M_g_s[:, n_comps + j] = t_s ** j
            
        try:
            p_s = solve_lin(M_g_s, y_g, sy, n_comps)
            r_s = (y_g - M_g_s @ p_s) / sy
            chi2_s.append(np.sum(r_s**2))
        except:
            chi2_s.append(chi2_min * 10) # Penalty
            
    chi2_s = np.array(chi2_s)
    
    # Нахождение минимума и сигмы
    min_idx = np.argmin(chi2_s)
    s_best = s_vals[min_idx]
    chi2_g_min = chi2_s[min_idx]
    
    nu_g = len(idx) - (n_comps + n_cont)
    if nu_g <= 0: nu_g = 1
    
    threshold = chi2_g_min * (1 + 1.0 / nu_g)
    # Полуширина интервала
    # Ищем точки пересечения с порогом слева и справа от минимума
    left_s = s_best
    right_s = s_best
    
    for i in range(min_idx, -1, -1):
        if chi2_s[i] > threshold:
            left_s = s_vals[i+1] if i+1 < len(s_vals) else s_vals[i]
            break
    else:
        left_s = s_vals[0]
        
    for i in range(min_idx, len(s_vals)):
        if chi2_s[i] > threshold:
            right_s = s_vals[i-1] if i-1 >= 0 else s_vals[i]
            break
    else:
        right_s = s_vals[-1]
        
    sigma_s = (right_s - left_s) / 2.0
    if sigma_s < 0.05: sigma_s = 0.05
    
    # E_node
    # Сумма(ампл_c * I_k * E_k) / Сумма(ампл_c * I_k)
    # amplitudes from base fit
    amps = p_base[:n_comps]
    
    num = 0.0
    den = 0.0
    for k, comp in enumerate(comps):
        amp = amps[k]
        if amp < 1e-6: continue # Пропуск незначимых
        for e_line, i_line in comp["lines"]:
            num += amp * i_line * e_line
            den += amp * i_line
            
    E_node = num / den if den > 0 else 0.0
    
    # Значимые компоненты
    Mw = M_g / sy[:, None]   # взвешенная матрица (как в решении)
    sigs = np.sqrt(np.abs(np.diag(np.linalg.pinv(Mw.T @ Mw))) * max(chi2_min / nu_g, 1.0))[:n_comps]
    # Простая оценка сигмы из ковариации линейной задачи сложна без полного J. 
    # Используем приближение: sig ~ amp / sqrt(N_eff)? Нет, лучше взять из res_base если доступно.
    # lsq_linear не возвращает ковариацию напрямую.
    # Для отчета используем отношение amp/sig >= 5. Если сигмы нет, считаем значимыми все с amp > 0.
    
    comps_info = {}
    for k, comp in enumerate(comps):
        comps_info[comp["name"]] = {"amp": float(amps[k]), "sig": float(sigs[k]) if k < len(sigs) else 0.0}
        
    return {
        "E_node": E_node,
        "s_best": s_best,
        "sigma_s": sigma_s,
        "chi2_nu": chi2_g_min / nu_g,
        "comps": comps_info,
        "nu": nu_g
    }

# --- Main ---
def main():
    global FWHM, KNOTS, NP
    FWHM = get_fwhm_func()
    
    tags = sys.argv[1:] if len(sys.argv) > 1 else ["bg", "bgw", "kcl", "sample"]
    
    synth_mode = os.environ.get("GS_SHAPE_SYNTH")
    is_synth = synth_mode is not None
    
    theta_true = None
    if is_synth:
        try:
            theta_true = [float(x) for x in synth_mode.split(",")]
            if len(theta_true) != NP:
                print(f"Ошибка: GS_SHAPE_SYNTH должен содержать {NP} чисел")
                return
        except ValueError:
            print("Ошибка разбора GS_SHAPE_SYNTH")
            return
            
    rng = np.random.default_rng(int(os.environ.get("GS_SHAPE_SEED", "1")))
    
    out_data = {}
    
    for tag in tags:
        print(f"\n== {tag} ==")
        KNOTS, NP = list(ALL_KNOTS), len(ALL_KNOTS) if ALL_KNOTS else DEG + 1   # полный набор узлов для каждого тега
        
        # Чтение данных или генерация
        if is_synth:
            d = C.read(tag) # Чтобы получить coeffs и n
            n_channels = len(d["counts"])
            coeffs = base_coeffs(d["coeffs"], n_channels)
            y = generate_synthetic(coeffs, n_channels, theta_true, rng)
        else:
            d = C.read(tag)
            y = np.asarray(d["counts"], float)
            n_channels = len(y)
            coeffs = base_coeffs(d["coeffs"], n_channels)
            
        # Подготовка групп
        active_groups = prepare_groups(coeffs, n_channels)
        
        if not active_groups:
            print("Нет активных групп для подгонки.")
            continue
            
        # Подгонка; #GS-37: при узлах — второй проход только с узлами значимых групп этого спектра
        try:
            fit = lambda: least_squares(objective, np.zeros(NP), args=(coeffs, y, active_groups), method="trf",
                                        x_scale=[1.0] * NP, diff_step=1e-3, max_nfev=200)
            res = fit()
            if KNOTS and AUTOKNOTS:
                keep, sig = knots_for(res.x, coeffs, y, active_groups)
                print("   Δχ² групп: " + ", ".join("%s %.0f" % kv for kv in sig.items()))
                print("   узлы: %s → %s" % (KNOTS, keep))
                if len(keep) < 2: raise RuntimeError("меньше двух значимых узлов: %s" % keep)
                if keep != KNOTS:
                    KNOTS, NP = keep, len(keep)
                    res = fit()
            theta_fit = res.x
            r_final, n_lin_params = build_residuals(theta_fit, coeffs, y, active_groups)
            
            chi2 = np.sum(r_final**2)
            nu = len(r_final) - n_lin_params - NP   # каналы ГРУПП, не всего спектра
            if nu <= 0: nu = 1
            
            chi2_nu = chi2 / nu
            
            # Ковариация theta
            J = res.jac
            cov_theta = np.linalg.pinv(J.T @ J) * (chi2 / nu)
            
            print(f"DEG={DEG}, θ = {theta_fit} ± {np.sqrt(np.diag(cov_theta))}, χ²/ν = {chi2_nu:.2f} (ν = {nu})")
            
            # Диагностика групп
            groups_info = {}
            for g in active_groups:
                info = analyze_group(theta_fit, coeffs, y, g, active_groups)
                if info:
                    groups_info[g["name"]] = info
                    
            # Вывод по группам
            print(f"   {'группа':<6} {'E_узла':>8} {'d(E)±σ':>12} {'остаток s±σ':>14} {'χ²/ν':>6} {'значимые компоненты'}")
            
            # d(E) и sigma на сетке для вывода
            E_file_grid = np.arange(30, 3001, 10)
            Bg = basis(E_file_grid)
            d_curve = Bg @ theta_fit
            sd_curve = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", Bg, cov_theta, Bg), 0.0))   # σ_d² = gᵀ·cov·g
                
            # Вывод групп
            for g_name, info in groups_info.items():
                # d(E) в центре группы? Или в E_node?
                # ТЗ: "d(E)±σ" в строке группы. Возьмем в E_node.
                g_vec_node = basis(info["E_node"])[0]
                d_node = float(g_vec_node @ theta_fit)
                sd_node = np.sqrt(max(g_vec_node @ cov_theta @ g_vec_node, 0.0))
                
                sig_comps = [k for k, v in info["comps"].items() if v["sig"] > 0 and v["amp"]/v["sig"] >= 5]
                
                print(f"   {g_name:<6} {info['E_node']:>8.1f} {d_node:>7.2f}±{sd_node:.2f} {info['s_best']:>6.2f}±{info['sigma_s']:.2f} {info['chi2_nu']:>6.2f} {sig_comps}")
                
            # d(E) в точках
            points = [80, 240, 350, 610, 910, 1460, 1765, 2615]
            print(f"   d(E) в точках:")
            for p in points:
                idx = np.argmin(np.abs(E_file_grid - p))
                print(f"      {p} кэВ: {d_curve[idx]:.3f} ± {sd_curve[idx]:.3f}")
                
            # Сравнение со своей шкалой
            E_own = C.energy_axis(tag, n_channels)
            E_file_ch = C.s.channel_to_energy(np.arange(n_channels), coeffs)
            diff_own = E_own - E_file_ch
            
            print(f"   Разность (Своя шкала - Файл):")
            for p in points:
                idx = np.argmin(np.abs(E_file_grid - p))
                # Найти ближайший канал к энергии p в файле
                ch_idx = np.argmin(np.abs(E_file_ch - p))
                print(f"      {p} кэВ: {diff_own[ch_idx]:.3f}")
                
            # Требует толкования
            issues = []
            for g_name, info in groups_info.items():
                if abs(info["s_best"]) > 3 * info["sigma_s"] and abs(info["s_best"]) > 1:
                    issues.append(f"Группа {g_name}: сдвиг {info['s_best']:.2f}")
                if info["chi2_nu"] > 3:
                    issues.append(f"Группа {g_name}: χ²/ν = {info['chi2_nu']:.2f}")
                    
            if chi2_nu > 3:
                issues.append(f"Спектр: χ²/ν = {chi2_nu:.2f}")
                
            print(f"   ТРЕБУЕТ ТОЛКОВАНИЯ: {'; '.join(issues) if issues else 'нет'}")
            
            # JSON сборка
            curve_data = {"E": E_file_grid.tolist(), "d": d_curve.tolist(), "sd": sd_curve.tolist()}
            
            groups_json = {}
            for g in active_groups:
                g_info = groups_info.get(g["name"], {})
                # Данные и модель для JSON
                lo, hi = g["lo"], g["hi"]
                ch_edges_file = C.s.channel_to_energy(np.arange(n_channels + 1) - 0.5, coeffs)
                d_edges = get_dE(ch_edges_file, theta_fit)
                edges_true = ch_edges_file + d_edges
                centers_true = 0.5 * (edges_true[:-1] + edges_true[1:])
                
                idx = g["idx"]

                y_g = y[idx].tolist()
                
                # Модель M*p
                comps = g["comps"]
                cdeg = g["cdeg"]
                n_comps = len(comps)
                n_cont = cdeg + 1
                
                if len(idx) > 0:
                    M_g = np.zeros((len(idx), n_comps + n_cont))
                    for k, comp in enumerate(comps):
                        cum_c = comp["cum"]
                        fine_edges = np.arange(0, g1s.FINE_E_MAX + FINE_STEP, FINE_STEP)
                        chan_edges_true = edges_true[idx[0] : idx[-1] + 2]
                        cum_vals = np.interp(chan_edges_true, fine_edges, cum_c)
                        M_g[:, k] = np.diff(cum_vals)
                    centers_g = centers_true[idx]
                    t = (centers_g - 0.5 * (lo + hi)) / (hi - lo)
                    for j in range(n_cont):
                        M_g[:, n_comps + j] = t ** j
                        
                    sy = np.sqrt(np.maximum(y[idx], 1.0))
                    lb = np.concatenate([np.zeros(n_comps), -np.inf * np.ones(n_cont)])
                    ub = np.inf * np.ones(n_comps + n_cont)
                    try:
                        model_g = (M_g @ solve_lin(M_g, y[idx], sy, n_comps)).tolist()
                    except:
                        model_g = [0.0]*len(idx)
                else:
                    model_g = []
                    
                groups_json[g["name"]] = {
                    "lo": lo, "hi": hi,
                    "E_node": g_info.get("E_node", 0),
                    "d": float(get_dE(g_info.get("E_node", 0), theta_fit)[0]) if g_info else 0,
                    "s": g_info.get("s_best", 0),
                    "sigma_s": g_info.get("sigma_s", 0),
                    "chi2_nu": g_info.get("chi2_nu", 0),
                    "comps": g_info.get("comps", {}),
                    "e": centers_true[idx].tolist(),
                    "data": y_g,
                    "model": model_g
                }
                
            # поля, которые подгонки и экспорт берут из своей шкалы (как в cal_own.json): refs — остаточный сдвиг s
            # значимых групп против порога #CAL-0 (0,25·ПШПВ); coeffs — E(канал) степени 4 для страницы
            chx = np.arange(n_channels, dtype=float); ef = C.s.channel_to_energy(chx, coeffs)
            e_true = ef + get_dE(ef, theta_fit)
            c_repr = np.polynomial.polynomial.polyfit(chx, e_true, C.REPR_ORDER)
            knot_groups = {min(active_groups, key=lambda q: abs(0.5 * (q["lo"] + q["hi"]) - k))["name"] for k in KNOTS}   # как knots_for
            refs = [{"group": gn, "E_lib": float(gi["E_node"]), "mu_file": float(gi["E_node"] - get_dE(gi["E_node"], theta_fit)[0]),
                     "resid_keV": float(gi["s_best"]), "limit_keV": float(0.25 * FWHM(gi["E_node"])),
                     "ok": bool(abs(gi["s_best"]) <= 0.25 * FWHM(gi["E_node"]))}
                    for gn, gi in groups_info.items() if gi["E_node"] > 0 and gn in knot_groups]
            out_data[tag] = {
                "refs": refs, "coeffs": c_repr.tolist(), "corr_nodes": [[k, float(t)] for k, t in zip(KNOTS, theta_fit)] if KNOTS else [],
                "repr_max_dev_keV": float(np.max(np.abs(np.polynomial.polynomial.polyval(chx, c_repr) - e_true))),
                "deg": DEG,
                "knots": KNOTS,   # пусто — полином DEG; иначе d(E) кусочно-линейная по этим узлам (кэВ шкалы файла)
                "theta": theta_fit.tolist(),
                "cov": cov_theta.tolist(),
                "chi2_nu": chi2_nu,
                "nu": nu,
                "file_coeffs": coeffs.tolist() if hasattr(coeffs, 'tolist') else list(coeffs),
                "curve": curve_data,
                "groups": groups_json
            }
            
        except Exception as e:
            print(f"Ошибка подгонки для {tag}: {e}")
            import traceback
            traceback.print_exc()

    # Сохранение JSON
    json_path = OUT_JSON
    if is_synth:
        base, ext = os.path.splitext(OUT_JSON)
        json_path = f"{base}_synth{ext}"
        
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(out_data, f, ensure_ascii=False, indent=1)
        
    if is_synth:
        print(f"\nСИНТЕТИКА: θ_true = {theta_true}, θ_fit = {theta_fit}")
        # d_fit - d_true в точках
        points = [80, 240, 350, 610, 910, 1460, 1765, 2615]
        print(f"   {'E':>6} | {'d_fit-d_true':>12} | {'σ_d':>8} | {'Отн. σ':>8}")
        for p in points:
            idx = np.argmin(np.abs(E_file_grid - p))
            d_true = np.interp(p, ALL_KNOTS, theta_true) if ALL_KNOTS else get_dE(np.array([p]), theta_true)[0]
            d_diff = d_curve[idx] - d_true
            sd = sd_curve[idx]
            rel = d_diff / sd if sd > 0 else 0
            print(f"   {p:>6} | {d_diff:>12.4f} | {sd:>8.4f} | {rel:>8.2f}")

if __name__ == "__main__":
    main()
