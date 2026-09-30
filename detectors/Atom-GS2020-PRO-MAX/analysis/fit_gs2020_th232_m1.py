# -*- coding: utf-8 -*-
r"""Метод 1 (шаблоны полного распада 8 звеньев) для КИ Th-232 в Маринелли на GS2020. Подгонка — импорт
mix_unfold_core.unfold (критерий A2 + E1, фон с фиксированным k). Спека: scripts\specs\SPEC-fit_gs2020_th232_m1.md"""

import sys
import os
import json
import math
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

os.environ.setdefault("SPECTRAVIBE_ROOT", r"D:\cloud-folder\Дозиметрия\ИИ\1 Скилы\0_Work\gamma-spectrum-analysis")
sys.path.insert(0, r"D:\repos-folder\repos\geant4-detector-models\common\py")
import becqmoni as bm

sys.path.insert(0, r"D:\repos-folder\repos\geant4-detector-models\detectors\Gamma-1S\analysis")
import mix_unfold_core as muc
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "cal"))
import gs2020_calib as CAL
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gs2020_extra_components as gx   # #GS-42: IB β⁻-звеньев в шаблоны М1 (β уже в распаде иона); GS_EXTRA=0 — без него
CAL_OWN = json.load(open(CAL.OUT_JSON, encoding="utf-8")) if os.path.exists(CAL.OUT_JSON) else None
print("ШКАЛА: %s" % CAL.OUT_JSON)   # #GS-28: какой файл своих шкал взят (база / гибрид) — артефакт в логе каждой подгонки
muc.g1s.FINE_E_MAX = 4800.0   # #SUM-1 (оператор 25.09): сетка размытия донора кончалась на 3300 кэВ и срезала сумм-пики Tl-208 3197/3475

REF = r"D:\cloud-folder\Дозиметрия\Спектры\Atom GS2020 PRO MAX\Референсы"
XML_SAMPLE = os.path.join(REF, "Калибровка Th-232 (без вычета фона).xml")
XML_BG = os.path.join(REF, "Фон лаба S31_18.xml")
BG_TAG = "bg"
# GS_BG_WATER=1 (оператор 26.09: «КИ с торием считать именно с фоном с водой»): фон — Маринелли 1 л с дист. водой,
# своя шкала по его реперам (тег "bgw" в cal/gs2020_calib.py); модельное ослабление GS_BG_T при этом запрещено.
BG_WATER = os.environ.get("GS_BG_WATER") == "1"
if BG_WATER:
    if os.environ.get("GS_BG_T"):
        raise SystemExit("ОТКАЗ: GS_BG_WATER=1 и GS_BG_T несовместимы — фон с водой уже ослаблен сосудом физически")
    XML_BG, BG_TAG = CAL.BKG_WATER_XML, "bgw"
OUT = os.environ.get("GS_OUT", r"C:\g4work\gs2020\run_marinelli\out_v5")  # #GS-3: GS_OUT — папка шаблонов (MgO 0,8 → out_v5_mgo08)
LO, HI = 150.0, 3600.0   # #SUM-1: верх окна 3600 захватывает сумм-пики Tl-208 3197 и 3475 кэВ.
# #XR-1 (25.09): попытка опустить LO до 25 кэВ ОТВЕРГНУТА — модель не описывает зону <150 кэВ:
# Ra224 уходит на 9303 Бк (×9,7 паспорта), Rn220 2533±1173 (ошибка ~= значению), χ²/ν подгонки 15,6
# (было 4,2 при LO=150). Рентген K/L при этом остаётся в библиотеке М2 и в шаблонах М1 (fluo=1) —
# физика в модели ЕСТЬ независимо от границы подгонки; ниже 150 кэВ модель просто не оценивается по χ².
PASSPORT_BQ = 910.0 * 1.052
PASSPORT_UNC = 0.06
CHAIN = [("Th232", 1.0), ("Ac228", 1.0), ("Th228", 1.0), ("Ra224", 1.0), ("Rn220", 1.0), ("Pb212", 1.0), ("Bi212", 1.0), ("Tl208", 0.3594)]
# #GS-42: Ra-228 — звено цепочки Th-232 в равновесии (BR 1), шаблона М1 у него нет (γ ниже окна); его β/IB входят
# только в столбец цепочки (главный результат), в поузловую подгонку — нет (нет амплитуды).
EXTRA_CHAIN_ONLY = [("Ra228", 1.0)]

PEAK_TABLE = [
    (238.197, 238.632, 24.421),
    (588.042, 583.187, 41.913),
    (915.581, 911.204, 57.092),
    (969.137, 964.766, 59.229),
    (973.341, 968.971, 59.395),
    (1453.398, 1460.822, 68.502),
    (2618.317, 2614.511, 121.205)
]


# #CAL-1 (оператор 25.09: «калибровка фона и образца не сделана независимо»): у фона СВОЯ шкала — реперы фона
# по его же спектру, центроид гауссиана+линейная подложка (scripts\cal\gs2020_cal_fit_both.py). 351/583/609/1764
# в фоне слиты в мультиплеты и в таблицу не входят.
BG_PEAK_TABLE = [(236.29, 238.632, None), (1454.27, 1460.822, None), (2609.97, 2614.511, None)]


def true_energy(e_file, table=None):
    table = PEAK_TABLE if table is None else table
    xs = [p[0] for p in table]
    ds = [p[1] - p[0] for p in table]
    return e_file + np.interp(e_file, xs, ds)


# GS_CAL_SL=<степень> — шкала образца по каналам пиков СпектраЛайн (оператор 26.09, скрин «Параметры пиков»): канал → E_lib
CAL_SL_DEG = 0 if CAL.SHAPE else int(os.environ.get("GS_CAL_SL", "0"))   # #GS-37: при шкале по форме SL не применяется
CAL_SL = CAL_SL_DEG > 0
SL_PEAKS = [(543.253, 238.632), (1297.795, 583.187), (1986.459, 911.204), (2098.076, 964.766),
            (2106.831, 968.971), (3129.695, 1460.822), (5537.873, 2614.511)]
# GS_CAL_SUM="6765:3187" — репер выше 2614: горб суммирования Tl-208 (канал — центр в измерении, энергия — центр
# того же горба в шаблоне Geant4 цепочки; оба не зависят от шкалы). Шкала файла ставит этот канал на ~3400 кэВ.
CAL_SUM = [tuple(float(x) for x in p.split(":")) for p in os.environ.get("GS_CAL_SUM", "").split(",") if p]


class Spec:
    def __init__(self, sp, tag="sample"):
        self.counts = list(sp.n)
        self.live_time = float(sp.live)
        self.n_channels = len(sp.n)
        self._e = CAL.energy_axis(tag, len(sp.n))   # #CAL-2: своя шкала КАЖДОГО спектра по его реперам (cal/gs2020_calib.py)
        if tag == "sample" and CAL_SL:              # #CAL-0: центроиды из таблицы пиков СпектраЛайн (README-референсы.md)
            pts = SL_PEAKS + CAL_SUM
            p = np.polyfit([q[0] for q in pts], [q[1] for q in pts], CAL_SL_DEG)
            self._e = np.polyval(p, np.arange(len(sp.n), dtype=float))
            print("#CAL-0 шкала образца по СпектраЛайн, степень %d: " % CAL_SL_DEG + "; ".join(
                "%.1f %+.2f" % (lib, np.polyval(p, c) - lib) for c, lib in pts)
                + "; E(кан 342)=%.1f E(кан 7016)=%.1f" % (np.polyval(p, 342.0), np.polyval(p, 7016.0)))
        if tag == "sample" and CAL.SHAPE and CAL_SUM:   # #GS-37: выше 2614 у шкалы по форме узла нет (за краем — константа),
            ch = np.arange(len(sp.n), dtype=float)       # шкала файла там ошибается на +230 кэВ → узел горба суммы, линейно от 2614
            c0 = float(np.interp(2614.511, self._e, ch))
            nodes = [(c0, 0.0)] + [(c, E - float(np.interp(c, ch, self._e))) for c, E in CAL_SUM]
            # парабола с нулевым наклоном в c0 (C¹): линейная стыковка давала излом наклона ~40 % на самом пике 2614
            # (форма 2614: χ²/ν 0,68 → 2,17); узел горба (последний) проходится точно
            a = nodes[-1][1] / (nodes[-1][0] - c0) ** 2
            self._e = self._e + np.where(ch <= c0, 0.0, a * (ch - c0) ** 2)
            print("#GS-37 шкала образца по форме + горб суммы: " + "; ".join("кан %.0f %+.1f кэВ" % q for q in nodes)
                  + "; E(кан 7016)=%.1f" % self._e[min(7016, len(sp.n) - 1)])

    def channel_to_energy(self, c):
        return float(self._e[int(c)])


# GS_BG_T="f,rhot": T(E) = 1 - f*(1 - exp(-mu/rho(E)*rhot)); mu/rho — вода NIST XCOM (см²/г) × 0,96 (эпоксид, Z/A)
BG_T = [float(x) for x in os.environ["GS_BG_T"].split(",")] if os.environ.get("GS_BG_T") else None
MU_E = [150, 200, 300, 400, 500, 600, 800, 1000, 1250, 1500, 2000, 3000]
MU_W = [0.1505, 0.1370, 0.1186, 0.1061, 0.0969, 0.0896, 0.0786, 0.0707, 0.0632, 0.0575, 0.0494, 0.0397]


def bg_transmission(e):
    mu = 0.96 * float(np.exp(np.interp(np.log(max(e, 150.0)), np.log(MU_E), np.log(MU_W))))
    return 1.0 - BG_T[0] * (1.0 - math.exp(-mu * BG_T[1]))


# Оператор 26.09 «с калибровкой по ПШПВ что-то не так»: фоновая K-40 (68,5 кэВ — уже закона ~79) шла в точки
# размытия, хотя страница заявляет её исключённой. GS_FWHM_K40=1 — прежнее поведение (для сверки); GS_BLUR — множитель ширины.
FWHM_K40 = os.environ.get("GS_FWHM_K40") == "1"
BLUR = float(os.environ.get("GS_BLUR", "1.0"))
# GS_TAIL — параметр левого хвоста ядра (донор Гамма-1С: TAIL_T=0.75, mix_unfold_g1s.py:24); не задан — донорский
TAIL = float(os.environ["GS_TAIL"]) if os.environ.get("GS_TAIL") else None


# ПШПВ по таблице СпектраЛайн (оператор 26.09 «точнее сделал»; README-референсы.md, скрин рядом) — по умолчанию.
# GS_FWHM_OLD=1 — прежние ПШПВ из PEAK_TABLE (таблица прибора 25.09; 2614: 121,205 против 107,420).
FWHM_SL = {238.632: 24.073, 583.187: 40.624, 911.204: 57.351, 964.766: 59.686, 968.971: 59.867,
           1460.822: 70.963, 2614.511: 107.420}
FWHM_OLD = os.environ.get("GS_FWHM_OLD") == "1"
# GS_FWHM_SCALE="238.632:1.05,583.187:1.12" — поточечные множители ПШПВ (#SHAPE-1: кривая из минимумов невязки формы)
FWHM_SCALE = {float(k): float(v) for k, v in (p.split(":") for p in os.environ.get("GS_FWHM_SCALE", "").split(",") if p)}


# GS_FWHM_CFW=1 — гладкий закон СпектраЛайн из Референсы\Calibr.cfw, секция [Calibration]: ПШПВ = Σ cᵢ·(√E)ⁱ
FWHM_CFW = os.environ.get("GS_FWHM_CFW") == "1"


def cfw_law():
    txt = open(os.path.join(REF, "Calibr.cfw"), encoding="cp1251").read()
    c = [float(x) for x in txt.split("[Calibration]")[1].split("Coeff=")[1].split()[0].split(",")]
    return lambda E: sum(ci * math.sqrt(E) ** i for i, ci in enumerate(c))


def write_fwhm_csv(path):
    with open(path, "w", encoding="utf-8") as f:
        f.write("E_keV,fwhm_keV\n")
        if FWHM_CFW:
            law = cfw_law()
            for E in np.arange(30.0, 4001.0, 10.0):
                f.write(f"{E},{law(E):.4f}\n")
            return
        for _, e_lib, fwhm_old in PEAK_TABLE:
            fwhm = (fwhm_old if FWHM_OLD else FWHM_SL[e_lib]) * FWHM_SCALE.get(e_lib, 1.0)
            if e_lib == 1460.822 and not FWHM_K40:
                continue
            f.write(f"{e_lib},{fwhm}\n")


SHAPE_PEAKS = [238.632, 338.32, 583.187, 727.33, 911.204, 968.971, 1592.5, 2614.511]


def shape_residual(e, net, var, model, fwhm_f):
    """#SHAPE-1: невязка ФОРМЫ пиков. В окне ±1,5 ПШПВ модель подгоняется к нетто как a·модель + b + c·(E−E0)
    (площадь и подложка свободны) — остаётся только форма. Возвращает ({E0: χ²/ν окна}, общий χ²/ν окон)."""
    out, X2, NU = {}, 0.0, 0
    for E0 in SHAPE_PEAKS:
        m = np.abs(e - E0) < 1.5 * fwhm_f(E0)
        w = 1.0 / np.sqrt(np.maximum(var[m], 1.0))
        A = np.vstack([model[m], np.ones(int(m.sum())), e[m] - E0]).T
        p = np.linalg.lstsq(A * w[:, None], net[m] * w, rcond=None)[0]
        x2 = float((((net[m] - A @ p) * w) ** 2).sum()); nu = int(m.sum()) - 3
        out[E0] = x2 / nu; X2 += x2; NU += nu
    return out, X2 / NU


def main():
    if not os.path.exists(XML_SAMPLE):
        raise SystemExit(f"ОТКАЗ: нет файла {XML_SAMPLE}")
    if not os.path.exists(XML_BG):
        raise SystemExit(f"ОТКАЗ: нет файла {XML_BG}")

    s = Spec(bm.read(XML_SAMPLE)[0])
    b = Spec(bm.read(XML_BG)[0], BG_TAG)
    for tg in ("sample", BG_TAG):
        print("#CAL-2 шкала %s — своя: " % tg + "; ".join("%.1f→%.1f (без него %+.2f)" % (q["mu_file"], q["E_lib"], q["resid_keV"]) for q in CAL_OWN[tg]["refs"]))
    if BG_T:  # экранирование фона сосудом: фон снят без сосуда (оператор 25.09 «делай»)
        b.counts = [c * bg_transmission(b.channel_to_energy(i)) for i, c in enumerate(b.counts)]
        print("ФОН ОСЛАБЛЕН сосудом: f=%g, ρt=%g г/см²; T(238)=%.3f T(583)=%.3f T(1461)=%.3f T(2615)=%.3f" % (
            BG_T[0], BG_T[1], *[bg_transmission(x) for x in (238.6, 583.2, 1460.8, 2614.5)]))

    print("Проверка калибровки:")
    for tg in ("sample", BG_TAG):
        for q in CAL_OWN[tg]["refs"]:
            print(f"  {tg}: File: {q['mu_file']:.3f}, Lib: {q['E_lib']:.3f}, Resid (leave-one-out): {q['resid_keV']:+.3f}")

    os.makedirs(OUT, exist_ok=True)
    fwhm_csv = os.path.join(OUT, "fwhm_points_gs2020.csv")
    write_fwhm_csv(fwhm_csv)

    templates = [(k, os.path.join(OUT, "mix_%s_npsmoff.csv" % k)) for k, _ in CHAIN]
    for _, path in templates:
        if not os.path.exists(path):
            raise SystemExit(f"ОТКАЗ: нет шаблона {path}")

    bg_e = np.array([b.channel_to_energy(i) for i in range(b.n_channels)])

    # #GS-42: IB β⁻-звеньев — в шаблон звена (на его распад) и в шаблон цепочки (× BR); β не добавляется (уже в распаде иона)
    comps = {k: gx.load(k, OUT, beta=False) for k, _ in CHAIN + EXTRA_CHAIN_ONLY} if gx.ENABLED else {}
    if not gx.ENABLED:
        print("ФИЗИКА (#GS-42): GS_EXTRA=0 — IB выключен (прежнее поведение)")
    plan_links = {p: [(comps[k], 1.0)] for k, p in templates if k in comps and comps[k]["ib"]}
    with gx.m1_ib(muc.g1s, plan_links):
        r = muc.unfold(s, b, templates, fwhm_csv, lo=LO, hi=HI, bg_energy_of_ch=bg_e, verbose=True, blur=BLUR, tail=TAIL)

    names = r["names"]
    activities_A2 = r["activities"]
    sd_counts = r["sd"]
    sel = r["sel"]
    net = r["net"]
    var = r["var"]
    e = r["e"]
    coef = r["coef"]
    cols = r["cols"]

    activities_E1 = r["e1"]["activities"]
    sd_counts_E1 = r["e1"]["sd"]

    sd_A2_Bq = sd_counts / s.live_time
    sd_E1_Bq = sd_counts_E1 / s.live_time

    model = (r["cols"] * r["coef"][:, None]).sum(axis=0)

    chi2 = float(np.sum((model[sel] - net[sel]) ** 2 / var[sel]))
    ndof = int(sel.sum()) - len(names)

    print("\nРезультаты подгонки:")
    print(f"{'Звено':<10} {'A2 Bq':>10} {'±sd Bq':>10} {'E1 Bq':>10} {'Chain-Eq A2/br':>15} {'Ratio Passport':>15}")
    print("-" * 75)

    chain_eq = []
    weights = []
    values = []
    used_links = []

    for i, (name, br) in enumerate(CHAIN):
        a2 = activities_A2[i]
        sd_bq = sd_A2_Bq[i]
        e1 = activities_E1[i]
        
        if np.isnan(sd_bq):
            sd_str = "nan"
        else:
            sd_str = f"{sd_bq:.3f}"

        chain_val = a2 / br
        ratio = chain_val / PASSPORT_BQ
        
        print(f"{name:<10} {a2:>10.3f} {sd_str:>10} {e1:>10.3f} {chain_val:>15.3f} {ratio:>15.4f}")

        chain_eq.append(chain_val)
        
        if not np.isnan(sd_bq) and sd_bq / a2 < 0.2:
            var_chain = (sd_bq / br) ** 2
            weights.append(1.0 / var_chain)
            values.append(chain_val)
            used_links.append(name)

    if weights:
        w_sum = sum(weights)
        weighted_mean = sum(w * v for w, v in zip(weights, values)) / w_sum
        weighted_sd = math.sqrt(1.0 / w_sum)
    else:
        weighted_mean = 0.0
        weighted_sd = 0.0

    print(f"\nВзвешенное среднее цепочки: {weighted_mean:.3f} ± {weighted_sd:.3f} Bq")
    print(f"Отношение к паспорту: {weighted_mean / PASSPORT_BQ:.4f}")
    print(f"Использованы звенья: {', '.join(used_links)}")

    print(f"\nChi2: {chi2:.3f}, ndof: {ndof}, Chi2/ndof: {chi2/ndof if ndof > 0 else 'N/A':.3f}")

    issues = []
    
    # (a) any PEAK_TABLE residual after correction |>0.5| keV
    for tg in ("sample", BG_TAG):
        for q in CAL_OWN[tg]["refs"]:
            if not q["ok"]:
                issues.append(f"{tg}: невязка калибровки {q['resid_keV']:+.2f} кэВ у репера {q['E_lib']} выше 0,25·ПШПВ")

    # (b) sign of (E_library − E_file) differs between rows
    signs = []
    for e_file, e_lib, _ in PEAK_TABLE:
        diff = e_lib - e_file
        if diff > 0:
            signs.append(1)
        elif diff < 0:
            signs.append(-1)
        else:
            signs.append(0)
    
    non_zero_signs = [s for s in signs if s != 0]
    if len(set(non_zero_signs)) > 1:
        rows_with_diff = []
        for idx, (e_file, e_lib, _) in enumerate(PEAK_TABLE):
            diff = e_lib - e_file
            if diff != 0:
                rows_with_diff.append(f"Row {idx} ({e_lib})")
        issues.append(f"Разные знаки разности E_lib - E_file в строках: {', '.join(rows_with_diff)}")

    # (c) any link with A2 = 0 (clamped)
    for i, name in enumerate(names):
        if activities_A2[i] == 0:
            issues.append(f"Звено {name} имеет активность 0 (зажато)")

    # (d) chi2/ndof > 2
    if ndof > 0 and chi2 / ndof > 2:
        issues.append(f"Chi2/ndof = {chi2/ndof:.3f} > 2")

    # (e) weighted chain activity differs from passport by more than 2 combined sigma
    if weighted_sd > 0:
        combined_sd = math.sqrt(weighted_sd ** 2 + (PASSPORT_BQ * PASSPORT_UNC) ** 2)
        diff = abs(weighted_mean - PASSPORT_BQ)
        if diff > 2 * combined_sd:
            issues.append(f"Взвешенная активность отличается от паспорта более чем на 2 сигмы")

    print("\nТРЕБУЕТ ТОЛКОВАНИЯ:")
    if issues:
        for issue in issues:
            print(issue)
    else:
        print("пусто")

    model_by_link = {}
    for i, name in enumerate(names):
        link_model = r["cols"][i] * r["coef"][i]
        model_by_link[name] = [round(float(x), 3) for x in link_model]

    # ГЛАВНЫЙ результат (оператор 25.09: «цепочка в равновесии»): один столбец — шаблон всей цепочки,
    # сумма сырых гистограмм звеньев при одинаковом n/BR (merge_templates_gs2020.py). Поузловая — диагностика.
    chain_tpl = os.path.join(OUT, "mix_Th232chain_npsmoff.csv")
    ib_items = [(k, br, comps[k]) for k, br in CHAIN + EXTRA_CHAIN_ONLY if k in comps and comps[k]["ib"]]
    with gx.m1_ib(muc.g1s, {chain_tpl: [(c, br) for _, br, c in ib_items]} if ib_items else {}):
        rc = muc.unfold(s, b, [("Th232chain", chain_tpl)], fwhm_csv,
                        lo=LO, hi=HI, bg_energy_of_ch=bg_e, verbose=False, blur=BLUR, tail=TAIL)
    mc = rc["cols"][0] * rc["coef"][0]
    chi2c = float(np.sum((mc[rc["sel"]] - rc["net"][rc["sel"]]) ** 2 / rc["var"][rc["sel"]]))
    ndofc = int(rc["sel"].sum()) - 1
    birgec = math.sqrt(max(chi2c / ndofc, 1.0))
    shp, shp_all = shape_residual(np.asarray(rc["e"]), np.asarray(rc["net"]), np.asarray(rc["var"]), mc,
                                  muc.g1s.make_fwhm(fwhm_csv))
    print("ФОРМА ПИКОВ (#SHAPE-1, χ²/ν в окнах ±1,5 ПШПВ): всего %.3f; " % shp_all
          + "; ".join("%.1f → %.2f" % (k, v) for k, v in shp.items()))
    Ac, dAc = float(rc["activities"][0]), float(rc["sd"][0] / s.live_time)
    Ae1 = float(rc["e1"]["activities"][0])
    print("\nЦЕПОЧКА В РАВНОВЕСИИ (метод 1): A2 %.1f ± %.1f (стат) ± %.1f (Бирге %.2f) Бк; E1 %.1f Бк; паспорт %.1f ± %.1f; "
          "отношение %.4f; χ²/ν %.3f" % (Ac, dAc, dAc * birgec, birgec, Ae1, PASSPORT_BQ, PASSPORT_BQ * PASSPORT_UNC,
                                           Ac / PASSPORT_BQ, chi2c / ndofc))
    extra = None
    if ib_items:   # #GS-42: IB звеньев в модели цепочки (амплитуда цепочки), окно и полосы <150, 150–400, >400 кэВ
        cols_ib = {k: gx.fold(c, muc.g1s.broaden, rc["ch_edges"], lambda E: BLUR * rc["fwhm"](E))["ib"][0] for k, _, c in ib_items}
        extra = {"amplitude": "цепочка, метод 1", "bands_keV": "lt150: E<150; 150_400: 150≤E<400; gt400: E≥400 (вся шкала)",
                 "var_excess_max": gx.m1_var_excess([(br, cols_ib[k], c["ib"]["n_eff"]) for k, br, c in ib_items],
                                                    rc["n_events"][0], rc["coef"][0], rc["var"], rc["sel"])}
        for k, br, c in ib_items:
            # ib_col: поканальный столбец IB звена k на сетке rc["e"]/model (та же величина, что фолдится в band_counts
            # ниже) — для слоя IB на странице (GS-42, п.2); БЕЗ округления, чтобы Σ(ib_col[sel]) == band_counts["window"].
            ib_col = br * cols_ib[k] * float(rc["coef"][0])
            extra[k] = {"br": br, "ib": gx.band_counts(ib_col, rc["e"], rc["sel"]),
                        "ib_col": [float(x) for x in ib_col],
                        "ib_file": c["ib"]["file"], "Y": c["ib"]["Y"], "drawn": c["ib"]["drawn"]}
        tot = sum(extra[k]["ib"]["window"] for k, _, _ in ib_items)
        print("IB (#GS-42, метод 1, цепочка): в окне %.0f отсчётов (%.2f %% модели в окне): %s; приближение дисперсии донора ≤ %.1e"
              % (tot, 100 * tot / mc[rc["sel"]].sum(), ", ".join("%s %.0f" % (k, extra[k]["ib"]["window"]) for k, _, _ in ib_items),
                 extra["var_excess_max"]))
    # Диагностика GS_PEAKWIN=k: та же цепочка, но подгонка только в окнах ±k·ПШПВ вокруг 20 линий библиотеки
    # (те же окна, что у метода 2) — сравнение методов без континуума между пиками.
    PW = float(os.environ.get("GS_PEAKWIN", "0"))
    if PW > 0:
        import yaml
        cfg = r"D:\repos-folder\repos\geant4-detector-models\detectors\Gamma-1S\web-th232\configs\th232.yaml"
        with open(cfg, encoding="utf-8") as fh:
            lines = [float(l["e_kev"]) for l in yaml.safe_load(fh)["library"]["lines"]]
        fw = muc.g1s.make_fwhm(fwhm_csv)
        w = np.zeros(len(rc["e"]), dtype=bool)
        for E in lines:
            w |= np.abs(rc["e"] - E) <= PW * fw(E)
        selw = rc["sel"] & w
        cw, sw, _ = muc.cb.fit_A2(rc["cols"], rc["counts"], rc["bg_scaled"], s.live_time / b.live_time, rc["n_events"], selw)
        print("ОКНА ПИКОВ ±%.2f·ПШПВ (каналов %d): цепочка М1 %.1f ± %.1f (стат) Бк, отношение к паспорту %.4f"
              % (PW, int(selw.sum()), cw[0] / s.live_time, sw[0] / s.live_time, cw[0] / s.live_time / PASSPORT_BQ))
    # вклад звена в модель цепочки: столбец звена (на распад) · BR · A цепочки · T; Σ по звеньям = модель цепочки
    stack = {nm: r["cols"][i] * br * Ac * s.live_time for i, (nm, br) in enumerate(CHAIN)}
    lost = float(np.abs(sum(stack.values())[rc["sel"]] - mc[rc["sel"]]).sum() / mc[rc["sel"]].sum())
    print("разложение цепочки по звеньям: относительная невязка Σ вкладов к модели %.2e" % lost)
    if lost > 1e-3:
        raise SystemExit("ОТКАЗ: Σ вкладов звеньев ≠ модели цепочки (%.3e)" % lost)
    chain_fit = {"A_Bq": Ac, "dA_stat_Bq": dAc, "dA_Bq": dAc * birgec, "birge": birgec, "E1_Bq": Ae1,
                 "chi2": chi2c, "ndof": ndofc, "model": [round(float(x), 3) for x in mc],
                 "stack": {k: [round(float(x), 3) for x in v] for k, v in stack.items()}}

    result_json = {
        "names": names,
        "activities_A2": [round(float(x), 3) for x in activities_A2],
        "sd_A2_Bq": [round(float(x), 3) for x in sd_A2_Bq],
        "activities_E1": [round(float(x), 3) for x in activities_E1],
        "sd_E1_Bq": [round(float(x), 3) for x in sd_E1_Bq],
        "chain_eq": [round(float(x), 3) for x in chain_eq],
        "weighted_chain_Bq": round(weighted_mean, 3),
        "weighted_chain_sd": round(weighted_sd, 3),
        "passport_Bq": PASSPORT_BQ,
        "chi2": round(chi2, 3),
        "ndof": ndof,
        "live_s": s.live_time,
        "k_bg": s.live_time / b.live_time,
        "lo": LO,
        "hi": HI,
        "e": [round(float(x), 3) for x in e],
        "net": [round(float(x), 3) for x in net],
        "model": [round(float(x), 3) for x in model],
        "model_by_link": model_by_link,
        "sel": [bool(x) for x in sel],
        "chain": chain_fit,
        "shape": {"all": shp_all, "peaks": {str(k): v for k, v in shp.items()}}
    }
    if extra:
        result_json["extra_components"] = extra

    json_path = os.path.join(OUT, "fit_m1%s%s%s.json" % ("_pw%g" % PW if PW > 0 else "", "_bgT%g_%g" % tuple(BG_T) if BG_T else ("_bgw" if BG_WATER else ""),
                                  ("_blur%g" % BLUR if BLUR != 1.0 else "") + ("_k40" if FWHM_K40 else "") + ("_tail%g" % TAIL if TAIL is not None else "") + ("_fwold" if FWHM_OLD else "") + ("_cfw" if FWHM_CFW else "") + ("_fwscale" if FWHM_SCALE else "") + ("_calsl%d" % CAL_SL_DEG if CAL_SL else "") + ("_calshape" if CAL.SHAPE else "") + ("_calsum" if CAL_SUM else "")))
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(result_json, f, ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
