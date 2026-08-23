# -*- coding: utf-8 -*-
"""НУКЛИДНОЕ шаблонное разложение фона (канон geant4-spectrum-pipeline, метод 1).

Отличие от fit_room_field.py: там ТРИ шаблона на цепочки K/Ra/Th с допущением
векового равновесия внутри каждой; здесь ВОСЕМЬ шаблонов на ЗВЕНЬЯ -
K-40, Ra-226, Pb-214, Bi-214, Pb-212(+Ra-224), Ac-228, Bi-212, Tl-208, -
и равновесие не постулируется, а ПРОВЕРЯЕТСЯ отношениями подобранных амплитуд.

Шаблоны: wallfield.exe nuc=<N> -> поле звена при 1 Бк/кг родителя ->
rc_curves (та же геометрия, что цепочечные базисы) -> bg_cyl_field_nuc_<N>.csv.
Нормировка - то же тождество Ф=4N/S, что в fit_room_field.fit().

ГЛАВНАЯ ЛОВУШКА (урок цепочечного фита, зафиксирован в fit_room_field.py):
базисы могут быть МАТЕМАТИЧЕСКИ НЕРАЗЛИЧИМЫ - комптоновский континуум
"забывает" родительский нуклид, и NNLS раскладывает вес между коррелирующими
столбцами произвольно. Для 8 звеньев риск ВЫШЕ, чем для 3 цепочек. Поэтому
диагностика вырожденности идёт ДО подгонки и печатается всегда: матрица
корреляций столбцов, число обусловленности, после подгонки - ошибки амплитуд
из ковариации. Амплитуда с ошибкой больше себя самой - не результат, а
артефакт вырожденности, и помечается явно.

Запуск: python fit_nuclides.py
"""
import glob
import json
import math
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "common", "py"))
import numpy as np
import paths
import rcspec
import read_rcxml
import fit_lines as fl

BUILD = str(paths.build("RadiaCode-103"))
RESULTS = os.path.abspath(os.path.join(HERE, "..", "results"))
# ДЕФОЛТ = bare (без сосуда): "Фон 7 дней без домика.xml" снят БЕЗ сосуда
# (оператор, 21.08) — измерению соответствует bare-геометрия, не m200.
# ДЕФОЛТ — МЕТОД 1 канона: полный распад ОДНОГО звена, nucleusLimits
# отсекает дочерние (D-001, 21.08). Таблица линий (старый метод)
# в архиве: results/_attic_table_method_20260821/, build/_attic_table_method_20260821/.
# Переопределяемо env-переменными — для точечного возврата к архиву без
# раздвоения скрипта на два файла (§33), не как рабочий режим.
BG_DIR = os.environ.get("G4MODELS_BG_DIR",
                        os.path.join(RESULTS, "bare", "background"))
BG_PREFIX = os.environ.get("G4MODELS_BG_PREFIX", "bg_bare_field_m1")
WF_PREFIX = os.environ.get("G4MODELS_WF_PREFIX", "wf_m1")

NUCS = ["K40", "Ra226", "Pb214", "Bi214", "Pb212", "Ac228", "Bi212", "Tl208"]
CYL = dict(r=45.0, z0=-45.0, z1=120.0)
MEASURED = os.path.join(r"C:\g4work\measured\RadiaCode-103",
                        "\u0424\u043e\u043d 7 \u0434\u043d\u0435\u0439 \u0431\u0435\u0437 \u0434\u043e\u043c\u0438\u043a\u0430.xml")
CAL_ROOM = [-3.711311, 2.444318, 0.000321]
MU_GLOB = "mu3_r500_"
MU_FLUX_PDG = 0.0167
MU_RDISK_MM = 500.0

# RAVNOVESIE Ra-226 i Th-232 ZADANO YAVNO (operator, 21.08: "torij i radij v
# ravnovesii dolzhny byt") - zvenya odnoj cepochki skladyvayutsya v ODIN
# prediktor DO podgonki (main(): merge_by_chain), a ne fitiruyutsya poodinochke
# s posleduyushchej proverkoj otnoshenij. K-40 - ne cepochka, svoya kolonka.
CHAIN_OF = {"K40": "K-40", "Ra226": "Ra-226", "Pb214": "Ra-226", "Bi214": "Ra-226",
            "Pb212": "Th-232", "Ac228": "Th-232", "Bi212": "Th-232", "Tl208": "Th-232"}

# Zvenja, stojashchie POSLE Rn-222 v cepochke Ra-226: tol'ko ih zadevaet
# emanacija. Ra-226 stoit DO radona i ostajotsja.
AFTER_RADON = {"Pb214", "Bi214"}
RN_LEAK = float(os.environ.get("G4MODELS_RN_LEAK", "0.0"))


def read_wallfield_total(path):
    s = 0.0
    for line in open(path, encoding="utf-8"):
        if line.startswith("#") or not line[:1].isdigit():
            continue
        s += float(line.split(",")[1])
    return s


def load_templates():
    r, hz = CYL["r"] / 10.0, 0.5 * (CYL["z1"] - CYL["z0"]) / 10.0
    area = 2 * math.pi * r * (r + 2 * hz)
    names, cols = [], []
    for n in NUCS:
        wf = os.path.join(BUILD, "%s_%s.csv" % (WF_PREFIX, n))
        bg = os.path.join(BG_DIR, "%s_%s.csv" % (BG_PREFIX, n))
        if not (os.path.exists(wf) and os.path.exists(bg)):
            print("[--] %s: shablon net" % n)
            continue
        flu = read_wallfield_total(wf)
        rate = flu * area / 4.0
        meta, hist = rcspec.read_spec(bg)
        t_run = float(meta["N_primaries"]) / rate
        cps = rcspec.fold(hist / t_run, "103")
        names.append(n)
        cols.append(cps)
        print("[%-6s] flu %.4e  rate %.2f 1/s  T %.3e s  cps/(Bq/kg) %.5f  hits %s"
              % (n, flu, rate, t_run, cps.sum(), meta.get("N_with_signal", "?")))
    return names, cols


def merge_by_chain(names, cols, rn_leak=None):
    """Skladyvaet zvenya odnoj cepochki v ODIN shablon (CHAIN_OF) - ravnovesie
    zadano do podgonki, ne proverjaetsya postfaktum (operator, 21.08).

    RN_LEAK - emanacija radona (#RN-1, operator 21.08). Rn-222 - GAZ, chast'
    ego uhodit iz materiala do raspada. Cepochka rvjotsja IMENNO na njom:
    zvenja POSLE radona (Pb-214, Bi-214) tereajut doljy f, a Ra-226 (186 keV)
    ostajotsja na meste. Ravnovesie po-prezhnemu ZADANO, no zadano NARUSHENNYM
    na izvestnuju velichinu - eto fizika, a ne podgonka po polosam.
    Th-232 NE zatragivaetsja: tam gaz Rn-220 zhivjot 55 s i ujti ne uspevaet.
    """
    f = RN_LEAK if rn_leak is None else float(rn_leak)
    merged = {}
    for n, c in zip(names, cols):
        w = (1.0 - f) if n in AFTER_RADON else 1.0
        key = CHAIN_OF.get(n, n)
        merged[key] = merged.get(key, 0.0) + w * c
    return list(merged.keys()), list(merged.values())


def load_muons():
    htot, ntot = None, 0.0
    for p in sorted(glob.glob(os.path.join(BUILD, MU_GLOB + "*.csv"))):
        if p.endswith(".sumw2.csv"):
            continue
        meta, h = rcspec.read_spec(p)
        ntot += float(meta["N_primaries"])
        htot = h if htot is None else htot + h
    if not ntot:
        return None, 0.0
    per = rcspec.fold(htot[:rcspec.NBINS] / ntot, "103")
    pdg = MU_FLUX_PDG * math.pi * (MU_RDISK_MM / 10.0) ** 2
    print("[mu    ] N %.3e  otklik_na_1 %.3e  PDG %.0f mu/s" % (ntot, per.sum(), pdg))
    return per, pdg


def main():
    smp = read_rcxml.read(MEASURED)[0]
    cnt = smp.counts[:-1].astype(float)
    ch = np.arange(len(cnt))
    e_meas = np.asarray(sum(c * ch ** i for i, c in enumerate(CAL_ROOM)))[:len(cnt)]
    cps_meas = cnt / smp.live
    print("izmerenie: %s, live %.2f sut" % (os.path.basename(MEASURED), smp.live / 86400))

    names, cols = load_templates()
    names, cols = merge_by_chain(names, cols)
    mu, pdg = load_muons()
    if mu is not None:
        names.append("mu")
        cols.append(mu)
    if len(names) < 2:
        raise SystemExit("shablonov malo - snachala progony")

    # NA SHKALU IZMERENIYA - CHEREZ fl.rebin_model_to_meas, NE np.interp.
    # np.interp beret modelnuyu PLOTNOST v odnoj tochke na kanal izmereniya i
    # nedoschityvaet model vo stolko raz, vo skolko realnyj kanal shire
    # modelnogo (u RC-103 ~2,77x). Eto uzhe bylo najdeno v proekte 12.08.2026
    # (sm. fit_lines.rebin_model_to_meas), i ya povtoril tu zhe oshibku 20.08:
    # K-40 poluchalos 910 Bq/kg vmesto ~330, polosy 0,775 vmesto ~1.
    A = np.zeros((len(e_meas), len(cols)))
    for k, c in enumerate(cols):
        A[:, k] = fl.rebin_model_to_meas(np.arange(len(c)) + 0.5, c, e_meas)

    print("")
    print("=== VYROZHDENNOST BAZISA (do podgonki) ===")
    norm = A / np.maximum(np.linalg.norm(A, axis=0), 1e-300)
    C = norm.T @ norm
    print("%-8s" % "" + "".join("%8s" % n for n in names))
    for i, n in enumerate(names):
        print("%-8s" % n + "".join("%8.3f" % C[i, j] for j in range(len(names))))
    sv = np.linalg.svd(norm, compute_uv=False)
    print("cond = %.1f  (>100 -> chast amplitud ne opredelitsya)" % (sv[0] / max(sv[-1], 1e-300)))
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            if C[i, j] > 0.995:
                print("  !! %s i %s korrelyaciya %.4f - razdelenie nedostoverno"
                      % (names[i], names[j], C[i, j]))

    counts = cps_meas * smp.live
    Ac = A * smp.live
    w = 1.0 / np.sqrt(np.maximum(counts, 1.0))
    from scipy.optimize import nnls
    amp, resid = nnls(Ac * w[:, None], counts * w)
    Aw = Ac * w[:, None]
    try:
        cov = np.linalg.inv(Aw.T @ Aw)
        sd = np.sqrt(np.maximum(np.diag(cov), 0.0))
    except np.linalg.LinAlgError:
        sd = np.full(len(amp), float("nan"))

    print("")
    print("=== AMPLITUDY (nuklidnoe razlozhenie) ===")
    print("%-8s %12s %12s %10s  %s" % ("zveno", "amplituda", "oshibka", "otn", "edinicy"))
    for n, a, s in zip(names, amp, sd):
        unit = "mu/s cherez disk" if n == "mu" else "Bq/kg"
        flag = ""
        if a > 0 and s == s and s > a:
            flag = "  <-- NE OPREDELENA"
        elif a == 0:
            flag = "  <-- obnulena NNLS"
        print("%-8s %12.2f %12.2f %10s  %s%s"
              % (n, a, s, ("%.0f%%" % (100 * s / a)) if a > 0 else "-", unit, flag))

    d = dict(zip(names, amp))
    ds = dict(zip(names, sd))

    def ratio(a, b):
        if a in d and b in d and d[b] > 0 and d[a] > 0:
            r = d[a] / d[b]
            er = r * math.sqrt((ds[a] / d[a]) ** 2 + (ds[b] / d[b]) ** 2)
            return r, er
        return float("nan"), float("nan")

    print("")
    print("=== RAVNOVESIE CEPOCHEK ===")
    for a, b, exp, what in (("Bi214", "Pb214", 1.00, "radon: pozdnie/rannie DPR"),
                            ("Pb214", "Ra226", 1.00, "radon: emanirovanie Ra-226"),
                            ("Tl208", "Ac228", 0.3594, "torij: vetvlenie Bi-212->Tl-208"),
                            ("Bi212", "Pb212", 1.00, "torij: Bi-212/Pb-212"),
                            ("Ac228", "Pb212", 1.00, "torij: Ac-228/Pb-212"),
                            ("Tl208", "Pb212", 1.00, "torij: Tl-208/Pb-212 (oba vidny)")):
        r, er = ratio(a, b)
        if r == r:
            print("  %-6s/%-6s = %6.2f +- %5.2f   ozhid %.4f   otn %.2f   (%s)"
                  % (a, b, r, er, exp, r / exp, what))

    pred = A @ amp
    print("")
    print("=== SVERKA PO POLOSAM ===")
    print("%-12s %10s %10s %8s" % ("polosa,keV", "izmereno", "model", "m/i"))
    for lo, hi in ((20, 100), (100, 300), (300, 700), (700, 1500),
                   (1500, 2000), (2000, 2400), (2400, 2830)):
        m = (e_meas >= lo) & (e_meas < hi)
        ym, pm = cps_meas[m].sum(), pred[m].sum()
        print("%5d-%-6d %10.5f %10.5f %8.3f" % (lo, hi, ym, pm, pm / ym if ym else float("nan")))
    m = (e_meas >= 20) & (e_meas < 2830)
    print("%-12s %10.5f %10.5f %8.3f" % ("polnyj", cps_meas[m].sum(), pred[m].sum(),
                                         pred[m].sum() / cps_meas[m].sum()))
    if "mu" in d and pdg:
        print("")
        print("muony: podobrano %.0f protiv PDG %.0f -> otn %.2f" % (d["mu"], pdg, d["mu"] / pdg))
    print("")
    print("!!! POKANALNYJ REZULTAT VYSHE - ARTEFAKT VYROZHDENNOSTI, NE FIZIKA.")
    print("    Rabochij metod - podgonka po linijam, nizhe.")
    mu_idx = names.index("mu") if "mu" in names else None
    return fit_by_lines(names, A, e_meas, cps_meas, smp.live, mu_idx, pdg)



# ---------------------------------------------------------------------------
# ПОДГОНКА ПО НЕТТО-ПЛОЩАДЯМ ЛИНИЙ (рабочий метод)
#
# Поканальный NNLS выше ОТВЕРГНУТ по факту прогона 20.08: корреляции столбцов
# 0.997-0.9995, cond=190, результат абсурден (Ra-226 4095 Бк/кг при физичных
# 10-40, K-40 обнулён при явной линии 1461, мюоны 8.8xPDG против честных
# 1.25xPDG из чистого окна). Это тот же класс, что зафиксирован в
# fit_room_field.py для трёх цепочек: комптоновский континуум в 20-300 кэВ
# несёт ~98% статистики и по форме одинаков у всех звеньев, поэтому по каналам
# звенья неразличимы, а NNLS раскладывает вес между ними произвольно.
#
# Уникальность звена - в ЕГО ЛИНИЯХ. Поэтому система строится на НЕТТО-площадях
# в окнах линий: подложка вычитается боковыми полосами ОДНИМ И ТЕМ ЖЕ
# алгоритмом в измерении и в каждом шаблоне, поэтому систематика вычитания
# сокращается, а матрица становится почти диагональной.
#
#   net_meas(окно k) = SUM_j a_j * net_model_j(окно k)
#
# Окна берутся +-0.75*FWHM вокруг линии, боковые - той же ширины по обе стороны.
# Мюонная амплитуда - из ЧИСТОГО окна 2700-2830 (выше всех гамма-линий).
# ---------------------------------------------------------------------------

# диагностические линии на звено (ENSDF/LNHB); только те, что реально видны
# на CsI при FWHM ~9% и не перекрыты соседом другого звена вплотную
DIAG = [
    ("K40",   1460.8),
    ("Pb214",  351.9),
    ("Bi214",  609.3), ("Bi214", 1120.3), ("Bi214", 1764.5),
    ("Pb212",  238.6),
    ("Ac228",  911.2), ("Ac228", 968.9),
    ("Bi212",  727.3),
    ("Tl208",  583.2), ("Tl208", 2614.5),
    ("Ra226",  186.2),
]


def net_window(spec, e_grid, e0, live=None):
    """NETTO v okne linii - CHEREZ GOTOVUYU fit_lines.line_net_area (§33).

    Sobstvennaya realizaciya OTVERGNUTA 20.08: okno +-0.75*FWHM bez korrekcii
    na pokrytie davalo K-40 = 733 Bq/kg protiv 324 u nezavisimogo fit_lines.py
    na teh zhe dannyh - dva svoih instrumenta rashodilis vdvoe. Zdes ta zhe
    konvenciya okon (nsig=2.5, gap=1.6, koefficient pokrytiya 0.9876), chto
    uzhe prinyata v proekte, poetomu rezultaty sravnimy naprjamuyu.
    """
    r = fl.line_net_area(e_grid, spec, e0)
    if r is None:
        return 0.0, float("nan")
    if live is None:
        return r["area"], float("nan")
    # sd u donora - v edinicah spektra (cps), no schitaetsya kak sqrt(gross+|cont|)
    # v teh zhe edinicah; perevodim v puassonovskuyu oshibku po OTSCHETAM
    gross_c = r["gross"] * live
    cont_c = abs(r["cont"]) * live
    return r["area"], math.sqrt(max(gross_c + cont_c, 1.0)) / live / 0.9876


def fit_by_lines(names, A, e_meas, cps_meas, live, mu_col_idx, pdg):
    """A - матрица [канал x звено] в cps на 1 Бк/кг (мюонный столбец - на 1 мюон/с)"""
    gamma_idx = [i for i in range(len(names)) if i != mu_col_idx]
    gnames = [names[i] for i in gamma_idx]

    rows, y, dy, labels, skipped, DIAG_E = [], [], [], [], [], []
    for nuc, e0 in DIAG:
        if CHAIN_OF.get(nuc, nuc) not in gnames:
            continue
        nm, dnm = net_window(cps_meas, e_meas, e0, live)
        row = [net_window(A[:, i], e_meas, e0)[0] for i in gamma_idx]
        if max(row) <= 0:
            skipped.append((nuc, e0, "shablon ne daet pika v okne"))
            continue
        # linija dolzhna byt VIDNA nad fonom: bez etogo v sistemu popadaet shum
        # (i otricatelnye netto), i matrica vyrozhdaetsya (cond ~1e16, 20.08)
        if not (dnm > 0 and nm > 3.0 * dnm):
            skipped.append((nuc, e0, "ne znachima: %.1f sigma" % (nm / dnm if dnm > 0 else 0)))
            continue
        rows.append(row)
        DIAG_E.append(e0)
        y.append(nm)
        dy.append(dnm)
        labels.append("%s %.1f" % (nuc, e0))
    # ZVENYA BEZ SOBSTVENNOJ ZNACHIMOJ LINII iz podgonki ISKLYUCHAYUTSYA.
    # Inache NNLS pripisyvaet im amplitudu iz ih VKLADA V CHUZHIE okna - eto ne
    # izmerenie zvena, a shum, i on portit sosednie amplitudy (proverено 20.08:
    # Bi-212 bez svoej linii poluchal 252 +- 362 Bq/kg i tyanul za soboj Pb-212).
    have_line = set(CHAIN_OF.get(nuc, nuc) for nuc, e0 in DIAG
                    if any(lab.startswith(nuc + " ") for lab in labels))
    keep = [k for k, n in enumerate(gnames) if n in have_line]
    dropped = [n for n in gnames if n not in have_line]
    for n in dropped:
        print("[-] %-6s - net ni odnoj znachimoj sobstvennoj linii, ISKLYUCHENO iz podgonki" % n)
    gnames = [gnames[k] for k in keep]
    gamma_idx = [gamma_idx[k] for k in keep]
    rows = [[r[k] for k in keep] for r in rows]

    # SOVMESTNYJ FIT MYUONOV I GAMMA (ispravleno 20.08 po zamechaniyu operatora
    # "ubav myuonov dobav toriya"). Bylo POSLEDOVATELNO: gamma po linijam ->
    # myuony iz OSTATKA v okne 2700-2828. Pri takoj sheme zanizhennyj Tl-208
    # avtomaticheski zavyshaet myuony: ves nepokrytyj gammoj schet v zhestkoj
    # oblasti uhodit v myuonnyj stolbec. Teper myuony - polnopravnyj stolbec
    # sistemy, a v sistemu dobavlena stroka "chistoe myuonnoe okno" (2700-2828,
    # vyshe vseh gamma-linij). Vnutri kazhdoj stroki izmerenie i model schitayutsya
    # ODINAKOVO: v strokah linij - netto nad podlozhkoj, v myuonnoj stroke -
    # polnyj integral okna (tam kontinuum, pika net).
    if mu_col_idx is not None:
        for k in range(len(rows)):
            rows[k].append(net_window(A[:, mu_col_idx], e_meas, DIAG_E[k])[0])
        msk = (e_meas >= 2700) & (e_meas < 2828)
        rows.append([A[:, i][msk].sum() for i in gamma_idx] + [A[:, mu_col_idx][msk].sum()])
        y.append(cps_meas[msk].sum())
        dy.append(math.sqrt(max(cps_meas[msk].sum() * live, 1.0)) / live)
        labels.append("myuonnoe okno 2700-2828")
        gnames = gnames + ["mu"]

    M = np.array(rows)
    y = np.array(y)
    dy = np.array(dy)

    print("")
    for nuc, e0, why in skipped:
        print("[-] %-6s %7.1f keV - ISKLYUCHENA: %s" % (nuc, e0, why))
    print("")
    print("=== NETTO V OKNAH ZNACHIMYH LINIJ ===")
    print("%-14s %12s %12s   %s" % ("liniya", "net izm cps", "sigma", "dominiruet v shablone"))
    for k, lab in enumerate(labels):
        j = int(np.argmax(M[k]))
        tot = sum(v for v in M[k] if v > 0)
        share = M[k][j] / tot if tot > 0 else float("nan")
        print("%-14s %12.3e %12.3e   %s (%.0f%% okna)"
              % (lab, y[k], dy[k], gnames[j], 100 * share))

    w = 1.0 / np.maximum(dy, 1e-300)
    from scipy.optimize import nnls
    amp, _ = nnls(M * w[:, None], y * w)
    Mw = M * w[:, None]
    try:
        cov = np.linalg.inv(Mw.T @ Mw)
        sd = np.sqrt(np.maximum(np.diag(cov), 0.0))
    except np.linalg.LinAlgError:
        sd = np.full(len(amp), float("nan"))

    sv = np.linalg.svd(M / np.maximum(np.linalg.norm(M, axis=0), 1e-300),
                       compute_uv=False)
    print("cond matricy linij = %.1f (protiv pokanalnoj - sm. vyshe)"
          % (sv[0] / max(sv[-1], 1e-300)))

    a_mu = d_mu = float("nan")

    print("")
    print("=== PROVERKA: POODINOCHNO PO KAZHDOJ LINII (bez NNLS) ===")
    print("%-14s %10s %10s   %s" % ("liniya", "a, Bq/kg", "sigma", "dolya svoego v okne"))
    for k, lab in enumerate(labels):
        nuc = lab.split()[0]
        g = CHAIN_OF.get(nuc, nuc)
        if g not in gnames:
            continue          # myuonnaya stroka - ne liniya, poodinochke ne schitaetsya
        j = gnames.index(g)
        own = M[k][j]
        tot_pos = sum(v for v in M[k] if v > 0)
        if own > 0:
            print("%-14s %10.1f %10.1f   %5.0f%%"
                  % (lab, y[k] / own, dy[k] / own, 100 * own / tot_pos if tot_pos else float("nan")))

    print("")
    print("=== AMPLITUDY PO LINIYAM ===")
    print("amplituda = aktivnost RODITELYA cepochki (Ra-226 / Th-232 / K-40),")
    print("vosstanovlennaya PO LINIJAM etogo zvena; pri ravnovesii vse ravny")
    print("%-8s %12s %12s %8s  %s" % ("zveno", "Bq/kg", "oshibka", "otn", ""))
    for n, a, s in zip(gnames, amp, sd):
        flag = "  <-- NE OPREDELENA" if (a > 0 and s == s and s > a) else (
               "  <-- obnulena" if a == 0 else "")
        print("%-8s %12.2f %12.2f %8s%s"
              % (n, a, s, ("%.0f%%" % (100 * s / a)) if a > 0 else "-", flag))
    if "mu" in gnames:
        a_mu = amp[gnames.index("mu")]
        d_mu = sd[gnames.index("mu")]

    # RAVNOVESIE Ra-226 i Th-232 UZHE ZADANO na urovne shablonov (main():
    # merge_by_chain) - kazhdaya cepochka fitiruetsya ODNOJ amplitudoj, zdes
    # nechego proverjat postfaktum i nechego "dostraivat": otdelnyh zvenev
    # v gnames bolshe net, kazhdaya cepochka uzhe polnaya odnoj kolonkoj.
    full = np.zeros(len(e_meas))
    for k, i in enumerate(gamma_idx):
        full += amp[k] * A[:, i]
    if mu_col_idx is not None and "mu" in gnames:
        full += amp[gnames.index("mu")] * A[:, mu_col_idx]

    print("")
    print("=== SVERKA PO POLOSAM (amplitudy iz linij, NE podognano po polosam) ===")
    print("%-12s %10s %10s %8s" % ("polosa,keV", "izmereno", "model", "m/i"))
    for lo, hi in ((20, 100), (100, 300), (300, 700), (700, 1500),
                   (1500, 2000), (2000, 2400), (2400, 2830)):
        m = (e_meas >= lo) & (e_meas < hi)
        ym, pm = cps_meas[m].sum(), full[m].sum()
        print("%5d-%-6d %10.5f %10.5f %8.3f" % (lo, hi, ym, pm, pm / ym if ym else float("nan")))
    m = (e_meas >= 20) & (e_meas < 2830)
    print("%-12s %10.5f %10.5f %8.3f" % ("polnyj", cps_meas[m].sum(), full[m].sum(),
                                         full[m].sum() / cps_meas[m].sum()))
    return dict(zip(gnames, amp)), a_mu

if __name__ == "__main__":
    main()
