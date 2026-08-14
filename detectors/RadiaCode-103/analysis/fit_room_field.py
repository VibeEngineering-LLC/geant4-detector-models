# -*- coding: utf-8 -*-
"""Подгонка трёх амплитуд поля ЕРН (K-40, ряд Ra-226, ряд Th-232) к измеренному
фону RC-103, вместо одного справочного (UNSCEAR) масштаба на бетон.

ПОЧЕМУ. Сверка модели с измерением поканально показала: до 1550 кэВ модель
завышает фон РОВНО одним множителем ~1,68, а полоса 1550-2700 кэВ (там живёт
Tl-208 2614,5) — недобирает почти вдвое. Один скаляр это не лечит: жёсткая
компонента (торий) лучше всех проходит свинец и определяет ответ на больших
толщинах, поэтому важно не среднее, а СОСТАВ поля. См. план, раздел
«Одним множителем не обойтись».

КОНВЕЙЕР (три стадии, каждая — воспроизводимый шаг):
  1. wallfield.exe <N> wf_<S>.csv <S>     — единичный (1 Бк/кг) отклик серии S
     в воздушной полости внутри бетона. Уже посчитано отдельно (см. лог).
  2. Спектр каждой серии -> GPS Arb макрос -> прогон rc_curves.exe (та же
     геометрия и конвенция, что run_bg.py: цилиндр Ф=4N/S, конфиг
     full/air/0.0012 — «пустой сосуд», штатный холостой набор). Даёт
     АБСОЛЮТНЫЙ спектр в кристалле, имп/с НА 1 Бк/кг активности серии.
  3. NNLS: measured_cps(E) = a_K*m_K(E) + a_Ra*m_Ra(E) + a_Th*m_Th(E),
     a >= 0 — амплитуды в Бк/кг, физический смысл прямой.

Запуск:
    python fit_room_field.py gen        — сгенерировать макросы GPS Arb
    python fit_room_field.py run [N]    — прогнать 3 серии через rc_curves
    python fit_room_field.py fit        — подогнать амплитуды, напечатать
"""
import math
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for _d in ("analysis", "drivers"):
    _p = os.path.join(HERE, "..", _d)
    if _p not in sys.path:
        sys.path.insert(0, _p)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "common", "py"))
import paths  # noqa: E402
import rcspec  # noqa: E402
import read_rcxml  # noqa: E402

SERIES = ("K", "Ra", "Th")
RESULTS = rcspec.RESULTS
BUILD = str(paths.build("RadiaCode-103"))
WF_EXE = os.path.join(BUILD, "wallfield.exe")
RC_EXE = os.path.join(BUILD, "rc_curves.exe")
COSMICMU_CSV = os.path.join(BUILD, "cosmicmu.csv")

# Тот же охватывающий цилиндр m200, что в run_bg.py — воспроизводить его
# здесь, а не импортировать, чтобы не тянуть в этот файл ThreadPoolExecutor
# драйвера ради одного словаря; числа СВЕРЕНЫ построчно с run_bg.py:66-67.
CYL_M200 = dict(r=45.0, z0=-45.0, z1=120.0)

# Личные измерения оператора — путь НЕ хардкодить (правило репозитория,
# common/py/paths.py). paths.measured() сам сообщит понятной ошибкой, если
# G4MODELS_MEASURED не задана и каталога измерений нет.
MEASURED_BG = str(paths.measured("RadiaCode-103") / "Фон 7 дней без домика.xml")


def field_mac_path(s):
    return os.path.join(RESULTS, "field_spectrum_%s.mac" % s)


def wf_csv_path(s):
    return os.path.join(BUILD, "wf_%s.csv" % s)


def bg_csv_path(s):
    d = rcspec.rdir("background")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, "bg_cyl_field_%s.csv" % s)


# --- шаг 1: чтение wallfield.csv (тот же формат, что analyze_field.py) ------
def read_wallfield(path):
    e, f = [], []
    for line in open(path, encoding="utf-8"):
        if line.startswith("#") or not line[:1].isdigit():
            continue
        a, b = line.split(",")
        e.append(float(a))
        f.append(float(b))
    return np.array(e), np.array(f)


# cosmicmu.csv шире, чем rcspec.NBINS=3201 (у мюона депозиты до 20 МэВ,
# основной пик 4-10 МэВ) — rcspec.read_spec() ТИХО обрезал бы всё за
# каналом 3200, поэтому свой читатель, без фиксированного размера массива.
def read_cosmicmu(path):
    meta = {}
    e, c = [], []
    for line in open(path, encoding="utf-8"):
        if line.startswith("#"):
            if "=" in line:
                k, v = line[1:].split("=", 1)
                meta[k.strip()] = v.strip()
        elif line[:1].isdigit():
            a, b = line.split(",")
            e.append(float(a))
            c.append(float(b))
    n = int(round(max(e))) + 2 if e else 1
    hist = np.zeros(n)
    for ei, ci in zip(e, c):
        hist[int(ei)] = ci
    return meta, hist


def gen_macros():
    for s in SERIES:
        p = wf_csv_path(s)
        if not os.path.exists(p):
            raise SystemExit("нет %s — сначала wallfield.exe N wf_%s.csv %s" % (p, s, s))
        e, flu = read_wallfield(p)
        out = field_mac_path(s)
        with open(out, "w", encoding="utf-8") as f:
            f.write("# Единичный (1 Бк/кг) отклик серии %s, wallfield.exe\n" % s)
            f.write("/gps/particle gamma\n/gps/ene/type Arb\n/gps/hist/type arb\n")
            for ei, fi in zip(e, flu):
                if fi > 0:
                    f.write("/gps/hist/point %.6f %.6e\n" % (ei / 1000.0, fi))
            f.write("/gps/hist/inter Lin\n")
        print("[gen]", s, "->", out, " (полный флюенс %.4e см^-2 с^-1)" % flu.sum())


# --- шаг 2: прогон через rc_curves, та же геометрия, что run_bg.py ---------
def run_series(n_events):
    c = CYL_M200
    for s in SERIES:
        out = bg_csv_path(s)
        if os.path.exists(out):
            print("[--]", s, "уже посчитано:", out)
            continue
        mac = os.path.join(BUILD, "field_run_%s.mac" % s)
        with open(mac, "w", encoding="utf-8") as f:
            f.write("\n".join([
                "/run/verbose 0", "/event/verbose 0", "/tracking/verbose 0",
                "/run/printProgress 0", "",
                "/gps/pos/type Surface",
                "/gps/pos/shape Cylinder",
                "/gps/pos/radius %.3f mm" % c["r"],
                "/gps/pos/halfz %.3f mm" % (0.5 * (c["z1"] - c["z0"])),
                "/gps/pos/centre 0 0 %.3f mm" % (0.5 * (c["z1"] + c["z0"])),
                "/gps/ang/type cos",
                "/control/execute %s" % field_mac_path(s).replace("\\", "/"),
                "/rc/outFile %s" % out.replace("\\", "/"),
                "/run/beamOn %d" % n_events, "",
            ]))
        t0 = time.time()
        log = out + ".log"
        with open(log, "w", encoding="utf-8") as lf:
            p = subprocess.run([RC_EXE, mac, "full", "air", "0.0012", "m200"],
                               cwd=BUILD, stdout=lf, stderr=subprocess.STDOUT)
        print("[%s] %s  %.1f мин" % ("ok" if p.returncode == 0 else "СБОЙ", s,
                                     (time.time() - t0) / 60), flush=True)


# --- шаг 3: NNLS-подгонка ----------------------------------------------------
# Нормировка — тождество Ф=4N/S (run_bg.py, wallfield.cc): темп первичных
# rate = fluence_total*area/4, реальное время прогона T = N/rate. Активность
# 1 Бк/кг уже "въедена" в флюенс на шаге 1 (wallfield с gSeries), поэтому
# здесь она не учитывается второй раз — только перевод розыгрыша в секунды.
def fit():
    # measured
    smp = read_rcxml.read(MEASURED_BG)[0]
    e_meas = smp.energy
    cps_meas = smp.counts / smp.live

    c = CYL_M200
    r, hz = c["r"] / 10, 0.5 * (c["z1"] - c["z0"]) / 10
    area = 2 * math.pi * r * (r + 2 * hz)

    bases = []
    for s in SERIES:
        e_flu, flu = read_wallfield(wf_csv_path(s))
        fluence_total = flu.sum()          # см^-2 с^-1, на 1 Бк/кг
        rate = fluence_total * area / 4.0  # первичных в секунду реального времени
        meta, hist = rcspec.read_spec(bg_csv_path(s))
        n = int(meta["N_primaries"])
        t_run = n / rate                   # секунд, эквивалентных этому прогону
        cps = hist / t_run                 # имп/с на канал, НА 1 Бк/кг серии s
        cps = rcspec.fold(cps, "103")       # свёртка с разрешением — измерение тоже размыто
        bases.append(cps)
        print("[%s] fluence_total=%.4e rate=%.2f 1/s t_run=%.1f с  cps_total=%.5f/(Бк/кг)"
              % (s, fluence_total, rate, t_run, cps.sum()))

    # 4-я колонка — космический мюонный континуум (см. cosmicmu.cc). ЕДИНИЦЫ
    # ДРУГИЕ: это не Бк/кг, а "мюонов/с через диск источника" — базис нормирован
    # НА ОДИН мюон (вероятность отклика), амплитуда a_mu подбирается NNLS так
    # же, как активности, а физический смысл (согласуется ли a_mu с реальным
    # потоком PDG ~1/см²/мин по площади диска) проверяется ПОСЛЕ подгонки,
    # не задаётся заранее.
    mu_amp_label = "mu (1/с через диск)"
    if os.path.exists(COSMICMU_CSV):
        meta_mu, hist_mu = read_cosmicmu(COSMICMU_CSV)
        n_mu = int(meta_mu["N_primaries"])
        cps_per_primary = hist_mu / n_mu    # отклик НА ОДИН мюон, безразмерно
        cps_per_primary = rcspec.fold(cps_per_primary, "103")
        bases.append(cps_per_primary)
        print("[mu] N_primaries=%d  hits=%s  отклик_на_1_мюон=%.3e"
              % (n_mu, meta_mu.get("N_with_signal", "?"), cps_per_primary.sum()))
        n_cols = 4
        col_names = list(SERIES) + ["mu"]
    else:
        print("[mu] cosmicmu.csv ещё не готов — подгонка БЕЗ мюонной компоненты "
              "(временно, до готовности прогона)")
        n_cols = 3
        col_names = list(SERIES)

    # приводим модель к энергетической шкале измерения; колонки разной длины
    # (K/Ra/Th — 3201 канал, mu — до ~20000), np.interp сам обрежет справа
    A_ch = np.zeros((len(e_meas), n_cols))
    for k, cps in enumerate(bases):
        e_model = np.arange(len(cps)) + 0.5
        A_ch[:, k] = np.interp(e_meas, e_model, cps, left=0, right=0)

    # ПОКАНАЛЬНЫЙ NNLS отвергнут: континуум K/Ra в низкой энергии на порядки
    # больше слабой линии Tl-208 2614,5, и обычный NNLS (минимизация суммы
    # квадратов без учёта структуры) обнулил Th-232 целиком — притом что
    # именно в полосе 1550-2700 кэВ измерение прямо требует ненулевого Th
    # (проверено: там pred=0, measured=0.00031). Правило контура —
    # chi2_by_energy_bands, скилл geant4-spectrum-pipeline v1.4.0 — разрезы
    # по группам нуклидов существуют ИМЕННО для такого случая. Подгонка по
    # ТЕМ ЖЕ 7 энергетическим полосам, что в поканальной сверке (план,
    # раздел «Одним множителем не обойтись») — не выбор круглых чисел, а
    # прямое переиспользование уже найденной диагностической разбивки.
    BANDS = [(20, 100), (100, 300), (300, 600), (600, 750),
             (750, 1400), (1400, 1550), (1550, 2700)]

    # ПОДГОНКА ПО 7 ШИРОКИМ ПОЛОСАМ ОТВЕРГНУТА (не то, что заявлено раньше в
    # этом комментарии) — проверка корреляции столбцов показала: K, Ra, Th
    # коррелируют между собой на 1,000 по этим 7 полосам. Комптоновский
    # континуум "забывает" исходный нуклид на низких энергиях, где сосредо-
    # точена вся статистика измерения (полосы 20-300 кэВ дают 98% сигнала) —
    # три серии математически НЕРАЗЛИЧИМЫ этой разбивкой, и NNLS произвольно
    # перекладывал вес между ними (305/127/0, затем 254/129/0, затем 0/145/0
    # при добавлении мюонов — псевдослучайные артефакты вырожденности, не
    # физика). Уникальность нуклида — в ЛИНИЯХ (1460,8 K-40; 609,3/1764,5
    # Ra-226; 583,2/911,2/2614,5 Th-232), не в интегралах широких полос.
    #
    # Возврат к ПОКАНАЛЬНОЙ подгонке — со СТАТИСТИЧЕСКИМ ВЕСОМ (Neyman chi2),
    # а не равным весом по SSE (первая отвергнутая попытка) и не в
    # пространстве СКОРОСТИ СЧЁТА (вторая попытка, тоже отвергнута): вес
    # 1/sqrt(cps) там срывался на редких каналах — один канал с 1 отсчётом
    # за 612250 с живого даёт cps=1,6е-6, вес ~770, и единственный шумный
    # канал забивал всю подгонку. Пол должен быть в пространстве ОТСЧЁТОВ
    # (считать целыми людьми, не долями в секунду): 1 отсчёт — осмысленный
    # минимум дисперсии по Пуассону, не зависит от времени набора.
    counts_meas = cps_meas * smp.live
    counts_pred_basis = A_ch * smp.live   # тот же базис, в отсчётах
    w = 1.0 / np.sqrt(np.maximum(counts_meas, 1.0))
    from scipy.optimize import nnls
    amp, resid = nnls(counts_pred_basis * w[:, None], counts_meas * w)

    print("\nПОДГОНКА ПО КАНАЛАМ, вес 1/sqrt(measured) (Neyman chi2):")
    print("%-12s %10s %10s %10s" % ("полоса,кэВ", "измерено", "модель", "модель/изм"))
    pred_check = A_ch @ amp
    for lo, hi in BANDS:
        m = (e_meas >= lo) & (e_meas < hi)
        ym, pm = cps_meas[m].sum(), pred_check[m].sum()
        print("%5d-%-6d %10.5f %10.5f %10.3f" % (lo, hi, ym, pm, pm / ym if ym > 0 else float("nan")))
    print("\nПОДОБРАННЫЕ АМПЛИТУДЫ:")
    for s, a in zip(col_names, amp):
        unit = "Бк/кг" if s != "mu" else mu_amp_label
        print("  %-3s %10.2f  %s" % (s, a, unit))
    print("невязка (норма, по полосам) = %.5f" % resid)

    if "mu" in col_names:
        # Сверка порядка величины: реальный поток PDG ~1 см^-2 мин^-1 через
        # площадь диска источника (см. cosmicmu.cc, R_DISK=70 мм) даёт
        # ожидаемый порядок a_mu. НЕ жёсткое ограничение подгонки — проверка
        # ПОСЛЕ, чтобы не подгонять модель под ожидание.
        r_disk_cm2 = math.pi * (7.0) ** 2   # R_DISK=70мм=7см, из cosmicmu.cc
        expect_mu = 1.0 / 60.0 * r_disk_cm2  # ~1/см^2/мин -> 1/с * площадь диска
        a_mu = amp[col_names.index("mu")]
        print("  сверка: поток PDG ~1/см²/мин * площадь диска %.0f см² = %.2f 1/с "
              "(подобрано %.2f, отношение %.2f)"
              % (r_disk_cm2, expect_mu, a_mu, a_mu / expect_mu if expect_mu else float("nan")))

    pred = A_ch @ amp
    chi2 = ((cps_meas - pred) ** 2 / np.maximum(cps_meas, cps_meas.mean())).sum()
    print("chi2 (по каналам, справочно, не критерий подгонки) = %.1f, каналов = %d"
          % (chi2, len(e_meas)))
    return dict(zip(col_names, amp))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "fit"
    if cmd == "gen":
        gen_macros()
    elif cmd == "run":
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 30_000_000
        run_series(n)
    elif cmd == "fit":
        fit()
    else:
        raise SystemExit("команды: gen | run [N] | fit")
