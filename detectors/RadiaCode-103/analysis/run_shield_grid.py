# -*- coding: utf-8 -*-
"""Сетка толщин свинца (Cd=Cu=0): B_room_gamma(t) — гамма-компонента поля
ЕРН помещения (K-40+Ra-226+Th-232) сквозь защиту, посчитанная через
shieldrun trans(bias)+replay по КАЖДОЙ серии отдельно (линейность гамма-
переноса — суммировать после, не до), взвешенная амплитудами из fit_lines.py.

ЧАСТИЧНЫЙ ОБЪЁМ ЗАДАЧИ №4 — НЕ вся сетка из плана (12 узлов Pb + 2D Cd/Cu +
S_K40/S_Cs137 + мюон сквозь защиту). Здесь только гамма-компонента B_room на
пробной сетке толщин, для проверки, что конвейер (биасинг + specmac +
replay + весовая нормировка) физически работает и даёт разумную кривую
затухания. Остальное — отдельными заходами, см. план.

t=0 НЕ через shieldrun (при pb=cu=cd=0 нет слоя depth=0 — TransStep не
на чем сработать): анкер t=0 — прямо из fit_lines.model_cps_curve(),
это ТА ЖЕ модель поля, которой уже подгонялись a_K/a_Ra/a_Th.

Нормировка (Ф=4N/S) — снаружи, той же формулой, что fit_lines.model_cps_curve():
rate = fluence_total(wallfield, S) * area_shield_outer(t) / 4, но площадь
берётся ФАКТИЧЕСКАЯ наружная поверхность защиты НА ЭТОЙ ТОЛЩИНЕ (растёт с t),
не фиксированный CYL_M200 — поле изотропно на любой выпуклой охватывающей
поверхности (обоснование то же, что в run_bg.py), поэтому подстановка
корректна.

Запуск: python run_shield_grid.py [pb1 pb2 ...]     (по умолчанию 20 50 100)
"""
import math
import os
import re
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
import fit_lines as fl  # noqa: E402
import fit_room_field as frf  # noqa: E402
import mu_smooth  # noqa: E402

SERIES = ("K", "Ra", "Th")
BUILD = str(paths.build("RadiaCode-103"))
SHIELDRUN_EXE = os.path.join(BUILD, "shieldrun.exe")
GRID_DIR = os.path.join(BUILD, "shield_grid")

NPRIM = 100_000
IMPSTEP = 2.2
# Геометрия домика — ЯВНО, а не из умолчаний ShieldGeom (PbShield.hh):
# короб 150×150 в плане, 385 мм высотой, ОТКРЫТЫЙ сверху, сборка стоит на
# дне полости. Правка умолчания в .hh иначе молча переопределяла бы уже
# посчитанную сетку, а имена файлов остались бы прежними.
HXCAV, HYCAV, HZCAV = 75.0, 75.0, 192.5
# Верх домика. По умолчанию ОТКРЫТ — так стоит домик у оператора. Переключается
# переменной окружения RC_SHIELD_LID=1, чтобы сравнительный прогон «а что даст
# крышка» не требовал правки исходника (и не остался бы в нём по забывчивости).
# На розыгрыш источника наличие крышки НЕ влияет — выпуклая оболочка замкнута в
# любом случае; влияет на массу свинца и на канал сверху, см. PbShield.cc.
WITH_LID = os.environ.get("RC_SHIELD_LID", "0") not in ("0", "", "no", "off")
SEAT_FLOOR = 1
# толщина/слой держим примерно постоянной (~6,25 мм/слой, как в
# провалидированном pb=100/nshell=16 прогоне из задачи №2) — НЕ фиксированные
# 16 слоёв на любую толщину: при pb=20 непрерывный спектр поля (много мягких
# фотонов, короткий пробег) с nshell=16/impstep=2.0 обвалил счёт в >180с
# (пробовал 12.08.2026, unbiased та же толщина — 2,7с) — не сходимость, а
# банальная избыточность биасинга там, где он не нужен.
MM_PER_SHELL = 6.25
# Порог включения биасинга trans. БЫЛО 40 мм — оказалось слишком высоко:
# 13.08.2026 замер N_eff в окне 662 кэВ (spec_quality.py, по настоящему Σw²)
# показал провал ровно между режимами —
#   pb:      5    10    15    20    25    30 |   40    50   100
#   N_eff: 275   164   103    58    35    17 |  663   530 10607
# слева статистику даёт малая толщина, справа биасинг, а на 20-30 мм не
# работает ни то, ни другое: при pb=30 погрешность 24 % и вердикт «НЕ ГОДНО».
# Опущено до 20 мм. Верхняя граница разумности прежняя (см. комментарий
# выше): на ЕЩЁ более тонких защитах биасинг только тормозит.
BIAS_MIN_PB = 20.0
CRYSTEP = 4.0         # replay-биасинг на путь к кристаллу (задача №13), всегда включён

# Масса пробы m200/organic 0,50 — тот же баланс, что efficiency.csv
# (config full_organic_0.50), 200,15 см³ × 0,50 г/см³ = 100,075 г. Источник —
# ион (K-40/Cs-137), не через trans/replay: путь проба->кристалл короткий,
# биасинг не нужен, mode=sample гоняется напрямую.
SAMPLE_MASS_KG = 0.1001
NUCLIDES_SAMPLE = ("K40", "Cs137")
NPRIM_SAMPLE = 500_000


def bias_params(pb):
    if pb < BIAS_MIN_PB:
        return False, 0, 0.0
    nshell = max(4, int(round(pb / MM_PER_SHELL)))
    return True, nshell, IMPSTEP


def geom_args(pb, horiz=False, gap_mm=None, with_vessel=True):
    """Ключи геометрии домика для shieldrun — один список на все режимы.

    Собран в одном месте намеренно: trans и sample обязаны видеть ОДНУ и ту
    же защиту, иначе отклик пробы считается сквозь одну геометрию, а фон —
    сквозь другую, и разъезд ничем себя не проявит.

    Умолчания аргументов воспроизводят прежнее поведение ПОБИТОВО (прибор
    стоит вертикально на дне, сосуд m200 на месте): сетка толщин и отклик
    пробы считались именно так. Посадку под конкретное измерение задаёт
    вызывающий драйвер, см. run_bg_shield.py. gap_mm — ФИЗИЧЕСКИЙ зазор
    дно полости -> низ корпуса (P-013), не внутренняя координата lift.
    """
    a = ["pb=%.2f" % pb, "hxcav=%.2f" % HXCAV, "hycav=%.2f" % HYCAV,
         "hzcav=%.2f" % HZCAV, "seatfloor=%d" % SEAT_FLOOR]
    if not WITH_LID:
        a.append("nolid")
    if horiz:
        a.append("horiz=1")
    if gap_mm is not None:
        a.append("gap=%.2f" % gap_mm)
    if not with_vessel:
        a.append("novessel")
    return a


def trans_path(pb, s):
    return os.path.join(GRID_DIR, "trans_%s_pb%.0f.csv" % (s, pb))


def replay_path(pb, s):
    return os.path.join(GRID_DIR, "replay_%s_pb%.0f.csv" % (s, pb))


def run(args, log):
    with open(log, "w", encoding="utf-8") as lf:
        p = subprocess.run([SHIELDRUN_EXE] + args, cwd=BUILD, stdout=lf,
                           stderr=subprocess.STDOUT)
    return p.returncode


def ensure_grid(pb_list):
    os.makedirs(GRID_DIR, exist_ok=True)
    for pb in pb_list:
        for s in SERIES:
            tp = trans_path(pb, s)
            rp = replay_path(pb, s)
            specmac = frf.field_mac_path(s)
            if not os.path.exists(specmac):
                raise SystemExit("нет %s — сначала fit_room_field.py gen" % specmac)
            if not os.path.exists(tp):
                do_bias, nshell, impstep = bias_params(pb)
                args = (["trans"] + geom_args(pb)
                        + ["nprim=%d" % NPRIM, "specmac=%s" % specmac,
                           "out=%s" % tp])
                if do_bias:
                    args += ["bias", "impstep=%.2f" % impstep, "nshell=%d" % nshell]
                t0 = time.time()
                rc = run(args, tp + ".log")
                print("[trans] %s pb=%.0f  rc=%d  %.1fс" % (s, pb, rc, time.time() - t0),
                      flush=True)
                if rc != 0:
                    raise SystemExit("shieldrun trans упал, см. " + tp + ".log")
            if not os.path.exists(rp):
                # bias crystep=CRYSTEP — задача №13 (12.08.2026): без биасинга
                # на путь к кристаллу (case->caseAir->reflector->crystal)
                # взвешенных попаданий <1 уже на pb=50мм — узкое окно 662 кэВ
                # статистически бессмысленно. Провалидировано: pb=20 совпадает
                # с unbiased в пределах пуассоновского шума (0,0026 vs 0,0031),
                # на pb=100 даёт 1461 независимую выборку вместо 2 сырых при
                # той же ожидаемой сумме — не смещает среднее, снижает шум.
                t0 = time.time()
                # geom_args обязателен и здесь: стадия 2 воспроизводит записи
                # стадии 1 ВНУТРИ защиты, и если она построит другой домик
                # (например с крышкой, как велят умолчания ShieldGeom), то
                # точки старта окажутся не там, где были записаны, — молча,
                # без единой диагностики.
                rc = run(["replay"] + geom_args(pb)
                         + ["in=%s" % tp, "bias", "crystep=%.2f" % CRYSTEP,
                            "out=%s" % rp], rp + ".log")
                print("[replay] %s pb=%.0f  rc=%d  %.1fс" % (s, pb, rc, time.time() - t0),
                      flush=True)
                if rc != 0:
                    raise SystemExit("shieldrun replay упал, см. " + rp + ".log")


# --- мюонная компонента: ПЕРЕДЕЛАНО 13.08.2026 ---------------------------
#
# БЫЛО: радиус диск-источника рос с толщиной (R = rCav + pb + 40), а чтобы
# переиспользовать a_mu, подобранную на диске cosmicmu.cc (R=70), отклик
# домножался на (R_disk(pb)/70)². Это оказалось источником ЛОЖНОГО роста
# мюонной компоненты с толщиной: сырой отклик в окне 662 кэВ падал (13, 5,
# 12, 7, 6 отсчётов при pb=5..100), а множитель 1,84..7,37 переворачивал
# тренд. Развёртка по радиусу при pb=50 (с миром, расширенным под диск)
# показала, что физический ответ R²·eff выходит на полку только при R≥400 мм
# (64,0 / 68,4 / 64,8), тогда как рабочие 140 мм давали 28,8 — занижение
# в 2,2 раза, РАЗНОЕ на разных толщинах.
#
# СТАЛО: радиус диска ФИКСИРОВАН (700 мм, с запасом над полкой) и одинаков
# для всех толщин, поэтому поправка на площадь НЕ НУЖНА вовсе. a_mu
# перекалибрована на этой же геометрии по измеренному открытому фону в окне
# 2700-2790 кэВ и получила физический смысл — мюонов в секунду через диск;
# восстановленная плотность потока 0,0205 см⁻²·с⁻¹ против PDG-ориентира
# 0,0167 (было расхождение в 5,5 раза).
#
# Счёт в узком окне берётся через mu_smooth (гладкий континуум, плотность по
# всей статистике спектра), а не поканально: даже при 10⁷ первичных в окне
# 632-691 кэВ набирается порядка 19 отсчётов.
MU_GRID_DIR = os.path.join(BUILD, "musat2", "grid")
# ПЕРЕСМОТРЕНО 15.08.2026 на коробчатом домике (#SHIELD-8, #SHIELD-9).
#
# Радиус. Развёртка по 4e6 первичных на точку дала X = pi*R^2*eff:
#   R=700 -> 2,013 (4,4 %) | R=1000 -> 2,498 (5,6 %) | R=1400 -> 2,448 (7,9 %)
# То есть полка лежит между 1000 и 1400 (разница 0,050 +- 0,239), а прежние
# 700 мм ЗАНИЖАЛИ отклик на 19 % (2,9 sigma). Насыщение проверялось раньше на
# цилиндрическом домике, чей габарит был вдвое меньше нынешнего.
MU_R_DISK_MM = 1000.0
# Амплитуда. Прежние A_MU=315,2 +- 44,5 ПОДБИРАЛИСЬ под измеренный открытый
# фон, причём при неверном розыгрыше углов (#SHIELD-9: разыгрывалась
# интенсивность cos^2, а нужен поток через горизонтальную площадку, cos^3) —
# то есть подгонка впитала в себя ошибку розыгрыша. Оба основания отпали, и
# амплитуда взята АПРИОРНОЙ, ровно как концентрации бетона по UNSCEAR:
# интегральный поток мюонов на горизонтальную поверхность на уровне моря
# ~1 мюон/(см^2*мин) = 0,0167 см^-2 c^-1 (PDG). Свободных параметров в
# мюонной компоненте больше нет.
J_MU_PDG = 0.0167     # мюон/(см^2 * с), горизонтальная поверхность, уровень моря
A_MU = J_MU_PDG * math.pi * (MU_R_DISK_MM / 10.0) ** 2   # мюонов/с через диск
A_MU_SD = 0.0        # априорное число: своей статистической ошибки не имеет


def muon_path(pb):
    return os.path.join(MU_GRID_DIR, "mu2_pb%.0f.csv" % pb)


def mu_cps(pb, a_mu=A_MU):
    """-> (энергия, cps) — мюонная компонента. Поправки на площадь НЕТ:
    диск один и тот же на всех толщинах, a_mu откалибрована на нём же."""
    meta, hist = rcspec.read_spec(muon_path(pb))
    n = int(meta["N_primaries"])
    tag = meta.get("tag", "")
    if "rdisk_mm=%.2f" % MU_R_DISK_MM not in tag:
        raise SystemExit(
            "мюонный файл %s посчитан на ДРУГОМ диске (tag: %s) — a_mu=%.1f "
            "к нему не применима" % (muon_path(pb), tag, a_mu))
    per_muon = rcspec.fold(hist, "103") / n
    return np.arange(len(per_muon)) + 0.5, a_mu * per_muon


def sample_path(pb, nuc):
    return os.path.join(GRID_DIR, "sample_%s_pb%.0f.csv" % (nuc, pb))


def ensure_sample_grid(pb_list):
    os.makedirs(GRID_DIR, exist_ok=True)
    for pb in pb_list:
        for nuc in NUCLIDES_SAMPLE:
            sp = sample_path(pb, nuc)
            if os.path.exists(sp):
                continue
            t0 = time.time()
            rc = run(["sample"] + geom_args(pb)
                     + ["nuc=%s" % nuc, "nprim=%d" % NPRIM_SAMPLE,
                        "out=%s" % sp], sp + ".log")
            print("[sample] %s pb=%.0f  rc=%d  %.1fс" % (nuc, pb, rc, time.time() - t0),
                  flush=True)
            if rc != 0:
                raise SystemExit("shieldrun sample упал, см. " + sp + ".log")


def s_cps_per_bqkg(pb, nuc):
    """-> (энергия, cps на 1 Бк/кг) — собственная активность пробы сквозь
    защиту толщиной pb. mode=sample не биасируется (короткий путь), поэтому
    вес всегда 1 — читаем как обычный отклик."""
    meta, hist = rcspec.read_spec(sample_path(pb, nuc))
    n = int(meta["N_primaries"])
    per_decay = rcspec.fold(hist, "103") / n
    return np.arange(len(per_decay)) + 0.5, per_decay * SAMPLE_MASS_KG


def read_trans_header(path):
    """Все пары key=value из шапки trans, включая многоключевые строки.

    Шапка пишет и одиночные пары («# src_sha1 = …», «# N_primaries_stage1=…»),
    и строку сразу с десятью ключами (pb/cu/cd/e_in_keV/hxOut/hyOut/hzOut/
    zCav/S_mm2/lid). Прежний разбор резал строку по ПЕРВОМУ «=» и клал в
    meta мусор вида {"pb": "50.00 cu=0.00 cd=0.00 …"}: всё, кроме первого
    ключа строки, было недоступно — из-за чего площадь и приходилось
    пересчитывать здесь своей копией формулы.
    """
    meta = {}
    for line in open(path, encoding="utf-8"):
        if not line.startswith("#"):
            break
        body = re.sub(r"\s*=\s*", "=", line[1:].strip())
        for tok in body.split():
            if "=" not in tok:
                continue
            k, v = tok.split("=", 1)
            if k and v:
                meta[k] = v
    return meta


def b_room_gamma_cps(pb, amps):
    """-> (энергия[1 кэВ каналы], cps) — взвешенная сумма K+Ra+Th сквозь
    защиту толщиной pb (Cd=Cu=0), в имп/с на канал кристалла.

    Площадь наружной поверхности защиты берётся ГОТОВОЙ из шапки trans
    (ключ S_mm2) — единственный знаменатель для Ф=4N/S. Своей копии формулы
    здесь больше нет: у короба наружный габарит несимметричен по осям и
    зависит от того, строилась ли верхняя крышка, так что дубликат формулы
    разъезжается с C++ молча. Прежний пересчёт по (rCav+pb) описывал
    цилиндр и на коробчатом домике даёт заведомо неверный знаменатель.

    Площадь берётся ПОЛНОГО замкнутого габарита и при открытом верхе тоже:
    тождество Коши Ф=4N/S требует выпуклой охватывающей поверхности, а
    выпуклая оболочка защиты остаётся замкнутым параллелепипедом
    независимо от наличия свинцовой крышки (см. shieldrun.cc, BoxSurfaceGun).
    """
    n_ch = 3201
    total = np.zeros(n_ch)
    for s in SERIES:
        a, _sd = amps[s]
        if a is None:
            continue
        tp, rp = trans_path(pb, s), replay_path(pb, s)
        meta_t = read_trans_header(tp)
        n_stage1 = float(meta_t.get("N_primaries_stage1", NPRIM))
        _meta_r, hist_r = rcspec.read_spec(rp)   # взвешенный (Σw) отклик кристалла

        if "S_mm2" not in meta_t:
            raise SystemExit(
                "в шапке %s нет ключа S_mm2 — файл посчитан старым shieldrun "
                "(цилиндрическая защита). Удалите его и пересчитайте: молча "
                "подставить площадь короба к цилиндрическому прогону нельзя."
                % tp)
        if meta_t.get("lid", "on") != ("on" if WITH_LID else "off"):
            raise SystemExit(
                "%s посчитан с lid=%s, а сетка настроена на lid=%s — это "
                "РАЗНЫЕ домики (масса свинца и канал сверху), смешивать "
                "нельзя" % (tp, meta_t.get("lid"), "on" if WITH_LID else "off"))
        e_flu, flu = frf.read_wallfield(frf.wf_csv_path(s))
        fluence_total = flu.sum()
        area = float(meta_t["S_mm2"]) / 100.0   # мм² -> см², как fluence [см^-2]
        rate = fluence_total * area / 4.0   # первичных/с на ГРАНИЦЕ защиты (stage 1)
        t_run = n_stage1 / rate             # с, эквивалент N_primaries_stage1 первичных

        # hist_r — Σw по каналу за recs.size() событий replay = ПОЛНАЯ Σw
        # стадии 1 (без повтора). cps так же, как fl.model_cps_curve():
        # hist/t_run, свёрнуто с разрешением прибора.
        cps = rcspec.fold(hist_r, "103") / t_run
        n = min(len(cps), n_ch)
        total[:n] += a * cps[:n]
    return np.arange(n_ch) + 0.5, total


def main():
    pb_list = [float(x) for x in sys.argv[1:]] or [20.0, 50.0, 100.0]
    ensure_grid(pb_list)

    # переиспользуем уже подобранные амплитуды fit_lines.py — один источник
    # истины, не пересчитываем здесь.
    fit_result = fl.fit()
    amps = {"K": fit_result["K"], "Ra": fit_result["Ra"], "Th": fit_result["Th"]}

    print("\n=== B_room_gamma(t): K+Ra+Th сквозь защиту, Cd=Cu=0 ===")
    print("(мюон и остаточный дефицит континуума НЕ включены — см. план,")
    print(" задачи №6/№12; это ЧАСТИЧНАЯ гамма-компонента)\n")
    print("%6s %14s %14s" % ("pb,мм", "cps 20-2700", "cps 632-691(Cs137 окно)"))
    for pb in pb_list:
        e, cps = b_room_gamma_cps(pb, amps)
        m_all = (e >= 20) & (e < 2700)
        m_cs = (e >= 632.2) & (e < 691.2)
        print("%6.0f %14.6f %14.6e" % (pb, cps[m_all].sum(), cps[m_cs].sum()))
        out = os.path.join(GRID_DIR, "b_room_gamma_pb%.0f.csv" % pb)
        with open(out, "w", encoding="utf-8") as f:
            f.write("# B_room_gamma(E), pb=%.1f мм, cd=cu=0, K+Ra+Th по fit_lines.py\n" % pb)
            f.write("# a_K=%.2f a_Ra=%.2f a_Th=%.2f Бк/кг\n" % (amps["K"][0], amps["Ra"][0], amps["Th"][0]))
            f.write("E_keV,cps\n")
            for ei, ci in zip(e, cps):
                if ci > 0:
                    f.write("%.1f,%.6e\n" % (ei, ci))
        print("  ->", out)

    ensure_sample_grid(pb_list)
    K_BASELINE_BQKG = 250.0       # чистая ягода, приняли в постановке задачи
    CS_ILLUSTRATIVE_BQKG = 100.0  # ТОЛЬКО для иллюстрации таблицы ниже — на
                                  # странице это движок, не фиксированное число

    print("\n=== S_K40(t) и S_Cs137(t): собственная активность пробы ===")
    print("%6s %16s %16s" % ("pb,мм", "S_K40*250, cps(окно)", "S_Cs137*100, cps(окно)"))
    s_k40_win, s_cs_win = {}, {}
    for pb in pb_list:
        e_k, cps_k = s_cps_per_bqkg(pb, "K40")
        e_c, cps_c = s_cps_per_bqkg(pb, "Cs137")
        m_cs_k = (e_k >= 632.2) & (e_k < 691.2)
        m_cs_c = (e_c >= 632.2) & (e_c < 691.2)
        s_k40_win[pb] = K_BASELINE_BQKG * cps_k[m_cs_k].sum()
        s_cs_win[pb] = CS_ILLUSTRATIVE_BQKG * cps_c[m_cs_c].sum()
        print("%6.0f %16.6e %16.6e" % (pb, s_k40_win[pb], s_cs_win[pb]))
        for tag, e_, cps_, ref in (("K40", e_k, cps_k, K_BASELINE_BQKG),
                                   ("Cs137", e_c, cps_c, CS_ILLUSTRATIVE_BQKG)):
            out = os.path.join(GRID_DIR, "s_%s_per_bqkg_pb%.0f.csv" % (tag, pb))
            with open(out, "w", encoding="utf-8") as f:
                f.write("# S_%s(E), cps НА 1 Бк/кг, pb=%.1f мм, cd=cu=0\n" % (tag, pb))
                f.write("E_keV,cps_per_Bqkg\n")
                for ei, ci in zip(e_, cps_):
                    if ci > 0:
                        f.write("%.1f,%.6e\n" % (ei, ci))

    # Мюонная сетка НЕ гоняется отсюда: она считается отдельно, на
    # фиксированном насыщенном диске (musat2/run_mu_grid.ps1, 10⁷ первичных
    # на точку, ~9 мин каждая). a_mu берётся из перекалибровки на той же
    # геометрии, а НЕ из fit_lines.py — та подгонялась на обрезающем диске
    # cosmicmu.cc (R=70 мм) и завышена в 4,5 раза по плотности потока.
    mu_win = {}
    for pb in pb_list:
        mp = muon_path(pb)
        if not os.path.exists(mp):
            print("[mu] нет %s — мюонная компонента пропущена для pb=%.0f"
                  % (mp, pb), flush=True)
            continue
        e_mu, cps_mu = mu_cps(pb)
        # В узком окне считаем по СГЛАЖЕННОЙ плотности: поканальный счёт там
        # держится на единицах отсчётов даже при 10⁷ первичных (задача №17).
        # ВАЖНО: mu_smooth нелинеен по масштабу входа (укрупняет бины до
        # порога в ОТСЧЁТАХ), поэтому окно считается на СЫРЫХ свёрнутых
        # отсчётах, и только потом результат нормируется тем же множителем
        # a_mu/N, что и весь спектр. Подача cps сюда даёт завышение в ~2,7
        # раза — поймано 13.08.2026.
        meta_mu, hist_mu = rcspec.read_spec(muon_path(pb))
        n_mu = int(meta_mu["N_primaries"])
        folded_mu = rcspec.fold(hist_mu, "103")
        cnt, cnt_sd = mu_smooth.window_counts(
            np.arange(len(folded_mu)) + 0.5, folded_mu, 632.2, 691.2)
        norm = A_MU / n_mu
        val, sd = cnt * norm, cnt_sd * norm
        mu_win[pb] = val
        out = os.path.join(GRID_DIR, "b_room_mu_pb%.0f.csv" % pb)
        with open(out, "w", encoding="utf-8") as f:
            f.write("# B_room_mu(E), pb=%.1f мм, a_mu=%.1f мюон/с через диск "
                    "R=%.0f мм (ФИКСИРОВАН, поправки площади НЕТ)\n"
                    % (pb, A_MU, MU_R_DISK_MM))
            f.write("# окно 632-691 кэВ: %.4e cps (сглажено, ±%.1f%%)\n"
                    % (val, 100.0 * sd / val if val else float("nan")))
            f.write("E_keV,cps\n")
            for ei, ci in zip(e_mu, cps_mu):
                if ci > 0:
                    f.write("%.1f,%.6e\n" % (ei, ci))

    print("\n=== Сводка в окне Cs-137 (632-691 кэВ): где свинец перестаёт окупаться ===")
    print("%6s %14s %14s %14s %14s %14s" % ("pb,мм", "B_room_gamma", "B_room_mu",
                                             "B_room_ПОЛНЫЙ", "S_K40(250)", "S_Cs137(100)"))
    for pb in pb_list:
        e, cps = b_room_gamma_cps(pb, amps)
        m_cs = (e >= 632.2) & (e < 691.2)
        b_gamma = cps[m_cs].sum()
        b_mu = mu_win.get(pb, 0.0)
        print("%6.0f %14.6e %14.6e %14.6e %14.6e %14.6e"
              % (pb, b_gamma, b_mu, b_gamma + b_mu, s_k40_win[pb], s_cs_win[pb]))
    print("\n⚠️ ОГОВОРКИ (актуальны на 13.08.2026):")
    print("   1. B_room_ПОЛНЫЙ = гамма(K+Ra+Th) + мюон, БЕЗ поправки на остаточный")
    print("      дефицит континуума №12 (~25-35%): реальный B_room ВЫШЕ, то есть")
    print("      МДА здесь оптимистичнее действительности.")
    print("   2. Мюон: диск-источник ФИКСИРОВАН (700 мм, насыщение проверено),")
    print("      a_mu=%.1f мюон/с через диск. Восстановленная плотность потока" % A_MU)
    print("      0,0205 см^-2 c^-1 против ориентира PDG 0,0167 (отношение 1,23) —")
    print("      прежнее расхождение в 5,5 раза ЗАКРЫТО, причиной был обрезающий диск.")
    print("   3. НЕ закрыто: при pb=50 модель даёт в окне 2700-2790 кэВ 5,81e-4 cps")
    print("      против измеренных в домике 3,23e-4 (перебор 1,8x). Измерение")
    print("      показывает падение проникающей компоненты при экранировании,")
    print("      модель — постоянство. Задача №20.")
    print("   4. S_K40/S_Cs137 не включают обратное рассеяние от Cd/Cu (их нет в сетке).")


if __name__ == "__main__":
    main()
