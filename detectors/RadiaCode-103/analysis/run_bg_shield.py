# -*- coding: utf-8 -*-
"""Фон внутри свинцового домика: стадия 1 (сквозь свинец) + массированная
стадия 2 (внутри полости) веером процессов, затем сборка спектра.

ЗАЧЕМ ОТДЕЛЬНЫЙ ДРАЙВЕР. run_shield_grid.py гоняет replay одним процессом и
ровно один раз по каждой записи. На прежней цилиндрической полости (r=50,
h=180 мм) этого хватало. Коробчатый домик оператора — 150x150x385 мм, полость
просторнее на два порядка, кристалл 10x10x10 мм лежит на её дне, и квант,
влетевший через открытый верх с высоты 385 мм, попадает в него с вероятностью
порядка 1e-6: одиночный прогон на 110 тыс. записей дал 0,25 взвешенных
попаданий (замер 15.08.2026). Спектр по 0,25 отсчёта не строится.

ЧТО С ЭТИМ ДЕЛАЮТ. Записи стадии 1 стоят дорого (перенос сквозь 50 мм свинца),
а перенос ВНУТРИ полости по той же записи можно разыгрывать сколько угодно раз
с новым зерном — это независимые выборки редкого события. Отсюда два ключа
shieldrun, добавленные 15.08.2026: repeat= (сколько раз пройти входной файл) и
seed= (без него параллельные процессы повторяли бы одну историю). Драйвер
раскладывает работу на NPROC процессов, каждому своё зерно, и складывает
взвешенные спектры.

НОРМИРОВКА. Тождество Коши: изотропное поле с флюенсом Ф на выпуклой
охватывающей поверхности S даёт rate = Ф*S/4 входящих частиц в секунду. S
берётся ГОТОВОЙ из шапки trans (ключ S_mm2), Ф — сумма по wf_<серия>.csv
(wallfield, 1 Бк/кг). Эквивалентное время прогона T = N_stage1*Sum(repeat)/rate,
и cps = Sum(w)/T. Ни одного подогнанного параметра.

Запуск:  python run_bg_shield.py [pb] [nprim_trans] [nproc] [repeat]
         python run_bg_shield.py 50 1000000 10 100"""
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "common", "py"))
import rcspec  # noqa: E402
import fit_room_field as frf  # noqa: E402
import run_shield_grid as g  # noqa: E402  — geom_args/пути/парсер шапки

SERIES = ("K", "Ra", "Th")

# --- ПОСАДКА ПРИБОРА ПОД «Фон домик 23 дня» (P-008/P-012/P-013) -------------
# Со слов оператора (17.08): «лежал горизонтально на картонной коробочке…
# кристаллом вниз», уточнение — «он не прижат. зазор 25 мм (пустая картонная
# коробочка)». Сосуда Маринелли при этом замере не было.
#
# Держится ЗДЕСЬ, а не в аргументах запуска, ровно потому, что забытый ключ
# уже стоил партии счёта: P-012 — прогоны шли вертикально и с lift=20, никакая
# печать этого не показывала, поймал оператор. Теперь конфигурация зашита в
# драйвер и продублирована в ИМЕНИ каталога результатов.
#
# Картон в модель не вводится: тонкий низкоплотный материал на пути 25 мм
# воздуха. Это ДОПУЩЕНИЕ, не измерение.
HORIZ = True
GAP_MM = 25.0          # физический зазор дно полости -> низ корпуса
WITH_VESSEL = False
# P-009: при горизонтальной посадке replay-биасинг на путь к кристаллу выходит
# за допустимый G4 диапазон отношений важности (GeomBias1001,
# ipre_over_ipost=0,0625) — две попытки починки провалились, причина не
# установлена. Смещение при этом не гипотетическое: отношение полос 20-60 к
# 700-1500 составило 1,71 без биасинга против 2,00 с ним. Считаем без него;
# статистика добирается числом процессов и repeat.
BIAS_REPLAY = False


def seat_tag():
    # RC_BG_TAG — явный override имени каталога. Нужен, чтобы СБОРЩИКИ
    # (bg_budget.py, plot_bg_compare.py) могли построить сводку по УЖЕ готовому
    # прогону другой посадки, пока считается новый: подменять сам драйвер или
    # копировать пути в сборщики значило бы завести вторую копию правила
    # именования. На сам расчёт не влияет — геометрию задают HORIZ/GAP_MM.
    t = os.environ.get("RC_BG_TAG")
    if t:
        return t
    return "%s_gap%.0f_%s" % ("horiz" if HORIZ else "vert", GAP_MM,
                              "vess" if WITH_VESSEL else "novess")


def geom(pb):
    return g.geom_args(pb, horiz=HORIZ, gap_mm=GAP_MM, with_vessel=WITH_VESSEL)
# Активности ЕРН в бетоне помещения оператора — ИЗМЕРЕННЫЕ, не типовые.
#
# Получены 15.08.2026 из площадей характеристических линий ОТКРЫТОГО фона
# («Фон 7 дней без домика.xml», analysis/fit_lines.py): каждая линия делится на
# модельный отклик той же серии при 1 Бк/кг, дальше взвешенное среднее по своим
# линиям. Это ИЗМЕРЕНИЕ активности источника по его же излучению — метрология, а
# не подгонка фона: ни один параметр переноса или геометрии под измерение не
# крутится, подбирается только концентрация в стене, и подбирается по ПИКАМ, а
# сверяется потом континуум.
#
# ПЕРЕПОДОБРАНЫ 18.08 (шаг П1 плана). Прежние значения 315,92 / 7,65 / 9,97
# получены на СТАРОМ разрешении (FWHM(662)=8,4 %, реально 9,83 %) и на ШТАТНОЙ
# энергошкале прибора, которая врёт (rms 6,10 кэВ на якорях). Оба входа с тех пор
# исправлены -> окна ROI сместились и расширились, площади линий пересчитаны.
#
#   K-40   523,02 +- 185,52 Бк/кг  (по линии 1460,8)   было 315,92
#   Ra-226   7,35 +-   1,70        (609,3; 1764,5 в этом фоне не видна, даёт
#                                   отрицательное нетто — в среднее не входит)
#   Th-232  12,61 +-   3,00        (583,2; 911,2; 2614,5 за краем спектра)
#
# НЕЗАВИСИМАЯ ПРОВЕРКА (не входила в подгонку): активности подобраны по ЛИНИЯМ,
# а сверяется КОНТИНУУМ по широким полосам открытого фона. Модель/измерение:
# 20-100 -> 0,884 · 100-300 -> 1,034 · 300-600 -> 1,004 · 600-750 -> 1,092 ·
# 750-1400 -> 1,065. На старых активностях те же полосы давали 0,476...0,612.
# То есть «недобор континуума открытого фона в 1,5-2 раза» (#SHIELD-16) был
# АРТЕФАКТОМ старых активностей, а не физикой.
#
# ПОЧЕМУ НЕ UNSCEAR. Раньше здесь стояли типовые медианы 400/40/30. С ними
# модель перебирала измеренный ОТКРЫТЫЙ фон в 1,4-1,7 раза во всём диапазоне до
# 1,5 МэВ — то есть бетон у оператора втрое-впятеро чище типового по радию и
# торию. Продолжать считать домик на UNSCEAR значило бы носить эту ошибку
# дальше, компенсируя её другими слагаемыми.
PRIOR_BQKG = {"K": 523.02, "Ra": 7.35, "Th": 12.61}

# P-006 (16.08): в наборе были только полный интеграл и узкие линии — по такой
# сетке нельзя отличить «модель занижает мягкую часть» от «нет источника в
# конкретной полосе», и нельзя сравнить линию с континуумом рядом с ней.
# Разложение измерения по полосам (23 дня, домик): 20-60 даёт 13,7 % полного
# счёта, 60-100 — 30,0 %, 100-300 — 35,1 %, то есть 79 % сигнала лежит ниже
# 300 кэВ, а окно «полный 20-3000» этим низом и определяется.
WINDOWS = (
    ("полный 20-3000", 20.0, 3000.0),
    # полосы: где именно расходится модель
    ("полоса 20-60 (Pb-210?)", 20.0, 60.0),
    ("полоса 60-100 (K-фл Pb)", 60.0, 100.0),
    ("полоса 100-300", 100.0, 300.0),
    ("полоса 300-700", 300.0, 700.0),
    ("полоса 700-1500", 700.0, 1500.0),
    ("полоса 1500-3000", 1500.0, 3000.0),
    # линии и континуум РЯДОМ с ними — иначе peak-to-total не разделить
    ("K-40 1461", 1400.0, 1520.0),
    ("континуум 1150-1380", 1150.0, 1380.0),
    ("Cs-137 632-691", 632.2, 691.2),
    ("Tl-208 2615", 2560.0, 2670.0),
    ("проникающая 2700-2790", 2700.0, 2790.0),
)

# --- Космические мюоны (#SHIELD-13, 16.08.2026) -----------------------------
# До этой правки b_room_prior.csv содержал ТОЛЬКО ЕРН (K/Ra/Th): cosmicmu.cc
# считался отдельным кодом и с ними никогда не складывался. Отсюда провал
# модели в окнах Tl-208 2615 и «проникающая 2700-2790», где мюонный континуум
# и есть основной вклад: сверка 16.08 давала там отношение модель/измерение
# 0,007 и 0,000 — то есть модели там не было вовсе, а не «модель занижает».
#
# Нормировка БЕЗ ПОДГОНКИ. Поток мюонов на горизонтальную поверхность на
# уровне моря — стандартная величина обзора Cosmic Rays (PDG): ~1 мюон на
# см² в минуту, то есть 1/60 = 0,0167 см^-2 с^-1. Умноженный на площадь
# диска-источника, он даёт темп в секундах реального времени; свободного
# параметра здесь нет. Прежняя амплитуда A_MU=315,2 была ПОДОГНАНА под
# измерение при неверном розыгрыше угла (P-003) и отменена — использовать
# её после исправления розыгрыша нельзя.
MU_FLUX_PDG = 0.0167          # мюон/(см²·с), горизонтальная поверхность
MU_RDISK_MM = 1000.0          # радиус диска в прогонах musat_box (#SHIELD-8)
# Список НЕ фиксирован: мюонные прогоны добавляются пачками по мере набора
# статистики (16.08.2026: было 2 файла и 161 попадание на 2e6 первичных —
# в узких окнах это единицы отсчётов; стало 16 файлов и 1215 на 1.6e7,
# погрешность 7,9 % -> 2,9 %). Жёсткий кортеж молча игнорировал бы новые
# файлы, и статистика росла бы «на диске», а не в результате.
MU_GLOB = "mu2_r1000_*.csv"


def muon_cps():
    """
    Возвращает скорость счёта космических мюонов, дисперсию и число первичных мюонов.
    
    Скорость счёта и дисперсия складываются по всем файлам, подходящим под MU_GLOB.
    Дисперсия учитывает веса из sidecar-файлов, если они существуют.
    Радиус диска взят из параметров прогона: при R=700 мм отклик был занижен на 19 %,
    так как диск обрезал наклонные треки, полка начинается около 1000 мм.
    """
    # Путь к каталогу данных
    data_dir = os.path.join(frf.BUILD, "musat_box")
    
    # Инициализация массивов и счётчика
    hist = np.zeros(rcspec.NBINS)
    hist2 = np.zeros(rcspec.NBINS)
    nprim = 0.0
    
    # Обработка каждого файла, подходящего под MU_GLOB (16.08.2026: было
    # MU_FILES — жёсткий кортеж на 2 имени, молча игнорировал доборные
    # прогоны; NameError не всплыл раньше, т.к. до сих пор не вызывалась
    # против свежих файлов).
    import glob
    for p in sorted(glob.glob(os.path.join(data_dir, MU_GLOB))):
        if p.endswith(".sumw2.csv"):
            continue
        meta, h = rcspec.read_spec(p)
        nprim += float(meta["N_primaries"])
        hist += h
        
        # Проверка наличия sidecar-файла с весами
        weights_path = p + ".sumw2.csv"
        if os.path.exists(weights_path):
            _, w = rcspec.read_spec(weights_path)
            hist2 += w
        else:
            # Если весов нет, используем h как верхнюю оценку
            hist2 += h
    
    # Если первичных мюонов нет — возвращаем None
    if nprim <= 0:
        return (None, None, 0.0)
    
    # Расчёт радиуса диска в см
    r_cm = MU_RDISK_MM / 10.0
    
    # Темп мюонов через диск
    rate = MU_FLUX_PDG * np.pi * r_cm * r_cm
    
    # Коэффициент масштабирования
    k = rate / nprim
    
    # Возвращаем результат: скорость счёта, дисперсия и число первичных мюонов
    return (rcspec.fold(hist, "103") * k, rcspec.fold(hist2, "103") * k * k, nprim)


def out_dir(pb):
    # Конфигурация — в ИМЕНИ каталога: прогон с крышкой и без неё нельзя
    # складывать в один спектр, а перепутать их, имея одинаковые имена файлов,
    # проще всего.
    # Посадка прибора — тоже в ИМЕНИ (P-012): прогон лёжа и стоя нельзя
    # складывать в один спектр, а прежние вертикальные файлы лежат в каталоге
    # со старым именем и молча подхватились бы как готовые.
    d = os.environ.get("RC_BG_DIR") or os.path.join(
        g.BUILD, "bg_shield",
        "pb%.0f_%s_%s" % (pb, "lid" if g.WITH_LID else "nolid", seat_tag()))
    os.makedirs(d, exist_ok=True)
    return d


def trans_file(pb, s):
    return os.path.join(out_dir(pb), "trans_%s.csv" % s)


def replay_file(pb, s, i):
    return os.path.join(out_dir(pb), "replay_%s_%02d.csv" % (s, i))


def spawn(args, log):
    lf = open(log, "w", encoding="utf-8")
    return subprocess.Popen([g.SHIELDRUN_EXE] + args, cwd=g.BUILD, stdout=lf,
                            stderr=subprocess.STDOUT), lf


def wait_all(jobs, what, strict):
    """strict=True — любой сбой останавливает расчёт; False — сбойные процессы
    просто теряются.

    Для стадии 2 нестрого СОЗНАТЕЛЬНО: процессы там независимы, каждый пишет
    свой файл целиком в конце прогона, и упавший не оставляет ни строки —
    потерять его значит потерять кусок статистики, а не испортить результат.
    Ронять из-за этого многочасовой веер нельзя, тем более что replay с
    биасингом уже показал ПЛАВАЮЩИЙ access violation (#SHIELD-5, 15.08.2026:
    один и тот же вызов на одном и том же файле дважды упал и трижды прошёл).
    Сколько процессов дошло — видно в печати, и в нормировку идёт фактическая
    Sum(repeat) по существующим файлам, а не ожидаемая.
    """
    bad = 0
    for (p, lf, tag) in jobs:
        rc = p.wait()
        lf.close()
        if rc != 0:
            bad += 1
            print("  [СБОЙ] %s rc=%d" % (tag, rc), flush=True)
    print("[%s] завершено, сбоев %d из %d" % (what, bad, len(jobs)), flush=True)
    if bad and strict:
        raise SystemExit("%s: часть процессов упала, см. логи рядом с выходом"
                         % what)


def stage1(pb, nprim):
    """trans по трём сериям разом — они независимы и делят машину поровну."""
    todo = [s for s in SERIES if not os.path.exists(trans_file(pb, s))]
    if not todo:
        print("[stage1] все три файла уже есть", flush=True)
        return
    do_bias, nshell, impstep = g.bias_params(pb)
    jobs = []
    t0 = time.time()
    for s in todo:
        tp = trans_file(pb, s)
        args = (["trans"] + geom(pb)
                + ["nprim=%d" % nprim, "specmac=%s" % frf.field_mac_path(s),
                   "out=%s" % tp])
        if do_bias:
            args += ["bias", "impstep=%.2f" % impstep, "nshell=%d" % nshell]
        p, lf = spawn(args, tp + ".log")
        jobs.append((p, lf, tp))
        print("[stage1] пуск %s -> %s" % (s, tp), flush=True)
    wait_all(jobs, "stage1", strict=True)
    print("[stage1] %.1f мин" % ((time.time() - t0) / 60), flush=True)


def stage2(pb, nproc, repeat):
    """
    веер replay-процессов; несколько попыток, потому что replay с importance biasing падает с access violation примерно на каждом шестом прогоне, файл при этом не создаётся, поэтому упавший процесс просто перезапускается с другим зерном; биасинг обязателен, так как даёт N_eff 16671 против 125 без него при согласующемся среднем.
    """
    t0 = time.time()
    for attempt in range(3):
        jobs = []
        all_done = True
        for si, s in enumerate(SERIES):
            for i in range(nproc):
                rp = replay_file(pb, s, i)
                if os.path.exists(rp):
                    continue
                all_done = False
                seed = 100000 * (attempt + 1) + 1000 * (si + 1) + i + 1
                args = ["replay"] + geom(pb) + [
                    "in=%s" % trans_file(pb, s),
                    "seed=%d" % seed,
                    "repeat=%d" % repeat,
                    "out=%s" % rp
                ]
                if BIAS_REPLAY:
                    args += ["bias", "crystep=%.2f" % g.CRYSTEP]
                p, lf = spawn(args, rp + ".log")
                jobs.append((p, lf, rp))
        if all_done:
            print("[stage2] попытка %d: все файлы уже есть" % (attempt + 1), flush=True)
            break
        else:
            print("[stage2] попытка %d, запущено процессов: %d" % (attempt + 1, len(jobs)), flush=True)
            wait_all(jobs, "stage2/попытка %d" % (attempt + 1), strict=False)
    print("[stage2] %.1f мин" % ((time.time() - t0) / 60), flush=True)


def combine(pb, s):
    """-> (cps на 1 Бк/кг серии, эффективное число попаданий) для серии."""
    meta_t = g.read_trans_header(trans_file(pb, s))
    if "S_mm2" not in meta_t:
        raise SystemExit("нет S_mm2 в шапке %s" % trans_file(pb, s))
    n_stage1 = float(meta_t["N_primaries_stage1"])
    area_cm2 = float(meta_t["S_mm2"]) / 100.0

    _e, flu = frf.read_wallfield(frf.wf_csv_path(s))
    rate = flu.sum() * area_cm2 / 4.0        # первичных/с на границе защиты

    hist = np.zeros(rcspec.NBINS)
    hist2 = np.zeros(rcspec.NBINS)
    total_repeat = 0.0
    nfile = 0
    for i in range(1000):
        rp = replay_file(pb, s, i)
        if not os.path.exists(rp):
            continue
        meta_r = g.read_trans_header(rp)
        rep = float(meta_r.get("repeat", 1))
        _m, h = rcspec.read_spec(rp)
        hist += h
        # Sum(w^2) пишется shieldrun отдельным файлом рядом; если его нет
        # (старый прогон) — считаем историю единичной, это верхняя оценка
        # точности, и её лучше видеть, чем молча остаться без ошибки.
        p2 = rp + ".sumw2.csv"
        if os.path.exists(p2):
            _m2, h2 = rcspec.read_spec(p2)
            hist2 += h2
        else:
            hist2 += h
        total_repeat += rep
        nfile += 1
    if not nfile:
        raise SystemExit("нет ни одного replay-файла для серии %s" % s)

    t_run = n_stage1 * total_repeat / rate
    cps = rcspec.fold(hist, "103") / t_run
    # Дисперсия — из Sum(w^2), который shieldrun пишет рядом со спектром. Без
    # неё нельзя сказать «сошлось» или «не сошлось»: при биасинге число
    # ОТСЧЁТОВ в канале ничего не говорит о точности, вес истории может
    # отличаться на порядки. Свёртка с разрешением перемешивает каналы, так
    # что fold(Sum w^2) — приближение дисперсии, годное для интеграла по
    # окну ШИРЕ разрешения (все окна ниже такие) и не годное поканально.
    var = rcspec.fold(hist2, "103") / (t_run * t_run)
    neff = (hist.sum() ** 2 / hist2.sum()) if hist2.sum() > 0 else 0.0
    print("  %-3s файлов %2d, Sum(repeat)=%.0f, событий %.3g, Sum(w)=%.4g, "
          "N_eff=%.0f, T=%.4g с"
          % (s, nfile, total_repeat, n_stage1 * total_repeat, hist.sum(),
             neff, t_run), flush=True)
    return cps, var


def main():
    pb = float(sys.argv[1]) if len(sys.argv) > 1 else 50.0
    nprim = int(sys.argv[2]) if len(sys.argv) > 2 else 1_000_000
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    repeat = int(sys.argv[4]) if len(sys.argv) > 4 else 100

    print("=== фон внутри домика, Pb %.0f мм, полость %.0fx%.0fx%.0f, крышка %s ==="
          % (pb, 2 * g.HXCAV, 2 * g.HYCAV, 2 * g.HZCAV,
             "есть" if g.WITH_LID else "НЕТ"))
    print("stage1 nprim=%d, stage2 %d проц. x repeat=%d на серию"
          % (nprim, nproc, repeat), flush=True)
    print("посадка: %s, зазор %.0f мм, сосуд %s, replay-биасинг %s"
          % ("ГОРИЗОНТАЛЬНО" if HORIZ else "вертикально", GAP_MM,
             "есть" if WITH_VESSEL else "НЕТ", "вкл" if BIAS_REPLAY else "ВЫКЛ"),
          flush=True)
    print("каталог: %s\n" % out_dir(pb), flush=True)

    stage1(pb, nprim)
    stage2(pb, nproc, repeat)

    print("\n[сборка]", flush=True)
    e = np.arange(rcspec.NBINS) + 0.5
    total = np.zeros(rcspec.NBINS)
    per_series = {}
    per_var = {}
    for s in SERIES:
        cps, var = combine(pb, s)
        per_series[s] = cps
        per_var[s] = var
        total += PRIOR_BQKG[s] * cps

    # #SHIELD-13/W-007: muon_cps() была написана, но НЕ подключалась к сумме —
    # b_room_prior.csv уходил без мюонной компоненты. Подключаю здесь.
    cps_mu, var_mu, nmu = muon_cps()
    if cps_mu is not None:
        total += cps_mu
        print("[мюоны] файлов учтено, N_primaries=%.4g, добавлено в total" % nmu, flush=True)
    else:
        print("!! [мюоны] нет данных musat_box/%s — total БЕЗ мюонов" % MU_GLOB, flush=True)

    print("\n=== B_room(E) по АПРИОРНЫМ концентрациям, без подгонки ===")
    print("бетон: " + ", ".join("%s=%.0f Бк/кг" % (k, v)
                                for k, v in PRIOR_BQKG.items()))
    # G2 (Codeaudit, 16.08): «сумма» ОБЯЗАНА браться из того же `total`, что
    # пишется в CSV, иначе таблица и файл расходятся — таблица показывала бы
    # спектр БЕЗ мюонов и читалась бы как «фикс не сработал». Мюоны выведены
    # отдельной колонкой, чтобы вклад был виден, а не спрятан в сумме.
    print("\n%-24s %13s %13s %13s %13s %13s"
          % ("окно", "K-40", "Ra-226", "Th-232", "мюоны", "сумма, cps"))
    for name, lo, hi in WINDOWS:
        m = (e >= lo) & (e < hi)
        parts = [PRIOR_BQKG[s] * per_series[s][m].sum() for s in SERIES]
        # Серии независимы (разные прогоны, разные зёрна) — дисперсии
        # складываются, амплитуда входит квадратом. G1 (Codeaudit): мюонная
        # дисперсия — независимое слагаемое; без неё погрешность занижена
        # именно в тех окнах, где мюоны и есть основной вклад.
        var_sum = sum(PRIOR_BQKG[s] ** 2 * per_var[s][m].sum() for s in SERIES)
        if var_mu is not None:
            var_sum += var_mu[m].sum()
        sd = np.sqrt(var_sum)
        tot = total[m].sum()
        # Codeaudit, 17.08: при отсутствии мюонных данных печатать ПРОЧЕРК, а не
        # 0.0 — ноль неотличим от «учтены и дали ноль в этом окне», если таблицу
        # скопируют отдельно от предупреждения выше.
        mu_str = ("%13.5e" % cps_mu[m].sum()) if cps_mu is not None else "%13s" % "—"
        print("%-24s %13.5e %13.5e %13.5e %s %13.5e  +- %.1f%%"
              % (name, parts[0], parts[1], parts[2], mu_str, tot,
                 100.0 * sd / tot if tot > 0 else float("nan")))

    out = os.path.join(out_dir(pb), "b_room_prior.csv")
    with open(out, "w", encoding="utf-8") as f:
        f.write("# B_room(E) внутри домика, БЕЗ подгонки\n")
        f.write("# pb=%.1f полость=%.0fx%.0fx%.0f lid=%s\n"
                % (pb, 2 * g.HXCAV, 2 * g.HYCAV, 2 * g.HZCAV,
                   "on" if g.WITH_LID else "off"))
        f.write("# " + " ".join("a_%s=%.0f" % (k, v)
                                for k, v in PRIOR_BQKG.items()) + " Bq/kg\n")
        # G3 (Codeaudit, 16.08): cps_total = ЕРН + мюоны, поэтому баланс
        # cps_total = cps_K+cps_Ra+cps_Th НЕ выполняется без колонки cps_mu.
        # Провенанс мюонной нормировки пишется здесь же, чтобы потребитель
        # (в т.ч. я сам через полгода) не гадал, дефект это или замысел.
        if cps_mu is not None:
            f.write("# muons: flux=%.4f cm^-2 s^-1 (PDG), r_disk=%.0f mm, "
                    "N_primaries=%.6g\n" % (MU_FLUX_PDG, MU_RDISK_MM, nmu))
        else:
            f.write("# muons: НЕТ ДАННЫХ, cps_mu=0, cps_total только ЕРН\n")
        f.write("E_keV,cps_total,cps_K,cps_Ra,cps_Th,cps_mu\n")
        for i in range(rcspec.NBINS):
            if total[i] > 0:
                f.write("%.1f,%.6e,%.6e,%.6e,%.6e,%.6e\n"
                        % (e[i], total[i], PRIOR_BQKG["K"] * per_series["K"][i],
                           PRIOR_BQKG["Ra"] * per_series["Ra"][i],
                           PRIOR_BQKG["Th"] * per_series["Th"][i],
                           cps_mu[i] if cps_mu is not None else 0.0))
    print("\n->", out)


if __name__ == "__main__":
    main()