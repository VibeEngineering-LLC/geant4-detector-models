import os
import sys
import numpy as np
from scipy.optimize import nnls
import json
import io
import math
import argparse

SPECTRAVIBE_ROOT = os.environ.get("SPECTRAVIBE_ROOT")
if not SPECTRAVIBE_ROOT:
    raise RuntimeError("Переменная окружения SPECTRAVIBE_ROOT не установлена")
sys.path.insert(0, os.path.join(SPECTRAVIBE_ROOT, "scripts"))
from gamma.io.lsrm_spe import read_lsrm_spe

PASSPORT = {"Am241": 4106.6, "Ti44chain": 1808.8,
            "Cs137chain": 1599.8, "Eu152": 2025.2}

# Параметр левого экспоненциального хвоста функции отклика (peak-image ЛСРМ).
# 0 — чистый гаусс, прежнее поведение. Значение 0,75 определено данными на
# линии RC-103 (профиль χ²/n.d.f. с минимумом, 28,2 → 5,57; PLAN-METHOD-1.md
# R.1). Для Гамма-1С своё значение НЕ измерялось — берётся как априор и
# проверяется прогоном; переопределяется ключом --tail.
TAIL_T = 0.75

# Сдвиг нумерации таблицы пиков LSRM относительно индекса массива отсчётов: position_ch ≈ индекс + 1,0 ± 0,15
# (самоаудит 11.09.2026 на 5 симметричных пиках, см. recalibrate_energy ниже). Полином энергии файла задан в
# нумерации таблицы. ЕДИНСТВЕННЫЙ источник значения (#CH-1, мера GEANT4:W-104) — потребители импортируют отсюда.
LSRM_CH_OFFSET = 1.0

# --- Ступень 2 §31.B: сгенерировано qwen3-coder:30b по спеке
# scripts/_spec_recalib.md. Перекалибровка шкалы по реперам самого
# спектра — требование #CAL-0: штатная шкала файла в зоне америция
# смещена на 0,375 ПШПВ (Am-241) и 0,485 ПШПВ (Eu-152 121,8), а порог
# подгонки заново — 0,25 ПШПВ. Выше 200 кэВ шкала файла точна.
RECAL_REFS = [59.541, 121.78, 244.70, 344.279, 511.0, 661.657, 778.90, 964.08, 1157.02, 1408.01]

def recalibrate_energy(spec, verbose=True, extra_refs=None, ch_offset=0.0):
    """extra_refs — внешние реперы [(ИНДЕКС массива этого спектра, E_кэВ)] ниже порога поиска пиков
    прибора (перенос K-рентгена доноров комплекта, analysis/low_scale_donors.py).
    ch_offset — сдвиг нумерации таблицы пиков LSRM относительно индекса массива отсчётов
    (11.09.2026, самоаудит А: position_ch ≈ индекс + 1,0 ± 0,15 на 5 симметричных пиках);
    из центров таблицы вычитается, многочлен файла вычисляется в нумерации таблицы."""
    table = (getattr(spec, "extras", None) or {}).get("lsrm_peaks_table") or []
    off = float(ch_offset)

    used_pairs = [(float(c), float(E)) for c, E in (extra_refs or [])]
    for e_ref in RECAL_REFS:
        threshold = max(6.0, 0.04 * e_ref)
        best_row = None
        min_diff = float('inf')
        for row in table:
            diff = abs(row["energy_keV"] - e_ref)
            if diff < min_diff and diff < threshold:
                min_diff = diff
                best_row = row
        if best_row is not None:
            used_pairs.append((best_row["position_ch"] - off, e_ref))
    
    if len(used_pairs) < 3:
        raise SystemExit(f"Найдено только {len(used_pairs)} реперных линий, требуется минимум 3")
    
    used_pairs.sort()
    channels = np.array([p[0] for p in used_pairs], dtype=float)
    energies = np.array([p[1] for p in used_pairs], dtype=float)
    if np.any(np.diff(channels) <= 0) or np.any(np.diff(energies) <= 0):
        # B-2: неверные внешние реперы давали немонотонную шкалу и −inf без отказа.
        raise SystemExit("ОТКАЗ: реперы шкалы не возрастают строго по каналу и энергии: %s"
                         % list(zip(channels.round(3), energies.round(3))))
    coef = np.polyfit(channels, energies, 2)

    if verbose:
        print(f"Перекалибровка шкалы: {len(used_pairs)} реперов, кусочная "
              f"интерполяция (невязка в реперах 0), вне "
              f"{channels[0]:.0f}…{channels[-1]:.0f} кан. — снизу "
              f"{'наклон первого отрезка' if extra_refs else 'квадратичная, сдвинутая к первому реперу'}, "
              f"сверху шкала файла, сдвинутая к реперу (D-018)")
    
    # ⚠ Был квадратичный polyfit — он сглаживал ошибку ГЛОБАЛЬНО и в зоне
    # америция оставлял сдвиг: максимум измерения приходился на 57,0 кэВ,
    # максимум модели — на 59,9, то есть модель правее на 2,9 кэВ. Отсюда
    # «недобор» слева от 60 вместе с избытком при 80–86: перекос шкалы, а не
    # эффективность. Реперы известны точно, между ними шкала гладкая — берём
    # кусочно-линейную интерполяцию, в реперах невязка ровно ноль.
    # D-018 (11.09.2026): выше последнего репера — не квадратичная
    # экстраполяция по реперам (к концу диапазона она уводила шкалу на
    # −72 кэВ), а форма шкалы самого файла, сдвинутая к последнему реперу.
    # Сумм-пики каскадов (1157+511) реперами не служат: их сдвиг — проверка
    # модели непропорциональности, а не калибровка.
    # Ниже первого репера (59,5 кэВ) опорных линий в файле нет — остаётся
    # квадратичная экстраполяция. Форма файла там проверена и отвергнута
    # измерением: χ²/ν 195,5 → 262,6, Am-241 −18 %, максимумы зоны
    # америция разошлись до 55,2 против 61,2 кэВ.
    f_file = lambda x: float(spec.channel_to_energy(x + off))

    def energy_from_channel(c):
        c = float(c)
        if channels[0] <= c <= channels[-1]:
            return float(np.interp(c, channels, energies))
        if c < channels[0] and extra_refs:
            # С внешними низкими реперами — продолжение наклона первого отрезка
            # (без скачка на первом репере; до порога окна подгонки 2–3 канала).
            k = (energies[1] - energies[0]) / (channels[1] - channels[0])
            return float(energies[0] + k * (c - channels[0]))
        if c < channels[0]:
            # P-029 (11.09.2026, самоаудит Б B-1): квадратичный многочлен по всем реперам через
            # первый репер НЕ проходит — был скачок −1,73 кэВ прямо на пике Am-241 и немонотонная
            # шкала. Та же форма, сдвинутая к первому реперу: непрерывно.
            return float(energies[0] + np.polyval(coef, c) - np.polyval(coef, channels[0]))
        return float(energies[-1] + f_file(c) - f_file(channels[-1]))

    energy_from_channel.pairs = list(zip(channels.tolist(), energies.tolist()))  # (индекс массива, кэВ) — для шкалы в свете
    return (energy_from_channel, len(used_pairs))


def read_template(path):
    """Читает шаблон из CSV-файла и возвращает гистограмму и параметры."""
    with open(path, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    
    # ⚠ Между шапкой и данными лежит таблица распределения по числу
    # комптоновских рассеяний, строки-комментарии с «#» и пустые строки.
    # Первая генерация обрывала разбор на первой строке без запятой.
    header = {}
    data_start = None
    for i, line in enumerate(lines):
        if line.startswith("bin_keV"):
            data_start = i + 1
            break
        if not line or line.startswith("#") or "," not in line:
            continue
        key, value = line.split(',', 1)
        if key and not key[0].isdigit():   # числовые ключи — таблица n_compt
            header[key] = value
    if data_start is None:
        raise SystemExit(f"ОТКАЗ: в {path} нет строки-маркера данных 'bin_keV'")
    for need in ('n_events_processed', 'npsm_enabled'):
        if need not in header:
            raise SystemExit(f"ОТКАЗ: в шапке {path} нет ключа '{need}'")

    n_events_processed = int(header['n_events_processed'])
    npsm_enabled = int(header['npsm_enabled'])
    
    hist = {}
    for line in lines[data_start:]:
        if not line or line.startswith("#"):
            continue
        energy, count_edep, count_light = map(float, line.split(','))
        key = energy
        value = count_light if npsm_enabled == 1 else count_edep
        if value > 0:
            hist[key] = value
    
    return hist, n_events_processed, npsm_enabled

def make_fwhm(points_csv):
    """Создаёт функцию для вычисления ПШПВ по энергии."""
    # Принимаем и путь, и уже открытый файловый объект: самопроверка подаёт
    # StringIO, рабочий режим — путь. Раньше функция умела только путь, и её
    # собственный тест падал на TypeError.
    points = []
    fh = points_csv if hasattr(points_csv, "read") else io.open(
        points_csv, "r", encoding="utf-8")
    try:
        for line in fh:
            if not line.strip() or line.startswith('E_keV'):
                continue
            parts = line.split(',')
            if len(parts) < 2:
                continue
            # Берём только первые ТРИ колонки: четвёртая (source) — текст,
            # и map(float, ...) по всей строке на ней падал.
            e = float(parts[0])
            fwhm = float(parts[1])
            points.append((e, fwhm))
    finally:
        if fh is not points_csv:
            fh.close()
    if len(points) < 2:
        raise SystemExit("ОТКАЗ: в файле точек ПШПВ меньше двух строк — "
                         "интерполировать не по чему")
    
    points.sort(key=lambda x: x[0])
    energies = np.array([p[0] for p in points])
    fwhms = np.array([p[1] for p in points])
    
    log_e = np.log(energies)
    log_fwhm = np.log(fwhms)
    
    def fwhm(E):
        if E < energies[0]:
            # Экстраполяция по логарифмическим точкам
            slope = (log_fwhm[1] - log_fwhm[0]) / (log_e[1] - log_e[0])
            return np.exp(log_fwhm[0] + slope * (np.log(E) - log_e[0]))
        elif E > energies[-1]:
            # Экстраполяция по логарифмическим точкам
            slope = (log_fwhm[-1] - log_fwhm[-2]) / (log_e[-1] - log_e[-2])
            return np.exp(log_fwhm[-1] + slope * (np.log(E) - log_e[-1]))
        else:
            # Линейная интерполяция по логарифмическим точкам
            idx = np.searchsorted(energies, E)
            if idx == 0:
                return fwhms[0]
            elif idx == len(energies):
                return fwhms[-1]
            else:
                # ⚠ Здесь была ошибка генерации: числитель брался в
                # логарифме, а знаменатель — в ЛИНЕЙНЫХ энергиях, отчего
                # интерполяция давала 8,93 вместо 30 в контрольной точке.
                # Обе величины обязаны быть в одной шкале.
                l1, l2 = log_e[idx-1], log_e[idx]
                f1, f2 = log_fwhm[idx-1], log_fwhm[idx]
                ratio = (np.log(E) - l1) / (l2 - l1)
                return float(np.exp(f1 + ratio * (f2 - f1)))
    
    return fwhm

def broaden_ch(hist, n_events, e_of_ch, fwhm, file_e, blur=1.0, sub=20):
    """Свёртка в КАНАЛАХ (11.09.2026, #AMT-3): отложение E → канал c(E) по шкале пробы e_of_ch,
    ширина — σ в каналах = blur·ПШПВ(E)/2,3548 / наклон шкалы ФАЙЛА file_e в этом канале (точки ПШПВ
    сняты в кэВ шкалы файла). Излом калибровочной шкалы не перекашивает пик, как у broaden()."""
    e = np.asarray(e_of_ch, dtype=float); n = len(e); ch = np.arange(n, dtype=float)
    fc = (np.arange(n * sub) + 0.5) / sub - 0.5            # центры под-бинов в каналах
    out = np.zeros(n * sub)
    for E, cnt in hist.items():
        if E < 1.0 or E < e[0] or E > e[-1]:
            continue                                        # без отложения или вне шкалы
        c = float(np.interp(E, e, ch))
        s = blur * fwhm(E) / 2.3548 / max(file_e(c + 0.5) - file_e(c - 0.5), 1e-9)
        lo, hi = np.searchsorted(fc, c - 12 * s), np.searchsorted(fc, c + 6 * s)
        z = (fc[lo:hi] - c) / s
        g = np.exp(-0.5 * z * z) if TAIL_T <= 0 else np.where(
            z >= -TAIL_T, np.exp(-0.5 * z * z), np.exp(np.clip(TAIL_T * z + 0.5 * TAIL_T ** 2, -700.0, 0.0)))
        if g.sum() > 0:
            out[lo:hi] += cnt * g / g.sum()
    return out.reshape(n, sub).sum(axis=1) / n_events


def broaden(hist, n_events, ch_edges, fwhm):
    """Раскладывает линии по каналам с учётом ширины."""
    # Мелкая сетка для интегрирования
    fine_grid = np.arange(0, 3300, 0.25)
    
    # Суммарная гистограмма
    broadened = np.zeros_like(fine_grid)
    
    for energy, count in hist.items():
        sigma = fwhm(energy) / 2.3548
        if sigma <= 0:
            continue
        
        # Ядро отклика: гаусс с ЛЕВЫМ экспоненциальным хвостом (peak-image
        # ЛСРМ). Форма взята у донора контура — `rcspec.fold`,
        # detectors/RadiaCode-103/analysis/rcspec.py:257-288: слева от z = −T
        # гаусс сшивается с экспонентой exp(T·z + T²/2), ядро нормируется
        # суммой, поэтому площадь сохраняется при любом T.
        # Зачем: реальный пик NaI не гауссов — измеритель отказывался мерить
        # Am-241 с χ²/dof = 81 «в окне не одна гауссова линия». В зоне
        # 57…78 кэВ три линии (Am-241 59,5; Ti-44 67,9 и 78,3) сливаются при
        # ширине 6–8 кэВ, и разделение их площадей зависит от ФОРМЫ крыла:
        # чистый гаусс отдаёт левому соседу меньше, чем есть на самом деле.
        z = (fine_grid - energy) / sigma
        if TAIL_T <= 0:
            g = np.exp(-0.5 * z * z)
        else:
            g = np.where(z >= -TAIL_T, np.exp(-0.5 * z * z),
                         np.exp(np.clip(TAIL_T * z + 0.5 * TAIL_T * TAIL_T,
                                        -700.0, 0.0)))
        s = np.sum(g)
        if s > 0:
            broadened += count * g / s
    
    # Интегрируем по каналам.
    # ⚠ Здесь было ДВА дефекта генерации. Первый: границы канала брались
    # searchsorted с side='left' и side='right', отчего точка на границе
    # попадала в оба соседних канала и сумма превышала единицу ровно на
    # одну точку сетки в канале (при шаге 0,25 и канале 3,0 — избыток
    # 1/12 = 8,3 %, измерено самопроверкой). Второй: кумулятивная сумма
    # пересчитывалась заново ВНУТРИ цикла по каналам.
    # Приём тот же, что в проверенном broaden_and_rebin контура: считаем
    # кумуляту один раз и снимаем её значения на границах интерполяцией.
    cum = np.cumsum(broadened)
    v = np.interp(ch_edges, fine_grid, cum, left=0.0, right=cum[-1])
    result = np.diff(v) / max(n_events, 1)
    
    return result

def selftest():
    """Выполняет тесты."""
    # Тест 1: make_fwhm
    points_csv = io.StringIO("E_keV,fwhm_keV,d_fwhm_keV,source\n100,10,0.5,measured\n1000,30,1.0,measured")
    fwhm_func = make_fwhm(points_csv)
    
    if abs(fwhm_func(100) - 10) > 0.01:
        print("SELFTEST FAILED")
        return 1
    if abs(fwhm_func(1000) - 30) > 0.01:
        print("SELFTEST FAILED")
        return 1
    if abs(fwhm_func(316.23) - 17.32) > 0.01:
        print("SELFTEST FAILED")
        return 1
    
    # Тест 2: broaden с одной линией и постоянной шириной
    hist = {662.0: 1000.0}
    n_events = 1000
    ch_edges = np.arange(0, 3000, 3.0)
    
    def fwhm_const(E):
        return 50.0
    
    col = broaden(hist, n_events, ch_edges, fwhm_const)
    total = np.sum(col)
    if abs(total - 1.0) > 0.01:
        print("SELFTEST FAILED")
        return 1
    
    # Тест 3: проверка влияния ширины
    def fwhm_double(E):
        return 100.0
    
    col2 = broaden(hist, n_events, ch_edges, fwhm_double)
    max1 = np.max(col)
    max2 = np.max(col2)
    if max1 <= max2:
        print("SELFTEST FAILED")
        return 1
    
    # Тест 3a: АБСОЛЮТНАЯ ширина. Тест 3 сравнивает две ширины между собой и
    # к общему масштабу слеп: домножение sigma на любой множитель проходило
    # незамеченным (поймано мутацией 11.09.2026), а это ровно тот дефект,
    # который делает пик в зоне америция вдвое шире или уже реального.
    # Меряем ПШПВ свёрнутой линии прямо по массиву — на половине высоты.
    # Сетка мелкая (1 кэВ): на канале 3 кэВ дискретность съедала разницу и
    # неверный делитель 2,0 вместо 2,3548 (ошибка на 18 %) проходил тест.
    fine = np.arange(0, 3000, 1.0)
    centers = 0.5 * (fine[:-1] + fine[1:])
    col_w = broaden({662.0: 1e6}, 1, fine, lambda E: 50.0)
    half = np.max(col_w) / 2.0
    above = centers[col_w >= half]
    if len(above) < 2:
        print("T3a FAIL: свёрнутая линия уже одного канала — ширина не задана")
        return 1
    got_fwhm = float(above[-1] - above[0])
    if abs(got_fwhm - 50.0) > 3.0:  # допуск втрое уже прежнего
        print(f"T3a FAIL: заданная ПШПВ 50 кэВ, получена {got_fwhm:.1f} кэВ")
        return 1

    # Тест 4: nnls
    A = np.array([[1, 0], [0, 1], [1, 0], [0, 1]])
    b = np.array([3, 7, 3, 7])
    coef, _ = nnls(A, b)
    if abs(coef[0] - 3) > 0.01 or abs(coef[1] - 7) > 0.01:
        print("SELFTEST FAILED")
        return 1
    
    print("SELFTEST OK")
    return 0

def main():
    parser = argparse.ArgumentParser()
    # Не required: с ключом --selftest файлы не нужны вовсе, а argparse
    # с required=True не даёт запустить самопроверку без них. Отсутствие
    # аргументов в рабочем режиме проверяется ниже явной ошибкой.
    parser.add_argument('--spe')
    parser.add_argument('--bg')
    parser.add_argument('--templates', nargs=4)
    parser.add_argument('--fwhm-points')
    parser.add_argument('--lo', type=float, default=40.0)
    parser.add_argument('--hi', type=float, default=1500.0)
    parser.add_argument('--json')
    parser.add_argument('--selftest', action='store_true')
    parser.add_argument('--tail', type=float, default=None,
                        help='параметр левого хвоста отклика (0 — чистый гаусс)')
    parser.add_argument('--recalibrate', action='store_true',
                        help='перекалибровать шкалу по реперам спектра (#CAL-0)')
    
    args = parser.parse_args()
    
    if args.selftest:
        return selftest()

    if args.tail is not None:
        globals()['TAIL_T'] = float(args.tail)
        print(f'Параметр хвоста отклика: T = {TAIL_T}')

    # Рабочий режим: аргументы обязательны. Проверка явная, потому что
    # required=True снят ради --selftest (см. выше).
    missing = [n for n, v in (("--spe", args.spe), ("--bg", args.bg),
                              ("--templates", args.templates),
                              ("--fwhm-points", args.fwhm_points)) if not v]
    if missing:
        raise SystemExit("ОТКАЗ: не заданы обязательные аргументы: "
                         + ", ".join(missing))

    # Расчёт ведёт ЯДРО, а не своя копия: здесь жила вторая реализация той же подгонки
    # (свой NNLS, свои веса, свой χ²), и после переноса критерия D-020 в ядро она молча
    # оставила командную строку на ПРЕЖНЕМ критерии — одни файлы давали разные активности
    # в CLI и в выгрузке. Импорт локальный: ядро импортирует этот модуль на верхнем уровне.
    import mix_unfold_core as core
    templates = []
    for t in args.templates:
        name, path = t.split('=', 1)
        templates.append((name, path))
    # Постановка сохранена дословно: фон канал в канал, свёртка в энергии, без шкалы света —
    # так CLI считала и раньше. Меняется ровно критерий подгонки, не расчёт.
    r = core.unfold(args.spe, args.bg, templates, args.fwhm_points, lo=args.lo, hi=args.hi,
                    recalibrate=args.recalibrate, tail=None, bg_energy_of_ch=None, verbose=True)
    spec, e, y, sel = r["spec"], r["e"], r["y"], r["sel"]
    bg_scaled, sigma, cols = r["bg_scaled"], r["sigma"], list(r["cols"])
    coef, activities, model = r["coef"], r["activities"], r["model"]
    chi2, ndof, n_refs = r["chi2"], r["ndof"], r["n_refs"]
    print("критерий %s: χ²/ν = %.3f (дисперсия критерия), χ²_ref/ν = %.3f (единая метрика)"
          % (r["crit"], chi2 / ndof, r["chi2_ref"] / ndof))
    print("вторая мера E1: невязка формы %.6f, активности %s"
          % (r["e1"]["tv"], np.array2string(r["e1"]["activities"], precision=3)))
    
    # Вклады в chi2 по полосам
    bands = [(40, 90), (90, 200), (200, 500), (500, 1000), (1000, 1500)]
    print("Вклады в chi2 по полосам:")
    for band in bands:
        mask = (e >= band[0]) & (e <= band[1])
        if np.sum(mask) == 0:
            continue
        chi2_band = np.sum(((model - y)[mask] / sigma[mask])**2)
        print(f"  {band[0]}-{band[1]} кэВ: {chi2_band / chi2 * 100:.1f}%")
    
    # Отношения расчёт/паспорт
    print("\nОтношения расчёт/паспорт:")
    for i, (name, _) in enumerate(templates):
        ratio = activities[i] / PASSPORT[name]
        print(f"  {name}: {ratio:.3f}")
    
    # Зона Америки (50-70 кэВ)
    mask_america = (e >= 50) & (e <= 70)
    if np.sum(mask_america) > 0:
        sum_meas = np.sum(y[mask_america])
        sum_model = np.sum(model[mask_america])
        ratio = sum_model / sum_meas if sum_meas != 0 else 0
        print(f"\nЗона Америки (50-70 кэВ):")
        print(f"  Сумма измерения: {sum_meas:.1f}")
        print(f"  Сумма модели: {sum_model:.1f}")
        print(f"  Отношение: {ratio:.3f}")
        
        contributions = []
        for i, (name, _) in enumerate(templates):
            contrib = np.sum(cols[i][mask_america])
            contributions.append(contrib)
        
        total = sum(contributions)
        if total > 0:
            print("  Вклады шаблонов:")
            for i, (name, _) in enumerate(templates):
                percent = contributions[i] / total * 100
                print(f"    {name}: {percent:.1f}%")
    
    # Запись в JSON
    if args.json:
        result = {
            "activities": dict(zip([t[0] for t in templates], activities)),
            "chi2": chi2,
            "ndof": ndof,
            "model": model.tolist(),
            "data": y.tolist()
        }
        with open(args.json, 'w', encoding='utf-8') as f:
            # numpy-типы (int64/float64) json не сериализует — приводим к
            # питоновским. allow_nan=False: NaN в результате означает, что
            # подгонка вырождена, и такой файл лучше не писать вовсе.
            def _plain(o):
                if isinstance(o, (np.integer,)):
                    return int(o)
                if isinstance(o, (np.floating,)):
                    return float(o)
                if isinstance(o, np.ndarray):
                    return o.tolist()
                raise TypeError(f"не сериализуется: {type(o).__name__}")
            json.dump(result, f, ensure_ascii=False, indent=2,
                      default=_plain, allow_nan=False)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
