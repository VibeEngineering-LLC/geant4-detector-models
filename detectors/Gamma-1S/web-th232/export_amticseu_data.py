# -*- coding: utf-8 -*-
"""Экспорт данных для страницы Mix AmTiCsEu -- Тест 3 скилла
geant4-spectrum-pipeline (11.08.2026), финальный тест конвейера после
Th-232 и Ra-226.

ГЛАВНОЕ СТРУКТУРНОЕ ОТЛИЧИЕ от export_data.py/export_ra226_data.py:
источник -- НЕ цепочка векового равновесия, а смесь ЧЕТЫРЁХ НЕЗАВИСИМЫХ
радионуклидов (искусственный калибровочный источник "№SRC-04" LSRM).
Поэтому:
  - НЕТ chain_<id>.csv (нечего физически "прогонять цепочкой" -- четыре
    нуклида не связаны распадом друг с другом). Каждая группа --
    ОТДЕЛЬНЫЙ Geant4-прогон (macros/decay_amticseu_isotopes.mac,
    configs/amticseu.yaml.geant4.isotope_runs), используемый НАПРЯМУЮ
    как столбец матрицы плана NNLS -- без промежуточного шага "доля от
    templ_total" (тот шаг в Th-232/Ra-226 нужен, чтобы получить
    физически точную ФОРМУ полной цепочки; здесь каждая группа УЖЕ
    физически полна сама по себе -- nucleusLimits Z-диапазон в макросе
    захватывает весь нужный каскад одной группы, включая дочерние
    состояния с короткими T1/2, см. комментарии в конфиге).
  - Метод 1: ЧЕТЫРЕ независимые свободные амплитуды + фон (5 столбцов
    NNLS), не одна амплитуда на "ветвь". Это то же самое обобщение,
    которое export_ra226_data.py уже применяло как ДИАГНОСТИЧЕСКИЙ
    разрез (radon_check: 2 амплитуды родитель/дочерние) -- здесь оно
    используется как ОСНОВНАЯ модель, потому что у источника и правда
    нет единой активности-родителя.
  - Метод 2 (F_B/TCS): библиотека и сумм-пики читаются ОБЫЧНЫМ
    механизмом ed.GAMMA_LIBRARY/ed.SUM_PEAKS (конфиг-драйв, тот же код,
    что и Th-232/Ra-226) -- geant4-spectrum-pipeline уже нуклид-
    агностичен на этом уровне. F_B данные добавлены в ОБЩИЙ файл
    data/conversion_coeff_sum_peak_levels.csv (211 строк, нуклиды этого
    источника дописаны 11.08.2026, тот же файл, что несёt обе цепочки
    Th-232/Ra-226 -- "conversion_coeff... Список нуклидов там расширен",
    см. докстринг export_data.py._sum_peaks_with_fb).
  - НЕТ варианта "full library" (полная непороговая библиотека) --
    упрощение, сознательно принятое по времени; методология описана в
    remarks, страница явно говорит "только отобранная библиотека".
  - НЕТ разложения К-рентгена дочерних отдельной сущностью (XRAY) --
    упрощение: шаблон iso_<key>.csv используется как есть (ПОЛНАЯ
    физика, включая рентген атомной релаксации и рентген-флуоресценцию
    защиты -- Geant4 их всё равно считает физически корректно в общей
    гистограмме; отдельная разбивка "сколько именно от рентгена" --
    только для наглядности отчёта Th-232/Ra-226, не физическая
    необходимость).
  - НЕТ radon_check/split3_check/radon_dpr_vs_ra_check -- диагностика
    утечки радона специфична для цепочки Ra-226, здесь нет аналога.

Запуск:
    python export_amticseu_data.py
Переменные окружения: G4MODELS_BUILD_GAMMA_1S (сетка отклика/шаблоны),
G4MODELS_AMTICSEU_BG_SPE (фон, приватный путь, в репозитории его нет),
SPECTRAVIBE_ROOT (ридер .spe).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault("G4MODELS_SOURCE_CONFIG",
                      os.path.join(HERE, "configs", "amticseu.yaml"))
os.environ.setdefault("PYTHONIOENCODING", "utf-8")

sys.path.insert(0, HERE)
import export_data as ed  # noqa: E402  -- ed._CFG уже наш конфиг (env выше)
import export_ra226_data as erd  # noqa: E402  -- переиспользуем run_method2/load_col/factory_fwhm_keV/read_pair (нуклид-агностичны, см. их докстринги)

# Тег сетки метода 2 -- risn379, ро=1,0 (см. drivers/run_grid.py risn379 1.0,
# build/Gamma-1S/grid/risn3791.00_E*.csv, 24 узла).
GRID_MATRIX = "risn379"
GRID_RHO = 1.0
GRID_PATTERN = "%s%.2f_E*.csv" % (GRID_MATRIX, GRID_RHO)


def read_pair():
    """Копия erd.read_pair() -- функция уже нуклид-агностична (читает
    ed._CFG["source"]["measured_sample_spe_rel"]/measured_background_
    spe_env"), но модуль export_ra226_data.py не выставляет наружу этот
    вызов иначе как через собственный main(); дублировать саму функцию
    дешевле, чем бороться с областью видимости."""
    return erd.read_pair()


def run_method1_independent(e, ch_edges, T, y_sel, bgm, sel, groups):
    """Метод 1 для источника из НЕЗАВИСИМЫХ нуклидов -- каждая группа
    получает СВОЙ свободный столбец NNLS напрямую из iso_<key>.csv, без
    промежуточного шага "доля от templ_total" (в Th-232/Ra-226 он нужен
    для физической точности формы полной цепочки; здесь каждая группа
    УЖЕ физически полна сама по себе, см. докстринг модуля).

    `groups` -- [(key, label_ru, label_en, color, br, note), ...] --
    тот же кортеж-формат, что ed.NUCS/erd.NUCS, br здесь уже применён
    ВНУТРИ Geant4-прогона там, где это возможно (Ti44chain/Eu152/Am241 --
    br=1,0, полная физика в одном прогоне); Cs137chain -- br=0,9457
    (доля распадов Cs-137, заселяющая Ba-137m, прогон спавнит сам
    изомер) домножается здесь программно, как и было задумано в конфиге.
    """
    hist_iso = {}
    missing = []
    for key, ru, en, color, br, note in groups:
        p = os.path.join(ed.BUILD, "iso_%s.csv" % key)
        if not os.path.isfile(p):
            missing.append(key)
            continue
        hist_iso[key] = ed.load_hist(p)
    if missing:
        raise SystemExit(
            "Нет МК-шаблонов групп AmTiCsEu: %s\n"
            "Запустить: cd %s && ./g1s.exe decay_amticseu_isotopes.mac vessel %.2f %s"
            % (missing, ed.BUILD, GRID_RHO, GRID_MATRIX))

    by_group_raw = {}
    for key, ru, en, color, br, note in groups:
        hist, N = hist_iso[key]
        by_group_raw[key] = ed.broaden_and_rebin(hist, N, ch_edges, True) * br

    neg = {k: float(v.min()) for k, v in by_group_raw.items()
          if float(v.min()) < -1e-12}
    if neg:
        raise SystemExit("отрицательные значения в шаблонах метода 1: %s" % neg)

    keys = [k for k, *_ in groups]
    # ОТКАЗ от отдельной сущности "XRAY" (директива оператора 11.08.2026:
    # "рентгеновские линии нуклидов всегда включай в шаблон нуклида") --
    # характеристический рентген (атомная релаксация, флуоресценция
    # защиты) физически сидит ВНУТРИ шаблона своей группы, iso_<key>.csv
    # несёт его целиком -- ничего вычитать/выделять не нужно, шаблон
    # используется как есть, той же группы цвет на графике.
    cols = [by_group_raw[k][sel] * T for k in keys] + [bgm]
    coef, dcoef, chi2, ndof, model = ed.fit_amplitudes(y_sel, cols, ed.SYS_FLOOR)

    result = {}
    stack = {}
    for i, key in enumerate(keys):
        result[key] = {"A_Bq": float(coef[i]), "dA_Bq": float(dcoef[i])}
        stack[key] = (by_group_raw[key] * float(coef[i]) * T).tolist()
    bg_amp, d_bg_amp = float(coef[-1]), float(dcoef[-1])

    A_rn = np.stack(cols, axis=1)
    sig = np.sqrt(np.maximum(y_sel, 1.0) + (ed.SYS_FLOOR * model) ** 2)
    cond = float(np.linalg.cond((A_rn.T * (1.0 / sig ** 2)) @ A_rn))

    return {
        "groups": result, "bg_amplitude": bg_amp, "d_bg_amplitude": d_bg_amp,
        "chi2": chi2, "ndof": ndof, "chi2_ndof": chi2 / ndof,
        "cond_number": cond, "n_channels_fit": int(sel.sum()),
    }, stack, {"template_decays": [{"nuclide": ru, "n": hist_iso[key][1]}
                                   for key, ru, en, color, br, note in groups]}


def main():
    meas, bg = read_pair()
    T = meas["live_s"]

    # ── энергетическая шкала: ЗАВОДСКАЯ + ЛОКАЛЬНАЯ поправка низа ────────
    # ИСПРАВЛЕНО 11.08.2026 (директива оператора «правь калибровку»).
    # Валидация ДО этой правки (argmax по 4 изолированным линиям 344,28/
    # 511/661,657/964,08, сверка с lsrm_peaks_table) дала ≤1,1 кэВ --
    # эти 4 точки ПОДТВЕРЖДЕНЫ и используются ниже как есть (заводская
    # верна там). Но узел 40-90 кэВ (Am-241, доминирующая линия 59,541
    # кэВ) в неё не входил -- и там нашёлся реальный локальный сдвиг,
    # ДВУМЯ независимыми методами, сошедшимися в КАНАЛЬНОМ пространстве
    # (не только в кэВ, где расхождение калибровок само по себе могло бы
    # объяснить сходство): (1) argmax по фон-вычтенным отсчётам, окно
    # ~0,55×ПШПВ, параболическое уточнение вершины -- канал 24,0-24,16
    # (расчёт воспроизводим: ../analysis/am241_anchor_check.py, добавлен
    # 11.08.2026 внешним аудитом -- прежде это утверждение не имело
    # сохранённого артефакта);
    # (2) скриншот оригинала SpectraLine (оператор), Am-241 подписан на
    # 54,32/54,63 кэВ -- при обратном пересчёте ЧЕРЕЗ ЭТУ ЖЕ заводскую
    # квадратику (коэффициенты из шапки .spe, общие с нашим e_of_ch) даёт
    # канал 24,49-24,60. Оба метода physически видят ОДИН И ТОТ ЖЕ канал
    # (разброс <0,6 канала ≈ 1,5 кэВ) -- сходимость в канальном
    # пространстве значит совпадение независимо от того, как каждый
    # метод переводит канал в энергию. Взят канал 24,3 (среднее).
    #
    # Поправка -- ТОТ ЖЕ рецепт, что и Ra-226 (export_ra226_data.py,
    # apply_xray_anchor_correction_sample, ra226-remarks.md §16):
    # квадратика (3 коэффициента) через 5 точек (канал=наблюдённый
    # максимум, кэВ=истинная энергия) -- 1 новая (24,3 -> 59,541) + 4
    # уже подтверждённые (заводская=истина в этих точках, поэтому канал
    # берётся ОБРАТНЫМ пересчётом заводской квадратики, не собственным
    # argmax -- у argmax на СЛАБЫХ линиях есть свой систематический
    # перекос от комптоновского континуума соседних сильных линий,
    # numerically проверено и ОТБРОШЕНО в этой же сессии: наивный argmax
    # ±10 кэВ на 964,057 кэВ давал уход −12,7 кэВ, хотя lsrm_peaks_table
    # для этой же линии подтверждает +1,10 -- разный метод, разная
    # надёжность, доверять узкому наивному окну на слабой линии нельзя).
    # Проверено (5 точек, невязка после рефита, ПЕРЕСЧИТАНО 11.08.2026
    # с финальным якорным каналом 24,3 -- прежняя версия комментария
    # считала с черновым 24,08, расхождение поймано внешним аудитом):
    # 59,541 −0,35; 344,279 +1,08; 511,0 −0,27; 661,657 −0,76; 964,057
    # +0,30 кэВ -- везде <1,1 кэВ, не хуже исходной валидации. Экстраполяция квадратики
    # ДАЛЬШЕ анкеров расходится с заводской (как и у Ra-226) -- поэтому
    # поправка ПЛАВНО сшивается с заводской (smoothstep) и уходит в ноль
    # к каналу 150 (~425 кэВ по заводской), выше -- чистая заводская, тот
    # же принцип "не трогать то, что не проверено".
    #
    # Фон в этой полосе НЕ поправлен: сознательно -- background-side
    # вклад в 40-90 кэВ пренебрежим (счёт фона ~100-160 против ~130000
    # net у образца, <0,15%), да и якоря (Am-241) в фоне физически нет
    # (источника не было при наборе фона) -- поправлять фон нечем и
    # незачем, в отличие от Ra-226 (там якорь -- К-рентген Pb защиты,
    # присутствует и в фоне тоже).
    _AM241_SAMPLE_ANCHOR_CH = np.array([24.3, 122.633, 178.884,
                                        229.605, 331.096])
    _AM241_SAMPLE_ANCHOR_KEV = np.array([59.541, 344.279, 511.0,
                                         661.657, 964.057])
    _AM241_BLEND_LO, _AM241_BLEND_HI = 60.0, 150.0

    def _apply_am241_anchor_correction(e_of_ch, ch):
        A = np.vstack([_AM241_SAMPLE_ANCHOR_CH ** k for k in range(3)]).T
        coef, *_ = np.linalg.lstsq(A, _AM241_SAMPLE_ANCHOR_KEV, rcond=None)
        e_new = coef[0] + coef[1] * ch + coef[2] * ch ** 2
        t = np.clip((ch - _AM241_BLEND_LO)
                    / (_AM241_BLEND_HI - _AM241_BLEND_LO), 0.0, 1.0)
        w_old = t * t * (3 - 2 * t)
        return (1.0 - w_old) * e_new + w_old * e_of_ch, coef

    e_factory = meas["e_of_ch"]
    e, _am241_coef = _apply_am241_anchor_correction(
        e_factory, np.arange(meas["n_channels"], dtype=float))
    energy_correction = {
        "applied": True,
        "method": "am241_anchor_local_quadratic_blend_to_factory",
        "anchor_kev": 59.541,
        "anchor_ch": 24.3,
        "sample_anchor_channels": _AM241_SAMPLE_ANCHOR_CH.tolist(),
        "sample_anchor_kev": _AM241_SAMPLE_ANCHOR_KEV.tolist(),
        "blend_lo_ch": _AM241_BLEND_LO, "blend_hi_ch": _AM241_BLEND_HI,
        "note": "Заводская шкала верна выше канала %d (~%.0f кэВ); ниже "
                "заменена локальным квадратиком по 5 якорям (1 новый -- "
                "Am-241 59,541 кэВ, канал 24,3 -- + 4 подтверждённых "
                "заводских), сшитым smoothstep-переходом. Подробности -- "
                "amticseu-remarks.md §11."
                % (int(_AM241_BLEND_HI),
                   float(np.polyval(meas["coefs"][::-1], _AM241_BLEND_HI))),
    }

    # ── энергетическая шкала ФОНА -- полный пересчёт по 7-линейной
    # методике ЛСРМ для ЕРН/фоновых спектров NaI (11.08.2026, замечание
    # оператора «фон не откалиброван», ИСПРАВЛЕНО повторно: первая
    # версия сдвигала только нуль по ОДНОМУ якорю (K-40) с наивным
    # широким окном 1350-1560 кэВ -- дала центроид 1444,46 (Δ=-16,36) и
    # была принята за чистый сдвиг нуля по прецеденту Ra-226. Оператор
    # указал прямо: "просто двигать - ошибка", применить методику
    # SpectraVibe (references/01_metadata_calibration.md, "Background-
    # only anchor heuristic" -- набор из 7 природных линий: Pb-212/214
    # 238,63, Pb-214 295,21/351,93, Tl-208+аннигиляция 511,0, Bi-214
    # 609,32/1120,29/1764,49, Ac-228 911,20, K-40 1460,822, Tl-208
    # 2614,511). С окном, отмасштабированным по ПШПВ прибора (0,6×ПШПВ,
    # не наивные широкие 100+ кэВ), K-40 САМ ПО СЕБЕ даёт другой,
    # гораздо более скромный центроид (1455,52, Δ=-5,30) -- прежнее
    # число было артефактом слишком широкого окна, та же ошибка, что уже
    # ловилась в этой сессии на образце (964,057 кэВ, наивный argmax).
    # По 10 природным линиям (7 основных + 511 отдельно от аннигиляции)
    # взвешенная (вес~sqrt(чистый счёт)) подгонка: линейная (degree=1)
    # RMS=1,81 кэВ, квадратичная (degree=2) RMS=1,86 кэВ (ХУЖЕ при
    # лишнем параметре -- кривизна шумовая, не физическая). Взят
    # линейный рефит ЦЕЛИКОМ (не локальная поправка -- якоря покрывают
    # весь диапазон 238-2614 кэВ). Ниже добавлен 11-й якорь (К-рентген
    # Pb защиты, 77,0 кэВ) по отдельной директиве оператора -- итоговый
    # RMS с ним 2,12 кэВ (хуже, чем без него, ожидаемо для статистически
    # слабого якоря, принят не по статистическим основаниям, а по прямой
    # инструкции). Воспроизводимо: analysis/bg_seven_line_anchor_check.py.
    #
    # 11-й якорь -- К-рентген Pb защиты (директива оператора "и ХРИ
    # свинца защиты на фоне есть... всегда проверяй и используй").
    # Физика подтверждена IAEA (Bi-207 EC->Pb-207, тот же механизм, что
    # даёт K-рентген Pb в спектре Ra-226 -- Kα1 74,970/36,5%, Kα2
    # 72,805/21,7%, Kβ1 84,986/12,5%, Kβ2 87,301/3,8%, интенсивностно-
    # взвешенное среднее кластера 76,66 кэВ). Центроид в ЭТОМ фоне НЕ
    # сходится с шириной окна: 75,88 (±5,8 кэВ) -> 72,39 (±10) -> 71,09
    # (±15) -> 69,46 (±20) -> 69,36 (±25 кэВ) -- дрейфует, в отличие от
    # всех остальных 10 якорей (слабый/размытый сигнал, низкая энергия,
    # возможна примесь K-рентгена йода NaI ~28-33 кэВ). По статистике
    # такой якорь не прошёл бы отбор -- ВКЛЮЧЁН по прямой директиве
    # оператора "просто прими для ХРИ свинца защиты 77" (тот же приём,
    # что якорь Ra-226 §16, округлённое табличное значение вместо
    # нестабильного измеренного). Канал взят НЕ по недостоверному
    # центроиду, а по argmax в узком (0,6×ПШПВ) окне вокруг
    # факторного предсказания канала для 77,0 кэВ (=34,4) -- канал=34,
    # netsum=77,2 (наименьший из 11 якорей, вес в МНК соответственно
    # мал). Воспроизводимо: analysis/bg_seven_line_anchor_check.py.
    _BG_ANCHORS_CH = np.array([34.0, 86.335, 106.048, 124.351, 177.959,
                               211.149, 313.183, 382.371, 495.287,
                               598.614, 882.950])
    _BG_ANCHORS_KEV = np.array([77.0, 238.63, 295.21, 351.93, 511.0,
                                609.32, 911.20, 1120.29, 1460.822,
                                1764.49, 2614.511])
    _BG_ANCHORS_NETSUM = np.array([77.2, 231.0, 124.4, 210.4, 286.6,
                                   278.4, 195.8, 214.8, 892.2, 255.3,
                                   194.7])
    _bg_w = np.sqrt(_BG_ANCHORS_NETSUM)
    _bg_A = np.vstack([np.ones_like(_BG_ANCHORS_CH), _BG_ANCHORS_CH]).T
    _bg_coef_new, *_ = np.linalg.lstsq(_bg_A * _bg_w[:, None],
                                       _BG_ANCHORS_KEV * _bg_w, rcond=None)
    _bg_c0_new, _bg_c1_new = float(_bg_coef_new[0]), float(_bg_coef_new[1])
    bg_e_corrected = (_bg_c0_new
                      + _bg_c1_new * np.arange(bg["n_channels"], dtype=float))

    bg_on_meas = np.interp(e, bg_e_corrected, bg["counts"].astype(float),
                           left=0.0, right=0.0)
    bg_scale_time = T / bg["live_s"]
    bg_scaled = bg_on_meas * bg_scale_time

    ch_edges = np.concatenate((
        [e[0] - 0.5 * (e[1] - e[0])],
        0.5 * (e[:-1] + e[1:]),
        [e[-1] + 0.5 * (e[-1] - e[-2])],
    ))

    # ── ПШПВ(E): заводская калибровка (та же директива "сам не
    # калибруй", что у Ra-226) ──────────────────────────────────────────
    if not meas.get("fwhm_coefs"):
        raise SystemExit(
            "в шапке .spe нет заводской ПШПВ-калибровки "
            "(stored_fwhm_calibration пуст) -- посчитать нечем.")
    fwhm_k, fwhm_p, fwhm_fit_rms_pct = erd.fit_power_law_to_factory_fwhm(
        meas["fwhm_coefs"], meas["fwhm_model"])
    ed.FWHM_LAW.update({"kind": "power", "k": fwhm_k, "p": fwhm_p})
    _ref_pts = (59.541, 121.7817, 344.2785, 511.0, 661.657, 964.057, 1157.022, 1408.013)
    fwhm_cal = {
        "source": "заводская (LSRM, полином в шапке .spe, z=sqrt(E))",
        "coefs": meas["fwhm_coefs"], "model": meas["fwhm_model"],
        "k": fwhm_k, "p": fwhm_p, "fit_rms_pct": fwhm_fit_rms_pct,
        "fwhm662_law": fwhm_k * 661.657 ** fwhm_p,
        "fwhm662_cs": ed.FWHM662,
        "res662_pct": 100.0 * fwhm_k * 661.657 ** fwhm_p / 661.657,
        "reference_points": [
            {"E_keV": E,
             "fwhm_factory_keV": erd.factory_fwhm_keV(meas["fwhm_coefs"],
                                                       meas["fwhm_model"], E),
             "fwhm_power_law_keV": fwhm_k * E ** fwhm_p}
            for E in _ref_pts],
    }
    print("ПШПВ: заводская калибровка, аппроксимация степенным законом "
          "k=%.4f p=%.4f (невязка аппроксимации %.2f%%), ПШПВ(662)=%.1f кэВ"
          % (fwhm_k, fwhm_p, fwhm_fit_rms_pct, fwhm_cal["fwhm662_law"]))

    eps_peak = ed.make_eps_peak_interp(os.path.join(ed.BUILD, "grid"),
                                       GRID_PATTERN)
    resp = ed.make_full_response(os.path.join(ed.BUILD, "grid"), ch_edges,
                                 True, eps_peak, GRID_PATTERN)

    sel = (e >= ed.E_FIT_LO) & (e <= ed.E_FIT_HI)
    y_sel = meas["counts"][sel].astype(float)
    bgm = bg_scaled[sel]

    groups = ed.NUCS  # [(key, ru, en, color, br, note), ...] из конфига
    keys = [k for k, *_ in groups]

    # ── метод 1: 4 независимые амплитуды + фон ──────────────────────────
    m1_result, m1_stack, m1_meta = run_method1_independent(
        e, ch_edges, T, y_sel, bgm, sel, groups)

    cfg_components = {c["key"]: c for c in ed._CFG["passport"]["components"]}
    mass_kg = ed._CFG["passport"]["mass_g"] / 1000.0
    days = ed._CFG["passport"]["days_pass_to_meas"]

    passport_by_group = {}
    for key in keys:
        c = cfg_components[key]
        decay_f = ed.decay_factor_years(c["half_life_years"], days)
        A_pass = c["bq_per_kg"] * mass_kg * decay_f
        dA_pass = A_pass * c["unc_pct"] / 100.0
        passport_by_group[key] = {
            "A_Bq": A_pass, "dA_Bq": dA_pass,
            "Bq_per_kg": c["bq_per_kg"], "unc_pct": c["unc_pct"],
            "half_life_years": c["half_life_years"], "decay_factor": decay_f,
        }
        m1_result["groups"][key]["A_over_passport"] = (
            m1_result["groups"][key]["A_Bq"] / A_pass)
    print("метод 1 (4 независимые амплитуды): chi2/ndof=%.2f" % m1_result["chi2_ndof"])
    for key in keys:
        g = m1_result["groups"][key]
        print("  %-12s A=%.1f+-%.1f Бк  ratio=%.3f"
              % (key, g["A_Bq"], g["dA_Bq"], g["A_over_passport"]))

    # ── метод 2: F_B/TCS-библиотека, ОДНА амплитуда НА ГРУППУ (не общая
    # для всех групп сразу -- та же логика, что в методе 1: run_method2
    # даёт by_nuc_w с формой НА ЕДИНИЦУ амплитуды каждого нуклида,
    # фитируем 4 независимых столбца вместо одного shape_total). ДВА
    # варианта библиотеки, как у Th-232/Ra-226 (замечание оператора
    # №2, 11.08.2026 -- вариант "все известные линии" был сознательно
    # пропущен как упрощение по времени, здесь восстановлен): "sel" --
    # отобранная (порог 1,3%, ed.GAMMA_LIBRARY), "full" -- ВСЕ известные
    # линии из ensdf_amticseu_chain_lines.csv (287 строк, тот же файл,
    # что уже использует _sum_peaks_with_fb, путь строится САМ по
    # source.id конфига -- ed.load_full_library без явного path). ─────
    def run_method2_variant(library):
        shape_total, by_nuc_w, lines_out, n_sum = erd.run_method2(
            library, ed.SUM_PEAKS, resp, e, ch_edges, keys)
        cols = [by_nuc_w[k][sel] * T for k in keys] + [bgm]
        coef, dcoef, chi2, ndof, model = ed.fit_amplitudes(
            y_sel, cols, ed.SYS_FLOOR)
        groups_out = {}
        for i, key in enumerate(keys):
            A_Bq, dA_Bq = float(coef[i]), float(dcoef[i])
            A_pass = passport_by_group[key]["A_Bq"]
            groups_out[key] = {"A_Bq": A_Bq, "dA_Bq": dA_Bq,
                               "A_over_passport": A_Bq / A_pass}
        bg_amp, d_bg_amp = float(coef[-1]), float(dcoef[-1])
        # predicted_net на строку -- сумма ПО ГРУППАМ (её собственная
        # амплитуда), не одна общая амплитуда, как в Th-232/Ra-226.
        for ln in lines_out:
            A_Bq = groups_out[ln["nuclide"]]["A_Bq"]
            ln["predicted_net"] = ln.get("weight_per_branch", 0.0) * A_Bq * T
        stack = {k: (by_nuc_w[k] * groups_out[k]["A_Bq"] * T).tolist() for k in keys}
        return {
            "groups": groups_out, "bg_amplitude": bg_amp, "d_bg_amplitude": d_bg_amp,
            "chi2": chi2, "ndof": ndof, "chi2_ndof": chi2 / ndof,
            "n_lines": len(library), "n_sum_peaks": n_sum,
            "n_sum_peaks_total": len(ed.SUM_PEAKS),
            "lines": lines_out, "stack": stack,
        }

    method2_sel = run_method2_variant(ed.GAMMA_LIBRARY)
    lib_full, full_skip = ed.load_full_library(nuc_keys=set(keys))
    method2_full = run_method2_variant(lib_full)
    print("библиотека 'все известные линии': %d строк (пропущено: рентген=%d, "
          "без интенсивности=%d, чужой нуклид=%d)"
          % (len(lib_full), full_skip["xray"], full_skip["no_intensity"],
             full_skip["other_nuclide"]))

    for tag, method2 in (("sel", method2_sel), ("full", method2_full)):
        print("метод 2 (%s, F_B/TCS, 4 независимые амплитуды): chi2/ndof=%.2f  "
              "линий=%d сумм=%d/%d"
              % (tag, method2["chi2_ndof"], method2["n_lines"],
                 method2["n_sum_peaks"], method2["n_sum_peaks_total"]))
        for key in keys:
            g = method2["groups"][key]
            print("  %-12s A=%.1f+-%.1f Бк  ratio=%.3f"
                  % (key, g["A_Bq"], g["dA_Bq"], g["A_over_passport"]))

    reference_lines = [[ln["E_keV"], ln["nuclide"]]
                      for ln in method2_sel["lines"] if ln["kind"] == "line"]

    palette = {n["key"]: n["color"] for n in ed._CFG["nuclides"]}
    label_ru = {n["key"]: n["label_ru"] for n in ed._CFG["nuclides"]}

    data = {
        "meta": {
            "detector": "Гамма-1С (УДС-ГЦ-63х63)",
            "vessel": "Маринелли 1 л, risn379 ρ=1,0 г/см³ (лёгкая "
                     "органо-минеральная основа с 20%% кальция, состав "
                     "по MATERIAL самого .spe; СВОЙ Geant4-прогон "
                     "11.08.2026)",
            "source_kind": "mix_independent_nuclides",
            "live_s": meas["live_s"], "real_s": meas["real_s"],
            "bg_live_s": bg["live_s"], "bg_real_s": bg["real_s"],
            "bg_scale_time": bg_scale_time,
            "start_time": meas["start"],
            "fwhm662_keV": ed.FWHM662,
            "e_fit_lo": ed.E_FIT_LO, "e_fit_hi": ed.E_FIT_HI,
            "cal_sample": {"coefs": meas["coefs"],
                          "order": len(meas["coefs"]) - 1,
                          "n_channels": meas["n_channels"]},
            "cal_bg": {"coefs": bg["coefs"], "order": len(bg["coefs"]) - 1,
                      "n_channels": bg["n_channels"]},
            "energy_correction": energy_correction,
            "bg_energy_correction": {
                "applied": True,
                "method": "eleven_anchor_weighted_linear_refit",
                "anchor_kev": _BG_ANCHORS_KEV.tolist(),
                "anchor_ch": _BG_ANCHORS_CH.tolist(),
                "coefs_old": bg["coefs"],
                "coefs_new": [_bg_c0_new, _bg_c1_new],
                "rms_residual_kev": float(np.sqrt(np.mean(
                    (_bg_A @ _bg_coef_new - _BG_ANCHORS_KEV) ** 2))),
                "note": "Полный пересчёт линейной калибровки фона (не "
                        "локальная поправка) по 11 якорям: 10 природных "
                        "линий (методика ЛСРМ для ЕРН/фона NaI, "
                        "SpectraVibe references/01_metadata_calibration.md) "
                        "+ К-рентген Pb защиты (77,0 кэВ, включён по прямой "
                        "директиве оператора, статистически слабее "
                        "остальных 10 -- центроид не сходится с шириной "
                        "окна, канал взят по argmax). Взвешенная МНК "
                        "(вес~sqrt(чистый счёт)). Первая версия (сдвиг нуля "
                        "по одному K-40 с наивным широким окном) отклонена "
                        "оператором как неверная -- "
                        "analysis/bg_seven_line_anchor_check.py. Подробности "
                        "-- amticseu-remarks.md §13.",
            },
            "level_note": "Библиотека линий -- IAEA Live Chart of "
                          "Nuclides (decay_rads), ✅ ПОЛНОСТЬЮ сверена "
                          "построчно с LNHB/DDEP (20/20 строк отобранной "
                          "библиотеки, макс. расхождение 1,94%%, "
                          "11.08.2026); полная (непороговая) библиотека "
                          "(%d линий, добавлена 11.08.2026 по замечанию "
                          "оператора, как на Th-232/Ra-226) LNHB-сверке "
                          "НЕ подвергалась -- там 📗, как и было. "
                          "К-рентген дочерних не выделен отдельной "
                          "сущностью (физика в шаблонах есть целиком, "
                          "отдельная разбивка -- только наглядность)."
                          % len(lib_full),
        },
        "fwhm_cal": fwhm_cal,
        "passport": passport_by_group,
        "nuclides": [{"key": k, "label_ru": label_ru[k], "color": palette[k]}
                    for k in keys],
        "spectrum": {
            "e_of_ch": e.tolist(),
            "counts": meas["counts"].tolist(),
            "bg_counts": bg_scaled.tolist(),
            "stack1": m1_stack,
        },
        "method1": m1_result,
        "method1_meta": m1_meta,
        "method2_sel": method2_sel,
        "method2_full": method2_full,
        "reference_lines": reference_lines,
    }

    out = os.path.join(HERE, "g1s_amticseu_data.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, separators=(",", ":"),
                  allow_nan=False)
    print("написано: %s (%d КБ)" % (out, os.path.getsize(out) // 1024))


if __name__ == "__main__":
    main()
