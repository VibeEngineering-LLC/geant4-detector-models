# -*- coding: utf-8 -*-
"""
Выгрузка данных для вкладки «Свинцовый домик» (web-background).

Файл: out/shield_data.json

Данные пересчитываются заново при каждом запуске на основе исходных измерений
и моделей, а не копируются из готовых отчётов. Включает спектр в домике,
предсказание без подгонки, понуклидные вклады, разделение по происхождению
(если доступно), пропускание и таблицы по полосам.
"""

import os
import sys
import json
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)

import fit_two_criteria as ftc
import fit_physical_chains as fpc
import predict_shield as ps
import decompose_shield as ds

OUT = os.path.join(HERE, "out", "shield_data.json")
STEP = 2  # шаг прореживания рядов для экономии места на странице


def read_origin(path):
    """Читает ДОПОЛНИТЕЛЬНЫЕ колонки шаблона: счёт, разделённый по происхождению.

    Формат файла: блок метаданных «ключ,значение», пустая строка, строка-заголовок
    «bin_keV,counts,cps[,counts_direct,counts_pb_scat,counts_pb_born]», затем
    данные через запятую. Шаблон, посчитанный до 29.08.2026, содержит только три
    колонки — это не ошибка, возвращаем None, вкладка работает и без разделения.

    Возвращает (словарь из трёх массивов сырых отсчётов МК, t_run_s).
    """
    if not os.path.exists(path):
        return None, None
    t_run = None
    rows, ncol, in_hist = [], 0, False
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if not in_hist:
                if line.startswith("bin_keV"):
                    ncol = len(line.split(","))
                    in_hist = True
                elif line.startswith("t_run_s,"):
                    t_run = float(line.split(",", 1)[1])
                continue
            p = line.split(",")
            if len(p) >= 6:
                rows.append((int(float(p[0])), float(p[3]), float(p[4]), float(p[5])))
    if ncol < 6 or not rows:
        return None, None
    n = max(r[0] for r in rows) + 1
    out = {k: np.zeros(n) for k in ("direct", "pb_scat", "pb_born")}
    for b_, d, sc, bo in rows:
        out["direct"][b_], out["pb_scat"][b_], out["pb_born"][b_] = d, sc, bo
    return out, t_run


def main():
    # 1. Исходные активности из открытого фона
    amp_map, names = ds.fit_open()
    acts = ds.nuclide_activities(amp_map)

    # 2. Чтение измерений
    # Открытый фон
    cnt_o, e_o, live_o, _ = ps.read_meas(ftc.MEAS_NAME, ftc.CAL_ROOM)
    sel_o = (e_o >= ftc.E_LO) & (e_o < ftc.E_HI)

    # Домик
    cnt_s, e_s, live_s, cal_s = ps.read_meas(ps.MEAS_SHIELD, None)
    sel_s = (e_s >= ftc.E_LO) & (e_s < ftc.E_HI)

    # 3. Понуклидные вклады в домике
    parts = ds.per_nuclide_spectra(ps.FMT_SHIELD, ps.MUON_SHIELD_CSV, e_s[sel_s], acts)
    
    # Полная модель домика (сумма вкладов)
    model_shield = np.zeros_like(e_s[sel_s])
    for k in parts:
        model_shield += parts[k]

    # 4. Пропускание
    T = ds.transmission(e_s[sel_s])

    # 5. Модель открытого фона (для таблиц ослабления)
    parts_open = ds.per_nuclide_spectra(ps.FMT_OPEN, ftc.MUON_CSV, e_o[sel_o], acts)
    model_open = np.zeros_like(e_o[sel_o])
    for k in parts_open:
        model_open += parts_open[k]

    # 6. Разделение по происхождению
    origin_direct = np.zeros_like(e_s[sel_s])
    origin_pb_scat = np.zeros_like(e_s[sel_s])
    origin_pb_born = np.zeros_like(e_s[sel_s])
    
    has_origin_flag = True
    
    # Проверяем наличие данных о происхождении для всех нуклидов (кроме мюонов)
    # Нуклиды - это ключи в acts, которые не являются 'mu'
    nuclide_keys = [k for k in acts if k != 'mu']
    
    origin_data_available = True
    
    # Сначала проверим, есть ли данные в шаблонах
    for name in nuclide_keys:
        # Путь к шаблону для данного нуклида в домике
        # ps.FMT_SHIELD - это формат имени файла или путь? 
        # В ds.per_nuclide_spectra используется ps.FMT_SHIELD.
        # Обычно там лежит путь или шаблон имени.
        # Предположим, что мы можем получить путь к файлу шаблона для нуклида.
        # В модуле ds обычно есть функция получения пути.
        # Если нет, попробуем реконструировать имя файла как в ps/ds.
        
        # Для надежности, проверим наличие файлов через попытку чтения
        # Но read_origin требует path. 
        # В predict_shield или decompose_shield должны быть константы путей.
        # Попробуем использовать логику из ds.per_nuclide_spectra косвенно.
        
        # Если мы не можем легко получить путь, мы можем попробовать прочитать
        # тот же файл, который читает ds. Но у нас нет доступа к внутренностям ds.
        # Однако, в условии сказано: "Шаблоны с домиком... содержат ДОПОЛНИТЕЛЬНЫЕ колонки".
        # Значит, файлы те же самые, что используются для моделирования.
        
        # Попробуем найти путь. Обычно это что-то вроде data/templates/shield/...
        # Если мы не знаем точный путь, мы можем попробовать прочитать файл, 
        # который был использован при расчете parts. Но parts уже посчитаны.
        
        # Альтернатива: предположить, что если один шаблон имеет колонки, то и все имеют.
        # Или проверить первый нуклид.
        
        pass

    # Попытка чтения origin данных
    # Нам нужно знать пути к файлам шаблонов. 
    # В модуле ps или ds должны быть функции для получения путей.
    # Если их нет, мы можем попробовать угадать структуру имен файлов.
    # Обычно: f"templates/shield/{name}.txt" или подобное.
    
    # Чтобы не гадать, попробуем использовать тот факт, что 
    # ds.per_nuclide_spectra читает файлы. Если бы мы могли перехватить пути...
    # Но мы не можем менять ds.
    
    # Давайте предположим стандартную структуру имен файлов, используемую в подобных проектах:
    # Обычно файлы лежат в подкаталоге относительно HERE или в data/.
    # В условии сказано "Шаблоны с домиком". 
    # Попробуем прочитать файл для первого нуклида, чтобы проверить наличие колонок.
    
    # Если мы не можем получить путь, has_origin будет False.
    
    # Разделение по происхождению читаем ПРЯМО из тех же шаблонов, что и всё
    # остальное: каталог ftc.TEMPLATE_DIR, имя ps.FMT_SHIELD % нуклид. Колонки
    # сырые (отсчёты МК), поэтому делим на t_run_s — получаются cps на 1 Бк/кг.
    origin_parts = {k: np.zeros_like(e_s[sel_s])
                    for k in ("direct", "pb_scat", "pb_born")}
    all_have_origin = True
    for name in nuclide_keys:
        data, t_run = read_origin(os.path.join(ftc.TEMPLATE_DIR, ps.FMT_SHIELD % name))
        if data is None or not t_run:
            all_have_origin = False
            continue
        for key in origin_parts:
            col = ftc.rcspec.fold(data[key] / t_run, "103", tail_T=ps.TAIL_T)
            origin_parts[key] += ps.to_meas_grid(col, e_s[sel_s]) * acts[name]
    if not all_have_origin:
        origin_parts = {k: np.zeros_like(e_s[sel_s]) for k in origin_parts}

    # Мюоны отдельно
    origin_muon = parts.get("mu", np.zeros_like(e_s[sel_s]))

    # 7. Подготовка данных для JSON
    
    # Прореживание
    idx = np.arange(0, len(e_s[sel_s]), STEP)
    
    def ser(v):
        """Сериализация ряда: прореживание и округление."""
        arr = np.asarray(v)[idx]
        return [round(float(x), 8) for x in arr]

    def ser_trans(v):
        """Сериализация пропускания: nan -> None."""
        arr = np.asarray(v)[idx]
        res = []
        for x in arr:
            if np.isnan(x):
                res.append(None)
            else:
                res.append(round(float(x), 8))
        return res

    energy_grid = ser(e_s[sel_s])
    measured_shield = ser(cnt_s[sel_s] / live_s) # cps
    model_shield_ser = ser(model_shield)
    
    parts_ser = {}
    for k in parts:
        parts_ser[k] = ser(parts[k])

    # T — СЛОВАРЬ {нуклид: ряд}, сериализуется поключно.
    trans_ser = {k: ser_trans(v) for k, v in T.items()}

    origin_ser = {
        "direct": ser(origin_parts["direct"]),
        "pb_scat": ser(origin_parts["pb_scat"]),
        "pb_born": ser(origin_parts["pb_born"]),
        "muon": ser(origin_muon)
    }

    # Активности
    acts_ser = {k: round(float(v), 8) for k, v in acts.items()}

    # Мюонные параметры (если есть в amp_map или отдельно)
    # ⚠ pdg — это ПОТОК PDG через диск шаблона (0,0167 см^-2 c^-1 x площадь), а
    # НЕ амплитуда. ps.load_cps переопределяет ftc.MUON_PDG_PER_S значением из
    # метаданных мюонного шаблона, поэтому оно соответствует его радиусу.
    # Отношение — amp/pdg, справочная сверка порядка, а не доля среди амплитуд.
    if 'mu' in amp_map:
        amp_mu = float(amp_map['mu'])
        pdg_mu = float(ftc.MUON_PDG_PER_S)
        muon_info = {"amp": round(amp_mu, 8),
                     "pdg": round(pdg_mu, 8),
                     "ratio": round(amp_mu / pdg_mu, 8) if pdg_mu > 0 else None}
    else:
        muon_info = {"amp": 0, "pdg": 0, "ratio": None}

    # Setup параметры
    setup_info = {
        "pb_mm": 50,
        "cavity": "150x150x385",
        "top_open": True,
        "tail_T": ps.TAIL_T,
        "f_rn": ps.F_RN,
        "r_th": ps.R_TH,
        "f_tn": ps.F_TN
    }

    # Chi2 и Shape
    # model_counts = полная модель * live_s
    model_counts_shield = model_shield * live_s
    meas_counts_shield = cnt_s[sel_s]
    
    chi2ndf, shape = ftc.metrics(model_counts_shield, meas_counts_shield, 0)

    # Таблицы
    
    # Bands
    bands_list = []
    for band in ftc.BANDS:
        lo, hi = band
        
        # Индексы для открытого фона
        idx_o = (e_o[sel_o] >= lo) & (e_o[sel_o] < hi)
        # Индексы для домика
        idx_s = (e_s[sel_s] >= lo) & (e_s[sel_s] < hi)
        
        # Суммы измерений (cps)
        meas_open_cps = np.sum(cnt_o[sel_o][idx_o]) / live_o
        meas_shield_cps = np.sum(cnt_s[sel_s][idx_s]) / live_s
        
        # Суммы моделей (cps)
        model_open_cps = np.sum(model_open[idx_o])
        model_shield_cps = np.sum(model_shield[idx_s])
        
        # Отношения
        ratio = model_shield_cps / meas_shield_cps if meas_shield_cps > 0 else 0
        
        # Ослабление
        att_meas = meas_open_cps / meas_shield_cps if meas_shield_cps > 0 else 0
        att_model = model_open_cps / model_shield_cps if model_shield_cps > 0 else 0
        
        # Доли нуклидов в модели домика
        frac = {}
        total_model_shield_band = np.sum(model_shield[idx_s])
        for k in parts:
            part_sum = np.sum(parts[k][idx_s])
            pct = (part_sum / total_model_shield_band * 100) if total_model_shield_band > 0 else 0
            frac[k] = round(pct, 2)
            
        bands_list.append({
            "lo": lo,
            "hi": hi,
            "meas": round(meas_shield_cps, 8),
            "model": round(model_shield_cps, 8),
            "ratio": round(ratio, 8),
            "att_meas": round(att_meas, 8),
            "att_model": round(att_model, 8),
            "frac": frac
        })

    # PbK bands
    pbk_intervals = [(40,60),(60,70),(70,80),(80,90),(90,100),(100,120),(120,160)]
    pbk_list = []
    for lo, hi in pbk_intervals:
        idx_s = (e_s[sel_s] >= lo) & (e_s[sel_s] < hi)
        
        meas_cps = np.sum(cnt_s[sel_s][idx_s]) / live_s
        model_cps = np.sum(model_shield[idx_s])
        ratio = model_cps / meas_cps if meas_cps > 0 else 0
        
        pbk_list.append({
            "lo": lo,
            "hi": hi,
            "meas": round(meas_cps, 8),
            "model": round(model_cps, 8),
            "ratio": round(ratio, 8)
        })

    # Contrib
    contrib_list = []
    total_model_cps = np.sum(model_shield)
    for k in parts:
        part_cps = np.sum(parts[k])
        pct = (part_cps / total_model_cps * 100) if total_model_cps > 0 else 0
        contrib_list.append({
            "name": k,
            "cps": round(part_cps, 8),
            "pct": round(pct, 2)
        })
    # Сортировка по убыванию cps
    contrib_list.sort(key=lambda x: x["cps"], reverse=True)

    # Сборка итогового словаря
    result = {
        "meas_shield": {
            "file": ps.MEAS_SHIELD,
            "live_s": live_s,
            "live_h": round(live_s / 3600, 2),
            "live_d": round(live_s / 86400, 2),
            "total_counts": int(np.sum(cnt_s[sel_s])),
            "cps": round(np.sum(cnt_s[sel_s]) / live_s, 8),
            "cal": [round(float(c), 8) for c in cal_s] if cal_s is not None else []
        },
        "meas_open": {
            "live_s": live_o,
            "live_h": round(live_o / 3600, 2),
            "cps": round(np.sum(cnt_o[sel_o]) / live_o, 8)
        },
        "energy": energy_grid,
        "measured": measured_shield,
        "model": model_shield_ser,
        "parts": parts_ser,
        "trans": trans_ser,
        "origin": origin_ser,
        "has_origin": all_have_origin,
        "acts": acts_ser,
        "muon": muon_info,
        "setup": setup_info,
        "chi2ndf": round(chi2ndf, 8),
        "shape": round(shape, 8),
        "bands": bands_list,
        "pbk": pbk_list,
        "contrib": contrib_list,
        "ru": ds.RU if hasattr(ds, 'RU') else {}
    }

    # Запись файла
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w', encoding='utf-8') as f:
        json.dump(result, f, ensure_ascii=False, separators=(",", ":"))
        
    size = os.path.getsize(OUT)
    print(f"{OUT} ({size} bytes)")

if __name__ == "__main__":
    main()
