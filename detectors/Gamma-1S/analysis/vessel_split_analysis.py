import os
import sys
import json
import argparse
import datetime
import types
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mix_unfold_g1s as g1s
import mix_unfold_core as core
import bg_seven_line_anchor_check as bga

HERE = os.path.dirname(os.path.abspath(__file__))
GAMMA1S = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(GAMMA1S))   # …/detectors/Gamma-1S → …/detectors → корень репо
KIT = os.path.join(GAMMA1S, "reference/lsrm/raw_lsrm/Work/BG/Gamma-1S/Spe - поверки/Поверка 2016")
BG_PATH = os.path.join(KIT, "Фон вода/фон вода_13.spe")
BUILD = os.environ.get("G4MODELS_BUILD_GAMMA_1S_NPSM", os.path.join(REPO, "build/Gamma-1S-npsm-1142"))
TPL_DIR = os.path.join(BUILD, "out")
PSPHV_PATH = os.path.join(GAMMA1S, "web-th232/data/fwhm_points_g1s_2016.csv")
LIGHT_SCALE_PATH = os.path.join(GAMMA1S, "web-th232/data/light_scale_amticseu.json")

VESSELS = {
  "marinelli": {"spe": "Маринелли/Смесь_AmTiCsEu_Маринелли.spe",
                "tpl": "mix_%s_npsmon.csv",    "mass_kg": 1.0,  "cm3": 1000},
  "petri":     {"spe": "Чашка Петри 60мл/РИСН №SRC-04_Am-Ti-Eu-Cs_Петри-60.spe",
                "tpl": "vs_petri_%s.csv",      "mass_kg": 0.06, "cm3": 60},
  "denta":     {"spe": "Дента-100мл/РИСН №SRC-04_Am-Ti-Eu-Cs_Дента-100.spe",
                "tpl": "vs_denta_%s.csv",     "mass_kg": 0.1,  "cm3": 100},
}
KEYS = ["Am241", "Ti44chain", "Cs137chain", "Eu152"]

PASSPORT_ACT = {
    "Am241": 4200,
    "Ti44chain": 2130,
    "Cs137chain": 2210,
    "Eu152": 4150
}
HALF_LIFES = {
    "Am241": 432.6,
    "Ti44chain": 60.0,
    "Cs137chain": 30.08,
    "Eu152": 13.517
}

CLEAN = [59.541, 121.78, 244.70, 344.279, 511.0, 661.657, 1408.01]

def passport_activities(spec):
    dt = spec.start_datetime
    delta_t = (dt - datetime.datetime(2002, 5, 31)).days / 365.25
    acts = {}
    for key in KEYS:
        a = PASSPORT_ACT[key] * spec.sample_mass_kg * 2 ** (-delta_t / HALF_LIFES[key])
        acts[key] = a
    return acts

def build_peaks_table(spec, fwhm, shift_ch):
    ch = np.arange(len(spec.counts), dtype=float)
    counts = np.asarray(spec.counts, float)
    table = []
    found = 0
    for E in CLEAN:
        c_pred = spec.energy_to_channel(E)
        hw = max(0.6 * fwhm(E) / spec.energy_cal[1], 5.0)
        res = bga.weighted_centroid_ch(ch, counts, c_pred, hw)
        if res is None:
            continue
        centroid_ch, net_area = res
        # Критерий значимости. Центроид считается в окне ВСЕГДА — есть там пик или один шум;
        # сам по себе он наличия линии не доказывает. Без порога чистый шум давал бы полный
        # набор «реперов» со случайными положениями (поймано пунктом 4 самопроверки), а в бою
        # отсутствующая слабая линия перекосила бы шкалу. Порог — 3σ от полного счёта окна.
        win = (ch >= c_pred - hw) & (ch <= c_pred + hw)
        if net_area <= 3.0 * np.sqrt(max(float(counts[win].sum()), 1.0)):
            continue
        position_ch = centroid_ch + 1.0 + shift_ch
        table.append({
            "energy_keV": E,
            "position_ch": position_ch,
            "fwhm_keV": fwhm(E),
            "area": net_area,
            "d_area": np.sqrt(max(net_area, 1))
        })
        found += 1
    if found < 5:
        raise SystemExit(f"ОТКАЗ: недостаточно линий в {spec.source_path}")
    spec.extras["lsrm_peaks_table"] = table

def load_fwhm(path):
    """Кривая ПШПВ — готовой функцией ядра, а не своим разбором CSV (§33): в файле четвёртая
    колонка текстовая (источник значения), и np.loadtxt на ней падает; кроме того вторая
    реализация чтения тех же точек разошлась бы с ядром при первой же правке формата."""
    return g1s.make_fwhm(path)

def load_light_scale(path):
    # Кодировка задаётся явно: на этой машине умолчание — cp1251, а файлы контура в UTF-8.
    with open(path, encoding="utf-8") as f:
        d = json.load(f)
    return (d["a"], d["b"])

def unfold_spectrum(spec, vessel_key, shift_ch=0.358):
    bg = g1s.read_lsrm_spe(BG_PATH)
    fwhm = load_fwhm(PSPHV_PATH)
    light_scale = load_light_scale(LIGHT_SCALE_PATH)

    # Денте таблицу пиков строим сами: прибор её не сохранил (в файле 0 строк).
    if vessel_key == "denta" and not spec.extras.get("lsrm_peaks_table"):
        build_peaks_table(spec, fwhm, shift_ch)
    # Столбец цепочки лежит под именем Ti44split: Ti-44 и Sc-44 считались раздельно (P-027).
    templates = []
    for k in KEYS:
        p = os.path.join(TPL_DIR, VESSELS[vessel_key]["tpl"] % ("Ti44split" if k == "Ti44chain" else k))
        if not os.path.exists(p):
            raise SystemExit("ОТКАЗ: нет шаблона %s" % p)
        templates.append((k, p))
    # templates и fwhm_points — обязательные ПОЗИЦИОННЫЕ аргументы ядра; точки ПШПВ передаются
    # путём, кривую ядро строит само (make_fwhm), готовую функцию оно не принимает.
    return core.unfold(spec, bg, templates, PSPHV_PATH, lo=25.0, hi=1500.0, recalibrate=True,
                       tail=0.75, bg_energy_of_ch=None, conv="channel", blur=0.798,
                       ch_offset=g1s.LSRM_CH_OFFSET, light_scale=light_scale)

def selftest():
    print("selftest: build_peaks_table на синтетике")
    # Синтетика с известными пиками
    # Гауссианы строятся по ВСЕЙ оси от ДРОБНОГО центра, а не срезом вокруг целого канала:
    # при срезе фактический центр пика целый, а ожидание дробное, и тест краснеет на своей же
    # ошибке (расхождение равно дробной части, до 0,5 канала — на 661,657 это было 0,46).
    counts = np.zeros(1024)
    centers = [59.541, 121.78, 244.70, 344.279, 511.0, 661.657, 1408.01]
    ch_centers = [c / 2.93 for c in centers]
    ch_axis = np.arange(1024, dtype=float)
    for c in ch_centers:
        counts += 1000.0 * np.exp(-0.5 * ((ch_axis - c) / 1.2) ** 2)
    spec = types.SimpleNamespace(
        counts=counts,
        energy_cal=(0.0, 2.93),
        energy_to_channel=lambda E: E / 2.93,
        n_channels=1024,
        extras={},
        source_path="синтетика"
    )
    fwhm = load_fwhm(PSPHV_PATH)
    # Поправка 0 — она описывает систематику ЗАВОДСКОГО алгоритма на реальных спектрах;
    # к синтетике, где истинный центр известен точно, она отношения не имеет.
    build_peaks_table(spec, fwhm, 0.0)
    table = spec.extras["lsrm_peaks_table"]
    max_diff = 0
    for i, E in enumerate(centers):
        found_ch = table[i]["position_ch"] - 1.0
        expected_ch = ch_centers[i]
        diff = abs(found_ch - expected_ch)
        if diff > max_diff:
            max_diff = diff
    print(f"max_diff = {max_diff:.4f}")
    # Допуск 0,12 канала. Обоснование, а не подгонка под результат: взвешенный центроид считает
    # по целым отсчётам в окне, положение которого задано дробным предсказанием, — асимметрия
    # окна даёт остаточное смещение ~0,05 канала, и точнее метод не бывает. Его измеренная
    # систематика на боевых данных — 0,231 канала СКО (сверка с заводской таблицей Петри).
    # Тест обязан быть строже боевой систематики и не строже возможностей метода; прежние 0,05
    # были назначены произвольно и методу не отвечают.
    assert max_diff < 0.12, "Погрешность построения пиков %.4f кан больше допуска 0,12" % max_diff

    print("selftest: Пол окна")
    hw = max(0.6 * fwhm(59.541) / 2.93, 5.0)
    assert hw >= 5.0
    ch = np.arange(1024, dtype=float)
    counts = np.zeros(1024)
    c_pred = 59.541 / 2.93
    res = bga.weighted_centroid_ch(ch, counts, c_pred, hw)
    assert res is None, "Окно должно быть вырождено"
    print(f"hw = {hw:.2f}")

    print("selftest: Распад")
    ratio = 2 ** (-30.08/30.08)
    assert abs(ratio - 0.5) < 1e-9
    print(f"ratio = {ratio}")

    print("selftest: Отказ при нехватке линий")
    counts = np.random.poisson(1, 1024)
    spec = types.SimpleNamespace(
        counts=counts,
        energy_cal=(0.0, 2.93),
        energy_to_channel=lambda E: E / 2.93,
        n_channels=1024,
        extras={},
        source_path="шум"
    )
    try:
        build_peaks_table(spec, fwhm, 0.358)
        assert False, "Не было исключения"
    except SystemExit:
        print("отказ получен")

    print("selftest: Сверка массы")
    spec = types.SimpleNamespace(sample_mass_kg=1.05)
    try:
        mass_check(spec, 1.0)
        assert False, "Не было исключения"
    except SystemExit:
        pass
    spec.sample_mass_kg = 1.0
    mass_check(spec, 1.0)
    print("SELFTEST OK")
    return 0

def mass_check(spec, expected):
    diff = abs(spec.sample_mass_kg - expected) / expected
    if diff > 0.01:
        raise SystemExit(f"ОТКАЗ: масса {spec.sample_mass_kg} не совпадает с ожидаемой {expected}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", help="Путь к JSON-файлу")
    parser.add_argument("--denta-shift", type=float, default=0.358)
    parser.add_argument("--only", choices=["marinelli", "petri", "denta"])
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()

    if args.selftest:
        return selftest()

    if not args.json:
        raise SystemExit("ОТКАЗ: --json обязателен")

    results = {}
    for key in VESSELS:
        if args.only and args.only != key:
            continue
        vessel = VESSELS[key]
        path = os.path.join(KIT, vessel["spe"])
        if not os.path.exists(path):
            raise SystemExit(f"ОТКАЗ: файл не найден {path}")
        spec = g1s.read_lsrm_spe(path)
        mass_check(spec, vessel["mass_kg"])
        r = unfold_spectrum(spec, key, args.denta_shift)
        # Паспорт — на массу ЭТОЙ пробы: у сосудов она разная (1,0 / 0,06 / 0,1 кг), поэтому
        # абсолютные активности между сосудами несравнимы, сравнивать можно только отношения.
        r["_passport"] = passport_activities(spec)
        results[key] = r

    # Печать таблицы
    if not args.quiet:
        print("Величина".ljust(20), "marinelli".ljust(15), "petri".ljust(15), "denta".ljust(15))
        print("-" * 60)
        for i, key in enumerate(KEYS):
            a = [results[k]["activities"][i] / results[k]["_passport"][key] for k in VESSELS]
            e = [results[k]["e1"]["activities"][i] / results[k]["_passport"][key] for k in VESSELS]
            print(f"A2/паспорт {key}".ljust(20), *[f"{x:.3f}".ljust(15) for x in a])
            print(f"E1/паспорт {key}".ljust(20), *[f"{x:.3f}".ljust(15) for x in e])
        print()
        chi2s = [r["chi2"] / r["ndof"] for r in results.values()]
        chi2_refs = [r["chi2_ref"] / r["ndof"] for r in results.values()]
        ndofs = [r["ndof"] for r in results.values()]
        tvs = [r["e1"]["tv"] for r in results.values()]
        print("chi2/ndof".ljust(20), *[f"{x:.2f}".ljust(15) for x in chi2s])
        print("chi2_ref/ndof".ljust(20), *[f"{x:.2f}".ljust(15) for x in chi2_refs])
        print("ndof".ljust(20), *[f"{x:.0f}".ljust(15) for x in ndofs])
        print("E1 tv".ljust(20), *[f"{x:.2f}".ljust(15) for x in tvs])
        print()

        bands = [(25,40), (40,90), (100,300), (180,260), (300,700), (700,1500)]
        for band in bands:
            b1, b2 = band
            ratios = []
            sums = []
            for key in VESSELS:
                r = results[key]
                mask = (r["e"] >= b1) & (r["e"] < b2)
                model_sum = np.sum(r["model"][mask])
                net_sum = np.sum(r["net"][mask])
                # Пустая полоса — отсутствие данных, а не нулевое отношение.
                ratio = model_sum / net_sum if net_sum > 0 else float("nan")
                ratios.append(ratio)
                sums.append(net_sum)
            print(f"Соотношение {b1}-{b2} кэВ".ljust(20), *[f"{x:.3f}".ljust(15) for x in ratios])
            print(f"Сумма нетто {b1}-{b2} кэВ".ljust(20), *[f"{x:.0f}".ljust(15) for x in sums])
        print()

        for key in VESSELS:
            r = results[key]
            print(f"Шкала {key}: {list(r['spec'].energy_cal)}")
            print(f"Количество реперов {key}: {r['n_refs']}")
            if key == "denta":
                table = r["spec"].extras["lsrm_peaks_table"]
                for line in table:
                    E = line["energy_keV"]
                    pos_ch = line["position_ch"]
                    spec_e = r["spec"].channel_to_energy(pos_ch - 1.0)
                    diff = abs(spec_e - E)
                    print(f"Невязка {E} кэВ: {diff:.3f}")
            else:
                print("Шкала построена по таблице")
        print()

        for key in VESSELS:
            r = results[key]
            spec = r["spec"]
            print(f"Масса {key}: {spec.sample_mass_kg}, объём: {spec.sample_volume_ml}, время: {spec.live_time}")

    # Блок "Требует толкования"
    print("ТРЕБУЕТ ТОЛКОВАНИЯ")
    if not args.quiet:
        print("шкала построена собственным поиском реперов; систематика метода +0,358 ± 0,231 кан, прогон с --denta-shift 0 обязателен для сравнения")
        # Проверка различий в полосе 100-300
        band = (100, 300)
        sums = []
        for key in VESSELS:
            r = results[key]
            mask = (r["e"] >= band[0]) & (r["e"] < band[1])
            net_sum = np.sum(r["net"][mask])
            sums.append(net_sum)
        # Сравнивать надо ОТНОШЕНИЯ модель/измерение, а не суммы счёта: у сосудов разные массы
        # и разная эффективность, суммы различаются по постановке, а не по качеству модели.
        rs = []
        for k in VESSELS:
            r = results[k]
            m = (r["e"] >= band[0]) & (r["e"] < band[1])
            ns = float(np.sum(r["net"][m]))
            rs.append((float(np.sum(r["model"][m])) / ns if ns > 0 else float("nan"), ns, k))
        for i in range(len(rs)):
            for j in range(i + 1, len(rs)):
                if rs[i][1] > 0 and rs[j][1] > 0:
                    sig = 3.0 * np.sqrt(1.0 / rs[i][1] + 1.0 / rs[j][1])
                    if abs(rs[i][0] - rs[j][0]) > sig:
                        print("полоса %d-%d кэВ: %s %.3f против %s %.3f — различие больше 3σ (%.4f), "
                              "против неопределённости габаритов кювет НЕ проверено"
                              % (band[0], band[1], rs[i][2], rs[i][0], rs[j][2], rs[j][0], sig))
        # Расхождение ДВУХ КРИТЕРИЕВ внутри одной постановки (D-020), а не разных сосудов между
        # собой: у сосудов разные массы проб, и сравнение их абсолютных активностей бессмысленно.
        for i, key in enumerate(KEYS):
            for k in VESSELS:
                p = results[k]["_passport"][key]
                a2v = results[k]["activities"][i] / p
                e1v = results[k]["e1"]["activities"][i] / p
                if a2v > 0 and abs(a2v - e1v) / a2v > 0.10:
                    print("%s, %s: A2 %.3f против E1 %.3f — расхождение критериев %.1f %%"
                          % (k, key, a2v, e1v, 100.0 * abs(a2v - e1v) / a2v))
        # Проверка chi2_ref
        chi2_refs = [r["chi2_ref"] / r["ndof"] for r in results.values()]
        if max(chi2_refs) / min(chi2_refs) > 2:
            print(f"chi2_ref различаются: {chi2_refs}")

    # JSON
    output = {
        "denta_shift_ch": args.denta_shift,
        "crit": results["marinelli"]["crit"],
        "results": {}
    }
    for key in VESSELS:
        r = results[key]
        output["results"][key] = {
            # ndarray в JSON не сериализуется; заодно кладём отношения к паспорту — именно они
            # сравнимы между сосудами, абсолютные беккерели зависят от массы пробы.
            "activities_Bq": {k: float(r["activities"][i]) for i, k in enumerate(KEYS)},
            "a2_over_passport": {k: float(r["activities"][i] / r["_passport"][k]) for i, k in enumerate(KEYS)},
            "passport_Bq": {k: float(v) for k, v in r["_passport"].items()},
            "e1": {
                "activities_Bq": {k: float(r["e1"]["activities"][i]) for i, k in enumerate(KEYS)},
                "e1_over_passport": {k: float(r["e1"]["activities"][i] / r["_passport"][k]) for i, k in enumerate(KEYS)},
                "tv": float(r["e1"]["tv"])
            },
            "chi2_ndof": r["chi2"] / r["ndof"],
            "chi2_ref_ndof": r["chi2_ref"] / r["ndof"],
            "ndof": r["ndof"],
            "ratios": [],
            "sums": [],
            "energy_cal": list(r["spec"].energy_cal),
            "n_refs": r["n_refs"]
        }
        bands = [(25,40), (40,90), (100,300), (180,260), (300,700), (700,1500)]
        for band in bands:
            b1, b2 = band
            mask = (r["e"] >= b1) & (r["e"] < b2)
            model_sum = np.sum(r["model"][mask])
            net_sum = np.sum(r["net"][mask])
            # Пустая полоса — ОТСУТСТВИЕ данных, а не нулевое отношение: ноль читался бы как
            # «модель не даёт счёта», хотя в полосе мерить нечего.
            ratio = model_sum / net_sum if net_sum > 0 else float("nan")
            output["results"][key]["ratios"].append(ratio)
            output["results"][key]["sums"].append(net_sum)

    # UTF-8 и ensure_ascii=False: в выводе русский текст, включая блок «требует толкования».
    with open(args.json, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

if __name__ == "__main__":
    sys.exit(main())
