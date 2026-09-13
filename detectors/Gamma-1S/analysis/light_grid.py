import sys
import os
import glob
import math
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

def build_nodes(build_out, prefix, read_template, light_to_energy, broaden_ch, e, file_e, fwhm, a, b, blur, paths=None):
    if paths is None:
        paths = glob.glob(os.path.join(build_out, prefix + "*.csv"))
    nodes = {}
    skipped = []
    empty = []
    for path in paths:
        name = os.path.basename(path)
        try:
            energy_str = name[len(prefix):-4]
            energy = float(energy_str)
        except ValueError:
            skipped.append(name)
            continue
        hist, n, npsm = read_template(path)
        if npsm != 1:
            raise SystemExit(f"Неправильный npsm в файле {name}")
        hist = {k: v for k, v in hist.items() if k >= 1}
        if sum(hist.values()) / n > 0.6:
            raise SystemExit(f"Слишком большой суммарный счёт после отбрасывания: {sum(hist.values()) / n} для энергии {energy}")
        e_of_ch = lambda c: np.interp(c, np.arange(len(e)), e)
        hist_e, st = light_to_energy(hist, a, b, e_of_ch, len(e))
        if st["sum_used"] <= 0:
            # Квант не регистрируется вовсе (напр. 11,89 кэВ в сосуде Маринелли — всё поглощают стенки
            # и матрица): это законный нулевой узел, а не отказ. Останавливать весь расчёт из-за линии
            # с физически нулевым вкладом нельзя, молчать о ней — тем более, поэтому она идёт в empty.
            nodes[energy] = (np.zeros(len(e)), 0.0, float("nan"), st["n_used"], st["n_in"])
            empty.append(energy)
            continue
        shape = broaden_ch(hist_e, n, e, fwhm, file_e, blur)
        E0 = energy
        # Окно расширено 12.09 по измерению: ниже 50 кэВ пик после перевода стоит на 4-5 кэВ НИЖЕ
        # номинала (26,34 → 22,05; 35,09 → 30,88), и с прежними 3,75 кэВ он не находился вовсе —
        # вклад K-рентгена цезия и L-рентгена нептуния молча обнулялся. Окно следует за светом.
        window = max(6.0, 1.5 * fwhm(E0))
        keys = [k for k in hist_e.keys() if abs(k - E0) <= window]
        if not keys:
            nodes[E0] = (shape, 0.0, float("nan"), st["n_used"], st["n_in"])
        else:
            E_peak = max(keys, key=lambda k: hist_e[k])
            eps = sum(hist_e[k] for k in hist_e.keys() if abs(k - E_peak) <= 1.0) / n
            nodes[E0] = (shape, eps, E_peak - E0, st["n_used"], st["n_in"])
    if not nodes:
        raise SystemExit(f"Не найдено файлов в {build_out} с префиксом {prefix}")
    return nodes, skipped, empty

def make_resp(nodes, tol=0.01):
    if not nodes:
        raise SystemExit("Пустой словарь узлов")
    def resp(E):
        closest = min(nodes.keys(), key=lambda x: abs(x - E))
        if abs(E - closest) > tol:
            raise SystemExit(f"Запрошенная энергия {E} не совпадает с ближайшей {closest}")
        shape, eps, _, _, _ = nodes[closest]
        return (shape, {}, eps)
    return resp

def summary(nodes, empty=None):
    if not nodes:
        return "Нет узлов"
    shifts = [shift for _, (_, _, shift, _, _) in nodes.items() if not math.isnan(shift)]
    lines = [f"Число узлов: {len(nodes)}"]
    if shifts:
        lines.append(f"Минимальное смещение: {min(shifts):.3f} кэВ")
        lines.append(f"Среднее смещение: {sum(shifts) / len(shifts):.3f} кэВ")
        lines.append(f"Максимальное смещение: {max(shifts):.3f} кэВ")
    else:
        lines.append("Нет узлов с конечным смещением")
    large_shifts = [f"{E}" for E, (_, _, shift, _, _) in nodes.items() if not math.isnan(shift) and abs(shift) > 2.0]
    no_peak = [f"{E}" for E, (_, _, shift, _, _) in nodes.items() if math.isnan(shift)]
    if large_shifts:
        lines.append(f"Энергии со смещением > 2 кэВ: {', '.join(large_shifts[:10])}")
    if no_peak:
        lines.append(f"Энергии без пика: {', '.join(no_peak)}")
    if empty:
        lines.append(f"Энергии без зарегистрированных событий: {', '.join(str(E) for E in empty)}")
    return "\n".join(lines)

def main():
    # Числа подобраны так, что проверка 1 краснеет на обеих порчах: eps вокруг E_peak=98 даёт 0,100; вокруг номинала E0=100 — 0,040; при окне 5,0 кэВ вместо 1,0 — 0,155.
    e = np.arange(200) * 1.0
    file_e = lambda c: 1.0 * c
    fwhm = lambda E: 5.0
    a, b, blur = 0.0, 1.0, 1.0
    base = {96: 10, 98: 100, 100: 40, 102: 5}
    rt = lambda p: (dict(base), 1000, 1)
    rt_zero = lambda p: ({**base, 0.5: 9000}, 1000, 1)
    rt_over = lambda p: ({0.5: 5000, 96: 10, 98: 700, 100: 40, 102: 5}, 1000, 1)
    rt_npsm0 = lambda p: (dict(base), 1000, 0)
    rt_empty = lambda p: ({0.5: 1000}, 1000, 1)
    rt_far = lambda p: ({94: 200, 95: 50, 99: 5}, 1000, 1)
    def light_to_energy(hist, a, b, e_of_ch, n_channels):
        new_hist = {float(e_of_ch(a + b * k)): v for k, v in hist.items()}
        s = sum(hist.values())
        return new_hist, {"n_in": s, "n_used": sum(new_hist.values()), "sum_in": s, "sum_used": sum(new_hist.values())}
    def broaden_ch(hist, n_events, e, fwhm, file_e, blur):
        return np.zeros(len(e))
    run = lambda read_tpl, paths: build_nodes(".", "gridon_E", read_tpl, light_to_energy, broaden_ch, e, file_e, fwhm, a, b, blur, paths=paths)
    # 1. Положение пика: максимум стоит на 98, а не на номинале 100.
    nodes, skipped, empty = run(rt, ["gridon_E100.0.csv"])
    if len(nodes) != 1 or abs(sorted(nodes)[0] - 100.0) > 1e-9:
        raise SystemExit(f"SELFTEST FAIL: положение пика — ожидался ровно один узел 100.0, получено {sorted(nodes)}")
    shape, eps, shift, n_used, n_in = nodes[100.0]
    if abs(shift + 2.0) > 1e-9:
        raise SystemExit(f"SELFTEST FAIL: положение пика — ожидалось смещение -2.0 кэВ, получено {shift}")
    if abs(eps - 0.1) > 1e-9:
        raise SystemExit(f"SELFTEST FAIL: пиковая эффективность — ожидалось 0.1, получено {eps}")
    # 2. Бин без отложения отброшен: счёт и эффективность не изменились, отказа по порогу нет.
    nodes2, _, _ = run(rt_zero, ["gridon_E100.0.csv"])
    shape2, eps2, shift2, n_used2, n_in2 = nodes2[100.0]
    if abs(eps2 - eps) > 1e-9 or abs(shift2 - shift) > 1e-9 or n_used2 != n_used or n_in2 != n_in:
        raise SystemExit(f"SELFTEST FAIL: отбрасывание бина без отложения — ожидалось eps={eps}, смещение={shift}, n_used={n_used}, n_in={n_in}; получено eps={eps2}, смещение={shift2}, n_used={n_used2}, n_in={n_in2}")
    # 3. Отказ по порогу 0,6 (ожидается SystemExit).
    failed = False
    try:
        run(rt_over, ["gridon_E100.0.csv"])
    except SystemExit:
        failed = True
    if not failed:
        raise SystemExit("SELFTEST FAIL: отказ по порогу 0,6 — при доле 0.755 отказа не было")
    # 4. Отказ на прогоне без непропорциональности, npsm=0 (ожидается SystemExit).
    failed = False
    try:
        run(rt_npsm0, ["gridon_E100.0.csv"])
    except SystemExit:
        failed = True
    if not failed:
        raise SystemExit("SELFTEST FAIL: npsm=0 — отказа не было")
    # 5. Разбор имён: дробная часть не отрезана, посторонний файл пропущен.
    nodes5, skipped5, _ = run(rt, ["gridon_E1157.022.csv", "gridon_eplusg_E1157.022.csv"])
    if len(nodes5) != 1 or abs(sorted(nodes5)[0] - 1157.022) > 1e-9:
        raise SystemExit(f"SELFTEST FAIL: разбор имён — ожидался ровно один узел 1157.022, получено {sorted(nodes5)}")
    if skipped5 != ["gridon_eplusg_E1157.022.csv"]:
        raise SystemExit(f"SELFTEST FAIL: разбор имён — ожидался пропуск ['gridon_eplusg_E1157.022.csv'], получено {skipped5}")
    # 6. Пустой узел: после отбрасывания значений <1 не осталось ничего — обязан быть узел, а не отказ.
    refused = None
    try:
        nodes6, _, empty6 = run(rt_empty, ["gridon_E11.89.csv"])
    except SystemExit as ex:
        refused = ex
    if refused is not None:
        raise SystemExit(f"SELFTEST FAIL: пустой узел — ожидался нулевой узел 11.89, получен отказ: {refused}")
    if sorted(nodes6) != [11.89] or empty6 != [11.89]:
        raise SystemExit(f"SELFTEST FAIL: пустой узел — ожидались узлы [11.89] и empty=[11.89], получено {sorted(nodes6)} и empty={empty6}")
    shape6, eps6, shift6, n_used6, n_in6 = nodes6[11.89]
    if eps6 != 0.0 or not math.isnan(shift6) or len(shape6) != len(e):
        raise SystemExit(f"SELFTEST FAIL: пустой узел — ожидалось eps=0.0, смещение=nan, длина формы {len(e)}; получено eps={eps6}, смещение={shift6}, длина {len(shape6)}")
    # 7. Дальний пик: 94,0 отстоит от номинала на 6,0 кэВ — вне старого окна 3,75, внутри нового 7,5.
    # Рядом с номиналом стоит мелкий ключ 99 (шум): со старым окном нашёлся бы он, смещение -1,0.
    nodes7, _, _ = run(rt_far, ["gridon_E100.0.csv"])
    shape7, eps7, shift7, n_used7, n_in7 = nodes7[100.0]
    if abs(shift7 + 6.0) > 1e-9:
        raise SystemExit(f"SELFTEST FAIL: дальний пик — ожидалось смещение -6.0 кэВ (пик 94.0), получено {shift7}")
    if abs(eps7 - 0.25) > 1e-9:
        raise SystemExit(f"SELFTEST FAIL: дальний пик — ожидалось eps=0.25 (ключи 94 и 95), получено {eps7}")
    print("SELFTEST OK")

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--selftest":
        main()
