# -*- coding: utf-8 -*-
"""Проверка энергетической калибровки записей комплекта и цена её дефекта.

ЗАЧЕМ. Площадь пика снималась окном, поставленным на ТАБЛИЧНУЮ энергию линии, а
калибровка бралась из файла как данность. Обе посылки надо проверять:

  1. Калибровка из файла может быть просто негодной — немонотонная шкала,
     вырожденная модель ПШПВ, полуширина вне физически возможной полосы.
  2. Даже при годной калибровке центроида пика может уехать от номинала. Тогда
     окно площади стоит не по центру: часть площади срезается, а полки фона
     заходят на склон пика. Оба эффекта смещают результат ВНИЗ и не видны ни по
     невязке подгонки, ни по статистике.

СВОИМИ СИЛАМИ ЭТО НЕ ДЕЛАЕТСЯ. У оператора есть SpectraVibe с готовыми
инструментами, и правило проекта — брать их, а не писать третью реализацию:

  gamma.calibration.calibration_gate.evaluate_calibration_gate — проверка самих
      коэффициентов калибровки, до всякого поиска пиков (монотонность шкалы,
      диапазон, монотонность и полоса ПШПВ);
  gamma.peaks.centroid_gost — центроида по ГОСТ 26874-86 §3.3.3 и вычитание
      пьедестала по §3.3.2, то есть канонические формулы, а не самодельные.

Нужна переменная SPECTRAVIBE_ROOT. Без неё скрипт скажет об этом и посчитает
центроиду своим грубым способом (becqmoni.peak_find), пометив результат как
несертифицированный — но проверку калибровочных коэффициентов сделать будет
нечем.

    python detectors/Gamma-1S/analysis/check_calibration.py

Порог внимания: сдвиг больше 0,15 ПШПВ либо изменение площади больше 1 %.
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import becqmoni as bm  # noqa: E402

SV = os.environ.get("SPECTRAVIBE_ROOT")
GOST = GATE = READ = None
if SV and os.path.isdir(os.path.join(SV, "scripts")):
    sys.path.insert(0, os.path.join(SV, "scripts"))
    try:
        from gamma.peaks import centroid_gost as GOST          # noqa: N812
        from gamma.calibration.calibration_gate import (        # noqa: N812
            evaluate_calibration_gate as GATE)
        from gamma.io.readers import read_spectrum as READ      # noqa: N812
    except Exception as exc:                                    # noqa: BLE001
        print("SpectraVibe найден, но модули не импортируются: %s" % exc)

# Аналитические линии по нуклиду записи — те же, по которым считаются активности
LINES = {
    "M_cs": [661.657], "M_k": [1460.822],
    "M_ra": [295.223, 351.932, 609.32, 1120.294],
    "Cs137": [661.657], "K40": [1460.822],
    "Ra226": [295.223, 351.932, 609.32, 1120.294, 1764.491],
    "Th232": [238.632, 583.187, 911.204, 2614.511],
    "Am241": [59.54], "Ba133": [80.997, 356.013], "Bi207": [569.7, 1063.66],
    "Cd109": [88.03], "Ce139": [165.86], "Co57": [122.06],
    "Co60": [1173.23, 1332.49], "Eu152": [121.78, 344.28, 1408.01],
    "Mn54": [834.85], "Na22": [1274.54], "Th228": [238.632, 583.187, 2614.511],
    "Y88": [898.04, 1836.06], "Zn65": [1115.54],
}
SHIFT_WARN, AREA_WARN = 0.15, 1.0


def lines_for(name):
    for key, ls in LINES.items():
        if key.lower() in name.lower():
            return ls
    return []


def gate_of(path):
    """Вердикт калибровочного шлюза SpectraVibe или None, если он недоступен."""
    if not (GATE and READ):
        return None
    try:
        return GATE(READ(path))
    except Exception as exc:                                    # noqa: BLE001
        return ("ошибка", str(exc))


def centroid_gost(sp, E0, fwhm_keV):
    """Центроида по ГОСТ 26874-86 §3.3.3.2, кэВ. None — если не вышло."""
    if GOST is None:
        return None
    try:
        ch0 = int(round(sp.channel(E0)))
        # ПШПВ в каналах: по производной калибровки в точке пика
        dE = sp.energy(ch0 + 1) - sp.energy(ch0)
        if dE <= 0:
            return None
        fw_ch = fwhm_keV / dE
        lo = max(0, int(ch0 - 3 * fw_ch))
        hi = min(len(sp.n), int(ch0 + 3 * fw_ch) + 1)
        if hi - lo < 7:
            return None
        counts = sp.n[lo:hi]
        ped = GOST.gost_select_pedestal_method(counts, ch0 - lo, fw_ch)
        net = ped.counts_net if ped.counts_net is not None else counts
        res = GOST.gost_centroid_weighted_mean(net, channel_offset=lo)
        if res is None or res.centroid_channel is None:
            return None
        return float(sp.energy(res.centroid_channel))
    except Exception:                                           # noqa: BLE001
        return None


def main():
    root = paths.ref("Gamma-1S")
    files = sorted(str(p) for p in root.rglob("*.xml"))
    if not files:
        print("Нет спектров комплекта в %s" % root)
        return 1
    if GOST is None:
        print("SPECTRAVIBE_ROOT не задан: канонической центроиды по ГОСТ и\n"
              "проверки калибровочных коэффициентов НЕ БУДЕТ, центроида\n"
              "посчитана грубой своей оценкой. Путь к SpectraVibe — в\n"
              "переменной SPECTRAVIBE_ROOT.\n")

    print("%-32s %8s %8s %7s %7s %8s %8s" %
          ("запись", "E ном.", "центр.", "сдвиг", "в ПШПВ", "ПШПВ %", "Δплощ."))
    rows, worst, gates = [], [], {}
    for p in files:
        name = os.path.basename(p)
        ls = lines_for(name)
        if not ls:
            continue
        g = gate_of(p)
        if g is not None and not isinstance(g, tuple):
            key = "прошла" if g.passed else "НЕ ПРОШЛА: " + str(g.reason)[:60]
            gates.setdefault(key, []).append(name)
        try:
            s, _ = bm.read(p)
        except Exception as exc:                                # noqa: BLE001
            print("%-32s не читается: %s" % (name[:32], exc))
            continue
        for E in ls:
            f = bm.peak_find(s, E)
            if f is None:
                continue
            cen_own, fw, _top = f
            if fw <= 0:
                continue
            cen = centroid_gost(s, E, fw)
            src = "ГОСТ"
            if cen is None:
                cen, src = cen_own, "своя"
            d = cen - E
            a_nom = bm.peak_area(s, E, fw, roi=1.0, side=1.0)
            a_cen = bm.peak_area(s, cen, fw, roi=1.0, side=1.0)
            if not a_nom or not a_cen or a_nom[0] <= 0:
                continue
            da = 100.0 * (a_cen[0] / a_nom[0] - 1.0)
            rows.append((name, E, cen, d, d / fw, 100 * fw / E, da, src))
            flag = "  <-" if (abs(d / fw) > SHIFT_WARN
                              or abs(da) > AREA_WARN) else ""
            print("%-32s %8.1f %8.1f %+7.2f %+7.3f %8.1f %+7.2f%s"
                  % (name[:32], E, cen, d, d / fw, 100 * fw / E, da, flag))
            if flag:
                worst.append((abs(da), name, E, d, da))

    if not rows:
        print("Ни одной линии не найдено.")
        return 1
    sh = [r[4] for r in rows]
    dk = [r[3] for r in rows]
    da = [r[6] for r in rows]
    src = {r[7] for r in rows}
    n = len(rows)
    mean_sh = sum(sh) / n
    sd = math.sqrt(sum((x - mean_sh) ** 2 for x in sh) / max(1, n - 1))
    print("\nлиний проверено: %d, центроида считалась: %s"
          % (n, ", ".join(sorted(src))))
    print("сдвиг центроида: средний %+.3f ПШПВ (СКО %.3f), размах %+.3f..%+.3f"
          % (mean_sh, sd, min(sh), max(sh)))
    print("в кэВ:           средний %+.2f, размах %+.2f..%+.2f"
          % (sum(dk) / n, min(dk), max(dk)))
    print("цена для площади: средняя %+.2f %%, размах %+.2f..%+.2f %%"
          % (sum(da) / n, min(da), max(da)))
    print("значимость среднего сдвига: %.1f сигмы среднего"
          % (abs(mean_sh) / (sd / math.sqrt(n)) if sd else 0.0))

    if gates:
        print("\nкалибровочный шлюз SpectraVibe:")
        for k, v in gates.items():
            print("   %-40s записей %d" % (k[:40], len(v)))
            if not k.startswith("прошла"):
                for x in v[:5]:
                    print("        ", x[:60])
    if worst:
        worst.sort(reverse=True)
        print("\nхудшие по цене:")
        for a, nm, E, d, dd in worst[:8]:
            print("   %-32s %.1f кэВ: сдвиг %+.2f кэВ, площадь %+.2f %%"
                  % (nm[:32], E, d, dd))
    print("\nОдносторонний средний сдвиг — системный дефект калибровки;\n"
          "разбросанный вокруг нуля — статистика поиска центроиды. Судить\n"
          "надо по среднему и его значимости, а не по отдельным записям.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
