"""Экспорт расчётных кривых в формат .efa ЛСРМ СпектраЛайн (Трек 1a).

Зачем. Оператор запускает СпектраЛайн с ДВУМЯ кривыми эффективности —
штатной поверочной и нашей расчётной — на одних и тех же спектрах. Вместе с
нашим конвейером это даёт дизайн 2×2 (кривая × алгоритм), раскладывающий
расхождение активностей на вклад кривой и вклад алгоритма обработки.

Как. Шаблоном служит НАСТОЯЩИЙ .efa той же геометрии из reference/lsrm
(версия 1.7.11918): у него берётся весь заголовок (детектор, геометрия,
объём, матрица, толщины) — так файл гарантированно читается программой.
Подменяются только строки узлов: энергия = наша eps_net, статистическая
погрешность, метка «G4MC» вместо нуклида, площадь = счёт узла сетки.

Секция Zones/Curve_* (полиномиальная аппроксимация зон) в экспорт НЕ
переносится: коэффициенты шаблона описывают ЧУЖУЮ кривую, а конвенция их
записи не документирована. Кривую по нашим узлам оператор перестраивает в
редакторе эффективности СпектраЛайн (штатная операция по точкам файла).

Конвенция eps. Экспортируется eps_net — площадь ±6 кэВ за вычетом левой
полки, на неразмытом спектре; ровно с этим числом весь проект сверялся с
ЛСРМ. Разница конвенций площади (наше окно против фита ±2,5 ПШПВ со
ступенькой) — предмет Трека 1c, и эксперимент 2×2 её и измеряет.

Обезличивание. Имя детектора берётся из шаблона (SN-01). Для запуска на
реальном приборе оператор задаёт переменную G1S_EFA_DETECTOR со своим
именем детектора — локально, в репозиторий такие файлы не коммитятся.

    python detectors/Gamma-1S/analysis/export_efa.py
    выход: results/efa_export/G4MC_<геометрия>.efa (CP-1251, CRLF)
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

RESULTS = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results"))

# (сетка -> файл кривой, геометрия ЛСРМ как в reference/lsrm)
EXPORTS = [
    ("rho1.60", "Маринелли"),
    ("denta1.60", "Дента"),
    ("petri1.60", "Петри"),
    ("p25cm", "Точечная-25см"),
]
# Точечная-5см: измеренного .efa в reference нет (только .efr) — шаблона
# заголовка для неё не существует, экспорт с шаблоном 25 см дал бы ложные
# метаданные расстояния. Отложено до появления шаблона.

DETECTOR = os.environ.get("G1S_EFA_DETECTOR")   # None -> из шаблона


def read_curve(grid):
    """[(E, eps_net, d_eps, net_counts)] из results/eff_<grid>.csv"""
    p = os.path.join(RESULTS, "eff_%s.csv" % grid)
    out = []
    for ln in open(p, encoding="utf-8"):
        if ln.startswith("#") or ln.startswith("E_keV") or not ln.strip():
            continue
        f = ln.strip().split(",")
        out.append((float(f[0]), float(f[1]), float(f[2]), float(f[5])))
    return out


def template_lines(geom):
    p = paths.efficiency_curve(geom, ext="efa")
    raw = open(str(p), "rb").read().decode("cp1251")
    return raw.replace("\r\n", "\n").split("\n")


if __name__ == "__main__":
    outdir = os.path.join(RESULTS, "efa_export")
    os.makedirs(outdir, exist_ok=True)
    for grid, geom in EXPORTS:
        tpl = template_lines(geom)
        head = []
        for ln in tpl:
            if not ln:
                continue
            key = ln.split("=", 1)[0]
            # конец заголовка: первая строка-узел (ключ — число) или зоны
            try:
                float(key)
                break
            except ValueError:
                pass
            if key in ("Zones",) or key.startswith(("Zone_", "Curve_")):
                break
            head.append(ln)
        if DETECTOR:
            newh = []
            det_old = None
            for ln in head:
                if ln.startswith("Detector="):
                    det_old = ln.split("=", 1)[1]
            for ln in head:
                if ln.startswith("[") and det_old:
                    ln = ln.replace(det_old, DETECTOR)
                if ln.startswith("Detector="):
                    ln = "Detector=" + DETECTOR
                newh.append(ln)
            head = newh

        rows = []
        for E, eps, d, cnt in read_curve(grid):
            dpct = 100.0 * d / eps if eps > 0 else 99.0
            rows.append("%g=%.6E,%.3f,G4MC,%.3f,%.3f,0"
                        % (E, eps, dpct, cnt, cnt * d / eps if eps else 0))

        body = "\r\n".join(head + rows) + "\r\n"
        fn = os.path.join(outdir, "G4MC_%s.efa" % geom)
        open(fn, "wb").write(body.encode("cp1251"))
        print("%-14s %2d узлов -> %s" % (geom, len(rows), fn))
    print("\nZones/Curve не экспортируются: кривую по узлам перестроить в"
          "\nредакторе эффективности СпектраЛайн (см. шапку скрипта).")
