"""Решающая проверка: чья эффективность верна на 662 кэВ в маринельке.

Возникло противоречие между двумя опорами.
  (а) Кривая ЛСРМ .efr для маринелли: eps(662) = 1,871e-2, а расчёт даёт
      2,234e-2 -> расчёт ВЫШЕ в 1,19 раза.
  (б) Разложение смеси РИСН-379 по паспортной активности: расчёт выше лишь
      в 1,05 раза (Cs-137 1519 против паспортных 1600 Бк/кг в маринельке,
      1678 против 1600 в Денте).
При этом точечные геометрии по кривым ЛСРМ дают 0,971 и 0,931, то есть
согласие. Выпадает именно маринелльная кривая.

Здесь противоречие решается напрямую и без всякой подгонки: берётся ТОТ САМЫЙ
источник, по которому ЛСРМ строил маринелльную кривую — Cs137_420-7-15, —
и его эффективность считается из площади пика и паспортной активности:
    eps_изм = R_пика / (A * p_gamma)
Никакого NNLS, никакой кривой ЛСРМ, одна линия.

Файлы читаются штатным ридером ЛСРМ-формата из проекта оператора
(scripts/gamma/io/lsrm_spe.py) — в .spe отсчёты упакованы двоично, своим
разбором их не взять.
"""
import math
import os
import re
import sys
from datetime import date

import numpy as np

# Корни путей — из переменных окружения (common/py/paths.py), чтобы в коде не
# было ни одного пути, привязанного к конкретной машине.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

RAG = str(paths.require_spectravibe("чтение сырых .spe ЛСРМ"))
sys.path.insert(0, os.path.join(RAG, "scripts"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gamma.io.lsrm_spe import read_lsrm_spe  # noqa: E402
import becqmoni as bm  # noqa: E402
from contam import dirty_shelves  # noqa: E402

VER = os.path.join(RAG, "detectors", "Gamma-1S", "raw_lsrm", "Work", "BG",
                   "Gamma-1S", "Spe - поверки", "Поверка 2016")
BUILD = str(paths.build("Gamma-1S"))

if not os.path.isdir(BUILD):
    raise SystemExit(
        "Нет каталога расчётных спектров %s.\n"
        "Они не коммитятся (сотни файлов), а воспроизводятся драйверами:\n"
        "    python detectors/Gamma-1S/drivers/run_grid.py\n"
        "    python detectors/Gamma-1S/drivers/run_all_grids.py\n"
        "Либо укажите G4MODELS_BUILD_GAMMA_1S на готовый каталог."
        % BUILD)
FWHM662 = 49.9

# Источники, по которым ЛСРМ строил маринелльную кривую (секции .efr).
# Активность и дата — из поля COMMENT самих файлов, читаются программно.
CASES = [
    ("Маринелли/Cs137_420-7-15_Маринелли.spe", 661.657, "Cs-137", 30.08, 0.8513),
    ("Маринелли/K40_420-7-21_Маринелли.spe", 1460.822, "K-40", 1.248e9, 0.1065),
    ("Маринелли/Th232_420-7-17_Маринелли.spe", 583.187, "Th-232", 1.405e10, 0.3089),
    ("Маринелли/Th232_420-7-17_Маринелли.spe", 2614.511, "Th-232", 1.405e10, 0.3585),
    ("Маринелли/Ra226_420-7-19_Маринелли.spe", 609.32, "Ra-226", 1600.0, 0.4601),
]
BG = "Фон вода/фон вода_13.spe"
# Списки линий и проверка полок — в общем модуле contam.py


class Wrap:
    """Переходник к becqmoni: те же поля, что у его Spectrum."""
    def __init__(self, sp):
        self.n = np.asarray(sp.counts, dtype=float)
        self.live = float(sp.live_time)
        self.real = float(sp.real_time)
        c = getattr(sp, "energy_cal", None)
        self.cal = list(c) if c else [0.0, 1.0]

    def energy(self, ch):
        ch = np.asarray(ch, dtype=float)
        return sum(c * ch ** k for k, c in enumerate(self.cal))

    def channel(self, E):
        ch = np.arange(len(self.n), dtype=float)
        return float(np.interp(E, self.energy(ch), ch))

    def counts_between(self, E0, E1):
        a, b = self.channel(E0), self.channel(E1)
        lo, hi = max(0, int(math.floor(a))), min(len(self.n), int(math.ceil(b)))
        return float(self.n[lo:hi].sum()), max(1, hi - lo)


def mc_eff(E, tag="rho1.60"):
    import glob
    import re
    for p in glob.glob(os.path.join(BUILD, "grid", tag + "_E*.csv")):
        m = re.search(r"_E(\d+\.\d)\.csv$", p)
        if not m or abs(float(m.group(1)) - E) > 1.0:
            continue
        hist, N = {}, None
        for line in open(p, encoding="utf-8"):
            if line.startswith("#"):
                if "N_primaries" in line:
                    N = int(line.split("=")[1])
                continue
            if line and line[0].isdigit():
                e, c = line.split(",")
                hist[float(e)] = int(c)
        peak = sum(c for e, c in hist.items() if abs(e - E) <= 6.0)
        return peak / N if N else None
    return None


if __name__ == "__main__":
    bgp = os.path.join(VER, BG.replace("/", os.sep))
    bgs = Wrap(read_lsrm_spe(bgp)) if os.path.exists(bgp) else None
    print("Решающая проверка на источниках, по которым построена кривая ЛСРМ\n")
    if bgs:
        print("фон: %s, живое %.0f с\n" % (BG, bgs.live))
    print("%-9s %8s %8s %8s %10s %11s %11s %7s %9s" %
          ("нуклид", "E, кэВ", "имп/с", "ро", "A, Бк", "eps изм", "eps МК",
           "МК/изм", ".efr"))
    EFR = {661.7: 1.8713e-2, 1460.8: 9.74254e-3, 583.2: 2.17122e-2,
           2614.5: 4.71386e-3, 609.3: 1.98283e-2}
    for rel, E, nuc, t12, pg in CASES:
        p = os.path.join(VER, rel.replace("/", os.sep))
        if not os.path.exists(p):
            print("%-9s нет файла %s" % (nuc, rel))
            continue
        sp = read_lsrm_spe(p)
        s = Wrap(sp)
        # COMMENT приходит с разрядкой (двухбайтовая запись) — сжимаем
        com = re.sub(r"\s+", "", " ".join(sp.comments or []))
        m_a = re.search(r"A=([\d.]+)", com)
        m_d = re.search(r"(\d\d)-(\d\d)-(\d{4})", com)
        mass = float(getattr(sp, "sample_mass_kg", 0) or 0)
        if not (m_a and mass):
            print("%-9s не разобран паспорт: %s" % (nuc, com[:50]))
            continue
        aspec = float(m_a.group(1))
        # ЛОВУШКА: у этих файлов end_datetime и file_created_datetime = None,
        # и без явной проверки поправка на распад молча получалась равной 1.
        # Для Cs-137 это ошибка в 18 %. Дата измерения — из MEASBEGIN
        # заголовка (extras) или из имени папки поверки.
        md = sp.end_datetime or sp.file_created_datetime
        if md is None:
            md = date(2016, 5, 30)          # «Поверка 2016», май
        else:
            md = md.date()
        k = 1.0
        if m_d:
            d0 = date(int(m_d.group(3)), int(m_d.group(2)), int(m_d.group(1)))
            k = 0.5 ** (((md - d0).days / 365.25) / t12)
        A = aspec * mass * k
        fw = FWHM662 * math.sqrt(E / 661.657)
        a = bm.peak_area(s, E, fw, roi=1.0, side=1.0)
        if a is None:
            print("%-9s ROI вне спектра" % nuc)
            continue
        rate = a[0] / s.live
        if bgs:
            ab = bm.peak_area(bgs, E, fw, roi=1.0, side=1.0)
            if ab:
                rate -= ab[0] / bgs.live
        # Проверка чистоты ФОНОВЫХ ПОЛОК: если в полку попадает сильная линия,
        # трапеция завышает подложку и площадь пика срезается. Так вышло на
        # 583,2 у тория — в левую полку (489–536) садится линия Tl-208 510,77.
        dirty = dirty_shelves(nuc, E, fw)
        eps = rate / (A * pg)
        mc = mc_eff(E)
        efr = EFR.get(round(E, 1))
        if dirty:
            print("%-9s %8.1f  ПОЛКА ЗАГРЯЗНЕНА линией %s кэВ — площадь "
                  "недостоверна, точка исключена"
                  % (nuc, E, ", ".join("%.1f" % x for x in dirty)))
            continue
        print("%-9s %8.1f %8.3f %8.2f %10.0f %11.4e %11s %7s %9s"
              % (nuc, E, rate, mass, A, eps,
                 "%.4e" % mc if mc else "-",
                 "%.3f" % (mc / eps) if mc else "-",
                 "%.3f" % (eps / efr) if efr else "-"))
    print("\nСтолбец «.efr» — отношение моего измерения к записанному ЛСРМ:")
    print("если он около 1,00, кривая ЛСРМ воспроизводима из сырых файлов и")
    print("паспортов, то есть верна, и расхождение принадлежит РАСЧЁТУ.")
