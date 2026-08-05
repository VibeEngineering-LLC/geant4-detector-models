# -*- coding: utf-8 -*-
"""Кривая эффективности по ППП из расчётных спектров сетки.

Наблюдаемая объявляется явно. Эффективность здесь — АБСОЛЮТНАЯ ПО ПИКУ ПОЛНОГО
ПОГЛОЩЕНИЯ на 4pi: доля испущенных источником квантов, дающих отсчёт в пике.
Прогон идёт конусом, поэтому счёт делится на N первичных и умножается на долю
телесного угла конуса, которую печатает сам exe (solid_angle_frac).

Площадь пика выводится ДВУМЯ способами, и обе величины пишутся в файл, потому
что они отвечают на разные вопросы:

  eps_peak     — СТРОГАЯ: в нерасплывшемся модельном спектре все события с
                 полным поглощением попадают ровно в канал E0, поэтому счёт
                 берётся из него напрямую. Конвенции нет вовсе: ни окна, ни
                 полок, ни исключений. Для модели это точное число, для
                 измерения такое невозможно;
  eps_peak_win — ОКОННАЯ: спектр размывается приборной ПШПВ и площадь берётся
                 `becqmoni.area_broadened` с полками, то есть той же
                 конвенцией, какой площадь снимается с ИЗМЕРЕННОГО спектра.
                 Сравнивать модель с кривой, построенной по измерениям, надо
                 именно этой величиной.

Из полок исключаются каналы пиков вылета характеристического рентгеновского
излучения кристалла (CsI: 28,6 и 30,97 кэВ) — иначе на мягком крае полка
садится на собственную особенность спектра и площадь занижается односторонне.

ВНИМАНИЕ к столбцу `shelf`. При исключении вылета левой полки может не
остаться, и `area_broadened` переходит на правую (режим `right`). Для
МОНОЭНЕРГЕТИЧЕСКОГО модельного спектра правая полка почти пуста, подложка
оценивается около нуля и оконная площадь завышается. На сетке 40–3000 кэВ это
даёт немонотонную ступень в области 122–200 кэВ — артефакт съёма, а не
свойство прибора. Строгая величина от него свободна, поэтому именно она
объявлена основной.

ПШПВ прибора взята из СОБСТВЕННОЙ записи Cs-137 (подгонка «гауссиана +
линейная подложка», `compare_cs137.py`): 41,60 кэВ на 661,657 кэВ, то есть
6,4 %. Закон ПШПВ(E) = ПШПВ(662)*sqrt(E/661,657) — приближение по одной
опорной точке; на краях диапазона оно даёт систематику формы пика, но не
площади, пока окно берётся в долях самой ПШПВ.

    python analysis/export_curve.py <каталог со спектрами> [<выход.csv>]
"""
import os
import sys
import math
import glob

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(_HERE, "..", "..", "..",
                                                 "common", "py")))
import becqmoni as bm   # noqa: E402

FWHM_662 = 41.60           # собственная запись Cs-137, подгонка пика
ESCAPES = (28.6, 30.97)    # K-вылет иода и цезия, кэВ (becqmoni.area_broadened)


def read_model(path):
    head, dic = {}, {}
    with open(path, encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("#"):
                if "=" in ln:
                    k, v = ln[1:].split("=", 1)
                    head[k.strip()] = v.strip()
            elif "," in ln and not ln.startswith("E_keV"):
                a, b = ln.split(",")
                dic[float(a)] = dic.get(float(a), 0.0) + float(b)
    return dic, head


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print("Нет каталога спектров. Сначала прогон:\n"
              "  asn16 curve_point_end10cm.mac\n"
              "он кладёт спектры в подкаталог spectra/ каталога сборки.")
        return 2
    src = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        _HERE, "..", "results", "eff_point_end10cm.csv")
    out = os.path.normpath(out)
    files = sorted(glob.glob(os.path.join(src, "*.csv")))
    files = [f for f in files if not f.endswith("_emit.csv")]
    if not files:
        print("В %s нет спектров." % src)
        return 2

    rows, stamps = [], set()
    for f in files:
        dic, head = read_model(f)
        e0 = float(head["E_prim_keV"])
        n = int(head["N_primaries"])
        frac = float(head["solid_angle_frac"])
        stamps.add(head.get("src_sha1", "?"))
        total = sum(dic.values())
        fwhm = FWHM_662 * math.sqrt(e0 / 661.657)
        # СТРОГАЯ площадь: канал полного поглощения нерасплывшегося спектра.
        # Гистограмма пишется как b = int(E/шаг), центр канала b + 0,5 —
        # поэтому канал ищется тем же floor, а не округлением к ближайшему:
        # при целом E0 (40, 50, 80 ...) ближайших каналов ДВА на равном
        # расстоянии, и выбор «первого попавшегося» уводил счёт в соседний
        # канал континуума. Расхождение доходило до трёх раз и выглядело как
        # физика: узлы с дробной энергией (59,541; 122,061; 661,657) давали
        # правильное значение, целые — заниженное.
        #
        # Берётся не один канал, а ±1,5 кэВ: сумма энерговыделений равна E0
        # лишь с точностью плавающей точки, и часть событий при целом E0
        # попадает в канал ниже. Континуум непосредственно под линией
        # пренебрежимо мал (рассеяние на малые углы), поэтому окно в три
        # канала конвенцией не является.
        kc = math.floor(e0) + 0.5
        fep = sum(c for e, c in dic.items() if abs(e - kc) <= 1.5)
        if fep <= 0:
            print("  ! %s: канал полного поглощения пуст при E0 = %.1f, "
                  "узел пропущен" % (os.path.basename(f), e0))
            continue
        wide = bm.broaden(dic, fwhm_at_662=FWHM_662)
        detail = {}
        area, _bgd = bm.area_broadened(wide, e0, fwhm, escapes=ESCAPES,
                                       detail=detail)
        rows.append(dict(E=e0, fep=fep, area=max(area, 0.0), n=n, frac=frac,
                         total=total, fwhm=fwhm,
                         shelf=detail.get("mode", "?")))

    if len(stamps) > 1:
        print("! В каталоге спектры от РАЗНЫХ сборок: %s\n"
              "  Кривая по смешанным прогонам недействительна." %
              ", ".join(sorted(stamps)))
        return 1

    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8", newline="") as g:
        g.write("# AtomSpectra Nano 16 PRO, точечный источник на оси кристалла,"
                " 10 см от наружной поверхности передней крышки (торец 18x15)\n")
        g.write("# наблюдаемая: абсолютная эффективность по ППП на 4pi\n")
        g.write("# src_sha1 = %s\n" % stamps.pop())
        g.write("# fwhm_662_keV = %.2f  (собственная запись Cs-137)\n"
                % FWHM_662)
        g.write("# escapes_keV = %s\n" % ", ".join("%.2f" % e for e in ESCAPES))
        g.write("# eps_peak — строгая (канал полного поглощения), основная\n")
        g.write("# eps_peak_win — оконная, конвенция измеренного спектра;"
                " при shelf=right завышена\n")
        g.write("E_keV,eps_peak,d_eps_peak,eps_peak_win,eps_total,fep_counts,"
                "area_counts,N_primaries,solid_angle_frac,fwhm_keV,shelf\n")
        for r in rows:
            eps = r["fep"] / r["n"] * r["frac"]
            d = eps / math.sqrt(max(r["fep"], 1.0))
            epsw = r["area"] / r["n"] * r["frac"]
            epst = r["total"] / r["n"] * r["frac"]
            g.write("%.3f,%.6e,%.3e,%.6e,%.6e,%.0f,%.0f,%d,%.8f,%.2f,%s\n"
                    % (r["E"], eps, d, epsw, epst, r["fep"], r["area"],
                       r["n"], r["frac"], r["fwhm"], r["shelf"]))
    print("записано: %s  (%d узлов)" % (out, len(rows)))
    print("%9s %12s %8s %12s %12s %7s %7s"
          % ("E, кэВ", "eps ППП", "стат,%", "eps окном", "eps полная",
             "пик/пол", "полка"))
    for r in rows:
        eps = r["fep"] / r["n"] * r["frac"]
        epsw = r["area"] / r["n"] * r["frac"]
        epst = r["total"] / r["n"] * r["frac"]
        print("%9.1f %12.4e %7.2f %12.4e %12.4e %7.3f %7s"
              % (r["E"], eps, 100.0 / math.sqrt(max(r["fep"], 1.0)),
                 epsw, epst, eps / epst, r["shelf"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
