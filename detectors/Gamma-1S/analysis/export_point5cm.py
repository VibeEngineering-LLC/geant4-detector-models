"""Кривая точечной 5 см в форматах .efa и .efr: наши eps в обвязке ЛСРМ.

Точечная 5 см раньше не экспортировалась: в committed reference/lsrm для
неё есть только .efr (точки), а .efa (точки + полином зон) нет — то есть
не было шаблона с секцией Zones, без которой СпектраЛайн файл отвергает.

Здесь шаблоном служит блок «Точечная-5см» РАБОЧЕГО мастера прибора: из
него берётся вся обвязка (заголовок, набор узловых энергий с нуклидами и
площадями, зоны), а значения eps подменяются нашими из сетки p5cm
(eps_net, лог-лог интерполяция между узлами сетки).

ПОЛИНОМ ЗОН ПЕРЕСЧИТЫВАЕТСЯ (efa_zones.py). Прежние версии оставляли
секцию Zones мастерской, и это обесценивало весь файл: программа считает
активность по полиному, а не по точкам. Проверено на Th-228 (точечная
5 см, четыре прогона 30.07.2026, в том числе после применения калибровки
по эффективности): фактически применённая эффективность равнялась узлам
ЛСРМ (2614,5 кэВ — 2,915E-03 против узла 2,895E-03), а не нашим
(3,697E-03), и активность совпадала со штатной до единицы. Выбор
детектора при мастерском полиноме ни на что не влияет.

Путь к мастеру задаётся переменной окружения G1S_LSRM_MASTER_EFA (файл
рабочего каталога прибора, в репозиторий не входит: несёт заводской
номер). Без неё скрипт сообщает, что делать, и выходит.

.efr — ПОЛНЫЙ набор: 12 блоков-источников, 24 линии от 59,5 до 2614,5
(Am-241, Cd-109, Ba-133, Co-57, Eu-152, Cs-137, Mn-54, Zn-65, Y-88,
Co-60, Na-22, Th-228). Шаблон берётся из committed reference — он уже
обезличен, и никаких подстановок для репозитория не нужно; для программы
пишется копия с настоящим именем детектора (G1S_LSRM_DETECTOR).
В блоках сохраняются исходные площади и опорные активности «0,0,1»:
редактор эффективности не пересчитывает eps из активности, а берёт
готовые значения и строит по ним зоны.
"""
import bisect
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import efa_zones  # noqa: E402

RESULTS = os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results"))
GEOM = "Точечная-5см"


def curve_nodes():
    pts = []
    with open(os.path.join(RESULTS, "eff_p5cm.csv"), encoding="utf-8") as fh:
        for ln in fh:
            if ln.startswith(("#", "E_keV")):
                continue
            f = ln.split(",")
            pts.append((float(f[0]), float(f[1]), float(f[3])))
    pts.sort()
    return pts


def make_interp(pts, col):
    xs = [math.log(e) for e, _n, _g in pts]
    ys = [math.log(v[col]) for v in pts]

    def f(E):
        x = math.log(E)
        i = max(1, min(bisect.bisect_left(xs, x), len(xs) - 1))
        t = (x - xs[i - 1]) / (xs[i] - xs[i - 1])
        return math.exp(ys[i - 1] + t * (ys[i] - ys[i - 1]))
    return f


def block_of(path, geom):
    raw = open(path, "rb").read().decode("cp1251")
    out, take = [], False
    for ln in raw.replace("\r\n", "\n").split("\n"):
        if ln.startswith("["):
            if take:
                break
            take = (";" + geom + "]") in ln
            if take:
                out = [ln]
                continue
        if take:
            out.append(ln)
    return out


if __name__ == "__main__":
    master = os.environ.get("G1S_LSRM_MASTER_EFA")
    if not master or not os.path.exists(master):
        raise SystemExit(
            "Не задан G1S_LSRM_MASTER_EFA — путь к рабочему .efa прибора\n"
            "(многоблочный файл с геометриями и секциями Zones). В"
            " репозитории его нет:\nон содержит заводской номер.")
    block = block_of(master, GEOM)
    if not block:
        raise SystemExit("в мастере нет блока %s" % GEOM)
    eps = make_interp(curve_nodes(), 1)          # eps_net

    out, report = efa_zones.rewrite_block(block, eps)
    n = sum(1 for ln in block
            if "=" in ln and "," in ln and ln[:1].isdigit())
    body = "\r\n".join(out)
    print("полином зон пересчитан под наши точки:")
    for z, emin, emax, k, rms, dmax in report:
        print("   зона %d: %.0f-%.0f кэВ, узлов %d, СКО %.2f %%,"
              " макс. отклонение %.2f %%" % (z, emin, emax, k, rms, dmax))

    # В репозиторий — только обезличенная копия: заводской номер прибора
    # из шапки мастера заменяется на псевдоним, как в committed reference.
    anon = re.sub(r"№\s*\d{3,4}-\d{2}", "№SN-01", body)
    dst = os.path.join(RESULTS, "efa_export", "G4MC_проба_Точечная-5см.efa")
    open(dst, "wb").write(anon.encode("cp1251"))
    print("узлов заменено %d -> %s (обезличено)" % (n, dst))

    # Для программы — с настоящей шапкой, иначе она не свяжет файл с
    # детектором. Этот файл живёт только в рабочем каталоге прибора.
    live = os.environ.get("G1S_LSRM_EFF_DIR")
    if live and os.path.isdir(live):
        p = os.path.join(live, "G4MC_проба_Точечная-5см.efa")
        open(p, "wb").write(body.encode("cp1251"))
        print("копия для программы: %s" % p)

    lsrm = {59.541: 3.976771e-2, 88.034: 4.32493e-2, 121.782: 4.096356e-2,
            238.632: 3.128426e-2, 661.657: 1.370872e-2, 1332.492: 6.465514e-3,
            2614.511: 2.894993e-3}
    print("\n%10s %12s %12s %9s" % ("E, кэВ", "наша", "ЛСРМ", "наша/ЛСРМ"))
    for E in sorted(lsrm):
        o = eps(E)
        print("%10.1f %12.4e %12.4e %9.3f" % (E, o, lsrm[E], o / lsrm[E]))

    # --- .efr: полный набор источников -------------------------------------
    # Шаблон — committed reference (уже обезличен). Заменяются значения eps
    # во ВСЕХ блоках; площади, опорные активности «0,0,1» и выходы линий
    # остаются исходными: редактор эффективности берёт eps готовыми.
    tpl = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "reference", "lsrm",
        "efficiency", "Gamma-1S_NaI_63x63_USB_SN-01",
        "УДС-ГЦ-63х63-USB__SN-01_-_Точечная-5см.efr")
    if not os.path.exists(tpl):
        print("\nшаблон .efr не найден: %s" % tpl)
        raise SystemExit(0)
    src = open(tpl, "rb").read().decode("cp1251")
    lines, m = [], 0
    for ln in src.replace("\r\n", "\n").split("\n"):
        if "=" in ln and "," in ln and ln[:1].isdigit():
            key, rest = ln.split("=", 1)
            parts = rest.split(",")
            parts[0] = "%.6E" % eps(float(key))
            ln = key + "=" + ",".join(parts)
            m += 1
        lines.append(ln)
    body_r = "\r\n".join(lines)
    dst_r = os.path.join(RESULTS, "efa_export", "G4MC_Точечная-5см.efr")
    open(dst_r, "wb").write(body_r.encode("cp1251"))
    print("\n.efr: заменено %d линий в %d блоках -> %s"
          % (m, body_r.count("["), dst_r))

    # Копия для программы — с настоящим именем детектора. Имя НАБОРА
    # (третье поле заголовка блока) программа берёт в метаданные результата;
    # у штатных наборов там стоит фамилия оператора, у наших — «GEANT4»
    # (указание оператора 29.07.2026), чтобы происхождение кривой было
    # видно в отчёте и наборы не путались со штатными.
    det = os.environ.get("G1S_LSRM_DETECTOR")
    if live and os.path.isdir(live):
        b = body_r
        if det:
            b = b.replace("УДС-ГЦ-63х63-USB №SN-01", det)

        def rename(mt):
            parts = mt.group(0)[1:-1].split(";")
            if len(parts) >= 3:
                nuc = re.search(r"[A-Z][a-z]?-\d+", parts[2])
                parts[2] = "GEANT4" + (" (%s)" % nuc.group(0) if nuc else "")
            return "[" + ";".join(parts) + "]"
        # без $ в шаблоне: строки заканчиваются CRLF, и «конец строки»
        # оказывается ПОСЛЕ \r — привязка по $ молча не срабатывает
        b = re.sub(r"^\[[^\r\n]*\]", rename, b, flags=re.M)
        p = os.path.join(live, "G4MC_Точечная-5см.efr")
        open(p, "wb").write(b.encode("cp1251"))
        print("копия для программы: %s" % p)
