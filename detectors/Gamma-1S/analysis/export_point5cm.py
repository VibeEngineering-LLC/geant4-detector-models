"""Кривая точечной 5 см в формате .efa: наши eps в родной обвязке ЛСРМ.

Точечная 5 см раньше не экспортировалась: в committed reference/lsrm для
неё есть только .efr (точки), а .efa (точки + полином зон) нет — то есть
не было шаблона с секцией Zones, без которой СпектраЛайн файл отвергает.

Здесь шаблоном служит блок «Точечная-5см» РАБОЧЕГО мастера прибора: из
него берётся вся обвязка (заголовок, набор узловых энергий с нуклидами и
площадями, зоны), а значения eps подменяются нашими из сетки p5cm
(eps_net, лог-лог интерполяция между узлами сетки).

Как и с маринелли: полином зон в шаблоне описывает ЧУЖУЮ кривую, поэтому
файл годится (а) как вход редактора эффективности, который перестроит
зоны по нашим точкам, (б) как проба «применяются ли точки вообще» — на
маринелли ответ был отрицательный, программа считает по полиному.

Путь к мастеру задаётся переменной окружения G1S_LSRM_MASTER_EFA (файл
рабочего каталога прибора, в репозиторий не входит: несёт заводской
номер). Без неё скрипт сообщает, что делать, и выходит.
"""
import bisect
import math
import os
import sys

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

    out, n = [], 0
    for ln in block:
        if "=" in ln and "," in ln and ln[:1].isdigit():
            key, rest = ln.split("=", 1)
            parts = rest.split(",")
            parts[0] = "%.6E" % eps(float(key))
            ln = key + "=" + ",".join(parts)
            n += 1
        out.append(ln)
    body = "\r\n".join(out)

    # В репозиторий — только обезличенная копия: заводской номер прибора
    # из шапки мастера заменяется на псевдоним, как в committed reference.
    import re
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
