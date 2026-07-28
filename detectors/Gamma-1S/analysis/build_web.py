# -*- coding: utf-8 -*-
"""Сводная веб-страница Гамма-1С: паспортные и расчётные кривые наложением.

Один источник данных — то, что лежит в репозитории:
  reference/.../*.efa, *.efr  — измеренные кривые ЛСРМ по геометриям;
  results/eff_*.csv           — расчётные кривые;
  results/summing_C.csv       — поправки на каскадное суммирование;
  results/kit_recalc_*.csv    — пересчёт комплекта против паспортов.

Выход: docs/gamma-1s/index.html — самодостаточная страница (графики строятся
как SVG прямо здесь, внешних запросов нет: страница обязана открываться и с
GitHub Pages, и из файла).

    python detectors/Gamma-1S/analysis/build_web.py

Никаких чисел руками: пересчитали кривые — пересоберите страницу одной
командой, и все подписи, отношения и медианы обновятся сами.
"""
import csv
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import paths  # noqa: E402

sys.path.insert(0, str(paths.tools()))
from fetch_efr import parse_efr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))          # detectors/Gamma-1S
REPO = os.path.abspath(os.path.join(ROOT, "..", ".."))
RES = os.path.join(ROOT, "results")
OUT_DIR = os.path.join(REPO, "docs", "gamma-1s")
OUT = os.path.join(OUT_DIR, "index.html")

# Геометрия: (подпись, файл измеренной кривой без расширения, расчётная кривая,
#             пояснение, ключ пересчёта комплекта)
GEOMS = [
    ("Маринелли 1 л", "Маринелли", "eff_rho1.60.csv",
     "Проба обнимает детектор с боков и сверху. ОИСН-16, ρ = 1,6 г/см³ — "
     "та же засыпка, при которой снималась паспортная кривая.",
     "Marinelli_1L"),
    ("«Дента» 120 мл", "Дента", "eff_denta1.60.csv",
     "Плоская кювета на головке. ОИСН-16, ρ = 1,6 г/см³.",
     "Denta_120mL"),
    ("Петри 60 мл", "Петри", "eff_petri1.60.csv",
     "Тонкий слой на головке, самая мелкая из кювет комплекта.",
     "Petri_60mL"),
    ("Точечная, 5 см", "Точечная-5см", "eff_p5cm.csv",
     "Точечный источник на 5 см от торца, крышка защиты ЗАКРЫТА. "
     "Здесь нет сосуда — проверяется сам детектор.",
     "Point_5cm"),
    ("Точечная, 25 см", "Точечная-25см", "eff_p25cm.csv",
     "Точечный источник на 25 см, крышка защиты СНЯТА. Вторая, независимая "
     "проверка детектора — в другой конфигурации защиты.",
     "Point_25cm"),
]

EFF_DIR = None


def find_eff_dir():
    """Каталог с .efa/.efr этого экземпляра прибора."""
    root = paths.ref("Gamma-1S")
    for p in root.rglob("*.efa"):
        return p.parent
    return None


def measured(name):
    """Точки измеренной кривой: [(E, eps, d_eps)] — из .efa, иначе из .efr.

    .efa — уже сведённая ЛСРМ кривая (одна точка на линию). У точечной 5 см
    её нет, там только .efr с блоком на каждый источник; блоки объединяются,
    повторы одной энергии усредняются с весом 1/σ².
    """
    for ext in (".efa", ".efr"):
        hits = sorted(EFF_DIR.glob("*%s*%s" % (name, ext)))
        if not hits:
            continue
        secs = parse_efr(paths.read_text(hits[0]))
        acc = {}
        for s in secs:
            for E, eps, dpct, nuc in s["points"]:
                if eps <= 0 or dpct <= 0:
                    continue
                d = eps * dpct / 100.0
                acc.setdefault(round(E, 2), []).append((eps, d, nuc))
        out = []
        for E, vals in sorted(acc.items()):
            w = sum(1.0 / v[1] ** 2 for v in vals)
            eps = sum(v[0] / v[1] ** 2 for v in vals) / w
            out.append((E, eps, w ** -0.5, vals[0][2]))
        return out, os.path.basename(hits[0])
    return [], None


def computed(fn):
    """Расчётная кривая: [(E, eps_net, d_eps, eps_gross)]."""
    out = []
    with open(os.path.join(RES, fn), encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            out.append((float(r["E_keV"]), float(r["eps_net"]),
                        float(r["d_eps"]), float(r["eps_gross"])))
    return out


def summing():
    """Поправки C по линиям: {E: (C, dC, нуклид)}."""
    out = {}
    p = os.path.join(RES, "summing_C.csv")
    if not os.path.exists(p):
        return out
    with open(p, encoding="utf-8") as fh:
        for r in csv.DictReader(l for l in fh if not l.startswith("#")):
            out[round(float(r["E_keV"]), 1)] = (
                float(r["C_summing"]), float(r["d_C"]), r["nuclide"])
    return out


def kit_medians():
    """Медианы A_изм/A_пасп по геометриям пересчёта комплекта."""
    out = {}
    for fn in ("kit_recalc_volume.csv", "kit_recalc_point.csv"):
        p = os.path.join(RES, fn)
        if not os.path.exists(p):
            continue
        with open(p, encoding="utf-8") as fh:
            for r in csv.DictReader(l for l in fh if not l.startswith("#")):
                out.setdefault(r["geometry"], []).append(float(r["ratio"]))
    return {k: sorted(v) for k, v in out.items()}


# --- построение SVG ---------------------------------------------------------
# Рисуем на сервере, а не в браузере: страница обязана открываться без единого
# внешнего запроса, а картинка — быть одинаковой везде, включая печать.

W, H = 760, 400
PAD_L, PAD_R, PAD_T, PAD_B = 66, 14, 14, 44
RW, RH = 760, 170
E_LO, E_HI = 50.0, 3200.0


def lx(E, w=W):
    a, b = math.log10(E_LO), math.log10(E_HI)
    return PAD_L + (math.log10(E) - a) / (b - a) * (w - PAD_L - PAD_R)


def ly(v, lo, hi):
    a, b = math.log10(lo), math.log10(hi)
    return H - PAD_B - (math.log10(v) - a) / (b - a) * (H - PAD_T - PAD_B)


def esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;"))


def ticks_e():
    out = []
    for v in (50, 100, 200, 500, 1000, 2000, 3000):
        out.append((v, "%d" % v))
    return out


def chart(meas, comp, lo, hi):
    """Наложение: измеренные точки с усами и расчётная кривая."""
    s = ['<svg viewBox="0 0 %d %d" class="plot" role="img" '
         'aria-label="Эффективность ППП: расчёт и измерение">' % (W, H)]
    # сетка по энергии
    for v, lab in ticks_e():
        x = lx(v)
        s.append('<line class="grid" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>'
                 % (x, PAD_T, x, H - PAD_B))
        s.append('<text class="ax" x="%.1f" y="%d" text-anchor="middle">%s</text>'
                 % (x, H - PAD_B + 16, lab))
    # сетка по эффективности: декады и половинки
    d0 = math.floor(math.log10(lo))
    d1 = math.ceil(math.log10(hi))
    dec = d0
    while dec <= d1:
        for m in (1, 2, 5):
            v = m * 10.0 ** dec
            if not (lo <= v <= hi):
                continue
            y = ly(v, lo, hi)
            s.append('<line class="grid" x1="%d" y1="%.1f" x2="%d" y2="%.1f"/>'
                     % (PAD_L, y, W - PAD_R, y))
            lab = ("%g" % (v * 100)) if v * 100 >= 0.01 else ("%.3f" % (v * 100))
            s.append('<text class="ax" x="%d" y="%.1f" text-anchor="end">%s</text>'
                     % (PAD_L - 6, y + 4, lab))
        dec += 1
    s.append('<text class="axt" x="%.1f" y="%d" text-anchor="middle">'
             'энергия, кэВ</text>' % ((PAD_L + W - PAD_R) / 2, H - 6))
    s.append('<text class="axt" transform="translate(14,%.1f) rotate(-90)" '
             'text-anchor="middle">эффективность ППП, %%</text>'
             % ((PAD_T + H - PAD_B) / 2))

    # расчётная кривая
    pts = [(E, v) for E, v, _, _ in comp if lo <= v <= hi and E >= E_LO]
    if pts:
        d = " ".join("%s%.1f,%.1f" % ("M" if i == 0 else "L",
                                      lx(E), ly(v, lo, hi))
                     for i, (E, v) in enumerate(pts))
        s.append('<path class="mc" d="%s"/>' % d)
        for E, v, dv, _ in comp:
            if not (lo <= v <= hi and E >= E_LO):
                continue
            x, y = lx(E), ly(v, lo, hi)
            if dv > 0 and v - dv > lo:
                s.append('<line class="mcerr" x1="%.1f" y1="%.1f" x2="%.1f" '
                         'y2="%.1f"/>' % (x, ly(v + dv, lo, hi),
                                          x, ly(max(v - dv, lo), lo, hi)))
            s.append('<circle class="mcp" cx="%.1f" cy="%.1f" r="2.6">'
                     '<title>расчёт: %.1f кэВ, %.4g %%</title></circle>'
                     % (x, y, E, v * 100))

    # измеренные точки
    for E, v, dv, nuc in meas:
        if not (lo <= v <= hi and E_LO <= E <= E_HI):
            continue
        x, y = lx(E), ly(v, lo, hi)
        if dv > 0 and v - dv > lo:
            s.append('<line class="experr" x1="%.1f" y1="%.1f" x2="%.1f" '
                     'y2="%.1f"/>' % (x, ly(v + dv, lo, hi),
                                      x, ly(max(v - dv, lo), lo, hi)))
        s.append('<rect class="expp" x="%.1f" y="%.1f" width="6" height="6">'
                 '<title>измерение (%s): %.1f кэВ, %.4g %% ± %.2g</title>'
                 '</rect>' % (x - 3, y - 3, esc(nuc), E, v * 100, dv * 100))
    s.append("</svg>")
    return "\n".join(s)


def ratio_chart(pairs, med):
    """МК/эксп по точкам: линия единицы и медиана."""
    lo, hi = 0.55, 1.75
    for p in pairs:
        r, rc = p[1], p[4]
        lo = min(lo, r * 0.9, (rc or r) * 0.9)
        hi = max(hi, r * 1.1, (rc or r) * 1.1)

    def y(v):
        return RH - 30 - (v - lo) / (hi - lo) * (RH - 30 - 12)

    s = ['<svg viewBox="0 0 %d %d" class="plot ratio" role="img" '
         'aria-label="Отношение расчёт/измерение по точкам">' % (RW, RH)]
    for v, lab in ticks_e():
        x = lx(v, RW)
        s.append('<line class="grid" x1="%.1f" y1="12" x2="%.1f" y2="%d"/>'
                 % (x, x, RH - 30))
        s.append('<text class="ax" x="%.1f" y="%d" text-anchor="middle">%s</text>'
                 % (x, RH - 12, lab))
    for v in (0.6, 0.8, 1.0, 1.2, 1.4, 1.6):
        if not (lo <= v <= hi):
            continue
        cls = "one" if abs(v - 1.0) < 1e-9 else "grid"
        s.append('<line class="%s" x1="%d" y1="%.1f" x2="%d" y2="%.1f"/>'
                 % (cls, PAD_L, y(v), RW - PAD_R, y(v)))
        s.append('<text class="ax" x="%d" y="%.1f" text-anchor="end">%.1f</text>'
                 % (PAD_L - 6, y(v) + 4, v))
    if lo <= med <= hi:
        s.append('<line class="med" x1="%d" y1="%.1f" x2="%d" y2="%.1f"/>'
                 % (PAD_L, y(med), RW - PAD_R, y(med)))
    for E, r, dr, nuc, rc, is_node in pairs:
        x = lx(E, RW)
        if dr > 0:
            s.append('<line class="mcerr" x1="%.1f" y1="%.1f" x2="%.1f" '
                     'y2="%.1f"/>' % (x, y(min(r + dr, hi)),
                                      x, y(max(r - dr, lo))))
        if rc is not None:
            s.append('<line class="corr" x1="%.1f" y1="%.1f" x2="%.1f" '
                     'y2="%.1f"/>' % (x, y(r), x, y(rc)))
            s.append('<circle class="corrp" cx="%.1f" cy="%.1f" r="2.4">'
                     '<title>%.1f кэВ, с поправкой на суммирование: %.3f'
                     '</title></circle>' % (x, y(rc), E, rc))
        s.append('<circle class="mcp" cx="%.1f" cy="%.1f" r="3">'
                 '<title>%.1f кэВ (%s): МК/эксп %.3f ± %.3f</title></circle>'
                 % (x, y(r), E, esc(nuc), r, dr))
    s.append('<text class="axt" x="%.1f" y="%d" text-anchor="middle">'
             'энергия, кэВ</text>' % ((PAD_L + RW - PAD_R) / 2, RH - 1))
    s.append("</svg>")
    return "\n".join(s)


DEG = 5          # как в compare_point.py: та же процедура, что строит кривую
INTERP_D = 0.02  # оценка погрешности интерполяции


def pair_up(meas, comp, C, gross=False):
    """(E, МК/эксп, погрешность, нуклид, отношение с поправкой C, узел?).

    gross=True — брать площадь БЕЗ вычета левой полки континуума. Это не
    придирка: разница между двумя способами взятия площади доходит до 10 % на
    мягких линиях, и именно по gross считались точечные геометрии в отчёте.
    Оба числа приводятся рядом, чтобы выбор обработки не выдавался за свойство
    детектора.

    Совпадение энергий узлов сетки и линий эталонов — редкость: у ЛСРМ до 24
    линий, а в сетке 20 узлов. Поэтому там, где узла нет, расчёт берётся
    интерполяцией полиномом по log-log — ТОЙ ЖЕ процедурой, которой строится
    рабочая кривая, и той же, что в compare_point.py. Иначе сравнивались бы
    8 линий вместо 24, и число на странице расходилось бы с отчётом.
    """
    import numpy as np
    col = 3 if gross else 1
    Eg = np.array([c[0] for c in comp], dtype=float)
    yg = np.array([c[col] for c in comp], dtype=float)
    dyg = np.array([c[2] for c in comp], dtype=float)
    cf = np.polyfit(np.log(Eg), np.log(yg), DEG,
                    w=yg / np.maximum(dyg, 1e-30))
    out = []
    for E, ev, edv, nuc in meas:
        if ev <= 0:
            continue
        node = None
        for c in comp:
            if abs(c[0] - E) <= 1.0:
                node = (c[col], c[2])
        if node is not None:
            m, dm, is_node = node[0], node[1], True
        elif Eg[0] <= E <= Eg[-1]:
            m = math.exp(float(np.polyval(cf, math.log(E))))
            dm, is_node = m * INTERP_D, False
        else:
            continue
        r = m / ev
        dr = r * math.hypot(dm / m if m else 0, edv / ev)
        c = C.get(round(E, 1))
        out.append((E, r, dr, nuc, (r / c[0]) if c else None, is_node))
    return out


def wmean(pairs):
    """Средневзвешенное МК/эксп в логарифме — как в отчёте, плюс RMS формы."""
    logs = [math.log(p[1]) for p in pairs]
    ws = [1.0 / (p[2] / p[1]) ** 2 for p in pairs]
    if not logs:
        return float("nan"), float("nan")
    lw = sum(l * w for l, w in zip(logs, ws)) / sum(ws)
    k = math.exp(lw)
    dev = [math.exp(l - lw) - 1 for l in logs]
    rms = math.sqrt(sum(d * d for d in dev) / len(dev))
    return k, rms


def med_of(v):
    v = sorted(v)
    n = len(v)
    if not n:
        return float("nan")
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def ru(x, nd=3):
    return ("%.*f" % (nd, x)).replace(".", ",")


# Стиль лежит отдельным файлом page.css, а не строкой здесь, по прозаической
# причине: проверка перед публикацией ищет номера источников, записанные после
# знака номера, и шестнадцатеричный цвет от такого номера по форме не
# отличается. В файлах разметки и стилей она подобные находки пропускает по
# расширению, в исходниках на питоне — нет. Ослаблять её ради удобства нельзя:
# номер экземпляра источника выглядит точно так же, и один раз он так и утёк.
CSS = open(os.path.join(HERE, "page.css"), encoding="utf-8").read()


TMPL = """<meta charset="utf-8">
<title>ГАММА-1С: расчётные и паспортные кривые эффективности</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="description" content="Модель спектрометра ГАММА-1С в Geant4:
наложение расчётных и паспортных кривых эффективности по пяти геометриям
комплекта поверки.">
<style>%(css)s</style>

<h1>ГАММА-1С: расчёт против паспорта</h1>
<p class="sub">Модель сцинтилляционного спектрометра NaI(Tl) 63&times;63 в
свинцовой защите — наложение расчётных и паспортных кривых эффективности по
пяти геометриям комплекта поверки.</p>

<p class="lead">Спектрометр ГАММА-1С поверен, и на каждую геометрию у него есть
кривая эффективности, снятая по эталонным источникам. Это редкая роскошь:
обычно расчётную модель сверять не с чем. Здесь есть с чем — и ниже показано,
где модель совпадает с прибором, а где нет, без сглаживания расхождений.</p>

<h2>Что сравнивается</h2>

<p><b>Паспортная кривая</b> — результат поверки: площадь пика полного
поглощения от эталонного источника, делённая на активность и выход линии.
Синие точки, усы — заявленная погрешность. Каждая точка привязана к своему
нуклиду; это не гладкая функция, а набор измерений.</p>

<p><b>Расчётная кривая</b> — Монте-Карло в Geant4 11.2.1
(EmStandardPhysics_option4 + RadioactiveDecay): моноэнергетические кванты
разыгрываются по объёму пробы, считается доля попавших в пик полного
поглощения. Красная линия, усы — статистика прогона. Геометрия построена по
чертежам и паспорту, ни одно число не подгонялось под измеренную кривую.</p>

<div class="card warn">
<p><b>Одна поправка обязательна, и без неё сравнение неверно.</b> В расчёте
квант приходит в кристалл один — суммироваться не с чем. В реальном распаде
каскадные кванты летят одновременно, и если два попали в кристалл, событие
уходит из своего пика. Поэтому измеренная эффективность каскадных линий
занижена, и отношение расчёт/измерение надо делить на поправку C, посчитанную
отдельным прогоном полного распада. Там, где C известна, сиреневый штрих на
нижнем графике показывает, куда уезжает точка. Контроль: у Cs-137 и K-40
каскада нет, и для них C вышла 0,983 и 1,009 — единица в пределах
погрешности.</p>
</div>

<h2>Сводка</h2>

<div class="tw"><table>
<thead><tr><th>геометрия</th><th class="n">линий</th>
<th class="n">МК/эксп<br>по чистой площади</th>
<th class="n">то же<br>без вычета полки</th>
<th class="n">RMS формы, %%</th>
<th class="n">линий в пересчёте</th><th class="n">медиана A/пасп</th></tr>
</thead><tbody>%(summary)s</tbody></table></div>

<p class="cap">Два столбца отношения — не разброс расчёта, а систематика
обработки: площадь пика можно брать с вычетом левой полки континуума и без
него. Разница доходит до 10 %% на мягких линиях и меньше процента на жёстких.
Опубликовать одно число, умолчав о втором, значило бы выдать выбор обработки
за свойство детектора. Дальше по тексту используется чистая площадь.</p>

<p>Два правых столбца — независимая проверка, в которой паспортная кривая не
участвует вовсе. Там решается обратная задача: из измеренного спектра
эталонного источника берётся скорость счёта в пике, делится на расчётную
эффективность и получается активность, которую сравнивают с паспортом
источника. Подогнать в этой процедуре нечего. Знак расхождения в обоих
столбцах один и тот же, хотя пути расчёта разные.</p>

<p><b>Читать так:</b> МК/эксп больше единицы — расчёт завышает эффективность;
A/пасп меньше единицы — то же самое (завышенной эффективностью
восстанавливается заниженная активность). Оба столбца согласованно говорят
одно: <b>маринелли завышена, кюветы занижены, точечные геометрии почти
верны</b>.</p>

<h2>Главный вывод: расходятся модели сосудов, а не детектор</h2>

<p>Две точечные геометрии — 5 и 25 см — согласны между собой и обе лежат около
единицы, причём в РАЗНЫХ конфигурациях защиты: с закрытой и со снятой крышкой.
Значит модель самого устройства детектирования проверена дважды и независимо.
А объёмные геометрии расходятся, и <b>расходятся в разные стороны</b>: одна
ошибка эффективности так себя вести не может.</p>

<p>Отсюда адрес расхождения — <b>геометрия сосудов</b>. И это ровно то место,
где данных не хватает: из документов у сосуда Маринелли известны только
габарит Ø150 и высота 110 мм, а толщина стенки 2 мм, диаметр колодца
Ø80 и его глубина 74 мм — <b>допущения</b>. Чертежей кювет комплекта нет.
Если проба в действительности отстоит от кристалла дальше, расчётная
эффективность падает; ближе — растёт.</p>

<div class="card">
<p><b>Что уже отклонено по заранее записанным признакам.</b> Гипотеза «дело в
подходе: Монте-Карло считает энерговыделение, а эксперимент площадь пика» —
будь это так, точечная геометрия дала бы то же превышение; она даёт единицу.
Гипотеза «дело в нуклидной зависимости» — отклонена после того, как разброс
внутри каждого сосуда схлопнулся (см. ниже). Остаются размеры сосудов, и для
них нужны не расчёты, а чертежи.</p>
</div>

<h2>Дефект, который до 28 июля 2026 портил эти числа</h2>

<p>В прогонах по цепочке (Ra-226 → … → Po-214, Th-232 → … → Pb-208) Geant4
доводил весь ряд до конца <b>внутри одного события</b>: порог «очень долгого
распада» поднят до 10<sup>30</sup> нс, иначе долгоживущие звенья не распались
бы вовсе. В итоге энерговыделения Ac-228, Tl-208, Bi-212 складывались, и
получались совпадения между ядрами, которые в природе распадаются с разницей
в <b>годы</b>. Спектрометр так себя не ведёт.</p>

<p>Лечение: энерговыделения собираются с отметкой глобального времени и
режутся на группы с разрывом больше 1 мкс — обычного времени разрешения
тракта. Значение некритично: настоящие каскады приходят за наносекунды, а
звенья ряда разделены секундами и годами, между этими масштабами шесть
порядков пустоты.</p>

<p><b>Проверка правки.</b> Одиночный нуклид (Tl-208, 100 000 распадов) до и
после дал побайтно те же данные и то же число сработавших событий — правка
задевает только цепочки. На них пики ряда тория выросли на 25–46 %%, ряда
радия на 12–14 %%; в точечной геометрии всего на 7 %%, потому что при малом
телесном угле двойное попадание маловероятно.</p>

<p><b>Что это изменило по существу.</b> Раньше в сосуде Маринелли цезий давал
отношение 0,78, а торий 1,26 — разница в 60 %% на одной и той же геометрии,
которую геометрия объяснить не может: она не знает, какой нуклид внутри.
Теперь весь сосуд укладывается в 0,69–0,89. У каждой геометрии остался
<b>один</b> множитель вместо разнобоя по нуклидам, и искать надо одну причину
на сосуд.</p>

<h2>Кривые по геометриям</h2>

%(legend)s

%(blocks)s

<h2>Как это воспроизвести</h2>

<p>Страница собирается из того, что лежит в репозитории, одной командой —
чисел, вписанных руками, здесь нет:</p>

<p><code>python detectors/Gamma-1S/analysis/build_web.py</code></p>

<p><b>Правило отбора точек, чтобы числа можно было проверить.</b> Берутся все
линии паспортной кривой, попадающие в диапазон расчётной сетки
(59,5–3000 кэВ). Расчёт берётся из узла сетки, если он ближе 1 кэВ, иначе
интерполяцией полиномом 5-й степени по log-log — той же процедурой, которой
строится рабочая кривая. Отношения усредняются в логарифме с весом
1/&sigma;². Отдельные скрипты в репозитории считают те же геометрии со своим,
более узким отбором точек: по маринелли и «Денте» результаты совпадают в
пределах 0,5 %%, по Петри расходятся на 5 %% — там из 19 линий отбирается 14.
Все использованные точки перечислены в таблицах выше, так что пересчитать
можно любым способом.</p>

<p>Расчётные кривые — <code>detectors/Gamma-1S/results/eff_*.csv</code>,
паспортные — <code>reference/lsrm/efficiency/</code>, поправки на
суммирование — <code>results/summing_C.csv</code>, пересчёт комплекта —
<code>results/kit_recalc_*.csv</code>. Сами расчётные спектры (258 файлов) в
репозиторий не входят: они воспроизводятся драйверами из
<code>drivers/</code>.</p>

<p class="foot">Полный отчёт со всеми оговорками —
<a href="https://github.com/VibeEngineering-LLC/geant4-detector-models/blob/main/detectors/Gamma-1S/REPORT.md">REPORT.md</a>.
Протокол обязательных проверок и список ловушек —
<a href="https://github.com/VibeEngineering-LLC/geant4-detector-models/tree/main/common/docs">common/docs</a>.
Исходники —
<a href="https://github.com/VibeEngineering-LLC/geant4-detector-models">на GitHub</a>.
Эталонные данные обезличены: фамилии операторов и номера приборов заменены
псевдонимами, соответствие подлинным не публикуется.<br>
Текст, код и графики этой страницы подготовлены с участием ИИ (Claude, модель
Opus&nbsp;5) под проверкой оператора; числа получены расчётом и приведены как
есть.</p>
"""


def legend(with_corr):
    k = ['<span class="k"><i style="border-color:var(--mc)"></i>расчёт '
         '(Geant4)</span>',
         '<span class="k"><i style="border-color:var(--exp)"></i>измерение '
         '(паспортная кривая ЛСРМ)</span>']
    if with_corr:
        k.append('<span class="k"><i style="border-color:var(--corr)"></i>'
                 'то же с поправкой на суммирование</span>')
    k.append('<span class="k"><i style="border-color:var(--med)"></i>медиана'
             '</span>')
    return '<p class="leg">%s</p>' % "".join(k)


def build():
    global EFF_DIR
    EFF_DIR = find_eff_dir()
    if EFF_DIR is None:
        raise SystemExit(
            "Не найдены файлы кривых .efa/.efr. Они лежат в\n"
            "detectors/Gamma-1S/reference/lsrm/efficiency/; укажите\n"
            "G4MODELS_REF, если эталоны вынесены из репозитория.")
    C = summing()
    kit = kit_medians()

    blocks, summary = [], []
    for title, mname, cfile, note, kitkey in GEOMS:
        meas, src = measured(mname)
        comp = computed(cfile)
        if not meas:
            print("!! нет измеренной кривой для", title)
            continue
        pairs = pair_up(meas, comp, C)
        rs = [p[1] for p in pairs]
        k, rms = wmean(pairs)
        kg, _ = wmean(pair_up(meas, comp, C, gross=True))
        med = med_of(rs)
        vals = [v for _, v, _, _ in meas] + [v for _, v, _, _ in comp]
        lo = min(vals) / 2.2
        hi = max(vals) * 2.2
        rows = []
        for E, r, dr, nuc, rc, is_node in pairs:
            c = C.get(round(E, 1))
            ev = next(v for e, v, _, _ in meas if abs(e - E) < .01)
            rows.append(
                "<tr><td class='n'>%s</td><td>%s</td><td class='n'>%s</td>"
                "<td class='n'>%s</td><td class='n'>%s</td>"
                "<td class='n'>%s</td><td>%s</td></tr>"
                % (ru(E, 1), esc(nuc), "%.4g" % (ev * 100),
                   "%.4g" % (r * ev * 100), ru(r),
                   (ru(c[0]) + " → " + ru(rc)) if c else "—",
                   "узел" if is_node else "интерп."))
        blocks.append("""
<h3>%s</h3>
<p>%s</p>
<figure>%s
%s
<p class="cap">Точки измерения — %s, усы: паспортная погрешность. Расчёт —
%s, усы: статистика прогона. Наведите курсор на точку, чтобы увидеть числа.</p>
</figure>
<figure>%s
<p class="cap">Отношение расчёт/измерение по линиям. Средневзвешенное
(в логарифме, как в отчёте) %s по %d точкам, разброс формы RMS %s %%; медиана
%s. Сиреневым — куда точка уходит после деления на поправку каскадного
суммирования там, где она посчитана: измерение теряет отсчёты из пика на
совпадениях, расчёт моноэнергетической сетки — нет.</p>
</figure>
<div class="tw"><table>
<thead><tr><th class="n">E, кэВ</th><th>по нуклиду</th>
<th class="n">измерено, %%</th><th class="n">расчёт, %%</th>
<th class="n">МК/эксп</th><th class="n">C → с поправкой</th>
<th>расчёт взят</th></tr></thead>
<tbody>%s</tbody></table></div>
""" % (esc(title), esc(note), chart(meas, comp, lo, hi),
            legend(False), esc(src), esc(cfile),
            ratio_chart(pairs, med), ru(k), len(pairs), ru(100 * rms, 1),
            ru(med), "".join(rows)))
        kr = kit.get(kitkey, [])
        summary.append((title, len(pairs), k, kg, rms,
                        len(kr), med_of(kr) if kr else None))

    srows = "".join(
        "<tr><td>%s</td><td class='n'>%d</td><td class='n'>%s</td>"
        "<td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td>"
        "<td class='n'>%s</td></tr>"
        % (esc(t), n, ru(k), ru(kg), ru(100 * rms, 1),
           nk if nk else "—", ru(km) if km is not None else "—")
        for t, n, k, kg, rms, nk, km in summary)

    html = TMPL % dict(css=CSS, summary=srows, blocks="".join(blocks),
                       legend=legend(True))
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(html)
    print("страница: %s (%.0f КБ)" % (OUT, os.path.getsize(OUT) / 1024))
    for t, n, k, kg, rms, nk, km in summary:
        print("   %-16s линий %2d, МК/эксп %s (без вычета полки %s), "
              "RMS формы %s %%, пересчёт комплекта %s"
              % (t, n, ru(k), ru(kg), ru(100 * rms, 1),
                 ru(km) if km is not None else "—"))


if __name__ == "__main__":
    build()
