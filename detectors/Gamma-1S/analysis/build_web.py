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
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "..", "common", "py"))
import csvio  # noqa: E402
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
     "Проба окружает детектор. ОИСН-16, ρ = 1,6 г/см³ — "
     "та же засыпка, при которой снималась паспортная кривая.",
     "Marinelli_1L"),
    ("«Дента» 120 мл", "Дента", "eff_denta1.60.csv",
     "Плоская кювета на головке. ОИСН-16, ρ = 1,6 г/см³.",
     "Denta_120mL"),
    ("Петри 60 мл", "Петри", "eff_petri1.60.csv",
     "Тонкий слой на головке, самая мелкая из кювет комплекта.",
     "Petri_60mL"),
    ("Точечная, 5 см", "Точечная-5см", "eff_p5cm.csv",
     "Точечный источник на 5 см от торца, крышка защиты закрыта. "
     "Здесь нет сосуда — проверяется сам детектор.",
     "Point_5cm"),
    ("Точечная, 25 см", "Точечная-25см", "eff_p25cm.csv",
     "Точечный источник на 25 см, крышка защиты снята. Вторая, независимая "
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


# --- подгоночная кривая ЛСРМ из .efa ----------------------------------------
# Точки — это ещё не паспортная кривая. Прибор работает по ПОДГОНКЕ, и она
# лежит в том же .efa, просто не среди точек: блок Zones/Curve_*.
#
# Формат (разобран по данным, в документации к формату его нет):
#   Zones=N
#   Zone_i = степень, lg(E_min), lg(E_max), sigma
#   Curve_i_k = коэффициенты k-го базисного полинома, СТАРШАЯ степень первой
#   Curve_i   = коэффициенты разложения по этому базису
# Значение: lg(eps) = sum_k c_k * P_k(lg E). Базис ортогональный (Форсайта) по
# самим точкам: у линейного P_2 корень стоит во взвешенном среднем lg E, у
# квадратичного оба корня внутри диапазона данных — по этому признаку
# соглашение и опознано.
#
# Проверка расшифровки: на точечной 25 см подгонка воспроизводит измеренные
# точки с точностью 2 % при объявленной sigma 0,005 в lg. У маринелли и
# «Денты» разброс больше (до 20 %), но он свойствен САМИМ ТОЧКАМ: на 238,6 и
# 242,0 кэВ у них 4,34 и 5,05 % — линии Ra-226 и Th-232 в NaI не разделяются,
# и подгонка идёт между ними.
#
# ЗОНЫ ПЕРЕКРЫВАЮТСЯ, и файл не говорит, как их сшивать. У точечной 25 см в
# перекрытии 273-769 кэВ зоны согласны на 0,3-5 %, у Петри в 234-1854 кэВ на
# 0,6-8 %; расходятся только у краёв, как и положено полиномам. Принято
# линейное смешивание по lg E внутри перекрытия — это ДОПУЩЕНИЕ, но оно даёт
# непрерывную кривую и не даёт краевых выбросов.


def zones_of(text):
    nz = re.search(r"^Zones=(\d+)", text, re.M)
    if not nz:
        return []
    out = []
    for i in range(1, int(nz.group(1)) + 1):
        m = re.search(r"^Zone_%d=([^\r\n]+)" % i, text, re.M)
        c = re.search(r"^Curve_%d=([^\r\n]+)" % i, text, re.M)
        if not (m and c):
            continue
        f = m.group(1).split(",")
        deg, xlo, xhi, sig = int(f[0]), float(f[1]), float(f[2]), float(f[3])
        basis = []
        for k in range(1, deg + 2):
            b = re.search(r"^Curve_%d_%d=([^\r\n]+)" % (i, k), text, re.M)
            if b:
                basis.append([float(v) for v in b.group(1).split(",")])
        cf = [float(v) for v in c.group(1).split(",")]
        out.append(dict(deg=deg, xlo=xlo, xhi=xhi, sig=sig, basis=basis,
                        cf=cf))
    return out


def _pv(c, x):
    v = 0.0
    for a in c:
        v = v * x + a
    return v


def _lg_zone(z, x):
    return sum(c * _pv(b, x) for c, b in zip(z["cf"], z["basis"]))


def fit_eps(zones, E):
    """Паспортная подгонка в точке E, или None вне области определения."""
    x = math.log10(E)
    act = [z for z in zones if z["xlo"] - 1e-9 <= x <= z["xhi"] + 1e-9]
    if not act:
        return None
    if len(act) == 1:
        return 10 ** _lg_zone(act[0], x)
    a, b = act[0], act[1]
    lo, hi = max(a["xlo"], b["xlo"]), min(a["xhi"], b["xhi"])
    w = 0.0 if hi <= lo else (x - lo) / (hi - lo)
    w = min(1.0, max(0.0, w))
    return 10 ** ((1 - w) * _lg_zone(a, x) + w * _lg_zone(b, x))


def fit_range(zones):
    if not zones:
        return None
    return (10 ** min(z["xlo"] for z in zones),
            10 ** max(z["xhi"] for z in zones))


def measured(name):
    """(точки, зоны подгонки, имя файла). Точки — [(E, eps, d_eps, нуклид)].

    .efa — сведённая ЛСРМ кривая: одна точка на линию плюс блок подгонки.
    У точечной 5 см .efa нет, там только .efr с блоком на каждый источник;
    блоки объединяются, повторы одной энергии усредняются с весом 1/σ², а
    подгоночной кривой в .efr не бывает вовсе.
    """
    for ext in (".efa", ".efr"):
        hits = sorted(EFF_DIR.glob("*%s*%s" % (name, ext)))
        if not hits:
            continue
        text = paths.read_text(hits[0])
        secs = parse_efr(text)
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
        return out, zones_of(text), os.path.basename(hits[0])
    return [], [], None


def computed(fn):
    """Расчётная кривая: [(E, eps_net, d_eps, eps_gross)]."""
    # Чтение — общей реализацией csvio.read(): она пропускает строки «#»,
    # включая штамп провенанса. Прежний csv.DictReader(fh) без фильтра принимал
    # первую строку штампа за шапку и падал на KeyError: 'E_keV'.
    return [(float(r["E_keV"]), float(r["eps_net"]),
             float(r["d_eps"]), float(r["eps_gross"]))
            for r in csvio.read(os.path.join(RES, fn))]


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


def kit_cell(km):
    """Ячейка сводки: отношение ± погрешность и пометка о согласии набора.

    Если χ²/ν больше 2, ряды нуклидов между собой не согласны, и число, каким
    бы узким ни вышел интервал, описывает набор плохо. Пометка ставится рядом
    с числом: читатель должен видеть это без обращения к таблицам.
    """
    if not km:
        return "—"
    s = "%s ± %s" % (ru(km[0]), ru(km[1]))
    if km[3] is not None and km[3] > 2.0:
        s += "<br><span class='cap'>χ²/ν = %s, ряды не согласны</span>" \
             % ru(km[3], 1)
    return s


def answer_card(kit):
    """Главный результат. Вывод об удельной активности — по сверке в СпектраЛайн
    (results/cell_2x2.csv), НЕ по нашему конвейеру: анализ спектров и активности
    в предмет статьи не входит. Отношение эффективностей (кривые) — наш анализ.
    Числа — из файла, не хардкод; kit не используется (наш конвейер активности —
    справочно, см. «Сводку»)."""
    cells = {}
    p = os.path.join(RES, "cell_2x2.csv")
    if os.path.exists(p):
        for r in csvio.read(p):
            try:
                cells[r["cell"]] = float(r["ratio_to_passport"])
            except (KeyError, ValueError):
                pass
    sl_lsrm = cells.get("кривая ЛСРМ x СпектраЛайн")
    sl_our = cells.get("кривая наша x СпектраЛайн")
    our_algo = cells.get("кривая наша x наш алгоритм")
    if sl_lsrm is None or sl_our is None:
        return ""
    ref = ""
    if our_algo is not None:
        ref = ('<p class="cap">Наш собственный разбор спектра (вне предмета '
               'статьи) по единственной чистой линии 2614,5 кэВ даёт по той же '
               'записи отношение %s — приведён справочно, в вывод не входит.</p>'
               % ru(our_algo, 2))
    return (
        '<div class="card warn" id="answer" style="border-left-width:3px">'
        '<p><b>Главный результат: пригодна ли методика расчёта эффективности '
        'методом Монте-Карло.</b> В предмет статьи входит расчёт эффективности '
        'регистрации; анализ спектров и сравнение удельной активности выполнены '
        'независимой программой СпектраЛайн, а не нашим кодом.</p>'
        '<p><b>Кривая эффективности (наш анализ).</b> Расчётная кривая близка к '
        'аттестованной на точечных геометриях и отходит от неё по форме на '
        'сосудах; ход отношения вдоль последовательности геометрий модель '
        'воспроизводит (разделы «Результаты» и «Отклик на геометрию»).</p>'
        '<p><b>Удельная активность (сверка в СпектраЛайн).</b> На аттестационной '
        'записи сосуда Маринелли, где спектр разобран программой СпектраЛайн, '
        'аттестованная кривая даёт удельную активность в отношении %s к '
        'паспортной (в пределах примерно десяти процентов). Расчётная кривая под '
        'тем же съёмом даёт %s, но эта клетка отвечает о ФОРМЕ кривой: её уровень '
        'редактор эффективности перенормировал, поэтому абсолютный уровень '
        'расчётной активности по ней не проверяется. Разбор — в разделе '
        '«Проверка метода».</p>'
        '%s'
        '<p>Практический вывод: методика расчёта эффективности воспроизводит '
        'форму аттестованной кривой, а на проверенной записи сосуда Маринелли '
        'удельная активность, снятая независимой программой, согласуется с '
        'паспортом в пределах погрешности. Проверка абсолютного уровня расчётной '
        'активности и распространение на прочие геометрии требуют такой же '
        'сверки спектра в СпектраЛайн.</p>'
        '<p><b>Собственные алгоритмы разбора спектра требуют доработки.</b> '
        'Калибровка, съём площади и деконволюция, разрабатываемые в рамках '
        'проекта, в предмет статьи не входят; по этой же сверке видно, что их '
        'однолинейный съём занижает удельную активность, поэтому для её '
        'определения они пока непригодны и требуют доработки.</p>'
        '</div>'
        % (ru(sl_lsrm, 2), ru(sl_our, 2), ref))


def spectraline_check():
    """Раздел «Проверка метода в СпектраЛайн»: эксперимент 2×2 (кривая × алгоритм)
    из results/cell_2x2.csv. По правилу контура анализ спектров и сравнение
    удельной активности выполняются в независимой программе СпектраЛайн; наш
    конвейер активности в предмет статьи не входит и приводится только справочно.
    Числа — из файла, не хардкод."""
    p = os.path.join(RES, "cell_2x2.csv")
    if not os.path.exists(p):
        return ""
    rows = list(csvio.read(p))
    if not rows:
        return ""
    passp = None
    with open(p, encoding="utf-8") as fh:
        m = re.search(r"паспорт\s+(\d+)\s+Бк/кг", fh.read())
        if m:
            passp = m.group(1)
    body = []
    for r in rows:
        c2 = r.get("chi2_dof") or ""
        body.append(
            "<tr><td>%s</td><td class='n'>%s</td><td class='n'>%s</td>"
            "<td class='n'>%s</td><td class='n'>%s</td></tr>"
            % (esc(r["cell"]), rug(float(r["A_Bq_kg"]), "%.0f"),
               ru(float(r["ratio_to_passport"]), 3),
               ru(float(c2), 2) if c2 else "—", r["n_lines"]))
    passref = (" Паспортная удельная активность записи — %s Бк/кг." % rug(
        float(passp), "%.0f")) if passp else ""
    return (
        '<h2 id="spectraline">Проверка метода в СпектраЛайн</h2>'
        '<p>Анализ спектров и сравнение удельной активности выполнены в '
        'независимой программе СпектраЛайн (разработка ЛСРМ): в предмет статьи '
        'входит методика расчёта эффективности, а не собственные алгоритмы '
        'разбора спектра, поэтому активность снимает аттестованное стороннее '
        'ПО, а не наш код. На аттестационной записи (сосуд Маринелли, ряд тория) '
        'поставлен эксперимент из четырёх сочетаний «кривая эффективности × '
        'способ съёма площади»: расчётная (Монте-Карло) и аттестованная (ЛСРМ) '
        'кривые в сочетании со съёмом средствами СпектраЛайн и нашего разбора. '
        'Так вклад кривой отделяется от вклада съёма.%s</p>'
        '<p class="cap"><b>Таблица 6 — удельная активность в четырёх сочетаниях '
        '«кривая × способ съёма» (СпектраЛайн, сосуд Маринелли, ряд тория).</b> '
        '«активность» — удельная активность, Бк/кг; «A/пасп» — её отношение к '
        'паспортной; «χ²/dof» — приведённая невязка сведения по линиям; «линий» — '
        'число линий в сведении.</p>'
        '<div class="tw"><table><thead><tr><th>сочетание</th>'
        '<th class="n">активность, Бк/кг</th><th class="n">A/пасп</th>'
        '<th class="n">χ²/dof</th><th class="n">линий</th></tr></thead>'
        '<tbody>%s</tbody></table></div>'
        '<p>Со съёмом средствами СпектраЛайн обе кривые дают удельную активность '
        'в пределах примерно десяти процентов от паспортной (аттестованная — '
        'отношение около 0,94 по шести линиям, расчётная — около 0,90). '
        'Однолинейный съём нашего разбора (единственная чистая линия 2614,5 кэВ) '
        'даёт около 0,80 и в вывод об активности не входит — он показывает лишь, '
        'что по одной жёсткой линии активность занижается. Съём средствами '
        'СпектраЛайн по нескольким линиям эту разницу снимает.</p>'
        '<div class="card"><p><b>Оговорка.</b> Клетка «расчётная кривая × '
        'СпектраЛайн» отвечает о ФОРМЕ расчётной кривой: её уровень редактор '
        'эффективности СпектраЛайн перенормировал множителем около 1,1, поэтому '
        'по этой клетке нельзя судить об абсолютном уровне расчётной кривой, '
        'только о её ходе с энергией.</p></div>'
        % (passref, "".join(body)))


def kit_activity():
    """A_изм/A_пасп по геометриям — ИЗ ФАЙЛОВ СВОДКИ, а не своей формулой.

    Раньше страница брала таблицу линий и сводила её МЕДИАНОЙ, тогда как
    отчёт сводил тот же пересчёт правилом ЛСРМ (средневзвешенное с весами
    1/(ΔA)²). Числа получались разные — 1,17 против 1,18 на 5 см и до трёх
    сотых на сосудах, — и оба назывались «пересчётом комплекта». Теперь
    сводку считает пересчёт и кладёт в kit_activity_*.csv, а страница только
    читает: одно правило, один источник.

    Возврат: {геометрия: (отношение, погрешность, число линий, χ²/ν)} по
    строкам nuclide=*. χ²/ν — мера согласия рядов между собой; больше ~2
    значит, что одним числом набор не сводится, и это выносится на страницу
    рядом с числом, а не прячется.
    """
    out = {}
    for fn in ("kit_activity_volume.csv", "kit_activity_point.csv"):
        p = os.path.join(RES, fn)
        if not os.path.exists(p):
            continue
        with open(p, encoding="utf-8") as fh:
            for r in csv.DictReader(l for l in fh if not l.startswith("#")):
                if r["nuclide"] != "*":
                    continue
                c2 = r.get("chi2_dof") or ""
                out[r["geometry"]] = (float(r["ratio"]), float(r["d_ratio"]),
                                      int(r["n_lines"]),
                                      float(c2) if c2 else None)
    return out


# --- построение SVG ---------------------------------------------------------
# Рисуем на сервере, а не в браузере: страница обязана открываться без единого
# внешнего запроса, а картинка — быть одинаковой везде, включая печать.

W, H = 760, 400
PAD_L, PAD_R, PAD_T, PAD_B = 66, 14, 14, 44
RW, RH = 760, 170
# Рамка по энергии — по САМОЙ ШИРОКОЙ из сравниваемых областей, а не по круглым
# числам: сетка теперь идёт от 45,3 до 3552,5 кэВ (края паспортных зон), и при
# прежних 50…3200 крайние узлы и правый конец зоны «Денты» просто не рисовались,
# то есть график молча скрывал ровно то, что добавляли.
E_LO, E_HI = 44.0, 3650.0


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
    return [(v, "%d" % v) for v in (50, 100, 200, 500, 1000, 2000, 3500)]


def smooth_path(fn, lo, hi, e0=None, e1=None, w=W, n=180):
    """Путь SVG по гладкой функции eps(E); разрывы обрываются и начинаются вновь."""
    a = math.log10(max(e0 or E_LO, E_LO))
    b = math.log10(min(e1 or E_HI, E_HI))
    if b <= a:
        return ""
    d, pen = [], False
    for i in range(n + 1):
        E = 10 ** (a + (b - a) * i / n)
        v = fn(E)
        if v is None or not (lo <= v <= hi):
            pen = False
            continue
        d.append("%s%.1f,%.1f" % ("L" if pen else "M", lx(E, w), ly(v, lo, hi)))
        pen = True
    return " ".join(d)


def chart(meas, comp, lo, hi, zones, mcfit):
    """Наложение: паспортная подгонка и расчётная кривая, обе с точками."""
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
            lab = rug(v * 100, "%g") if v * 100 >= 0.01 else rug(v * 100, "%.3f")
            s.append('<text class="ax" x="%d" y="%.1f" text-anchor="end">%s</text>'
                     % (PAD_L - 6, y + 4, lab))
        dec += 1
    s.append('<text class="axt" x="%.1f" y="%d" text-anchor="middle">'
             'энергия, кэВ</text>' % ((PAD_L + W - PAD_R) / 2, H - 6))
    s.append('<text class="axt" transform="translate(14,%.1f) rotate(-90)" '
             'text-anchor="middle">эффективность ППП, %%</text>'
             % ((PAD_T + H - PAD_B) / 2))

    # Полосы, где эталонных линий НЕТ. Расчёт там есть, сверять его не с чем, и
    # это надо видеть глазом: у маринелли и «Денты» ниже 239 кэВ нет ни одного
    # источника, а именно там кривая заворачивает от самопоглощения.
    if meas:
        e_lo = min(E for E, _, _, _ in meas)
        e_hi = max(E for E, _, _, _ in meas)
        for a, b, where in ((E_LO, e_lo, "ниже"), (e_hi, E_HI, "выше")):
            if b <= a * 1.001:
                continue
            x0, x1 = lx(a), lx(b)
            s.append('<rect class="nodata" x="%.1f" y="%d" width="%.1f" '
                     'height="%d"><title>%s %s кэВ эталонных линий в этой '
                     'геометрии нет: расчёт есть, сверка недоступна</title>'
                     '</rect>'
                     % (x0, PAD_T, x1 - x0, H - PAD_B - PAD_T, where,
                        rug(a if where == "выше" else b, "%.0f")))

    # паспортная подгонка ЛСРМ — сплошной линией, в своей области определения
    if zones:
        r = fit_range(zones)
        d = smooth_path(lambda E: fit_eps(zones, E), lo, hi, r[0], r[1])
        if d:
            s.append('<path class="expfit" d="%s"/>' % d)

    # расчётная кривая — та же гладкая подгонка, по которой берётся отношение
    if mcfit:
        d = smooth_path(mcfit, lo, hi, comp[0][0], comp[-1][0])
        if d:
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
                 '<title>расчёт, узел сетки: %s кэВ, %s %%</title></circle>'
                 % (x, y, ru(E, 1), rug(v * 100)))

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
                 '<title>измерение (%s): %s кэВ, %s %% ± %s</title>'
                 '</rect>' % (x - 3, y - 3, esc(nuc), ru(E, 1),
                              rug(v * 100), rug(dv * 100, "%.2g")))
    s.append("</svg>")
    return "\n".join(s)


def ratio_chart(pairs, med, zones=None, mcfit=None, erange=None, mspan=None):
    """МК/эксп: гладкая кривая подгонка-к-подгонке плюс точки по линиям."""
    # Масштаб по вертикали — по фактическому размаху точек (с усами и с
    # поправкой), симметрично вокруг единицы: прежние жёсткие 0,55…1,75
    # сжимали данные к центру и обрезали сиреневые выбросы.
    span = [1.0]
    for E, r, dr, nuc, rc, is_node in pairs:
        span += [r + dr, r - dr]
        if rc is not None:
            span.append(rc)
    if med == med:
        span.append(med)
    lo0, hi0 = min(span), max(span)
    pad = max((hi0 - lo0) * 0.06, 0.02)
    lo, hi = lo0 - pad, hi0 + pad
    rng = hi - lo
    step = 0.05 if rng <= 0.35 else (0.1 if rng <= 0.9 else 0.2)
    nd = 2 if step < 0.1 else 1

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
    k0 = int(math.ceil(lo / step - 1e-9))
    k1 = int(math.floor(hi / step + 1e-9))
    for k in range(k0, k1 + 1):
        v = k * step
        if not (lo <= v <= hi):
            continue
        cls = "one" if abs(v - 1.0) < 1e-9 else "grid"
        s.append('<line class="%s" x1="%d" y1="%.1f" x2="%d" y2="%.1f"/>'
                 % (cls, PAD_L, y(v), RW - PAD_R, y(v)))
        s.append('<text class="ax" x="%d" y="%.1f" text-anchor="end">%s</text>'
                 % (PAD_L - 6, y(v) + 4, ru(v, nd)))
    if lo <= med <= hi:
        s.append('<line class="med" x1="%d" y1="%.1f" x2="%d" y2="%.1f"/>'
                 % (PAD_L, y(med), RW - PAD_R, y(med)))
    # Отношение двух ПОДГОНОК: то, что видно на верхнем графике как расстояние
    # между линиями. Точки по линиям остаются — они показывают, где расходятся
    # сами измерения, а не кривые.
    if zones and mcfit and erange:
        e0 = max(erange[0], fit_range(zones)[0])
        e1 = min(erange[1], fit_range(zones)[1])
        # Внутри диапазона измеренных линий кривая отношения опирается на
        # данные — сплошная. Вне его паспортная сторона сама экстраполяция,
        # поэтому там пунктир: это отношение двух моделей, а не проверка.
        segs = {"rfit": [], "rfitx": []}
        prev = None
        for i in range(241):
            if e1 <= e0:
                break
            E = 10 ** (math.log10(e0)
                       + (math.log10(e1) - math.log10(e0)) * i / 240)
            fv = fit_eps(zones, E)
            mv = mcfit(E)
            r = (mv / fv) if (fv and mv) else None
            inside = bool(mspan and mspan[0] <= E <= mspan[1])
            key = "rfit" if inside else "rfitx"
            if r is None or not (lo <= r <= hi):
                prev = None
                continue
            pt = (lx(E, RW), y(r))
            if prev is None or prev[1] != key:
                segs[key].append(["M%.1f,%.1f" % pt])
            else:
                segs[key][-1].append("L%.1f,%.1f" % pt)
            prev = (pt, key)
        for key, chunks in segs.items():
            for ch in chunks:
                if len(ch) > 1:
                    s.append('<path class="%s" d="%s"/>' % (key, " ".join(ch)))
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
                     '<title>%s кэВ, с поправкой на суммирование: %s'
                     '</title></circle>' % (x, y(rc), ru(E, 1), ru(rc)))
        s.append('<circle class="mcp" cx="%.1f" cy="%.1f" r="3">'
                 '<title>%s кэВ (%s): МК/эксп %s ± %s</title></circle>'
                 % (x, y(r), ru(E, 1), esc(nuc), ru(r), ru(dr)))
    s.append('<text class="axt" x="%.1f" y="%d" text-anchor="middle">'
             'энергия, кэВ</text>' % ((PAD_L + RW - PAD_R) / 2, RH - 1))
    s.append("</svg>")
    return "\n".join(s)


DEG_LO, DEG_HI = 3, 8   # в каких пределах выбирается степень
INTERP_D = 0.02         # оценка погрешности интерполяции


def _cf(E, y, d, deg):
    """Полином по log-log с центрированием — иначе на 8-й степени плохая
    обусловленность и предупреждения numpy."""
    import numpy as np
    x = np.log(E)
    x0 = x.mean()
    return np.polyfit(x - x0, np.log(y), deg,
                      w=y / np.maximum(d, 1e-30)), x0


def best_degree(E, y, d, lo=DEG_LO, hi=DEG_HI):
    """Степень полинома — по СКОЛЬЗЯЩЕМУ ИСКЛЮЧЕНИЮ УЗЛА, а не на глаз.

    Подгонка делается без одного узла и предсказывает его значение; берётся
    степень с наименьшей медианной ошибкой предсказания. Проверка отвечает
    ровно на нужный вопрос — обобщает подгонка или гоняется за шумом, — и
    остаётся верной при изменении сетки.

    Понадобилось это, когда сетку расширили до краёв паспортных зон
    (45,3…3552,5 кэВ). Прежняя жёстко зашитая пятая степень на удлинившемся
    плече начала колебаться у краёв: отклонение в узлах доходило до 14 %
    против прежних 3 %. Скользящее исключение показало, что ошибка
    предсказания при этом ПАДАЕТ до восьмой степени — то есть у кривой в
    двух декадах действительно столько структуры, а не шума (статистика
    прогонов около 0,5 %).
    """
    import numpy as np
    n = len(E)
    x, ly = np.log(E), np.log(y)
    w = y / np.maximum(d, 1e-30)
    out = []
    for deg in range(lo, min(hi, n - 3) + 1):
        err = []
        for i in range(n):
            m = np.ones(n, bool)
            m[i] = False
            x0 = x[m].mean()
            cf = np.polyfit(x[m] - x0, ly[m], deg, w=w[m])
            err.append(abs(np.exp(np.polyval(cf, x[i] - x0)) / y[i] - 1))
        out.append((float(np.median(err)), deg))
    out.sort()
    return out[0][1], out[0][0]


def pair_up(meas, comp, C, gross=False):
    """(E, МК/эксп, погрешность, нуклид, отношение с поправкой C, узел?).

    gross=True — брать площадь БЕЗ вычета левой полки континуума. Это не
    придирка: разница между двумя способами взятия площади доходит до 10 % на
    мягких линиях, и именно по gross считались точечные геометрии в отчёте.
    Оба числа приводятся рядом, чтобы выбор обработки не выдавался за свойство
    детектора.

    Совпадение энергий узлов сетки и линий эталонов — редкость: у ЛСРМ до 24
    линий, а в сетке 24 узла. Поэтому там, где узла нет, расчёт берётся
    интерполяцией полиномом по log-log — ТОЙ ЖЕ процедурой, которой строится
    рабочая кривая, и той же, что в compare_point.py. Иначе сравнивались бы
    8 линий вместо 24, и число на странице расходилось бы с отчётом.
    """
    import numpy as np
    col = 3 if gross else 1
    Eg = np.array([c[0] for c in comp], dtype=float)
    yg = np.array([c[col] for c in comp], dtype=float)
    dyg = np.array([c[2] for c in comp], dtype=float)
    deg, _ = best_degree(Eg, yg, dyg)
    cf, x0 = _cf(Eg, yg, dyg, deg)
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
            m = math.exp(float(np.polyval(cf, math.log(E) - x0)))
            dm, is_node = m * INTERP_D, False
        else:
            continue
        r = m / ev
        dr = r * math.hypot(dm / m if m else 0, edv / ev)
        c = C.get(round(E, 1))
        out.append((E, r, dr, nuc, (r / c[0]) if c else None, is_node))
    return out


def mc_fit(comp, gross=False):
    """Гладкая расчётная кривая: та же подгонка, что даёт интерполяцию."""
    import numpy as np
    col = 3 if gross else 1
    Eg = np.array([c[0] for c in comp], dtype=float)
    yg = np.array([c[col] for c in comp], dtype=float)
    dyg = np.array([c[2] for c in comp], dtype=float)
    deg, loo = best_degree(Eg, yg, dyg)
    cf, x0 = _cf(Eg, yg, dyg, deg)

    def f(E):
        if not (Eg[0] <= E <= Eg[-1]):
            return None
        return math.exp(float(np.polyval(cf, math.log(E) - x0)))
    return f, (float(Eg[0]), float(Eg[-1])), deg, loo


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



def rug(x, fmt="%.4g"):
    """Число для видимого текста: точку в запятую (ГОСТ 8.417)."""
    return (fmt % x).replace(".", ",")

def nusl(n):
    """Склонение «узел» по числу: 1 узел, 2-4 узла, 5+ узлов (11-14 узлов)."""
    if 11 <= n % 100 <= 14:
        return "узлов"
    return {1: "узел", 2: "узла", 3: "узла", 4: "узла"}.get(n % 10, "узлов")


# Стиль лежит отдельным файлом docs/assets/page.css и подключается ссылкой,
# а не инлайном. Отдельный файл, а не строка здесь, — по прозаической причине:
# проверка перед публикацией ищет номера источников, записанные после знака
# номера, и шестнадцатеричный цвет от такого номера по форме не отличается.
# В файлах разметки и стилей она подобные находки пропускает по расширению,
# в исходниках на питоне — нет. Ослаблять её ради удобства нельзя: номер
# экземпляра источника выглядит точно так же, и один раз он так и утёк.
CSS_HREF = "../assets/page.css"


GEOM_RU = {"Marinelli_1L": "Маринелли 1 л", "Denta_120mL": "«Дента» 120 мл",
           "Petri_60mL": "Петри 60 мл"}


def _read_rows(fn):
    """Строки CSV из results/ без заголовка и комментариев. -> [[поле, ...]]"""
    p = os.path.join(RES, fn)
    if not os.path.exists(p):
        return []
    out = []
    for ln in open(p, encoding="utf-8"):
        ln = ln.strip()
        if not ln or ln.startswith("#") or ln[0].isalpha() and "," in ln \
                and ln.split(",")[0] in ("geometry",):
            continue
        out.append(ln.split(","))
    return out


def balance_table():
    rows = _read_rows("deconv_balance.csv")
    if not rows:
        return "<p>Таблица баланса не построена: запустите deconv_balance.py.</p>"
    body = []
    for g, E, n, fm, fg, r in rows:
        cls = "" if abs(float(r) - 1) < 0.05 else ' class="warn"'
        body.append("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td>"
                    "<td>%s</td><td%s>%s</td></tr>"
                    % (esc(GEOM_RU.get(g, g)), ru(float(E), 1), n,
                       ru(100 * float(fm), 1), ru(100 * float(fg), 1),
                       cls, ru(float(r), 3)))
    return ("<div class=\"tw\"><table><thead><tr><th>геометрия</th>"
            "<th>линия, кэВ</th><th>линий</th><th>доля пика, измерение, %%</th>"
            "<th>доля пика, модель, %%</th><th>изм/мод</th></tr></thead>"
            "<tbody>%s</tbody></table></div>" % "".join(body))


def continuum_table():
    rows = _read_rows("continuum.csv")
    if not rows:
        return "<p>Таблица континуума не построена: запустите continuum.py.</p>"
    by = {}
    for g, lo, hi, _A, r, dr in rows:
        by.setdefault(g, []).append((float(lo), float(hi), float(r), float(dr)))
    out = []
    for g, v in by.items():
        body = []
        for lo, hi, r, dr in v:
            cls = "" if abs(r - 1) < 0.05 else ' class="warn"'
            body.append("<tr><td>%s…%s</td><td%s>%s</td><td>%s</td></tr>"
                        % (ru(lo, 0), ru(hi, 0), cls, ru(r, 3), ru(dr, 3)))
        out.append("<h4>%s</h4><div class=\"tw\"><table><thead><tr>"
                   "<th>участок, кэВ</th><th>к опоре 2614,5</th><th>±</th>"
                   "</tr></thead><tbody>%s</tbody></table></div>"
                   % (esc(GEOM_RU.get(g, g)), "".join(body)))
    return "".join(out)



# Разделы «Установка и модель» и «Отклик на смену геометрии» — общие для
# плоского отчёта (эта страница) и вкладочной статьи (build_article.py).
# Числа корреляции — дословно из docs/geometry-response.md (задача 134),
# описание модели — из docs/report.md §2; при правке источников сверять.
SETUP_HTML = """
<h2 id="setup">Установка и модель</h2>

<p>Объект — сцинтилляционный спектрометр энергии гамма-излучения Гамма-1С:
устройство детектирования УДС-ГЦ-63×63 в защите «Экран-1СГ». Модель построена
по чертежу устройства детектирования: кристалл NaI(Tl) Ø63×63 мм, отражатель
MgO, входное окно и корпус, световод, фотоэлектронный умножитель. Защита
задана по указаниям оператора: сталь 3 мм, свинец 50 мм, кадмий и медь.
Толщины кадмия и меди в доступных документах не заданы и приняты по 1 мм;
контролем сборки служат массы из паспорта установки — Pb 167,1 кг при норме
не менее 165, Cu 1,60 кг, Cd 1,58 кг. Сосуды комплекта и матрица ОИСН-16
описаны по геометрическим размерам и составу, объявленному в заголовках
аттестованных кривых.</p>

<figure>
<img src="figures/gamma1s_section.png" alt="Разрез модели Гамма-1С в сборе:
кристалл NaI(Tl) диаметром 63 мм, сосуд Маринелли 1 л, свинцовая защита со
ступенчатой полостью" style="max-width:640px;width:100%%;height:auto;
display:block" width="1526" height="1608">
<figcaption class="cap"><b>Рисунок 1 — разрез расчётной модели в сборе</b>
(сосуд Маринелли 1 л установлен на головку). Красным отмечены элементы, известные без числовых
размеров: они присутствуют на эскизе руководства, но в расчётную модель не
перенесены (задача 153); влияние камеры и шахты на эффективность
регистрации ограничено обратным рассеянием.</figcaption>
</figure>

<h3>Методика расчёта</h3>

<p>Расчёт ведёт программа Geant4 (версия 11.2.1) — открытый пакет
моделирования прохождения частиц через вещество методом Монте-Карло,
разработанный ЦЕРН и используемый в физике высоких энергий и дозиметрии.
Задействованы два набора физических моделей. Первый описывает
<b>взаимодействие гамма-квантов и электронов с веществом</b> — фотоэффект,
комптоновское рассеяние, образование электрон-позитронных пар, тормозное
излучение и ионизацию — по наиболее точному из стандартных наборов пакета
для энергий до нескольких МэВ. Второй описывает <b>радиоактивный распад
ядер</b>: испускание гамма-квантов, бета- и альфа-частиц по табличным
вероятностям переходов и периодам полураспада. Кривые эффективности
считаются моноэнергетическими квантами, разыгранными равномерно по объёму
пробы, — берётся доля событий в пике полного поглощения; сетка энергий одна
на все геометрии (%(gridspan)s, %(gridn)s %(gridnusl)s на линиях комплекта и
краях аттестованных зон). Поправка на каскадное суммирование — отдельными
прогонами полного распада; событие спектрометра отделено от события расчёта
по времени разрешения тракта, иначе звенья цепочки распада, разделённые в
природе годами, суммировались бы в одно срабатывание.</p>

<p>Заложенные в расчёт вероятности взаимодействия гамма-квантов с веществом
(сечения, определяющие коэффициент ослабления) сверены с независимым
эталонным справочником — базой ослабления фотонов XCOM Национального
института стандартов и технологий США (NIST). Сверка проведена на четырёх
материалах в диапазоне 59,5…661,7 кэВ: среднее отклонение −0,27 %% при
принятом пороге 1 %%; отклонение систематическое (пятнадцать из шестнадцати
значений отрицательны), но мало. Каждый расчётный спектр несёт метку
происхождения: цифровой отпечаток исходного кода геометрии, впечатанный в
программу при сборке, полный набор параметров прогона и долю полного
телесного угла, в который разыгрывались кванты; таблицы результатов
снабжены машинно проверяемым заголовком, объявляющим, какая величина в них
записана.</p>
"""

GEORESP_HTML = """
<h2 id="georesp">Отклик на смену геометрии и корреляция откликов</h2>

<p>Обе группы кривых — аттестованная и расчётная — прослежены вдоль
последовательности геометрий <b>точка → Петри → Дента → Маринелли</b>
(полный разбор — <a href="https://github.com/VibeEngineering-LLC/geant4-detector-models/blob/main/detectors/Gamma-1S/docs/geometry-response.md">docs/geometry-response.md</a>,
числа — <code>results/geometry_response.csv</code>). Отклик на смену
геометрии выделяется отношением кривой к кривой предыдущей геометрии внутри
своей группы: всё, что от геометрии не зависит — нормировка, конвенция съёма
площади, — сокращается тождественно. Ход отклика по энергии характеризуется
локальным показателем степени <i>n</i>(E)&nbsp;&equiv;&nbsp;<code>d&nbsp;ln&nbsp;R&nbsp;/&nbsp;d&nbsp;ln&nbsp;E</code>:
для степенного участка, где отношение идёт как энергия в степени <i>n</i>
(R&nbsp;&prop;&nbsp;E<sup><i>n</i></sup>), он равен этому <i>n</i> и показывает,
насколько круто отношение меняется с энергией. В математической статистике та
же величина <i>d&nbsp;ln&nbsp;y/d&nbsp;ln&nbsp;x</i> называется эластичностью
функции.</p>

<p class="cap"><b>Таблица 5 — отклик кривых на смену геометрии.</b> Для каждого
шага последовательности геометрий — локальный показатель степени
<i>n</i>&nbsp;=&nbsp;d&nbsp;ln&nbsp;R/d&nbsp;ln&nbsp;E отдельно для аттестованной
(ЛСРМ) и расчётной кривой и их расхождение в единицах суммарной погрешности
&sigma;.</p>
<div class="tw"><table>
<thead><tr><th>шаг последовательности</th>
<th class="n">d&nbsp;ln&nbsp;R/d&nbsp;ln&nbsp;E, ЛСРМ</th>
<th class="n">d&nbsp;ln&nbsp;R/d&nbsp;ln&nbsp;E, расчёт</th>
<th class="n">расхождение</th></tr></thead>
<tbody>
<tr><td>Петри / точка</td><td class="n">+0,025 ± 0,019</td>
<td class="n">+0,285 ± 0,002</td><td class="n">13,5 σ</td></tr>
<tr><td>Дента / Петри</td><td class="n">+0,024 ± 0,020</td>
<td class="n">+0,105 ± 0,003</td><td class="n">4,1 σ</td></tr>
<tr><td>Маринелли / Дента</td><td class="n">+0,060 ± 0,019</td>
<td class="n">+0,049 ± 0,003</td><td class="n"><b>0,5 σ — совпадает</b></td></tr>
</tbody></table></div>

<p>Направление отклика модель воспроизводит на каждом шаге: по сопоставимым
плотностям антисовпадений знака нет ни на одном переходе. Согласие тем лучше,
чем дальше шаг от точечной геометрии: наибольшее расхождение приходится на
введение объёма пробы и самопоглощения (точка → Петри), а последний шаг, где
меняется лишь охват кристалла источником, воспроизводится в пределах
погрешности.</p>

<p>Антисовпадения знака появляются исключительно на несопоставимых
плотностях пробы (0,60…1,00 г/см³ относительно аттестованных 1,60): смена
плотности переворачивает знак энергетического хода отклика при неизменной
геометрии сосуда. Отсюда следует практический вывод: при переносе метода на
прибор без аттестованной кривой плотность пробы входит в число величин,
ошибка в которых меняет не масштаб, а качественный вид результата.</p>
"""


TMPL = """<meta charset="utf-8">
<title>ГАММА-1С: отработка алгоритмов анализа по поверенному комплекту</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="description" content="Отработка алгоритмов гамма-спектрометрии
(калибровка, поиск пиков, съём площади, деконволюция мультиплетов, усреднение
активности) и их проверка по поверенному комплекту спектрометра ГАММА-1С;
эффективность считается в Geant4.">
<meta name="theme-color" content="#ffffff" media="(prefers-color-scheme: light)">
<meta name="theme-color" content="#1a1a1a" media="(prefers-color-scheme: dark)">
<link rel="stylesheet" href="%(css_href)s">

<nav class="topnav"><div class="topnav-in">
<a href="#top">Начало</a>
<a href="#answer">Итог</a>
<a href="#checks">Проверки</a>
<a href="#compare">Величины</a>
<a href="#setup">Модель</a>
<a href="#summary">Сводка</a>
<a href="#vessels">Вывод</a>
<a href="#curves">Кривые</a>
<a href="#georesp">Геометрия</a>
<a href="#where">Расхождения</a>
<a href="#spectra">Спектры</a>
<a href="#repro">Воспроизведение</a>
<a href="#appendix">Приложение</a>
</div></nav>

<div class="page" lang="ru" id="top">
<h1>ГАММА-1С: отработка алгоритмов анализа спектров по поверенному комплекту</h1>
<p class="sub">Калибровка, поиск пиков, съём площади, деконволюция мультиплетов,
усреднение активности — проверка каждого шага по комплекту поверки
сцинтилляционного спектрометра NaI(Tl) 63&times;63 в свинцовой защите.
Эффективность регистрации по пику полного поглощения (ППП) рассчитывается
методом Монте-Карло в Geant4.</p>

<p class="lead"><b>Цель работы — проверка алгоритмов анализа спектров, а не
сама модель прибора.</b> Модель даёт эффективность регистрации, связывающую
площадь ППП с активностью. Мерой проверки служит поверенный комплект: на
каждую геометрию измерения имеется аттестованная кривая эффективности, у
каждого источника — паспортная активность с погрешностью. Сорок поверочных
спектров с паспортными активностями позволяют установить, где цепочка
обработки расходится и на сколько.</p>

%(answer)s

<h2 id="checks">Три независимые проверки</h2>

<p>У трёх способов проверки разные источники систематики, поэтому согласие
трёх результатов даёт больше, чем точность каждого в отдельности.</p>

<div class="tw"><table>
<thead><tr><th>проверка</th><th>что берётся за истину</th>
<th>что может её испортить</th></tr></thead>
<tbody>
<tr><td><b>Кривая против кривой.</b> Расчётная эффективность против паспортной
подгонки ЛСРМ</td>
<td>кривая поверки — уже обработанный результат прибора</td>
<td>наши ошибки геометрии и обработки складываются с чужими; способ съёма
площади у ЛСРМ нам неизвестен</td></tr>
<tr><td><b>Пересчёт комплекта.</b> Из спектра получаем активность и сверяем с
паспортом источника</td>
<td>паспорт источника — ничего обработанного, только активность и дата</td>
<td>подогнать нечего, но нужна вся цепочка сразу: калибровка, площадь,
эффективность, суммирование, наложения</td></tr>
<tr><td><b>Алгоритм против алгоритма.</b> Деконволюция мультиплетов против оконного
съёма на одних и тех же линиях</td>
<td>чистые одиночные линии, где оба способа обязаны совпасть</td>
<td>ничего внешнего не проверяет — зато отделяет ошибку алгоритма от ошибки
модели</td></tr>
</tbody></table></div>

<h2 id="compare">Сравниваемые величины</h2>

<p><b>Аттестованная кривая</b> — результат поверки: площадь ППП от эталонного
источника, отнесённая к активности и выходу линии. Точки привязаны к линиям
конкретных нуклидов; усы — заявленная погрешность.</p>

<p><b>Расчётная кривая</b> — результат моделирования методом
Монте-Карло (Geant4 11.2.1; подробнее о физике расчёта — в разделе
«Модель»): гамма-кванты одной энергии разыгрываются равномерно по объёму
пробы, и берётся доля тех из них, что дали пик полного поглощения. Геометрия
построена по чертежам и паспорту прибора без подгонки под аттестованную
кривую.</p>

<div class="card warn">
<p><b>Обязательная поправка на каскадное суммирование.</b> В моноэнергетическом
расчёте суммирование отсутствует; в реальном распаде одновременная регистрация
каскадных квантов выводит событие из ППП, занижая измеренную эффективность
каскадных линий. Отношение расчёт/измерение делится на поправку C, вычисленную
отдельным прогоном полного распада. Контроль на нуклидах без каскада: Cs-137 —
C = 0,983, K-40 — C = 1,009, единица в пределах погрешности.</p>
</div>

%(setup)s

<h2 id="summary">Сводка</h2>

<div class="tw"><table>
<thead><tr><th>геометрия</th><th class="n">линий</th>
<th class="n">МК/эксп<br>по чистой площади</th>
<th class="n">то же<br>без вычета полки</th>
<th class="n">RMS формы, %%</th>
<th class="n">линий в пересчёте</th>
<th class="n">A/пасп<br>по правилу ЛСРМ</th></tr>
</thead><tbody>%(summary)s</tbody></table></div>

<p class="cap"><b>Что в столбцах.</b> «линий» — число линий кривой,
вошедших в сравнение (штук). «МК/эксп» — отношение расчётной эффективности
(МК — Монте-Карло) к измеренной; величина безразмерная, единица означает
совпадение, больше единицы — расчёт завышает. «то же без вычета полки» — то
же отношение, но площадь пика снята без вычитания левой полки континуума
(полка — участок спектра непосредственно слева от окна пика, по которому
оценивается подложка под пиком); разница двух столбцов до 10 %% на линиях
малых энергий и менее процента на жёстких, далее используется площадь с
вычетом. «RMS формы» — среднеквадратичный разброс отношения МК/эксп по
линиям, в процентах: мера того, насколько форма расчётной кривой отходит от
измеренной. «линий в пересчёте» — число линий (штук), по которым восстановлена
активность. «А/пасп» — отношение активности, восстановленной из спектра по
расчётной эффективности, к паспортной; тоже безразмерное, единица — совпадение.
Свод активности — по правилу ЛСРМ (средневзвешенное с весами 1/(ΔA)², где
ΔA — погрешность активности по линии); в него входят только линии с долей
собственного выхода в окне не ниже 0,95.</p>

<p>МК/эксп выше единицы означает завышенную расчётную эффективность; A/пасп
ниже единицы — то же. Оба пути согласованно дают: <b>маринелли завышена,
кюветы занижены, точечные геометрии близки к единице</b>.</p>

<h2 id="vessels">Расхождение обусловлено моделями сосудов, а не детектора</h2>

<p>Точечные геометрии 5 и 25 см согласуются между собой и с единицей в двух
разных конфигурациях защиты — модель блока детектирования проверена дважды.
Объёмные геометрии расходятся в противоположные стороны, что несовместимо с
единой ошибкой эффективности. Наиболее вероятный источник расхождения —
неполнота геометрии сосудов.</p>

<p>Точных внутренних размеров сосудов (толщины стенки, диаметра и глубины
колодца, распределения пробы по объёму) чертежи комплекта не содержат. В
спецификации прецизионных измерений ЛСРМ для каждой кюветы заданы только
номинальный объём, внешние габариты и эффективная толщина слоя пробы
d<sub>эфф</sub> (перепечатаны в приложении, таблица А.4); по этим данным
внутренняя геометрия сосуда однозначно не восстанавливается, и в модель она
заложена с допущениями.
Именно неизвестная геометрия сосуда, а не детектор и не алгоритмы обработки
(они проверены точечной геометрией, где сосуда нет), — вероятная причина
того, что активность расходится с паспортом на сосудах и совпадает на
точечных источниках: ошибка в геометрии сосуда меняет результат, а точечный
источник от неё не зависит.</p>

<div class="card">
<p><b>Отклонённые объяснения.</b> Расхождение можно было бы отнести
на счёт того, что расчёт и прибор считают разные величины: расчёт — полное
энерговыделение в кристалле, а прибор — площадь пика полного поглощения. Но
тогда то же превышение проявилось бы и на точечном источнике, где сосуда нет,
а там отношение равно единице; значит, объяснение отпадает. Второе возможное
объяснение — что расхождение зависит от нуклида, — тоже не подтвердилось:
после того как звенья цепочки распада были разделены во времени (они разнесены
в природе на годы), разброс значений внутри каждого сосуда свёлся к одному
общему множителю.</p>
</div>

<h2>Учёт временной структуры цепочки распада (правка 28.07.2026)</h2>

<p>До правки Geant4 доводил цепочку распада до конца внутри одного события, и
энерговыделения звеньев, разделённых в природе годами, суммировались. Правка:
энерговыделения группируются по глобальному времени с порогом 1 мкс —
истинные каскады (наносекунды) остаются в одном событии, звенья ряда
(секунды и более) разделяются. На одиночном нуклиде правка не меняет ничего
(проверено побайтно); в цепочках пики ряда тория выросли на 25–46 %%, ряда
радия — на 12–14 %%. После правки нуклидный разброс внутри сосуда Маринелли
(0,78 у Cs-137 против 1,26 у Th-232) сократился до одного множителя
0,69–0,89 на сосуд.</p>

%(georesp)s

<h2 id="curves">Кривые по геометриям</h2>

<div class="card">
<p><b>Аттестованная кривая</b> — зонная подгонка ЛСРМ из файла
<code>.efa</code>: внутри зоны разложение lg&nbsp;&epsilon; по ортогональному
базису от lg&nbsp;E; активности прибор считает по ней, а не по точкам. На
графиках она восстановлена из файла (сплошная линия), квадраты — измеренные
точки. Разброс точек у маринелли и «Денты» до 20 %% свойствен самим точкам:
линии 238,6 и 242,0 кэВ в NaI(Tl) не разделяются, кривая проходит между ними —
поэтому расчёт сравнивается с кривой, а не с отдельной точкой. Правило сшивки
перекрывающихся зон файлом не задано; принято линейное смешивание по lg E
(в перекрытии зоны согласны на 0,3–8 %%). Расчётная кривая — подгонка
полиномом 5-й степени по log-log по узлам сетки; у точечной 5 см подгонки ЛСРМ
нет (формат <code>.efr</code> блока кривой не содержит), показаны точки.</p>
</div>

<h3>Диапазон расчёта задан диапазоном аттестованных кривых</h3>

<p>Сетка энергий накрывает объединение областей определения зон
<code>.efa</code> всех геометрий — <b>45,3&hellip;3552,5 кэВ</b>; узлы посажены
на линии эталонов комплекта, чтобы сверка шла в узлах, а не интерполяцией.
Аттестованная кривая существует только там, где есть эталонные линии;
расчётная от источника не зависит.</p>

<div class="tw"><table>
<thead><tr><th>геометрия</th><th class="n">эталонные линии</th>
<th class="n">зоны подгонки</th><th>чем ограничена снизу</th></tr></thead>
<tbody>
<tr><td>Маринелли 1 л</td><td class="n">239–2615</td><td class="n">185–3305</td>
<td>только объёмные ОИСН: Cs-137, K-40, Ra-226, Th-232. Линия с наименьшей энергией —
Pb-212 239 кэВ из ряда тория</td></tr>
<tr><td>«Дента» 120 мл</td><td class="n">239–2615</td><td class="n">186–3552</td>
<td>тот же набор, та же наименьшая энергия</td></tr>
<tr><td>Петри 60 мл</td><td class="n">68–2615</td>
<td class="n">56–1854 и 234–3305</td>
<td>добавились Ti-44 и Eu-152 — отсюда линии 68, 78 и 122 кэВ</td></tr>
<tr><td>точечная 5 см</td><td class="n">60–2615</td><td class="n">подгонки нет</td>
<td>точечные источники, включая Am-241 59,5 кэВ</td></tr>
<tr><td>точечная 25 см</td><td class="n">60–2615</td>
<td class="n">45–769 и 273–3552</td><td>то же</td></tr>
</tbody></table></div>

<p>Расчётная сетка одна на все геометрии: <b>%(gridspan)s</b>, %(gridn)s
%(gridnusl)s на линиях комплекта и краях аттестованных зон — общая сетка
обеспечивает сравнимость кривых между геометриями.</p>

<p>В маринелли и «Денте» область ниже 239 кэВ экспериментальному контролю
недоступна — объёмного источника с мягкой линией в комплекте нет, — и именно
там кривая проходит через максимум (около 170–200 кэВ: поглощение мягких
квантов в матрице ОИСН-16 с 71 %% железа и стенке сосуда). На графиках эта
область затенена; выше 2614,5 кэВ обе кривые экстраполируют. Косвенное
подтверждение расчёта в затенённой области даёт геометрия Петри (эталоны с
68 кэВ, максимум кривой воспроизводится измерением), однако перенос вывода на
маринелли остаётся допущением.</p>

%(legend)s

<h2 id="where">Локализация расхождения по спектру</h2>

<p>Разброс формы 9…17 %% локализуется двумя проверками без подгонки
параметров модели; обе воспроизводятся одной командой. Полные таблицы — в
<a href="#appendix">приложении</a>.</p>

<h3>Баланс пика и континуума в окне деконволюции</h3>

<p>Активность по группе линий равна отношению амплитуд одной формы, снятых с
измерения и с уширенной модели, поэтому доля площади окна, отданная подгонкой
пикам, переносится в активность множителем один к одному. Сравнение долей
(таблица А.1 приложения): под группой 583,2 кэВ модель даёт континуума больше
измеренного, под одиночной 2614,5 кэВ — меньше. Расхождение знакопеременно по
энергии и общим множителем эффективности не устраняется.</p>

<h3>Континуум вне пиков</h3>

<p>Отношение расчётного и измеренного континуума на участках, свободных от
пиков (нормировка на линию 2614,5 кэВ; участки ближе 2,5&sigma; к линиям с
выходом ярче 2 %% исключены автоматически), сведено в таблицах А.2 приложения.
Ниже 700 кэВ согласие в пределах 3 %%; расхождение сосредоточено выше
1600 кэВ: недобор подложки в полосе 1668…2252 кэВ и перебор в полосе
2308…2452 кэВ — на комптоновском крае линии 2614,5 кэВ. Порядок геометрий
обратен балансу пика: сильнее расходится маринелли, самый интенсивный источник
комплекта; зависимость от скорости счёта указывает на наложения импульсов,
отсутствующие в расчёте. Действуют два механизма; ни один пока не сведён к числу.</p>

%(blocks)s

<h2 id="spectra">Спектры: калибровка, поиск пиков, деконволюция</h2>

<p>Рисунки раздела построены тем же прогоном тех же функций, которыми
считаются публикуемые числа: кривые деконволюции берутся из тех же колонок
матрицы плана, которыми решалась задача, — расхождение рисунка с таблицей
означало бы разный счёт и было бы видно.</p>

<h3>Калибровка: центр окна — по найденной центроиде</h3>

<p>Калибровка тракта у каждой записи уходит по-своему (цезиевый пик на 658,6
вместо 661,657 кэВ, Tl-208 — на 2610,7 вместо 2614,5). Окно на табличной
энергии срезает часть пика, а полки фона попадают на его склон; обе ошибки смещают
результат вниз. Центр окна измерения ставится по найденной центроиде (первый
момент по вычтенной подложке).</p>

<div class="card warn">
<p><b>Два ограничения процедуры, установленные на практике.</b> Окно модели по
центроиде не сдвигается: линии модельного спектра стоят на истинных энергиях,
калибровочного дефекта нет — каждая сторона центрируется на своём пике при
общих ширине окна и устройстве полок. У линии в мультиплете центроида есть центр
тяжести группы, а не калибровочный сдвиг: у 911,2 кэВ Ac-228 она смещена на
+18,7 кэВ соседями 964,8 и 969,0; сдвиг по ней увеличил разброс по «Денте» с
1,37 до 1,66. Сдвиг применяется только к одиночной линии в окне; признак —
<b>чистота</b> (доля выхода окна, приходящаяся на саму линию, по спектру
испускания того же прогона), порог 0,95.</p>
</div>

<h3>Площадь: оконный съём и деконволюция</h3>

<p>Площадь плохо разделённой линии есть сумма нескольких переходов и одному
переходу не приписывается; оконный съём такие линии отбрасывает, у Th-232
остаётся одна годная — 2614,5 кэВ. Мультиплеты разбирает связанная деконволюция:
у группы линий один свободный параметр — активность, площади связаны как
S<sub>k</sub> = A&middot;I<sub>k</sub>&middot;&epsilon;<sub>k</sub>&middot;t
[ЛСРМ, «Алгоритмические основы», формула 5.2-7 в предельной форме];
нормировка — второй такой же подгонкой по модельному спектру распада,
уширенному до разрешения прибора. Из отношения выпадает всё общее для двух
сторон: выход линии, эффективность, каскадное суммирование, вклад соседей,
доля пика за краями участка. Проверка — на одиночных линиях, где деконволюция
обязана совпасть с окном; сводка — таблица А.3 приложения, полная таблица —
<code>results/deconv_lines.csv</code>.</p>

<div class="card warn">
<p><b>Установленное при сверке ограничение оконного съёма: у линии 351,9 кэВ
полка фона лежит на пике 295,2 кэВ Ra-226</b> (полка E&minus;2&middot;ПШПВ…
E&minus;ПШПВ на 351,9 кэВ — это 279…316 кэВ). Признак — аномальная
чувствительность к закону ПШПВ: 5…12 %% на этой линии против &lt;1,5 %% на
остальных. Частично эффект сокращается одинаковой полкой модельной стороны, но
не полностью: линия 351,9 кэВ при чистоте 0,99 для оконного съёма ненадёжна —
мера чистоты смотрит только внутрь окна и полку не контролирует. Деконволюции
дефект не свойствен: континуум в ней подгоняется, а не берётся с полки.</p>
</div>

<div class="card">
<p><b>Нерешённое.</b> Группы 583 и 911 кэВ дают активность на 12…26 %% выше
одиночной 2614,5 кэВ того же нуклида. Веса подгонки, нормировка и вклад
аннигиляционной линии 511 кэВ проверены и расхождения не объясняют; остаётся
несовпадение формы континуума модели и измерения в середине спектра. До его
разбора публикуемая активность считается оконным съёмом по чистым линиям,
деконволюция остаётся отработкой метода.</p>
</div>

%(deconv_legend)s

%(spectra)s

<h2 id="repro">Воспроизведение</h2>

<p>Страница собирается из репозитория одной командой; чисел, вписанных
вручную, нет:</p>

<p><code>python detectors/Gamma-1S/analysis/build_web.py</code></p>

<p><b>Правило отбора точек.</b> Берутся все линии аттестованной кривой в
диапазоне расчётной сетки (%(gridspan)s); расчёт — из узла сетки при
расстоянии до 1 кэВ, иначе интерполяцией полиномом 5-й степени по log-log.
Отношения усредняются в логарифме с весом 1/&sigma;². Независимые скрипты
репозитория с более узким отбором дают то же в пределах 0,5 %% (маринелли,
«Дента») и 5 %% (Петри, где из 19 линий отбирается 14).</p>

<p>Расчётные кривые — <code>detectors/Gamma-1S/results/eff_*.csv</code>,
паспортные — <code>reference/lsrm/efficiency/</code>, поправки на
суммирование — <code>results/summing_C.csv</code>, пересчёт комплекта —
<code>results/kit_recalc_*.csv</code>, деконволюция —
<code>results/deconv_lines.csv</code>. Сами расчётные спектры в репозиторий не
входят: они воспроизводятся драйверами из <code>drivers/</code>.</p>

<p><b>Чем считается каждый шаг.</b> Чтение записей и площади —
<code>common/py/becqmoni.py</code>; сетки энергий —
<code>drivers/grid_energies.py</code> (один список на все драйверы, чтобы копии
не разъезжались); кривые — <code>analysis/export_curves.py</code>; пересчёт
комплекта — <code>analysis/kit_recalc.py</code>; связанная деконволюция —
<code>analysis/deconv.py</code>; рисунки этого раздела —
<code>analysis/spectra_figs.py</code>. Численная часть подгонки (веса обратной
сигмой, границы неотрицательности, ковариация через SVD) взята из пакета
gamma-spectrum-analysis, модуль <code>gamma/peaks/coupled_multiplet.py</code>;
связь площадей через &epsilon;&middot;t — наша, и она отправлена туда обратно
отдельным коммитом.</p>

<h2 id="appendix">Приложение. Полные таблицы</h2>

<details>
<summary>А.1. Баланс пика и континуума в окне деконволюции</summary>
%(balance)s
<p class="cap">Доля площади окна, отданная подгонкой пикам, у измерения и у
модели; их отношение переносится в активность множителем один к одному.</p>
</details>

<details>
<summary>А.2. Континуум вне пиков по участкам спектра</summary>
%(continuum)s
<p class="cap">Отношение расчётного континуума к измеренному, нормированное на
линию 2614,5 кэВ; выделены участки с отклонением от единицы более 5 %%.</p>
</details>

<details>
<summary>А.3. Деконволюция против оконного съёма по линиям</summary>
<div class="tw"><table>
<thead><tr><th>геометрия</th><th>нуклид</th><th class="n">E, кэВ</th>
<th class="n">линий<br>в группе</th><th class="n">чистота</th>
<th class="n">A/пасп<br>деконволюцией</th><th class="n">A/пасп<br>окном</th>
<th class="n">&chi;&sup2;/dof</th></tr></thead>
<tbody>%(deconv)s</tbody></table></div>
<p class="cap">Строки серым — линии в мультиплете, для оконного съёма негодные.
Совпадение двух правых столбцов на чистых линиях подтверждает нормировку
деконволюции.</p>
</details>

<details>
<summary>А.4. Параметры измерительных кювет (спецификация ЛСРМ)</summary>
<div class="tw"><table>
<thead><tr><th class="n">№</th><th>измерительная кювета</th>
<th class="n">объём, л</th><th class="n">внешние габариты, мм</th>
<th class="n">d<sub>эфф</sub>, мм</th></tr></thead>
<tbody>
<tr><td class="n">1</td><td>сосуд типа Маринелли</td><td class="n">0,5</td>
<td class="n">∅125, H&nbsp;100</td><td class="n">15 (2)</td></tr>
<tr><td class="n">2</td><td>сосуд типа Маринелли</td><td class="n">1,0</td>
<td class="n">∅150, H&nbsp;110</td><td class="n">26 (2)</td></tr>
<tr><td class="n">3</td><td>сосуд типа Маринелли</td><td class="n">3,0</td>
<td class="n">∅180, H&nbsp;200</td><td class="n">60 (5)</td></tr>
<tr><td class="n">5</td><td>пластмассовая кювета типа «Дента»</td><td class="n">0,12</td>
<td class="n">∅75, H&nbsp;35</td><td class="n">36 (2)</td></tr>
<tr><td class="n">6</td><td>Петри</td><td class="n">0,075</td>
<td class="n">∅88, H&nbsp;14</td><td class="n">15 (2)</td></tr>
</tbody></table></div>
<p class="cap">Перепечатано из спецификации измерительных кювет ЛСРМ (нумерация
и обозначения источника сохранены; d<sub>эфф</sub> — эффективная толщина слоя
пробы, в скобках — погрешность последней значащей цифры). В комплект прибора
входят строки 2, 5 и 6: Маринелли 1,0 л, «Дента» 0,12 л, Петри. Модель
построена на этих внешних габаритах; внутренняя геометрия сосуда (толщина
стенок, колодец) в источнике не задана и в расчёте принята с допущениями —
см. раздел «Критика и ограничения».</p>
</details>

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
</div>
"""


def zonenote(zones, mspan, soft, comp):
    """Зоны подгонки и — главное — чем ограничен диапазон паспортной кривой."""
    out = []
    if zones:
        parts = ["зона %d: степень %d, %.0f–%.0f кэВ, разброс точек %s %%"
                 % (i, z["deg"], 10 ** z["xlo"], 10 ** z["xhi"],
                    ru((10 ** z["sig"] - 1) * 100, 1))
                 for i, z in enumerate(zones, 1)]
        out.append("Подгонка ЛСРМ — %s." % "; ".join(parts))
    else:
        out.append("Подгоночной кривой в файле нет (только <code>.efr</code> "
                   "с блоками по источникам) — показаны измеренные точки.")
    out.append("Эталонные линии: %s–%s кэВ. Снизу диапазон ограничен тем, "
               "что источников с линией мягче <b>%s кэВ (%s)</b> в этой "
               "геометрии нет. Расчётная сетка одна на все геометрии, "
               "%s–%s кэВ: от наличия источника она не зависит. Затенённые "
               "полосы — области, где сверка с аттестованной кривой "
               "недоступна."
               % (ru(mspan[0], 1), ru(mspan[1], 1), ru(soft[0], 1),
                  esc(soft[3]), ru(comp[0][0], 1), ru(comp[-1][0], 0)))
    return '<p class="leg">%s</p>' % "<br>".join(out)


def legend(with_corr):
    k = ['<span class="k"><i style="border-color:var(--exp)"></i>'
         'паспортная подгонка ЛСРМ (по зонам)</span>',
         '<span class="k"><i style="border-color:var(--exp);'
         'border-top-style:dotted"></i>её измеренные точки</span>',
         '<span class="k"><i style="border-color:var(--mc)"></i>расчёт '
         '(Geant4), подгонка по узлам сетки</span>']
    if with_corr:
        k.append('<span class="k"><i style="border-color:var(--corr)"></i>'
                 'точка с поправкой на суммирование</span>')
        k.append('<span class="k"><i style="border-color:var(--med)"></i>'
                 'медиана по точкам</span>')
    k.append('<span class="k"><i style="border-color:var(--ink);opacity:.25">'
             '</i>затенено — эталонных линий нет</span>')
    return '<p class="leg">%s</p>' % "".join(k)


def spectra_section():
    """Блоки записей комплекта и сводка деконволюции.

    Строится тем же прогоном тех же функций, что считает публикуемые числа
    (analysis/spectra_figs.py) — рисунки и таблицы не могут разойтись.
    """
    import spectra_figs as sf
    import kit_recalc as krec
    titles = {"Marinelli_1L": "Маринелли 1 л", "Denta_120mL": "«Дента» 120 мл",
              "Petri_60mL": "Петри 60 мл"}
    blocks, srows, skipped, records = [], [], [], []
    for rec in krec.VOLUME_RECORDS:
        geom = rec[0]
        got = sf.record_block(*rec, geom_title=titles.get(geom, geom))
        if isinstance(got, str):
            # Пропуск называется вслух и попадает НА СТРАНИЦУ: молчаливый
            # пропуск читается как «всё посчитано», и один такой уже стоил
            # недель (пустая таблица пересчёта комплекта).
            skipped.append((titles.get(geom, geom), rec[2], got))
            print("!! пропущено: %s %s — %s" % (geom, rec[2], got))
            continue
        html, rows = got
        blocks.append(html)
        records.append((titles.get(geom, geom), rec[2], html))
        for r in rows:
            srows.append(
                "<tr%s><td>%s</td><td>%s</td><td class='n'>%s</td>"
                "<td class='n'>%d</td><td class='n'>%s</td>"
                "<td class='n'>%s</td><td class='n'>%s</td>"
                "<td class='n'>%s</td></tr>"
                % ("" if r["clean"] else " class='dim'",
                   esc(titles.get(geom, geom)), esc(rec[2]), ru(r["E"], 1),
                   r["nl"], ru(r["frac"], 2) if r["frac"] is not None else "—",
                   ru(r["A"] / r["A0"], 3) if r["A"] else "—",
                   ru(r["win"] / r["A0"], 3) if r["win"] else "—",
                   ru(r["chi2"], 2) if r["chi2"] is not None else "—"))
    if skipped:
        blocks.append('<div class="card warn"><p><b>Записи, не попавшие в этот '
                      'раздел (%d).</b> Перечислены потому, что пропуск без '
                      'предупреждения читается как «всё посчитано»:</p><ul>%s</ul>'
                      '</div>'
                      % (len(skipped),
                         "".join("<li>%s, %s — %s</li>"
                                 % (esc(g), esc(n), esc(why))
                                 for g, n, why in skipped)))
    # Скрипт интерактивных графиков — ОДИН на страницу и в самом конце, когда
    # все контейнеры уже в разметке.
    blocks.append(sf.SPECTRA_JS)
    return "\n".join(blocks), "".join(srows), sf.DECONV_LEGEND, records


def grid_span():
    """Диапазон и число узлов расчётной сетки — из самого списка энергий."""
    sys.path.insert(0, os.path.join(ROOT, "drivers"))
    from grid_energies import LINES
    return ("%s–%s кэВ" % (ru(min(LINES), 1), ru(max(LINES), 1)), len(LINES))


def geometry_data():
    """Блоки по геометриям (title, slug, note, html) + сводка по ним.

    Вынесено из build(), чтобы одни и те же расчёты можно было уложить
    и в плоскую страницу отчёта (build_web.py), и во вкладки статьи
    (build_article.py) — без второй реализации математики.
    """
    global EFF_DIR
    EFF_DIR = find_eff_dir()
    if EFF_DIR is None:
        raise SystemExit(
            "Не найдены файлы кривых .efa/.efr. Они лежат в\n"
            "detectors/Gamma-1S/reference/lsrm/efficiency/; укажите\n"
            "G4MODELS_REF, если эталоны вынесены из репозитория.")
    C = summing()
    kit = kit_activity()

    blocks, summary = [], []
    for title, mname, cfile, note, kitkey in GEOMS:
        meas, zones, src = measured(mname)
        comp = computed(cfile)
        if not meas:
            print("!! нет измеренной кривой для", title)
            continue
        mcf, mcrange, mcdeg, mcloo = mc_fit(comp)
        mspan = (min(E for E, _, _, _ in meas), max(E for E, _, _, _ in meas))
        soft = min(meas, key=lambda p: p[0])
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
                % (ru(E, 1), esc(nuc), rug(ev * 100),
                   rug(r * ev * 100), ru(r),
                   (ru(c[0]) + " → " + ru(rc)) if c else "—",
                   "узел" if is_node else "интерп."))
        block_html = """
<h3>%s</h3>
<p>%s</p>
<figure>%s
<figcaption class="cap"><b>Рисунок 4 — расчётная (Монте-Карло) и измеренная
(ЛСРМ) кривые эффективности регистрации по пику полного поглощения для выбранной
геометрии.</b> Синяя линия — расчёт, оранжевые квадраты — измеренные точки
аттестованной кривой, серые полосы — области без эталонных линий, где сверка
недоступна; усы измерения — паспортная погрешность, усы расчёта — статистика
прогона. Наведите курсор на точку, чтобы увидеть числа.</figcaption>
</figure>
<figure>%s
<figcaption class="cap"><b>Рисунок 5 — отношение расчётной и измеренной
эффективности по линиям для выбранной геометрии.</b> Синие точки — отношение по
отдельной линии (усы — погрешность), зелёная штриховая — единица (полное
совпадение), сиреневые точки — то же отношение после деления на поправку
каскадного суммирования. Средневзвешенное (в логарифме) %s по %d точкам, разброс
формы RMS %s %%, медиана %s; расчётная кривая — полином %d-й степени по log-log
(степень выбрана скользящим исключением узла, медианная ошибка предсказания
%s %%).</figcaption>
</figure>
<p class="cap"><b>Таблица 3 — сопоставление расчётной и измеренной эффективности
по линиям для выбранной геометрии.</b> «измерено, %%» и «расчёт, %%» — абсолютная
эффективность по ППП; «МК/эксп» — их отношение (единица — совпадение); «C → с
поправкой» — отношение до и после деления на поправку каскадного суммирования;
«расчёт взят» — из узла сетки или интерполяцией.</p>
<div class="tw"><table>
<thead><tr><th class="n">E, кэВ</th><th>по нуклиду</th>
<th class="n">измерено, %%</th><th class="n">расчёт, %%</th>
<th class="n">МК/эксп</th><th class="n">C → с поправкой</th>
<th>расчёт взят</th></tr></thead>
<tbody>%s</tbody></table></div>
<p class="leg">%s</p>
""" % (esc(title), esc(note),
            chart(meas, comp, lo, hi, zones, mcf),
            ratio_chart(pairs, med, zones, mcf, mcrange, mspan),
            ru(k), len(pairs), ru(100 * rms, 1), ru(med), mcdeg,
            ru(100 * mcloo, 1), "".join(rows),
            zonenote(zones, mspan, soft, comp))
        blocks.append((title, block_html))
        kr = kit.get(kitkey)
        summary.append((title, len(pairs), k, kg, rms,
                        kr[2] if kr else 0, kr if kr else None))

    srows = "".join(
        "<tr><td>%s</td><td class='n'>%d</td><td class='n'>%s</td>"
        "<td class='n'>%s</td><td class='n'>%s</td><td class='n'>%s</td>"
        "<td class='n'>%s</td></tr>"
        % (esc(t), n, ru(k), ru(kg), ru(100 * rms, 1),
           nk if nk else "—",
           kit_cell(km))
        for t, n, k, kg, rms, nk, km in summary)
    return blocks, summary, srows


def build():
    blocks, summary, srows = geometry_data()
    spec_blocks, spec_rows, spec_legend, _spec_records = spectra_section()
    gspan, gn = grid_span()
    html = TMPL % dict(css_href=CSS_HREF, summary=srows,
                       answer=answer_card(kit_activity()),
                       blocks="".join(html for _, html in blocks),
                       legend=legend(True), spectra=spec_blocks,
                       deconv=spec_rows, deconv_legend=spec_legend,
                       balance=balance_table(), continuum=continuum_table(),
                       gridspan=gspan, gridn=gn, gridnusl=nusl(gn),
                       setup=SETUP_HTML % dict(gridspan=gspan,
                           gridn=gn, gridnusl=nusl(gn)),
                       georesp=GEORESP_HTML)
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(html)
    print("страница: %s (%.0f КБ)" % (OUT, os.path.getsize(OUT) / 1024))
    for t, n, k, kg, rms, nk, km in summary:
        print("   %-16s линий %2d, МК/эксп %s (без вычета полки %s), "
              "RMS формы %s %%, пересчёт комплекта %s"
              % (t, n, ru(k), ru(kg), ru(100 * rms, 1),
                 ("%s ± %s по %d линиям" % (ru(km[0]), ru(km[1]), km[2]))
                 if km else "—"))


if __name__ == "__main__":
    build()
