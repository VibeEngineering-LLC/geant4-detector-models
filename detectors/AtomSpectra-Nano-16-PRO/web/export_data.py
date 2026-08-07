# -*- coding: utf-8 -*-
"""Выгрузка расчёта Geant4 в JSON для интерактивной страницы.

Свёртка с приборным разрешением, ширина канала отображения и нормировка на
квант в 4pi скопированы ДОСЛОВНО из штатных иллюстраций проекта
(`analysis/draw_channels.py`, `analysis/response_matrix.py`), чтобы страница
и рисунки не могли разойтись.

Проверки при выгрузке — обычные `raise`, а не `assert`: под `python -O`
инструкция `assert` вырезается, и сторож молча исчезает.

    python export_data.py [<каталог спектров>] [<каталог вывода>]
"""
import io
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_CHAN = r"C:/g4work/asn16/build/spectra_resp50"
SRC_RES = (r"C:/g4work/geant4-detector-models/detectors/"
           r"AtomSpectra-Nano-16-PRO/results")
RUN_LOG = r"C:/g4work/asn16/build/run_resp.log"

FWHM_662 = 41.60
E_MAX = 3200.0
MEC2 = 510.99895          # кэВ, энергия покоя электрона
WANTED = [180.0, 1480.0, 3000.0]
PEAK_HALF = 1.0           # полуширина окна пика полного поглощения, кэВ

# Порядок и подписи — по физическому смыслу. Цвета сюда НЕ входят: они
# заданы в таблице стилей отдельно для светлой и тёмной темы (контраст
# заливки к фону не может быть одинаковым на белом и на почти чёрном).
# Каждая подпись — пара (ru, en); переключается на странице кнопкой.
ORDER = [
    ("photo",      "фотоэффект, ничего не вылетело",
                   "photoelectric, nothing escaped"),
    ("compt_full", "комптон и поглощение, без вылета",
                   "Compton + absorption, no escape"),
    ("pair_full",  "пары, оба 511 поглощены",
                   "pair, both 511 keV absorbed"),
    ("compt_esc1", "однократный комптон, квант ушёл",
                   "single Compton, γ escaped"),
    ("compt_escN", "многократный комптон, квант ушёл",
                   "multiple Compton, γ escaped"),
    ("pair_esc1",  "пары, вылетел один 511",
                   "pair, one 511 keV escaped"),
    ("pair_esc2",  "пары, вылетели оба",
                   "pair, both 511 keV escaped"),
    ("brems_esc",  "вылет тормозного",
                   "bremsstrahlung escape"),
    ("xray_esc",   "вылет характеристического рентгена",
                   "K X-ray escape"),
    ("external",   "вторичные из корпуса и обёртки",
                   "secondaries from housing and wrapping"),
    ("other",      "остаточный канал (сторож)",
                   "residual channel (guard)"),
]


class Bad(SystemExit):
    """Отказ выгрузки. Отдельный класс, чтобы сборка не путала его с ошибкой."""


def need(cond, msg):
    if not cond:
        raise Bad("выгрузка отклонена: " + msg)


def g4(x):
    """4 значащие цифры; ноль остаётся нулём."""
    return 0.0 if x == 0 else float("%.4g" % x)


def fwhm(e):
    return FWHM_662 * math.sqrt(max(e, 1.0) / 661.657)


def step_for(e0):
    """Ширина канала отображения: четыре канала на полуширину."""
    return max(5.0, round(fwhm(e0) / 4.0 / 5.0) * 5.0)


def read_chan(path):
    head, names, rows = {}, None, []
    for ln in io.open(path, encoding="utf-8"):
        ln = ln.rstrip("\n")
        if ln.startswith("#"):
            if "=" in ln:
                k, v = ln.lstrip("# ").split("=", 1)
                head[k.strip()] = v.strip()
            continue
        p = ln.split(",")
        if names is None:
            names = p[1:]
            continue
        rows.append((float(p[0]), [float(x) for x in p[1:]]))
    return head, names, rows


def read_spectrum(path):
    head, rows = {}, []
    for ln in io.open(path, encoding="utf-8"):
        ln = ln.rstrip("\n")
        if ln.startswith("#"):
            if "=" in ln:
                k, v = ln.lstrip("# ").split("=", 1)
                head[k.strip()] = v.strip()
            elif "CsI" in ln:
                # Первая строка шапки: «# ATOMSPECTRA …, CsI(Tl) 18x15x60 mm».
                # Пишется прогоном ИЗ ГЕОМЕТРИИ, поэтому переживает правку
                # размеров кристалла — в отличие от строки в шаблоне.
                m = re.search(r"(CsI\(Tl\)\s*[\d\s x×]*mm)", ln)
                need(m, "первая строка шапки %s не разбирается: %r"
                     % (os.path.basename(path), ln))
                s = re.sub(r"(\d)\s*x\s*(?=\d)", r"\1 × ", m.group(1))
                head["__crystal__"] = s.replace(" mm", " мм")
            continue
        if not ln or ln.startswith("E_keV"):
            continue
        e, c = ln.split(",")
        rows.append((float(e), float(c)))
    return head, rows


def broaden(pairs, w, nch, step):
    """Свёртка с приборным разрешением; вес w переводит отсчёты в 4pi."""
    out = [0.0] * nch
    for e, c in pairs:
        if c <= 0:
            continue
        s = fwhm(e) / 2.3548
        lo, hi = max(0.0, e - 4 * s), min(E_MAX, e + 4 * s)
        acc, norm = [], 0.0
        for k in range(int(lo / step), min(int(hi / step) + 1, nch)):
            g = math.exp(-0.5 * (((k + 0.5) * step - e) / s) ** 2)
            acc.append((k, g))
            norm += g
        if norm <= 0:
            out[min(int(e / step), nch - 1)] += c * w
            continue
        for k, g in acc:
            out[k] += c * w * g / norm
    return out


def geant4_version(log_path):
    """Версия из журнала ИМЕННО ЭТОГО прогона.

    Заголовок установки на диске тоже её содержит, но с прогоном никак не
    связан: установок может быть несколько, а журнал печатает ту, что
    отработала.
    """
    need(os.path.exists(log_path), "нет журнала прогона %s" % log_path)
    txt = io.open(log_path, encoding="utf-8-sig", errors="replace").read(4000)
    m = re.search(r"Geant4 version Name:\s*geant4-(\d+)-(\d+)-patch-(\d+)", txt)
    if m:
        return "%d.%d.%d" % (int(m.group(1)), int(m.group(2)), int(m.group(3)))
    m = re.search(r"Geant4 version Name:\s*geant4-(\d+)-(\d+)\b", txt)
    need(m, "в журнале %s не нашлась строка версии Geant4" % log_path)
    return "%d.%d" % (int(m.group(1)), int(m.group(2)))


def export_tabs(src):
    files = [f for f in sorted(os.listdir(src)) if f.endswith("_chan.csv")]
    need(files, "в %s нет файлов *_chan.csv" % src)
    index = {}
    for f in files:
        head, names, rows = read_chan(os.path.join(src, f))
        index[float(head["E_prim_keV"])] = (f, head, names, rows)

    tabs, stamp, git_describe, crystal, bin_kev = [], None, None, None, None
    for e_want in WANTED:
        e0 = min(index, key=lambda x: abs(x - e_want))
        need(abs(e0 - e_want) < 0.01,
             "нет узла %s, ближайший %s" % (e_want, e0))
        fname, head, names, rows = index[e0]
        stamp = head.get("src_sha1", "?")
        git_describe = head.get("git_describe", "?")
        bin_kev = float(head["bin_keV"].split()[0])
        w = float(head["solid_angle_frac"]) / float(head["N_primaries"])
        step = step_for(e0)
        nch = int(E_MAX / step)
        cols = [(k + 0.5) * step for k in range(nch)]

        shead, spec = read_spectrum(os.path.join(
            src, fname.replace("_chan.csv", ".csv")))
        crystal = shead.get("__crystal__", crystal)

        # СТОРОЖ: сумма каналов обязана совпасть с полным спектром до отсчёта.
        # Прогон это печатает, но страница собирается отдельно от прогона —
        # проверка повторяется здесь, на тех же файлах.
        sum_chan = sum(sum(v) for _, v in rows)
        sum_spec = sum(c for _, c in spec)
        need(sum_chan == sum_spec, "узел %.0f: сумма каналов %.0f, спектр %.0f"
             % (e0, sum_chan, sum_spec))
        j_other = names.index("other")
        need(sum(v[j_other] for _, v in rows) == 0,
             "узел %.0f: остаточный канал населён" % e0)

        # Доля событий, попавших В ПИК, — не то же, что доля событий без
        # вылета. Часть первичных теряет энергию в крышке и корпусе ДО входа
        # в кристалл: из кристалла ничего не вылетает, но энерговыделение
        # уже меньше полного. Поэтому обе величины считаются отдельно и
        # называются по-разному.
        in_peak = lambda e: abs(e - (e0 + 0.5)) <= PEAK_HALF
        n_peak = sum(c for e, c in spec if in_peak(e))
        n_nofly = sum(v[names.index(k)] for _, v in rows
                      for k in ("photo", "compt_full", "pair_full")
                      if k in names)
        # Состав ОКНА ПИКА по каналам. Считать долю канала в пике как его
        # долю среди событий без вылета нельзя: множества не вложены и
        # различаются в полтора раза. В окне пика на 180 кэВ есть даже
        # событие канала вылета рентгена — вылетевший квант унёс меньше
        # ширины окна.
        peak_ch = {}
        for k in names:
            j = names.index(k)
            c = sum(v[j] for e, v in rows if in_peak(e))
            if c:
                peak_ch[k] = int(c)
        need(abs(sum(peak_ch.values()) - n_peak) < 0.5,
             "узел %.0f: состав окна пика %d не сходится со счётом %d"
             % (e0, sum(peak_ch.values()), n_peak))

        # У маркера две подписи: короткая — на графике, полная — в
        # заголовке всплывающего пояснения. Полная не помещается вдоль линии
        # и налезает на кривые.
        markers = [{"e": e0, "short": {"ru": "пик", "en": "peak"},
                    "kind": "peak", "id": "peak"},
                   {"e": 2.0 * e0 * e0 / (MEC2 + 2.0 * e0),
                    "short": {"ru": "комптоновский край",
                              "en": "Compton edge"},
                    "kind": "edge", "id": "edge"}]
        feat = {"e_compton": 2.0 * e0 * e0 / (MEC2 + 2.0 * e0)}

        def argmax_of(nm):
            if nm not in names:
                return None, 0.0
            j = names.index(nm)
            pts = [(e, v[j]) for e, v in rows if v[j] > 0]
            if not pts:
                return None, 0.0
            # Ничья разрешается в пользу МЕНЬШЕЙ энергии — детерминированно и
            # не зависит от того, как отсортирован список.
            best = max(pts, key=lambda p: (p[1], -p[0]))
            return best[0], sum(c for _, c in pts)

        # Положения пиков вылета берутся ИЗ ФОРМУЛЫ (E0 − 511, E0 − 1022), а
        # данные служат проверкой: argmax сырой гистограммы обязан лечь в тот
        # же килоэлектронвольт. Обратный порядок — «маркер туда, где повыше» —
        # при десятках отсчётов ставит маркер на шум.
        for nm, dE, lab_ru, lab_en, mid in (
                ("pair_esc1", 511.0, "вылет 511",  "511 keV escape",  "esc511"),
                ("pair_esc2", 1022.0, "вылет 1022", "1022 keV escape", "esc1022")):
            e_th = e0 - dE
            if e_th <= 0:
                continue
            e_ms, n = argmax_of(nm)
            if n < 10:
                continue
            need(abs(e_ms - e_th) <= 1.0,
                 "узел %.0f: пик %s по данным %.1f, по формуле %.1f"
                 % (e0, nm, e_ms, e_th))
            markers.append({"e": e_th,
                            "short": {"ru": lab_ru, "en": lab_en},
                            "kind": "escape", "id": mid})
            feat["e_" + mid] = e_th

        # Вылет K-рентгена формулы не имеет: уходит любая линия K-серии иода
        # или цезия. Положение берётся по данным, и это оговаривается.
        e_x, n_x = argmax_of("xray_esc")
        if n_x >= 100:
            j = names.index("xray_esc")
            # Ничья — в пользу меньшей энергии, тем же правилом, что и
            # argmax_of: иначе «первый» и «второй» максимумы могли бы
            # указывать на один и тот же бин, разрешённый по-разному.
            peaks = sorted((p for p in ((e, v[j]) for e, v in rows) if p[1] > 0),
                           key=lambda p: (-p[1], p[0]))[:2]
            markers.append({"e": e_x,
                            "short": {"ru": "вылет K-рентгена",
                                      "en": "K X-ray escape"},
                            "kind": "escape", "id": "xray"})
            feat["e_xray"] = e_x
            feat["d_xray"] = g4(e0 + 0.5 - e_x)
            need(not peaks or peaks[0][0] == e_x,
                 "узел %.0f: два правила разрешения ничьи разошлись "
                 "(%.1f и %.1f)" % (e0, peaks[0][0], e_x))
            if len(peaks) > 1 and abs(peaks[1][0] - e_x) > 1.0:
                feat["e_xray2"] = peaks[1][0]
                feat["d_xray2"] = g4(e0 + 0.5 - peaks[1][0])

        # Обрыв континуума в числах: полосы одинаковой ширины под краем и
        # над ним. Ширина полосы (50 кэВ) назначена: она должна быть заметно
        # шире канала отображения и заметно уже расстояния до пика.
        if "compt_esc1" in names:
            j1 = names.index("compt_esc1")
            ec, band = feat["e_compton"], 50.0
            feat["edge_band"] = band
            feat["edge_below"] = int(sum(v[j1] for e, v in rows
                                         if ec - band <= e < ec))
            feat["edge_above"] = int(sum(v[j1] for e, v in rows
                                         if ec <= e < ec + band))

        # Зоны спектра — из формул, а не «на глаз». Комптоновский континуум
        # от нуля до кинематического края; зазор между краем и пиком
        # заполняется многократным комптоновским рассеянием; пик — окно
        # ±0,5 ПШПВ вокруг E0; выше пика тянется хвост от свёртки с
        # приборным разрешением и от каналов, где энерговыделение может
        # превысить E0 (нет таких на моноисточнике, но окно шкалы шире).
        xmax = min(E_MAX, e0 * 1.15)
        fw = fwhm(e0)
        e_ed = 2.0 * e0 * e0 / (MEC2 + 2.0 * e0)
        peak_lo = max(e_ed, e0 - fw / 2.0)
        peak_hi = min(xmax, e0 + fw / 2.0)
        zones = [{"lo": 0.0, "hi": e_ed, "id": "cont",
                  "label": {"ru": "комптоновский континуум",
                            "en": "Compton continuum"}}]
        gap_lo, gap_hi = e_ed, peak_lo
        if gap_hi - gap_lo > fw * 0.6:
            zones.append({"lo": gap_lo, "hi": gap_hi, "id": "gap",
                          "label": {"ru": "между краем и пиком",
                                    "en": "between edge and peak"}})
        zones.append({"lo": peak_lo, "hi": peak_hi, "id": "peak",
                      "label": {"ru": "пик полного поглощения",
                                "en": "full-energy peak"}})
        if xmax - peak_hi > fw * 0.4:
            zones.append({"lo": peak_hi, "hi": xmax, "id": "over",
                          "label": {"ru": "хвост выше пика",
                                    "en": "tail above peak"}})
        # Каждая зона обязана быть непустой и лежать внутри поля графика:
        # рисуется в ту же координатную сетку и обязана с ней сходиться.
        for z in zones:
            need(0 <= z["lo"] < z["hi"] <= xmax,
                 "узел %.0f: зона %s = [%.1f, %.1f] вне поля" %
                 (e0, z["id"], z["lo"], z["hi"]))

        total = broaden([(e, sum(v)) for e, v in rows], w, nch, step)
        channels = []
        for name, label_ru, label_en in ORDER:
            if name not in names:
                continue
            j = names.index(name)
            raw = sum(v[j] for _, v in rows)
            if raw == 0:
                continue
            channels.append({
                "key": name, "label": {"ru": label_ru, "en": label_en},
                "pct": g4(100.0 * raw / sum_chan), "n": int(raw),
                "ys": [g4(v) for v in broaden([(e, v[j]) for e, v in rows],
                                              w, nch, step)],
            })

        tabs.append({
            "e0": e0, "step": step,
            "n_primaries": int(head["N_primaries"]),
            "n_signal": int(sum_chan),
            "n_peak": int(n_peak),
            "peak_pct": g4(100.0 * n_peak / sum_chan),
            "peak_ch": peak_ch,
            "nofly_pct": g4(100.0 * n_nofly / sum_chan),
            "peak_half": PEAK_HALF,
            "solid_angle_frac": float(head["solid_angle_frac"]),
            "feat": feat, "markers": markers, "zones": zones, "xmax": xmax,
            "xs": cols, "total": [g4(v) for v in total],
            "channels": channels,
        })
    return tabs, stamp, git_describe, crystal, bin_kev


def export_composition(src):
    """Доли каналов на ВСЕХ узлах сетки, без спектров.

    Три вкладки — три среза; состав отклика меняется по всей шкале, и видно
    это только на всех узлах разом. Данных нужно немного: по одному числу на
    канал и узел.
    """
    es, rows_out, keys = [], [], None
    for fn in sorted(f for f in os.listdir(src) if f.endswith("_chan.csv")):
        head, names, rows = read_chan(os.path.join(src, fn))
        if keys is None:
            keys = [k for k, _, _ in ORDER if k in names]
        tot = sum(sum(v) for _, v in rows)
        if tot == 0:
            continue
        # Сторожа замыкания стоят на ВСЕХ узлах, а не на трёх показанных:
        # полосу состава рисуют все 61, и населённый остаточный канал
        # уехал бы в неё незамеченным.
        need(sum(v[names.index("other")] for _, v in rows) == 0,
             "узел %s: остаточный канал населён" % head["E_prim_keV"])
        es.append(float(head["E_prim_keV"]))
        rows_out.append([g4(sum(v[names.index(k)] for _, v in rows) / tot)
                         for k in keys])
    for e, r in zip(es, rows_out):
        # Порог по разрядности округления выгрузки (4 значащие цифры на
        # долю, одиннадцать каналов), а не «на глаз»: прежние 2e-3 были в
        # двадцать раз слабее достижимого и пропустили бы реальный перекос.
        need(abs(sum(r) - 1.0) < 1e-3,
             "узел %.0f: сумма долей %.6f" % (e, sum(r)))

    # Зоны состава: переход от преобладания фотоэффекта к комптоновскому
    # рассеянию — на пороге 0,5, канал `photo`; порог рождения пар — 1022
    # кэВ по физике (2·m_e·c²), а не подстройкой под данные.
    jp = keys.index("photo")
    soft_hi = es[0]
    for e, r in zip(es, rows_out):
        if r[jp] >= 0.5:
            soft_hi = e
    pair_lo = 2.0 * MEC2
    zones = [{"lo": es[0], "hi": soft_hi, "id": "soft",
              "label": {"ru": "фотоэффектная полоса",
                        "en": "photoelectric-dominated"}},
             {"lo": soft_hi, "hi": pair_lo, "id": "compton",
              "label": {"ru": "комптоновская полоса",
                        "en": "Compton-dominated"}},
             {"lo": pair_lo, "hi": es[-1], "id": "pair",
              "label": {"ru": "полоса рождения пар",
                        "en": "pair-production band"}}]
    for z in zones:
        need(es[0] <= z["lo"] < z["hi"] <= es[-1] + 1e-6,
             "зона состава %s = [%.1f, %.1f] вне сетки" %
             (z["id"], z["lo"], z["hi"]))
    return {"es": es, "keys": keys, "f": rows_out, "zones": zones}


def read_matrix(path):
    head, cols, rows, es = {}, None, [], []
    for ln in io.open(path, encoding="utf-8"):
        ln = ln.rstrip("\n")
        if ln.startswith("#"):
            if "=" in ln:
                k, v = ln.lstrip("#@ ").split("=", 1)
                head[k.strip()] = v.strip()
            continue
        p = ln.split(",")
        if cols is None:
            cols = [float(x) for x in p[1:]]
            continue
        es.append(float(p[0]))
        rows.append([g4(float(x)) for x in p[1:]])
    return head, cols, es, rows


def export_matrix():
    h1, c1, e1, raw = read_matrix(os.path.join(
        SRC_RES, "response_matrix_raw_10keV.csv"))
    h2, c2, e2, bro = read_matrix(os.path.join(
        SRC_RES, "response_matrix_10keV.csv"))
    need(c1 == c2 and e1 == e2, "сетки raw/broadened разошлись")
    return {"cols": c1, "es": e1, "raw": raw, "broadened": bro,
            "stamp": h2.get("src.spectra_sha1", "?")}


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else SRC_CHAN
    out = sys.argv[2] if len(sys.argv) > 2 else HERE
    tabs, stamp, git_describe, crystal, bin_kev = export_tabs(src)
    matrix = export_matrix()
    need(matrix["stamp"] == stamp, "штампы матрицы и каналов разошлись: %s / %s"
         % (matrix["stamp"], stamp))
    comp = export_composition(src)
    labels = {k: {"ru": ru, "en": en} for k, ru, en in ORDER}
    order = [k for k, _, _ in ORDER
             if any(c["key"] == k for t in tabs for c in t["channels"])]

    es = matrix["es"]
    # Сетка объявлена равномерной. Последний интервал короче (2980 -> 3000),
    # и это ЕДИНСТВЕННОЕ исключение — оно проверяется, а не замалчивается.
    steps = sorted({round(es[i + 1] - es[i], 6) for i in range(len(es) - 2)})
    need(len(steps) == 1, "сетка неравномерна: %s" % steps)
    last = round(es[-1] - es[-2], 6)

    data = {
        "stamp": stamp, "git_describe": git_describe,
        "mec2": MEC2, "fwhm_662": FWHM_662, "order": order,
        "labels": labels, "comp": comp,
        "run": {
            "geant4": geant4_version(RUN_LOG),
            "physics": "EmStandard_option4",
            "crystal": crystal,
            "bin_keV": bin_kev,
            "n_nodes": len(es), "e_lo": es[0], "e_hi": es[-1],
            "e_step": steps[0], "e_step_last": last,
            "n_primaries": tabs[0]["n_primaries"],
            "cell_keV": matrix["cols"][1] - matrix["cols"][0],
            "cone_deg": 35.0,
            "solid_angle_frac": tabs[0]["solid_angle_frac"],
        },
        "tabs": tabs, "matrix": matrix,
    }
    # Доля телесного угла в шапке прогона обязана отвечать объявленному
    # конусу: именно на неё делится счёт при переходе к 4pi.
    want = (1.0 - math.cos(math.radians(data["run"]["cone_deg"]))) / 2.0
    need(abs(want - data["run"]["solid_angle_frac"]) < 1e-6,
         "доля телесного угла %.8f не отвечает конусу %.0f° (%.8f)"
         % (data["run"]["solid_angle_frac"], data["run"]["cone_deg"], want))

    path = os.path.join(out, "asn16_data.json")
    with io.open(path, "w", encoding="utf-8") as g:
        json.dump(data, g, ensure_ascii=False, separators=(",", ":"))
    print("записано %s, %.0f КБ" % (path, os.path.getsize(path) / 1024.0))
    for t in tabs:
        print("  E0=%6.0f шаг %4.1f  сигнал %6d  в пике %5.2f %%  "
              "без вылета %5.2f %%  каналов %d"
              % (t["e0"], t["step"], t["n_signal"], t["peak_pct"],
                 t["nofly_pct"], len(t["channels"])))
    print("  матрица %d × %d, состав на %d узлах, Geant4 %s, штамп %s"
          % (len(es), len(matrix["cols"]), len(comp["es"]),
             data["run"]["geant4"], stamp))
    return 0


if __name__ == "__main__":
    sys.exit(main())
