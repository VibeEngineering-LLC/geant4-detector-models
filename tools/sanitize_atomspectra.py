# -*- coding: utf-8 -*-
"""Подготовка XML-спектра AtomSpectra к публикации в открытом репозитории.

`tools/anonymize.py` рассчитан на фамилии и заводские номера в приборных
форматах ЛСРМ. В XML AtomSpectra опознаётся другое, и его надо убирать
прицельно, не трогая ни одного отсчёта:

  * `DeviceConfigReference/Guid` и `ROIConfigReference/Guid` — идентификаторы
    ПРОФИЛЕЙ в установке программы у конкретного оператора. Сами по себе это
    случайные UUID, но они одинаковы во всех его файлах и связывают их между
    собой и с ним;
  * `DeviceConfigReference/Name` — имя профиля, набранное человеком; здесь в нём
    оказалось лишнее («1.» и название другого прибора из соседнего профиля);
  * `BackgroundSpectrumFile` — имя файла фона; здесь «Фон дом …», то есть место
    измерения.

Что НЕ трогается: отсчёты, калибровки, времена, число каналов, параметры формы
линии — то есть всё, на чём стоит расчёт. Скрипт проверяет это сам: после замен
сверяет, что сумма отсчётов и все числовые узлы спектра совпадают с исходником.

    python tools/sanitize_atomspectra.py <исходник.xml> <выход.xml>
"""
import io
import os
import re
import sys

# что и на что меняется; значения нарочно неинформативные и постоянные,
# чтобы файлы разных дней не связывались между собой
GUID_STUB = "00000000-0000-0000-0000-000000000000"
DEVICE_NAME = "AtomSpectra Nano 16 PRO"
BG_NAME = "background.xml"


def sanitize(text):
    changes = []

    def sub_tag(pattern, repl, label, t):
        new, n = re.subn(pattern, repl, t, flags=re.S)
        if n:
            changes.append("%s: %d замен" % (label, n))
        return new

    t = text
    t = sub_tag(r"(<DeviceConfigReference>.*?<Name>)(.*?)(</Name>)",
                lambda m: m.group(1) + DEVICE_NAME + m.group(3),
                "имя профиля прибора", t)
    t = sub_tag(r"(<DeviceConfigReference>.*?<Guid>)(.*?)(</Guid>)",
                lambda m: m.group(1) + GUID_STUB + m.group(3),
                "GUID профиля прибора", t)
    t = sub_tag(r"(<ROIConfigReference>.*?<Guid>)(.*?)(</Guid>)",
                lambda m: m.group(1) + GUID_STUB + m.group(3),
                "GUID профиля областей интереса", t)
    t = sub_tag(r"(<BackgroundSpectrumFile>)(.*?)(</BackgroundSpectrumFile>)",
                lambda m: m.group(1) + BG_NAME + m.group(3),
                "имя файла фона", t)
    return t, changes


def numeric_fingerprint(text):
    """Все числа внутри спектральных узлов — для сверки «данные не тронуты»."""
    body = re.sub(r"<(DeviceConfigReference|ROIConfigReference|SampleInfo|"
                  r"BackgroundSpectrumFile)>.*?</\1>", "", text, flags=re.S)
    body = re.sub(r"<BackgroundSpectrumFile>.*?</BackgroundSpectrumFile>", "",
                  body, flags=re.S)
    return re.findall(r"-?\d+\.?\d*(?:[eE][+-]?\d+)?", body)


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    src, dst = sys.argv[1], sys.argv[2]
    raw = io.open(src, encoding="utf-8-sig").read()
    out, changes = sanitize(raw)

    a, b = numeric_fingerprint(raw), numeric_fingerprint(out)
    if a != b:
        raise SystemExit("ОТКАЗ: изменились числовые данные (%d против %d узлов)"
                         % (len(a), len(b)))

    left = []
    for pat, why in [(r"[Дд]ом\b", "упоминание дома"),
                     (r"RadiaScan", "имя другого прибора"),
                     (r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-(?!0{12})[0-9a-f]{12}",
                      "непустой GUID")]:
        if re.search(pat, out):
            left.append(why)
    if left:
        raise SystemExit("ОТКАЗ: осталось после чистки — " + "; ".join(left))

    os.makedirs(os.path.dirname(os.path.abspath(dst)), exist_ok=True)
    io.open(dst, "w", encoding="utf-8", newline="").write(out)
    print("исходник: %s (%d КБ)" % (src, len(raw) // 1024))
    for c in changes:
        print("  " + c)
    print("числовые данные совпали: %d узлов" % len(a))
    print("записано:", dst)


if __name__ == "__main__":
    main()
