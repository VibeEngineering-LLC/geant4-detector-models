# -*- coding: utf-8 -*-
"""#GS-78: правила вычитки русского текста GS2020 (patches/ru_rules.json): js() - донорский JS, data() - строки данных."""
import io, json, os, re
_R = json.load(io.open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "patches", "ru_rules.json"), encoding="utf-8"))
_D = [(re.compile(p), r) for p, r in _R["data"]]

def js(text):
    text = re.sub(r"\\u([0-9a-fA-F]{4})", lambda m: chr(int(m.group(1), 16)) if int(m.group(1), 16) >= 0x80 else m.group(0), text)   # часть строк донорского JS записана \uXXXX
    for old, new in _R["js"]:
        text = text.replace(old, new)
    return text

_TH = re.compile(r"(?<![\d,.])(\d{4,})(?=(?:,\d+)?[  ]+(?:г|кг|с|Бк)\b)")

def fill(f):
    """#GS-78 r6: значения меток {{...}}: целые от 4 цифр перед г/кг/с/Бк — с неразрывным пробелом тысяч (1 085 г)."""
    w = lambda g: (lambda: _TH.sub(lambda m: format(int(m.group(1)), ",").replace(",", " "), g()))
    return {k: (w(v) if callable(v) else v) for k, v in f.items()}

def data(text):
    for rx, rep in _D:
        text = rx.sub(rep, text)
    return text

def clean_json(path):
    t = io.open(path, encoding="utf-8", newline="").read()
    n = data(t)
    if n != t:
        json.loads(n)
        io.open(path, "w", encoding="utf-8", newline="").write(n)
    return abs(len(n) - len(t))
