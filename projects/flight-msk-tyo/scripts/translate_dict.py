"""Словарный перевод строк реакций/материалов RU->EN. Используется make_page_data_en.py."""
import json, re, os

_dir = os.path.dirname(os.path.abspath(__file__))
REACTIONS = json.load(open(os.path.join(_dir, "ru_en_reactions.json"), encoding="utf-8"))
TERMS = json.load(open(os.path.join(_dir, "ru_en_terms.json"), encoding="utf-8"))
_TERMS_SORTED = sorted(TERMS.items(), key=lambda kv: -len(kv[0]))

def has_cyrillic(s):
    return bool(re.search(r"[А-Яа-яЁё]", s))

def translate_reaction(s):   # непокрытая строка с кириллицей -> None (сигнал доделать словарь)
    key = s.replace("**", "")   # словарь построен без markdown-разметки (упрощение — жирный не переносим)
    if key in REACTIONS:
        return REACTIONS[key]
    return key if not has_cyrillic(key) else None

def translate_place(s):
    out = s.replace("**", "")
    for ru, en in _TERMS_SORTED:
        out = out.replace(ru, en)
    return out if not has_cyrillic(out) else None

def translate_yield(s):   # СНАЧАЛА терм-подстановка (ключи словаря — с запятой), ПОТОМ запятая->точка
    out = s.replace("**", "")
    for ru, en in _TERMS_SORTED:
        out = out.replace(ru, en)
    out = re.sub(r"(?<=\d),(?=\d)", ".", out)
    return out if not has_cyrillic(out) else None
