import sys
import re
import json
from pathlib import Path
from argparse import ArgumentParser

sys.stdout.reconfigure(encoding="utf-8")

# Модульный уровень: маркеры производных чисел на русском
DERIVED_MARKERS = (
    "расчёт", "расчет", "вычисл", "оцен", "наш", "мы ", "получ", "размах",
    "отношение", "во сколько", "×", "раз", "%", "процент", "производн",
    "следств", "гипотез", "порядок", "приблиз", "≈", "~"
)

# Модульный уровень: шумовые токены для пропуска
NOISE = {
    "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "100",
    "2009", "2011", "2012", "2014", "2015", "2022", "2023", "2026"
}

def extract_numbers(text: str) -> set[str]:
    """Извлекает числа из текста и возвращает их нормализованные строки."""
    text = text.replace("−", "-")
    matches = re.findall(r"\d+(?:[.,]\d+)?", text)
    numbers = set()
    for match in matches:
        norm = match.replace(",", ".")
        numbers.add(norm)
    return numbers

def variants(num: str) -> set[str]:
    """Возвращает множество вариантов написания числа."""
    result = {num}
    if "." in num:
        # Заменить точку на запятую
        result.add(num.replace(".", ","))
        # Убрать ведущие нули и точку, если есть
        stripped = num.rstrip("0").rstrip(".")
        if stripped:
            result.add(stripped)
        # Целая часть
        integer_part = num.split(".")[0]
        if len(integer_part) > 1:
            result.add(f"{integer_part[0]}.{integer_part[1:]}")
            result.add(f"{integer_part[0]},{integer_part[1:]}")
    # Удалить точку и проверить
    no_dot = num.replace(".", "")
    if len(no_dot) >= 2:
        result.add(no_dot)
    return result

def read_pdf_text(path: Path) -> str:
    """Читает текст из PDF файла.

    Нормализация десятичного разделителя. В математическом шрифте этого класса
    статей десятичная точка извлекается как ДВОЕТОЧИЕ: «5:96» вместо «5.96»
    (проверено grep-ом по извлечённому тексту Nat Commun 14:7790 — «5.96»
    встречается 0 раз, «5:96» 1 раз). Без этой замены инструмент систематически
    не находит дробные числа из формул и таблиц и выдаёт ложную тревогу.
    Побочный эффект: пострадает запись времени вида «12:30», в научном тексте
    это приемлемо.
    """
    import fitz
    doc = fitz.open(path)
    text = ""
    for page in doc:
        text += page.get_text() + "\n"
    doc.close()
    return re.sub(r"(?<=\d):(?=\d)", ".", text)

def main():
    parser = ArgumentParser()
    parser.add_argument("--doc", type=Path, required=True, help="Путь к Markdown файлу")
    parser.add_argument("--src", action="append", type=Path, required=True, help="Путь к PDF файлу (повторяется)")
    parser.add_argument("--json", type=Path, help="Путь для записи JSON отчета")

    args = parser.parse_args()

    # Читаем документ
    with open(args.doc, encoding="utf-8") as f:
        doc_lines = [line.rstrip() for line in f.readlines()]

    # Собираем текст из всех PDF
    full_text = ""
    for src_path in args.src:
        full_text += read_pdf_text(src_path) + "\n"

    # Извлекаем числа из источников
    source_numbers = extract_numbers(full_text)

    found = []
    derived = []
    orphan = []

    for line_num, line in enumerate(doc_lines, 1):
        if not line.strip() or line.startswith(("```", "|---", "---")):
            continue

        is_derived_line = any(marker in line.lower() for marker in DERIVED_MARKERS)
        numbers_in_line = extract_numbers(line)

        for num in sorted(numbers_in_line):
            # Пропускаем шум
            if num in NOISE or len(num.replace(".", "")) < 2:
                continue

            record = {
                "line": line_num,
                "num": num,
                "text": line.strip()[:110]
            }

            if variants(num) & source_numbers:
                found.append(record)
            elif is_derived_line:
                derived.append(record)
            else:
                orphan.append(record)

    # Выводим отчет
    print("=" * 78)
    print(f"ДОКУМЕНТ: {args.doc}")
    for src_path in args.src:
        print(f"ИСТОЧНИК: {src_path.name}")
    print("=" * 78)
    print(f"чисел найдено в источнике дословно     : {len(found)}")
    print(f"НЕТ В ИСТОЧНИКЕ, строка помечена       : {len(derived)}")
    print(f"НЕТ В ИСТОЧНИКЕ, без метки             : {len(orphan)}")
    print(f"ВСЕГО БЕЗ ОПОРЫ В ИСТОЧНИКЕ            : {len(derived) + len(orphan)}")
    print()

    # #SA-3: метка «производное» НЕ выводит число из-под проверки, она лишь
    # смягчает категорию. Обе корзины предъявляются, иначе подменённое число
    # прячется за первым попавшимся символом «≈» или «%» в строке.
    if orphan or derived:
        print("-" * 78)
        print("ТРЕБУЕТ ТОЛКОВАНИЯ — числа, которых НЕТ в тексте источника.")
        print("Каждое обязано получить письменный ответ: чужой источник /")
        print("вычислено из чего / наша величина / ошибка.")
        print("-" * 78)
        for rec in orphan:
            print(f"  БЕЗ МЕТКИ  стр.{rec['line']:>4}  [{rec['num']}]  {rec['text']}")
        for rec in derived:
            print(f"  помечено   стр.{rec['line']:>4}  [{rec['num']}]  {rec['text']}")
    else:
        print("(все числа документа найдены в источнике)")

    # JSON отчет
    if args.json:
        report = {
            "found": len(found),
            "derived": derived,
            "orphan": orphan
        }
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(report, f, ensure_ascii=False, indent=1)
        print(f"\nмашинный отчёт: {args.json}")

    return 1 if (orphan or derived) else 0

if __name__ == "__main__":
    sys.exit(main())
