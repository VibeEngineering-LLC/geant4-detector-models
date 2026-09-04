# -*- coding: utf-8 -*-
"""Стерильная вычитка документа локальной моделью (ступень 2 лестницы IRON MODE).

Модель получает ТОЛЬКО текст части документа и критерии проверки. Ни выводов
автора, ни ожидаемых ответов в промпте нет — иначе это не проверка, а
подтверждение названного ответа.
"""
import os
import sys
import json
import pathlib
import argparse

sys.stdout.reconfigure(encoding="utf-8")
# Путь к helper'у guarded_generate. Профиль пользователя НЕ зашивается: имя
# профиля не должно попадать в публикуемый репозиторий, и на другой машине
# жёсткий путь всё равно не работает. Порядок поиска: явная переменная
# окружения -> домашний каталог текущего пользователя.
_helper_dir = os.environ.get("CLAUDE_WORKFLOW_SCRIPTS")
if not _helper_dir:
    _helper_dir = str(pathlib.Path.home() / ".claude" / "skills" / "workflow" / "scripts")
if not pathlib.Path(_helper_dir, "vram_guard_reference.py").exists():
    print(f"Не найден vram_guard_reference.py в: {_helper_dir}\n"
          f"Укажите каталог явно: set CLAUDE_WORKFLOW_SCRIPTS=<путь>")
    sys.exit(1)
sys.path.insert(0, _helper_dir)
from vram_guard_reference import guarded_generate

CRIT = """Ты — придирчивый научный рецензент. Перед тобой ФРАГМЕНТ справочника по
радиоактивным рядам урана и тория. Найди в нём дефекты. Отвечай по-русски.

Ищи строго это:
1. ВНУТРЕННИЕ ПРОТИВОРЕЧИЯ: два места фрагмента утверждают несовместимое.
2. ФИЗИЧЕСКИЕ ОШИБКИ: неверная формула, неверный порядок величины, перепутанные
   единицы, невозможное соотношение активностей.
3. УТВЕРЖДЕНИЯ БЕЗ ИСТОЧНИКА: число или факт подан как установленный, но ссылки
   на источник рядом нет.
4. ПОДМЕНА: расчётное или номинальное значение подано как измеренное.
5. НЕОДНОЗНАЧНОСТЬ: формулировка допускает два прочтения, из которых одно неверно.

НЕ придирайся к стилю, оформлению, длине и порядку разделов.
Если дефектов нет — верни пустой список. Не выдумывай дефекты ради заполнения ответа.

Формат ответа — JSON: {"findings":[{"severity":"high","quote":"точная цитата","problem":"в чём дефект"}]}

ФРАГМЕНТ:
"""


def main():
    parser = argparse.ArgumentParser(description="Стерильная вычитка markdown-документа")
    parser.add_argument("doc", type=pathlib.Path, help="Путь к markdown-файлу")
    parser.add_argument("--model", default="qwen3.6:27b", help="Имя модели Ollama")
    parser.add_argument("--chunk", type=int, default=14000, help="Максимальная длина части в символах")
    parser.add_argument("--crit", type=pathlib.Path, default=None,
                        help="Файл с критериями проверки. Без него берутся встроенные "
                             "(заточены под справочник по рядам U/Th). ВНИМАНИЕ: критерии "
                             "не должны содержать ожидаемых ответов — иначе это не "
                             "проверка, а подтверждение (#SA-7).")
    args = parser.parse_args()

    global CRIT
    if args.crit:
        if not args.crit.exists():
            print(f"Файл критериев не найден: {args.crit}")
            sys.exit(1)
        CRIT = args.crit.read_text(encoding="utf-8")

    if not args.doc.exists():
        print(f"Файл не найден: {args.doc}")
        sys.exit(1)

    text = args.doc.read_text(encoding="utf-8")

    # Разрезание на части по границам разделов
    parts = []
    current_part_lines = []
    current_length = 0

    for line in text.splitlines():
        # Если строка начинается с "## ", это начало нового раздела
        if line.startswith("## "):
            # Если текущая часть уже превышает лимит, закрываем её
            if current_length > args.chunk and current_part_lines:
                parts.append("\n".join(current_part_lines))
                current_part_lines = []
                current_length = 0
        
        current_part_lines.append(line)
        current_length += len(line) + 1  # +1 для символа новой строки

    # Добавляем последнюю часть, если она не пуста
    if current_part_lines:
        parts.append("\n".join(current_part_lines))

    print(f"частей: {len(parts)}, модель: {args.model}")

    for i, part in enumerate(parts):
        print("=" * 90)
        print(f"ЧАСТЬ {i + 1}/{len(parts)} ({len(part)} символов)")
        print("=" * 90)

        try:
            prompt = CRIT + part
            result = guarded_generate(
                model=args.model,
                prompt=prompt,
                fmt="json",
                want_gpu=True,
                priority=50,
                max_wait_s=900,
                temperature=0,
                num_ctx=32768,
                think=False,
                extra_options={"num_predict": 4000}
            )

            # Извлечение ответа
            if isinstance(result, dict):
                response = result.get("response", "")
            else:
                response = str(result)

            if not response:
                print("(модель ничего не вернула)")
                continue

            try:
                data = json.loads(response)
                findings = data.get("findings", [])
                if not findings:
                    print("(дефектов не найдено)")
                else:
                    for finding in findings:
                        severity = finding.get("severity", "unknown")
                        problem = finding.get("problem", "")
                        quote = finding.get("quote", "")
                        # Обрезаем цитату до 300 символов
                        if len(quote) > 300:
                            quote = quote[:300] + "..."
                        print(f"[{severity}] {problem}")
                        print(f"    цитата: {quote}")
            except json.JSONDecodeError:
                print("(не JSON, сырой ответ)")
                if len(response) > 3000:
                    print(response[:3000] + "...")
                else:
                    print(response)

        except Exception as e:
            print(f"ОТКАЗ Ollama: {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
