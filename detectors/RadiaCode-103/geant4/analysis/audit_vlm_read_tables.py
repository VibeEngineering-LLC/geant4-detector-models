# -*- coding: utf-8 -*-
"""СТЕРИЛЬНЫЙ проход зрячей моделью по растровым таблицам статей.

Назначение: независимое прочтение таблицы с картинки локальной VLM, чтобы сверить
её с прочтением агента. Промпт НЕ содержит ожидаемых значений (#SA-7): модель
получает задачу «перепиши таблицу», а не вопрос «правильно ли там 556».

Использование:
    python audit_vlm_read_tables.py <картинка.png> [<картинка2.png> ...] [--model M]

Выводит сырой ответ модели по каждой картинке. Сравнение с нашими числами —
ОТДЕЛЬНЫЙ шаг, механический, не здесь.
"""
import sys
import base64
import json
import argparse
import urllib.request

sys.stdout.reconfigure(encoding="utf-8")

PROMPT = (
    "Ниже изображение таблицы из научной статьи. Перепиши её содержимое как есть, "
    "строка за строкой, включая заголовки колонок и единицы измерения. "
    "Числа переписывай ТОЧНО как напечатано, ничего не округляй, ничего не пересчитывай. "
    "Если значение не читается, напиши '?' вместо него. Не добавляй пояснений и выводов."
)


def ask(model, image_path, prompt, host="http://127.0.0.1:11434"):
    with open(image_path, "rb") as f:
        img = base64.b64encode(f.read()).decode("ascii")
    payload = {
        "model": model,
        "prompt": prompt,
        "images": [img],
        "stream": False,
        "options": {"temperature": 0.0, "num_ctx": 8192},
    }
    req = urllib.request.Request(
        host + "/api/generate",
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=900) as r:
        data = json.loads(r.read().decode("utf-8"))
    # Ярус 0 по #EVAL-1: механика проверяет причину остановки и непустоту ответа.
    reason = data.get("done_reason")
    text = (data.get("response") or "").strip()
    if reason != "stop":
        raise RuntimeError(f"done_reason={reason!r}, ответ не завершён штатно")
    if not text:
        raise RuntimeError("пустой ответ модели")
    return text


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("images", nargs="+")
    ap.add_argument("--model", default="qwen2.5vl:7b")
    ap.add_argument("--prompt", default=PROMPT)
    a = ap.parse_args(argv)
    rc = 0
    for img in a.images:
        print(f"\n===== {img} | {a.model} =====")
        try:
            print(ask(a.model, img, a.prompt))
        except Exception as exc:
            print(f"ОТКАЗ: {exc!r}")
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))