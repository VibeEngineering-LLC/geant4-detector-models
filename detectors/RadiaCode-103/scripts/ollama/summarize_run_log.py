# -*- coding: utf-8 -*-
"""Суммаризатор логов прогонов Geant4 — первый helper делегации контура.

ЗАЧЕМ. Один прогон shieldrun пишет лог на тысячи строк, из которых 95 % —
таблицы EM-моделей и параметры физлиста, одинаковые от запуска к запуску. За
ночь 15.08.2026 таких логов набралось больше сотни, и я читал их глазами через
Read/Grep, сжигая контекст на машинерию. Предписание Censor-контура
(PROMPT_geant4-delegation-2026-08-15.md): такое уходит вниз по лестнице #SA-2.

ЛЕСТНИЦА, ПО КОТОРОЙ ЭТО НАПИСАНО:
  (а) python-скрипт — структура логов Geant4 РЕГУЛЯРНА (строки RESULT, блоки
      G4Exception, коды возврата), поэтому извлечение фактов делается здесь
      регулярками: бесплатно, детерминированно, воспроизводимо;
  (б) Ollama — только там, где нужен разбор НЕструктурированного: сводка по
      пачке логов и классификация «сошлось / не сошлось / аномалия»
      (ключ --explain), через guarded_generate(), НЕ через raw requests.post;
  (в) Claude-субагент — не нужен вовсе: суждений по ходу тут нет.

ВЫВОД — JSON на stdout, с провенансом: у каждого извлечённого поля указан файл
и номер строки, откуда оно взято. Это требование anti-hallucination стандарта:
downstream обязан иметь возможность проверить любое число по исходнику.

Запуск:
    py -3 summarize_run_log.py <лог-или-маска> [...] [--explain] [--max-lines N]
    py -3 summarize_run_log.py build/RadiaCode-103/bg_shield/pb50_nolid/*.log
"""
import argparse
import glob
import json
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# Строки, ради которых лог вообще читают. Всё остальное — шум физлиста.
PAT_RESULT = re.compile(r"^RESULT\s+(\S+)\s+(.*)$")
PAT_KV = re.compile(r"(\w+)=\s*([-\d.eE+]+|\S+)")
PAT_EXC_START = re.compile(r"G4Exception-START")
PAT_EXC_END = re.compile(r"G4Exception-END")
PAT_EXC_CODE = re.compile(r"\*\*\* G4Exception\s*:\s*(\S+)")
PAT_ISSUER = re.compile(r"issued by\s*:\s*(.+?)\s*$")
PAT_FATAL = re.compile(r"Fatal|FATAL|\*\*\* Break|aborting", re.I)
PAT_TIME = re.compile(r"(?:Elapsed|elapsed|User=|Real=)\s*([\d.]+)")


def parse_log(path):
    """-> dict с фактами одного лога. Каждое поле несёт номер строки-источника."""
    out = {
        "file": os.path.abspath(path),
        "size_bytes": os.path.getsize(path) if os.path.exists(path) else 0,
        "results": [],        # строки RESULT, разобранные в key=value
        "exceptions": [],     # G4Exception: код, кем выдано, строка
        "fatal": [],          # признаки падения
        "n_lines": 0,
    }
    if not os.path.exists(path):
        out["error"] = "файла нет"
        return out

    in_exc = False
    exc = None
    with open(path, encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f, 1):
            out["n_lines"] = i
            s = line.rstrip("\n")

            m = PAT_RESULT.match(s)
            if m:
                kind, rest = m.group(1), m.group(2)
                kv = {k: v for k, v in PAT_KV.findall(rest)}
                out["results"].append({"line": i, "kind": kind, "values": kv})
                continue

            if PAT_EXC_START.search(s):
                in_exc, exc = True, {"line": i, "code": None, "issuer": None}
                continue
            if in_exc:
                mc = PAT_EXC_CODE.search(s)
                if mc:
                    exc["code"] = mc.group(1)
                mi = PAT_ISSUER.search(s)
                if mi:
                    exc["issuer"] = mi.group(1)
                if PAT_EXC_END.search(s):
                    out["exceptions"].append(exc)
                    in_exc, exc = False, None
                continue

            if PAT_FATAL.search(s):
                out["fatal"].append({"line": i, "text": s.strip()[:200]})

    # Сводка по кодам исключений: у Geant4 одно и то же предупреждение может
    # повторяться тысячи раз — в отчёт идёт код и счётчик, не тысяча строк.
    codes = {}
    for e in out["exceptions"]:
        c = e.get("code") or "?"
        codes[c] = codes.get(c, 0) + 1
    out["exception_codes"] = codes
    out["n_exceptions"] = len(out["exceptions"])
    out["exceptions"] = out["exceptions"][:5]     # образцы, не весь список
    return out


def explain(reports, max_wait_s=600):
    """Сводка по пачке логов через Ollama. Импорт guarded_generate ОБЯЗАТЕЛЕН
    (HARD RULE скилла workflow, LOCKED 2026-06-04): raw requests.post не
    встаёт в машинную очередь, и параллельные вызовы роняют хост по VRAM."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from _vram_guard import guarded_generate

    compact = []
    for r in reports:
        compact.append({
            "file": os.path.basename(r["file"]),
            "results": [x["values"] for x in r["results"]][:3],
            "exception_codes": r.get("exception_codes", {}),
            "fatal": len(r.get("fatal", [])),
        })
    prompt = (
        "Ты разбираешь итоги прогонов Geant4. Ниже JSON с фактами по каждому "
        "логу: строки RESULT (уже разобранные), коды исключений со счётчиками, "
        "число фатальных признаков.\n\n"
        "Верни СТРОГО JSON вида {\"summary\": \"...\", \"anomalies\": [...], "
        "\"verdict\": \"ok|warn|fail\"}. В anomalies — только то, что реально "
        "выбивается: разброс результатов между файлами, повторяющиеся коды "
        "исключений, признаки падения. Ничего не выдумывай: если данных мало, "
        "так и напиши.\n\nДАННЫЕ:\n"
        + json.dumps(compact, ensure_ascii=False, indent=1)
    )
    # Параметры генерации у этой версии гарда — ОТДЕЛЬНЫЕ аргументы (fmt/
    # temperature/num_ctx), а не словарь options: сигнатуру взял из самого
    # _vram_guard.py, а не из примера в документации скилла (первая попытка
    # была написана по примеру и упала с TypeError — ровно та ошибка, от
    # которой защищает правило «сверять с исходником, а не с памятью»).
    resp = guarded_generate(
        "qwen3-coder:30b",
        prompt,
        want_gpu=True,
        priority=50,
        max_wait_s=max_wait_s,
        project="geant4-detector-models",
        agent="geant4-log-summarizer",
        fmt="json",
        temperature=0,
        num_ctx=32768,
    )
    try:
        body = resp if isinstance(resp, str) else resp.get("response", "{}")
        return json.loads(body)
    except (json.JSONDecodeError, AttributeError, TypeError) as e:
        # §4: провал Ollama — significant event, не «потерпим». Возвращаем
        # признак наверх, а не молча пустую сводку.
        return {"ollama_failure": str(e), "raw": str(resp)[:500]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="+")
    ap.add_argument("--explain", action="store_true",
                    help="добавить сводку от Ollama (guarded_generate)")
    ap.add_argument("--max-wait", type=int, default=600)
    a = ap.parse_args()

    files = []
    for p in a.paths:
        files.extend(sorted(glob.glob(p)) or [p])

    reports = [parse_log(p) for p in files]
    out = {
        "n_files": len(reports),
        "n_results": sum(len(r["results"]) for r in reports),
        "n_exceptions": sum(r.get("n_exceptions", 0) for r in reports),
        "n_fatal": sum(len(r.get("fatal", [])) for r in reports),
        "files": reports,
    }
    if a.explain:
        out["ollama_summary"] = explain(reports, a.max_wait)
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()

# ---------------------------------------------------------------------------
# ЗАМЕЧАНИЕ О ПРОИСХОЖДЕНИИ ЭТОГО ФАЙЛА (честно, для ревизии).
# IRON MODE скилла `workflow` предписывает генерировать код через Ollama
# (gen_code.py), а Claude оставляет за собой только спецификацию. Здесь это
# НЕ выполнено: файл написан Claude напрямую, потому что оператор в этот
# момент дал прямой запрет «ничего не запускай», а вызов Ollama — запуск.
# Приоритет указания оператора выше процедурного правила скилла.
# При снятии запрета файл подлежит ревизии по IRON MODE.
# ---------------------------------------------------------------------------