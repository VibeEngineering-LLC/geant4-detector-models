"""
regression_check.py — объективный шлюз §34-петли «регрессия физики на эталонных задачах»
контурa GEANT4.

Скрипт запускает бинарник rc_curves.exe на эталонных макросах и сравнивает число
зарегистрированных событий (hits) с сохранённым эталоном. Допуск вычисляется как
3σ по Пуассону от эталонного значения, однако он является конфигурационным ожиданием
«с запасом», а не измеренной дисперсией: при фиксированном коде и отсутствии явного
RNG-seed Geant4 использует детерминированный дефолтный движок, обеспечивая полное
воспроизведение результата. Допуск существует для защиты от случайной рандомизации
seed в будущем (правило #AH-1: конфигурационное ожидание — не то же самое, что измерение).

Часть §34-петли контура GEANT4, карточка audit/loops.md.
"""

import sys
import os
import re
import json
import subprocess
import time
from pathlib import Path

# Настройка кодировки вывода для корректной работы с кириллицей в консоли Windows
sys.stdout.reconfigure(encoding="utf-8")

GEANT4_ENV = {
    "GEANT4_ROOT": r"C:\geant4",
    "G4NEUTRONHPDATA": r"C:\geant4\share\data\G4NDL4.7",
    "G4LEDATA": r"C:\geant4\share\data\G4EMLOW8.5",
    "G4LEVELGAMMADATA": r"C:\geant4\share\data\PhotonEvaporation5.7",
    "G4RADIOACTIVEDATA": r"C:\geant4\share\data\RadioactiveDecay5.6",
    "G4PARTICLEXSDATA": r"C:\geant4\share\data\G4PARTICLEXS4.0",
    "G4SAIDXSDATA": r"C:\geant4\share\data\G4SAIDDATA2.0",
    "G4ABLADATA": r"C:\geant4\share\data\G4ABLA3.3",
    "G4INCLDATA": r"C:\geant4\share\data\G4INCL1.2",
    "G4ENSDFSTATEDATA": r"C:\geant4\share\data\G4ENSDFSTATE2.3",
}

EXE = r"D:\Claude_files\repos\geant4-detector-models\build\RadiaCode-103\rc_curves.exe"
HERE = Path(os.path.dirname(os.path.abspath(__file__)))
BASELINE_PATH = HERE / "baseline.json"
STATE_PATH = HERE / "state.json"

REFS = [
    {"name": "cs137_662", "macro": "ref_cs137.mac"},
    {"name": "k40_1461", "macro": "ref_k40.mac"},
]

SIGMA_K = 3  # множитель допуска (3 сигма Пуассона)


def run_macro(macro_name):
    """
    Запускает бинарник с указанным макросом и возвращает количество hits.
    Поднимает RuntimeError при ошибке выполнения или отсутствии результата.
    """
    env = os.environ.copy()
    # Добавляем путь к Geant4 в начало PATH
    geant4_bin = r"C:\geant4\bin"
    if "PATH" in env:
        env["PATH"] = f"{geant4_bin};{env['PATH']}"
    else:
        env["PATH"] = geant4_bin
    
    # Устанавливаем переменные окружения Geant4
    for key, value in GEANT4_ENV.items():
        env[key] = value

    macro_path = os.path.join(HERE, macro_name)
    
    try:
        result = subprocess.run(
            [EXE, macro_path, "bare"],
            env=env,
            capture_output=True,
            text=True,
            timeout=120
        )
    except Exception as e:
        raise RuntimeError(f"Ошибка запуска процесса: {e}")

    if result.returncode != 0:
        raise RuntimeError(
            f"Процесс завершился с кодом {result.returncode}.\n"
            f"STDOUT (последние 2000 символов):\n{result.stdout[-2000:]}\n"
            f"STDERR (последние 2000 символов):\n{result.stderr[-2000:]}"
        )

    # Поиск строки RESULT и парсинг hits
    match = re.search(r"hits=\s*(\d+)", result.stdout)
    if not match:
        raise RuntimeError(
            f"Строка RESULT не найдена в выводе.\n"
            f"STDOUT (последние 2000 символов):\n{result.stdout[-2000:]}\n"
            f"STDERR (последние 2000 символов):\n{result.stderr[-2000:]}"
        )

    return int(match.group(1))


def load_baseline():
    """
    Читает файл baseline.json. Возвращает словарь или None, если файл отсутствует.
    """
    if not BASELINE_PATH.exists():
        return None
    
    try:
        with open(BASELINE_PATH, "r", encoding="utf-8") as f:
            return json.load(f)
    except (json.JSONDecodeError, IOError) as e:
        print(f"Ошибка чтения эталона: {e}")
        return None


def save_state(results, timestamp_str):
    """
    Сохраняет состояние проверки в state.json.
    results: список словарей с результатами по каждому эталону.
    timestamp_str: строка времени проверки.
    """
    overall_status = "OK" if all(r["status"] == "OK" for r in results) else "FAIL"
    
    state_data = {
        "checked_at": timestamp_str,
        "overall": overall_status,
        "details": results
    }

    try:
        with open(STATE_PATH, "w", encoding="utf-8") as f:
            json.dump(state_data, f, ensure_ascii=False, indent=2)
    except IOError as e:
        print(f"Ошибка сохранения состояния: {e}")


def main():
    # 1. Разбор аргументов командной строки
    update_baseline_mode = "--update-baseline" in sys.argv

    # 2. Загрузка эталона
    baseline = load_baseline()
    
    if baseline is None and not update_baseline_mode:
        print("БАЗОВЫЙ ЭТАЛОН НЕ НАЙДЕН — запустите с --update-baseline один раз")
        sys.exit(2)

    # 3. Запуск макросов и сбор результатов
    current_hits = {}
    for ref in REFS:
        try:
            hits = run_macro(ref["macro"])
            current_hits[ref["name"]] = hits
            print("[%s] hits=%d" % (ref["name"], hits))
        except RuntimeError as e:
            print(f"КРИТИЧЕСКАЯ ОШИБКА при проверке {ref['name']}: {e}")
            sys.exit(1)

    # 4. Обычный режим проверки
    if not update_baseline_mode:
        results_list = []
        has_fail = False
        
        for ref in REFS:
            name = ref["name"]
            base_hits = baseline[name]["hits"]
            new_hits = current_hits[name]
            
            # Вычисление допуска: 3 сигма Пуассона, минимум 1
            tol = max(1, round(SIGMA_K * (base_hits ** 0.5)))
            
            if abs(new_hits - base_hits) <= tol:
                status = "OK"
            else:
                status = "FAIL"
                has_fail = True
            
            print("  %-12s baseline=%-8d new=%-8d допуск=±%-6d -> %s" % (
                name, base_hits, new_hits, tol, status
            ))
            
            results_list.append({
                "name": name,
                "baseline_hits": base_hits,
                "new_hits": new_hits,
                "tolerance": tol,
                "status": status
            })

        # Сохранение состояния
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        save_state(results_list, timestamp)

        if has_fail:
            print("РЕГРЕССИЯ: физика изменилась за пределы допуска")
            sys.exit(1)
        else:
            print("OK: все эталоны в допуске")
            sys.exit(0)

    # 5. Режим обновления эталона
    else:
        new_baseline = {}
        today_date = time.strftime("%Y-%m-%d")
        
        for ref in REFS:
            name = ref["name"]
            new_baseline[name] = {
                "hits": current_hits[name],
                "date": today_date
            }
        
        try:
            with open(BASELINE_PATH, "w", encoding="utf-8") as f:
                json.dump(new_baseline, f, ensure_ascii=False, indent=2)
            print(f"Эталон обновлён: {BASELINE_PATH}")
            sys.exit(0)
        except IOError as e:
            print(f"Ошибка записи эталона: {e}")
            sys.exit(1)


if __name__ == "__main__":
    main()
