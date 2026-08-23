"""
geom_sanity_check.py: Скрипт-шлюз §34-петли «response-geometry-sanity-radiacode103».

В отличие от regression_check.py (который ловит изменения физики на неизменной геометрии),
этот скрипт проверяет, что после НАМЕРЕННОЙ правки геометрии детектора кристалл всё ещё
физически вменяемо откликается на гамма-кванты. Сравнение с прошлым baseline не производится,
так как изменение геометрии закономерно меняет отклик.

См. audit/loops.md для описания контура GEANT4.
"""

import sys
import os
import json
import hashlib
import time
from pathlib import Path

# Ступень -1a: переиспользование логики запуска из соседнего скрипта
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from regression_check import run_macro, EXE, GEANT4_ENV

# Настройка кодировки вывода для корректного отображения кириллицы на Windows
sys.stdout.reconfigure(encoding="utf-8")

# --- Константы и конфигурация ---

SANITY_REFS = [
    {"name": "e0060", "macro": "sanity_060.mac", "e_kev": 60.0},
    {"name": "e0662", "macro": "ref_cs137.mac", "e_kev": 661.7},
    {"name": "e1461", "macro": "ref_k40.mac", "e_kev": 1460.8},
    {"name": "e2614", "macro": "sanity_2614.mac", "e_kev": 2614.5},
]

N_EVENTS = 30000
MIN_HITS = 30
MAX_HITS = N_EVENTS

GEOM_FILES = [
    r"D:\Claude_files\repos\geant4-detector-models\detectors\RadiaCode-103\geometry\RCDetector.cc",
    r"D:\Claude_files\repos\geant4-detector-models\detectors\RadiaCode-103\geometry\RCDetector.hh",
]

STATE_PATH = Path(__file__).parent / "geom_sanity_state.json"


def geometry_hash():
    """
    Вычисляет SHA1 хеш от конкатенации байтового содержимого файлов геометрии.
    Возвращает первые 12 символов hex-строки для краткости.
    """
    sha1 = hashlib.sha1()
    for file_path in GEOM_FILES:
        try:
            with open(file_path, "rb") as f:
                sha1.update(f.read())
        except FileNotFoundError:
            raise RuntimeError(f"Файл геометрии не найден: {file_path}")
        except IOError as e:
            raise RuntimeError(f"Ошибка чтения файла геометрии {file_path}: {e}")
    
    return sha1.hexdigest()[:12]


def main():
    # 3. Разбор аргументов
    force_run = "--force" in sys.argv

    # 1. Вычисление текущего хеша геометрии
    current_hash = geometry_hash()
    print(f"Текущий хеш геометрии: {current_hash}")

    # 4. Чтение прежнего состояния
    prev_state = None
    if STATE_PATH.exists():
        try:
            with open(STATE_PATH, "r", encoding="utf-8") as f:
                prev_state = json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            print(f"Предупреждение: ошибка чтения состояния {STATE_PATH}: {e}. Считаем состояние отсутствующим.")
            prev_state = None

    # 5. Проверка на необходимость пропуска
    if not force_run and prev_state is not None:
        saved_hash = prev_state.get("geometry_hash")
        saved_overall = prev_state.get("overall")
        
        if saved_hash == current_hash and saved_overall == "OK":
            print(f"АКТУАЛЬНО: геометрия не менялась с последней успешной проверки ({current_hash}), sanity-прогон пропущен")
            sys.exit(0)

    # 6. Прогон всех контрольных точек
    print("Запуск sanity-проверок...")
    details = []
    fail_count = 0
    
    for ref in SANITY_REFS:
        name = ref["name"]
        e_kev = ref["e_kev"]
        macro_name = ref["macro"]
        
        try:
            hits = run_macro(macro_name)
            
            # Проверка диапазона hits
            if MIN_HITS <= hits <= MAX_HITS:
                status = "OK"
                error_msg = None
            else:
                status = "FAIL"
                error_msg = f"hits={hits} вне диапазона [{MIN_HITS}, {MAX_HITS}]"
                fail_count += 1
                
            print(f"  %-8s E=%8.1f кэВ  hits=%-8s -> %s" % (name, e_kev, hits, status))
            
            details.append({
                "name": name,
                "e_kev": e_kev,
                "hits": hits,
                "status": status,
                "error": error_msg
            })
            
        except RuntimeError as e:
            # Ошибка запуска или парсинга
            status = "FAIL"
            error_msg = str(e)
            fail_count += 1
            
            print(f"  %-8s E=%8.1f кэВ  hits=%-8s -> %s" % (name, e_kev, "-", status))
            
            details.append({
                "name": name,
                "e_kev": e_kev,
                "hits": None,
                "status": status,
                "error": error_msg
            })

    # 7. Определение общего статуса
    overall = "FAIL" if fail_count > 0 else "OK"
    
    # 8. Запись нового состояния
    new_state = {
        "geometry_hash": current_hash,
        "overall": overall,
        "checked_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "details": details
    }
    
    try:
        with open(STATE_PATH, "w", encoding="utf-8") as f:
            json.dump(new_state, f, ensure_ascii=False, indent=2)
    except IOError as e:
        print(f"Ошибка записи состояния в {STATE_PATH}: {e}")

    # 9 и 10. Вывод результата и выход
    if overall == "FAIL":
        print(f"SANITY FAIL: геометрия правлена (RCDetector.cc/.hh), но отклик детектора на {fail_count} из {len(SANITY_REFS)} контрольных энергий не проходит проверку вменяемости — см. geom_sanity_state.json")
        sys.exit(1)
    else:
        print(f"SANITY OK: геометрия физически вменяема на всех {len(SANITY_REFS)} контрольных энергиях (hash={current_hash}). НАПОМИНАНИЕ: это НЕ заменяет полный производственный пересчёт results/bare/response/resp_*.csv (28 точек, ручное действие, см. DECISIONS.md D-002) — та сетка используется в подгонке (anchor_lines.py) и должна быть пересчитана вручную при значимой правке геометрии.")
        sys.exit(0)


if __name__ == "__main__":
    main()
