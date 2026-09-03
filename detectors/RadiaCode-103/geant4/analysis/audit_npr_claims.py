# -*- coding: utf-8 -*-
import sys
import os
import re
import numpy as np
from collections import defaultdict

sys.stdout.reconfigure(encoding="utf-8")

def parse_md_table(lines):
    """Разбор таблицы из markdown. Возвращает список кортежей."""
    result = []
    for line in lines:
        if not line.startswith("|"):
            continue
        parts = [p.strip() for p in line.split("|")[1:-1]]
        if len(parts) < 4:
            continue
        try:
            e = float(parts[0].replace("−", "-").replace(",", "."))
            v1 = float(parts[1].replace(",", "."))
            v2 = float(parts[2].replace(",", "."))
            v3 = float(parts[3].replace(",", "."))
            result.append((e, v1, v2, v3))
        except ValueError:
            continue
    return result

def load_curve(path):
    """Загрузка CSV-файла с кривой. Возвращает два массива."""
    xs = []
    ys = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split(",")
            if len(parts) < 2:
                continue
            try:
                x = float(parts[0])
                y = float(parts[1])
                xs.append(x)
                ys.append(y)
            except ValueError:
                continue
    return np.array(xs), np.array(ys)

def interp_log(E, xs, ys):
    """Линейная интерполяция по log10(E)."""
    if E <= 0:
        return 0.0
    log_x = np.log10(xs)
    log_E = np.log10(E)
    idx = np.searchsorted(log_x, log_E)
    if idx == 0:
        return ys[0]
    elif idx >= len(ys):
        return ys[-1]
    else:
        x0, x1 = log_x[idx-1], log_x[idx]
        y0, y1 = ys[idx-1], ys[idx]
        t = (log_E - x0) / (x1 - x0)
        return y0 + t * (y1 - y0)

def check(name, got, want, tol):
    """Сравнение с допуском."""
    diff = abs(got - want)
    verdict = "OK" if diff <= tol else "FAIL"
    print(f"{name}: {got:.6f} vs {want:.6f}, diff={diff:.6f}, {verdict}")
    return verdict == "OK"

def p1_csitable():
    """Проверка таблицы CsI(Tl)."""
    doc_path = r"D:\GoogleDrive\Рабочая папка ИИ\GEANT4\references\nonproportionality-csi-tl-concept-2026-09-03.md"
    with open(doc_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    
    table_start = -1
    for i, line in enumerate(lines):
        if line.startswith("## 4."):
            table_start = i + 1
            break
    
    if table_start == -1:
        print("P1: FAIL — не найден раздел ## 4.")
        return False

    table_lines = []
    for i in range(table_start, len(lines)):
        line = lines[i]
        if line.startswith("## "):
            break
        if line.startswith("|"):
            table_lines.append(line)
    
    table_data = parse_md_table(table_lines)
    csv_paths = [
        os.path.join(os.path.dirname(__file__), "out", "payne3_fig2_CsI_Tl_T_m40.csv"),
        os.path.join(os.path.dirname(__file__), "out", "payne3_fig2_CsI_Tl_T_0.csv"),
        os.path.join(os.path.dirname(__file__), "out", "payne3_fig2_CsI_Tl_T_p40.csv")
    ]
    
    if not all(os.path.exists(p) for p in csv_paths):
        print("P1: FAIL — не все CSV-файлы найдены.")
        return False

    csv_data = [load_curve(p) for p in csv_paths]
    ok_count = 0
    fail_count = 0
    
    for e, d1, d2, d3 in table_data:
        v1 = interp_log(e, csv_data[0][0], csv_data[0][1])
        v2 = interp_log(e, csv_data[1][0], csv_data[1][1])
        v3 = interp_log(e, csv_data[2][0], csv_data[2][1])
        
        diff1 = abs(d1 - v1)
        diff2 = abs(d2 - v2)
        diff3 = abs(d3 - v3)
        
        verdict1 = "OK" if diff1 <= 0.0006 else "FAIL"
        verdict2 = "OK" if diff2 <= 0.0006 else "FAIL"
        verdict3 = "OK" if diff3 <= 0.0006 else "FAIL"
        
        print(f"{e:.1f}: doc={d1:.3f}, csv={v1:.6f}, diff={diff1:.6f}, {verdict1}")
        print(f"{e:.1f}: doc={d2:.3f}, csv={v2:.6f}, diff={diff2:.6f}, {verdict2}")
        print(f"{e:.1f}: doc={d3:.3f}, csv={v3:.6f}, diff={diff3:.6f}, {verdict3}")
        
        if verdict1 == "OK":
            ok_count += 1
        else:
            fail_count += 1
        if verdict2 == "OK":
            ok_count += 1
        else:
            fail_count += 1
        if verdict3 == "OK":
            ok_count += 1
        else:
            fail_count += 1
    
    print(f"P1: {ok_count} OK, {fail_count} FAIL")
    return fail_count == 0

def p2_nai():
    """Проверка таблицы NaI(Tl)."""
    doc_path = r"D:\GoogleDrive\Рабочая папка ИИ\GEANT4\references\nonproportionality-csi-tl-concept-2026-09-03.md"
    with open(doc_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    
    table_start = -1
    for i, line in enumerate(lines):
        if line.startswith("### 8.2"):
            table_start = i + 1
            break
    
    if table_start == -1:
        print("P2: FAIL — не найден раздел ### 8.2.")
        return False

    table_lines = []
    for i in range(table_start, len(lines)):
        line = lines[i]
        if line.startswith("## "):
            break
        if line.startswith("|"):
            table_lines.append(line)
    
    table_data = parse_md_table(table_lines)
    csv_paths = [
        os.path.join(os.path.dirname(__file__), "out", "payne3_fig2_NaI_Tl_T_m40.csv"),
        os.path.join(os.path.dirname(__file__), "out", "payne3_fig2_NaI_Tl_T_0.csv"),
        os.path.join(os.path.dirname(__file__), "out", "payne3_fig2_NaI_Tl_T_p40.csv")
    ]
    
    if not all(os.path.exists(p) for p in csv_paths):
        print("P2: FAIL — не все CSV-файлы найдены.")
        return False

    csv_data = [load_curve(p) for p in csv_paths]
    ok_count = 0
    fail_count = 0
    
    for e, d1, d2, d3 in table_data:
        v1 = interp_log(e, csv_data[0][0], csv_data[0][1])
        v2 = interp_log(e, csv_data[1][0], csv_data[1][1])
        v3 = interp_log(e, csv_data[2][0], csv_data[2][1])
        
        diff1 = abs(d1 - v1)
        diff2 = abs(d2 - v2)
        diff3 = abs(d3 - v3)
        
        verdict1 = "OK" if diff1 <= 0.0006 else "FAIL"
        verdict2 = "OK" if diff2 <= 0.0006 else "FAIL"
        verdict3 = "OK" if diff3 <= 0.0006 else "FAIL"
        
        print(f"{e:.1f}: doc={d1:.3f}, csv={v1:.6f}, diff={diff1:.6f}, {verdict1}")
        print(f"{e:.1f}: doc={d2:.3f}, csv={v2:.6f}, diff={diff2:.6f}, {verdict2}")
        print(f"{e:.1f}: doc={d3:.3f}, csv={v3:.6f}, diff={diff3:.6f}, {verdict3}")
        
        if verdict1 == "OK":
            ok_count += 1
        else:
            fail_count += 1
        if verdict2 == "OK":
            ok_count += 1
        else:
            fail_count += 1
        if verdict3 == "OK":
            ok_count += 1
        else:
            fail_count += 1
    
    print(f"P2: {ok_count} OK, {fail_count} FAIL")
    return fail_count == 0

def p3_maxima():
    """Проверка максимумов кривых."""
    doc_path = r"D:\GoogleDrive\Рабочая папка ИИ\GEANT4\references\nonproportionality-csi-tl-concept-2026-09-03.md"
    
    # CsI(Tl) максимум
    cs_e = 8.6
    cs_val = 1.218
    
    # NaI(Tl) максимум
    na_e = 12.7
    na_val = 1.162
    
    csv_paths = [
        os.path.join(os.path.dirname(__file__), "out", "payne3_fig2_CsI_Tl_T_0.csv"),
        os.path.join(os.path.dirname(__file__), "out", "payne3_fig2_NaI_Tl_T_0.csv")
    ]
    
    if not all(os.path.exists(p) for p in csv_paths):
        print("P3: FAIL — не все CSV-файлы найдены.")
        return False

    cs_data = load_curve(csv_paths[0])
    na_data = load_curve(csv_paths[1])
    
    max_idx_cs = np.argmax(cs_data[1])
    max_e_cs = cs_data[0][max_idx_cs]
    max_val_cs = cs_data[1][max_idx_cs]
    
    max_idx_na = np.argmax(na_data[1])
    max_e_na = na_data[0][max_idx_na]
    max_val_na = na_data[1][max_idx_na]
    
    ok1 = check("P3 CsI(Tl) E", max_e_cs, cs_e, 0.05)
    ok2 = check("P3 CsI(Tl) val", max_val_cs, cs_val, 0.0006)
    ok3 = check("P3 NaI(Tl) E", max_e_na, na_e, 0.05)
    ok4 = check("P3 NaI(Tl) val", max_val_na, na_val, 0.0006)
    
    return ok1 and ok2 and ok3 and ok4

def p4_geant_materials():
    """Проверка материалов в Geant4."""
    path = r"C:\g4work\thirdparty\geant4-11.2.1\source\materials\src\G4NistMaterialBuilder.cc"
    
    if not os.path.exists(path):
        print("P4: SKIP — файл не найден.")
        return True
    
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    
    found1 = False
    found2 = False
    line_num1 = -1
    line_num2 = -1
    
    for i, line in enumerate(lines):
        if 'AddMaterial("G4_CESIUM_IODIDE", 4.51, 0, 553.1, 2);' in line:
            found1 = True
            line_num1 = i + 1
        if 'AddMaterial("G4_SODIUM_IODIDE", 3.667, 0, 452., 2);' in line:
            found2 = True
            line_num2 = i + 1
    
    if not found1 or not found2:
        print("P4: FAIL — не найдены нужные строки.")
        return False
    
    print(f"P4 CsI: строка {line_num1}")
    print(f"P4 NaI: строка {line_num2}")
    return True

def p5_thresholds():
    """П5. Пороги продукции из таблицы couples в логе прогона.

    Формат Geant4: гамма и e- стоят на ОДНОЙ строке `Energy thresholds`, а не на
    двух следующих (ошибка первой редакции, снята при аудите 03.09.2026), и блок
    материала может стоять как выше, так и ниже любого другого — искать по всему
    файлу, а не вперёд от найденного.
    """
    import re
    log_path = os.path.join(os.path.dirname(__file__), "..", "run_field", "logs",
                            "_chk_origin_on.log")
    if not os.path.exists(log_path):
        print("P5: SKIP — файл лога не найден:", log_path)
        return True
    with open(log_path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    unit = {"eV": 1e-3, "keV": 1.0, "MeV": 1e3, "GeV": 1e6}

    def thresholds_for(material):
        """Вернуть (gamma_keV, e_keV) для материала или None."""
        for i, ln in enumerate(lines):
            if f"Material : {material}" in ln:
                for j in range(i + 1, min(i + 6, len(lines))):
                    if "Energy thresholds" in lines[j]:
                        txt = lines[j]
                        g = re.search(r"gamma\s+([\d.eE+-]+)\s*(eV|keV|MeV|GeV)", txt)
                        e = re.search(r"(?<![A-Za-z])e-\s+([\d.eE+-]+)\s*(eV|keV|MeV|GeV)", txt)
                        if not g or not e:
                            return None
                        return (float(g.group(1)) * unit[g.group(2)],
                                float(e.group(1)) * unit[e.group(2)])
        return None

    results = []
    for material, want_g, want_e in (("G4_Pb", 18.2, 154.4),
                                     ("G4_CESIUM_IODIDE", 6.03, 94.5)):
        got = thresholds_for(material)
        if got is None:
            print(f"P5: FAIL — не найден блок порогов для {material}")
            results.append(False)
            continue
        results.append(check(f"P5 {material} gamma", got[0], want_g, 0.1))
        results.append(check(f"P5 {material} e-", got[1], want_e, 0.1))
    ok = all(results)
    print(f"P5: {sum(results)} OK, {len(results) - sum(results)} FAIL")
    return ok


def p6_deex():
    """Проверка результатов deex."""
    base_dir = os.path.join(os.path.dirname(__file__), "..", "run_field", "output")
    
    files = [f for f in os.listdir(base_dir) if f.startswith("_deexchk_Bi214_") and f.endswith(".csv")]
    
    if len(files) != 3:
        print("P6: FAIL — не найдены все файлы.")
        return False
    
    data = {}
    for f in files:
        path = os.path.join(base_dir, f)
        with open(path, "r", encoding="utf-8") as file:
            lines = file.readlines()
        
        header = {}
        table_start = -1
        for i, line in enumerate(lines):
            if not line.strip():
                table_start = i + 1
                break
            if "," in line:
                k, v = line.split(",", 1)
                header[k.strip()] = v.strip()
        
        if table_start == -1:
            print("P6: FAIL — не найдена таблица.")
            return False
        
        table_lines = lines[table_start:]
        table_data = []
        for line in table_lines:
            parts = line.strip().split(",")
            if len(parts) >= 2:
                try:
                    bin_kev = float(parts[0])
                    counts = float(parts[1])
                    table_data.append((bin_kev, counts))
                except ValueError:
                    continue
        
        data[f] = {
            "header": header,
            "table": table_data
        }
    
    # Проверить n_events и em_deex
    n_events = set()
    em_deex = set()
    for name, info in data.items():
        n_events.add(info["header"].get("n_events"))
        em_deex.add(info["header"].get("em_deex"))
    
    if len(n_events) != 1 or list(n_events)[0] != "100000000":
        print("P6: FAIL — n_events не совпадают или не равны 100000000.")
        return False
    
    if len(em_deex) != 3:
        print("P6: FAIL — не все значения em_deex найдены.")
        return False
    
    # Проверить em_cut_mm
    for name, info in data.items():
        cut = info["header"].get("em_cut_mm")
        if cut != "0.05":
            print(f"P6: FAIL — em_cut_mm не равен 0.05 в {name}.")
            return False
    
    # Пересчитать отношения
    std_data = data["_deexchk_Bi214_std.csv"]["table"]
    deex_data = data["_deexchk_Bi214_deex.csv"]["table"]
    max_data = data["_deexchk_Bi214_max.csv"]["table"]
    
    def get_counts(data, start, end):
        total = 0
        for bin_kev, counts in data:
            if start <= bin_kev < end:
                total += counts
        return total
    
    # 20–100 keV
    std_20_100 = get_counts(std_data, 20, 100)
    deex_20_100 = get_counts(deex_data, 20, 100)
    max_20_100 = get_counts(max_data, 20, 100)
    
    ratio_deex_std_20_100 = deex_20_100 / std_20_100 if std_20_100 > 0 else 0
    ratio_max_std_20_100 = max_20_100 / std_20_100 if std_20_100 > 0 else 0
    
    # 100–300 keV
    std_100_300 = get_counts(std_data, 100, 300)
    deex_100_300 = get_counts(deex_data, 100, 300)
    max_100_300 = get_counts(max_data, 100, 300)
    
    ratio_deex_std_100_300 = deex_100_300 / std_100_300 if std_100_300 > 0 else 0
    ratio_max_std_100_300 = max_100_300 / std_100_300 if std_100_300 > 0 else 0
    
    # Проверить значения
    ok1 = check("P6 deex/std 20–100", ratio_deex_std_20_100, 1.001, 0.002)
    ok2 = check("P6 max/std 20–100", ratio_max_std_20_100, 1.027, 0.002)
    ok3 = check("P6 deex/std 100–300", ratio_deex_std_100_300, 1.006, 0.002)
    ok4 = check("P6 max/std 100–300", ratio_max_std_100_300, 0.996, 0.002)
    
    return ok1 and ok2 and ok3 and ok4

def p7_files():
    """Проверка наличия файлов."""
    # PDF
    pdf_dir = r"D:\GoogleDrive\Рабочая папка ИИ\GEANT4\references\pdf"
    if not os.path.exists(pdf_dir):
        print("P7: FAIL — директория PDF не найдена.")
        return False
    
    pdf_files = [f for f in os.listdir(pdf_dir) if f.endswith(".pdf")]
    print(f"P7 PDF: {len(pdf_files)} файлов")
    
    # Книги
    books_dir = r"D:\GoogleDrive\Дозиметрия\Книги\Гамма-спектрометрия"
    if not os.path.exists(books_dir):
        print("P7: FAIL — директория книг не найдена.")
        return False
    
    book_files = [f for f in os.listdir(books_dir) if any(kw in f for kw in ["Nonproportionality", "Khodyuk", "COSINE", "BATSE"])]
    print(f"P7 Книги: {len(book_files)} файлов")
    
    # Скрипты
    script_dir = os.path.dirname(__file__)
    scripts = ["digitize_payne_fig2.py", "compare_physics_runs.py", "calib_nonlinear.py"]
    found_scripts = [s for s in scripts if os.path.exists(os.path.join(script_dir, s))]
    print(f"P7 Скрипты: {len(found_scripts)} файлов")
    
    return True

def p8_links():
    """Проверка ссылок вида `файл:строка`."""
    doc_path = r"D:\GoogleDrive\Рабочая папка ИИ\GEANT4\references\nonproportionality-csi-tl-concept-2026-09-03.md"
    
    with open(doc_path, "r", encoding="utf-8") as f:
        content = f.read()
    
    pattern = r"`([^`]+)\.([ch][ch]):(\d+)`"
    matches = re.findall(pattern, content)
    
    root1 = r"D:\Claude_files\repos\geant4-detector-models"
    root2 = r"C:\g4work\thirdparty\geant4-11.2.1\source"
    
    ok_count = 0
    fail_count = 0
    
    for name, ext, line_num in matches:
        full_name = f"{name}.{ext}"
        
        found_file = None
        for root in [root1, root2]:
            for dirpath, _, filenames in os.walk(root):
                if full_name in filenames:
                    found_file = os.path.join(dirpath, full_name)
                    break
            if found_file:
                break
        
        if not found_file:
            print(f"P8: SKIP — файл {full_name} не найден.")
            continue
        
        with open(found_file, "r", encoding="utf-8") as f:
            lines = f.readlines()
        
        line_num = int(line_num)
        if len(lines) < line_num:
            print(f"P8: FAIL — файл {full_name}, строка {line_num} не существует.")
            fail_count += 1
        else:
            print(f"P8: OK — {full_name}:{line_num}")
            print(f"  {lines[line_num-1].strip()}")
            ok_count += 1
    
    print(f"P8: {ok_count} OK, {fail_count} FAIL")
    return fail_count == 0

def selftest():
    """Самопроверка."""
    # Тест 1
    test_line = "| 10 | 1,296 | 1,216 | 1,136 |"
    lines = [test_line]
    result = parse_md_table(lines)
    if not (len(result) == 1 and abs(result[0][0] - 10.0) < 0.001 and abs(result[0][1] - 1.296) < 0.001):
        print("SELFTEST FAIL: разбор таблицы не работает")
        return 1
    
    # Тест 2
    test_line = "| 10 | 1,296 | 1,216 | 1,136 |"
    lines = [test_line]
    result = parse_md_table(lines)
    if not (len(result) == 1 and abs(result[0][0] - 10.0) < 0.001 and abs(result[0][1] - 1.296) < 0.001):
        print("SELFTEST FAIL: разбор таблицы с минусом не работает")
        return 1
    
    # Тест 3
    if check("selftest", 1.226, 1.216, 0.0006):
        print("SELFTEST FAIL: сравнение не умеет краснеть")
        return 1
    
    # Тест 4
    xs = np.array([10, 100])
    ys = np.array([1.0, 2.0])
    val = interp_log(31.62, xs, ys)
    if abs(val - 1.5) > 0.01:
        print("SELFTEST FAIL: интерполяция не логарифмическая")
        return 1
    
    print("SELFTEST OK")
    return 0

def main(argv):
    """Основная функция."""
    doc_path = r"D:\GoogleDrive\Рабочая папка ИИ\GEANT4\references\nonproportionality-csi-tl-concept-2026-09-03.md"
    selftest_mode = False
    
    if "--selftest" in argv:
        return selftest()
    
    if "--doc" in argv:
        idx = argv.index("--doc")
        if idx + 1 < len(argv):
            doc_path = argv[idx+1]
    
    results = []
    
    try:
        results.append(("P1", p1_csitable()))
    except Exception as e:
        print(f"P1: ERROR — {e}")
        results.append(("P1", False))
    
    try:
        results.append(("P2", p2_nai()))
    except Exception as e:
        print(f"P2: ERROR — {e}")
        results.append(("P2", False))
    
    try:
        results.append(("P3", p3_maxima()))
    except Exception as e:
        print(f"P3: ERROR — {e}")
        results.append(("P3", False))
    
    try:
        results.append(("P4", p4_geant_materials()))
    except Exception as e:
        print(f"P4: ERROR — {e}")
        results.append(("P4", False))
    
    try:
        results.append(("P5", p5_thresholds()))
    except Exception as e:
        print(f"P5: ERROR — {e}")
        results.append(("P5", False))
    
    try:
        results.append(("P6", p6_deex()))
    except Exception as e:
        print(f"P6: ERROR — {e}")
        results.append(("P6", False))
    
    try:
        results.append(("P7", p7_files()))
    except Exception as e:
        print(f"P7: ERROR — {e}")
        results.append(("P7", False))
    
    try:
        results.append(("P8", p8_links()))
    except Exception as e:
        print(f"P8: ERROR — {e}")
        results.append(("P8", False))
    
    # Сводка
    print("\nСводка:")
    ok = 0
    fail = 0
    skip = 0
    
    for name, res in results:
        if res is True:
            ok += 1
        elif res is False:
            fail += 1
        else:
            skip += 1
    
    print(f"Проверка | OK | FAIL | SKIP")
    print("-" * 30)
    for name, res in results:
        status = "OK" if res is True else ("FAIL" if res is False else "SKIP")
        print(f"{name:6} | {status}")
    
    print(f"\nИтого: {ok} OK, {fail} FAIL, {skip} SKIP")
    
    return 1 if fail > 0 else 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
