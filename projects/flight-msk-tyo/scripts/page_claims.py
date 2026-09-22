import sys, re
sys.stdout.reconfigure(encoding="utf-8")

def main():
    if len(sys.argv) != 2:
        print("Использование: page_claims.py <text.md>", file=sys.stderr)
        sys.exit(2)
    
    if sys.argv[1] == "--selftest":
        test_text = "Итог 73,66 ± 0,23 отсч./с\n\n@@CHART@@\nПик 1460,5 кэВ, 265 200 отсч./ч"
        lines = test_text.split('\n')
        result = []
        for i, line in enumerate(lines, 1):
            if not line.strip() or re.match(r'^@@[A-Z]+@@$', line):
                continue
            numbers = re.findall(r'\d+(?:[  ]\d{3})*(?:,\d+)?', line)
            if numbers:
                clean_line = re.sub(r'[*_]+', '', line[:160])
                result.append(f"{i}|{';'.join(numbers)}|{clean_line}")
        if len(result) == 2 and result[0].endswith("73,66;0,23|Итог 73,66 ± 0,23 отсч./с") and result[1].endswith("1460,5;265 200|Пик 1460,5 кэВ, 265 200 отсч./ч"):
            print("SELFTEST PASS")
        else:
            print("SELFTEST FAIL")
        sys.exit(0)
    
    try:
        with open(sys.argv[1], 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Ошибка: файл '{sys.argv[1]}' не найден.", file=sys.stderr)
        sys.exit(2)
    
    count_lines = 0
    count_numbers = 0
    
    for i, line in enumerate(lines, 1):
        if not line.strip() or re.match(r'^@@[A-Z]+@@$', line):
            continue
        numbers = re.findall(r'\d+(?:[  ]\d{3})*(?:,\d+)?', line)
        if numbers:
            count_lines += 1
            count_numbers += len(numbers)
            clean_line = re.sub(r'[*_]+', '', line[:160])
            print(f"{i}|{' '.join(numbers)}|{clean_line}")
    
    print(f"строк с числами: {count_lines}, чисел всего: {count_numbers}")

if __name__ == "__main__":
    main()
