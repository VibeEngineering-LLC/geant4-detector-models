import re

def parse_lines(path):
    with open(path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    start = None
    for i, line in enumerate(lines):
        if line.startswith('## 9. БОЛЬШАЯ СВОДНАЯ ТАБЛИЦА ЛИНИЙ'):
            start = i + 1
            break
    else:
        raise ValueError("table not found")
    
    table_lines = []
    for i in range(start, len(lines)):
        line = lines[i]
        if line.startswith('## '):
            break
        if line.startswith('|') and not line.startswith('| E, кэВ |') and not line.startswith('|---|'):
            table_lines.append(line)
    
    result = []
    for line in table_lines:
        cells = [cell.strip().replace('**', '') for cell in line.split('|')[1:-1]]
        if len(cells) < 5:
            cells.extend([''] * (5 - len(cells)))
        
        energy_str = cells[0]
        energies = re.split(r'\s*/\s*', energy_str)
        energies = [re.sub(r'\s+', '', e).replace(',', '.') for e in energies]   # \s ловит и неразрывные пробелы
        
        numbers = []
        for e in energies:
            match = re.search(r'[0-9]+(?:\.[0-9]+)?', e)
            if match:
                num = float(match.group())
                if 1 <= num <= 20000:
                    numbers.append(num)
        
        if not numbers:
            continue
            
        reaction, yield_, place, cls = cells[1], cells[2], cells[3], cells[4]
        text = f"{reaction} | {place} | {cls}"
        
        for energy in numbers:
            result.append({
                "E_keV": energy,
                "text": text,
                "reaction": reaction,
                "yield": yield_,
                "place": place,
                "cls": cls
            })
    
    return sorted(result, key=lambda x: x["E_keV"])
