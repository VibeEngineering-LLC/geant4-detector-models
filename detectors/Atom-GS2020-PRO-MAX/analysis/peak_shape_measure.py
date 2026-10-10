import sys
import json
import argparse
import numpy as np

GAUSS_RATIO = {
    0.5: 2.3548,
    0.2: 3.5885,
    0.1: 4.2919
}

def subtract_continuum(E, y, peak_keV, side_keV=30.0, half_win_keV=60.0):
    left_band = (E < peak_keV - half_win_keV) & (E >= peak_keV - half_win_keV - side_keV)
    right_band = (E > peak_keV + half_win_keV) & (E <= peak_keV + half_win_keV + side_keV)
    
    if np.sum(left_band) < 3 or np.sum(right_band) < 3:
        raise ValueError("боковая полоса пуста")
        
    band_indices = np.where(left_band | right_band)[0]
    coeffs = np.polyfit(E[band_indices], y[band_indices], 1)
    continuum = np.polyval(coeffs, E)
    return y - continuum

def cross_level(E, y, i_top, level, step):
    i = i_top
    while 0 <= i < len(y) and y[i] >= level:
        i += step
    if i < 0 or i >= len(y):
        return float("nan")
    # Linear interpolation
    if step > 0:
        i_prev = i - 1
        if i_prev < 0:
            return float("nan")
        x1, x2 = E[i_prev], E[i]
        y1, y2 = y[i_prev], y[i]
    else:
        i_next = i + 1
        if i_next >= len(y):
            return float("nan")
        x1, x2 = E[i], E[i_next]
        y1, y2 = y[i], y[i_next]
    if y2 == y1:
        return (x1 + x2) / 2
    t = (level - y1) / (y2 - y1)
    return x1 + t * (x2 - x1)

def measure_shape(E, y, peak_keV, side_keV=30.0, half_win_keV=60.0):
    net = subtract_continuum(E, y, peak_keV, side_keV, half_win_keV)
    search_mask = np.abs(E - peak_keV) <= half_win_keV
    search_indices = np.where(search_mask)[0]
    if len(search_indices) == 0:
        raise ValueError("пик не найден: высота <= 0")
    i_top = search_indices[np.argmax(net[search_indices])]
    h = net[i_top]
    if h <= 0:
        raise ValueError("пик не найден: высота <= 0")
    
    result = {
        "E_вершины": E[i_top],
        "высота": h,
        "ПШПВ": 0.0,
        "FW1/5": 0.0,
        "FW1/10": 0.0,
        "лево_1/2": 0.0,
        "право_1/2": 0.0,
        "лево_1/10": 0.0,
        "право_1/10": 0.0,
        "FW1/5_к_ПШПВ": 0.0,
        "FW1/10_к_ПШПВ": 0.0,
        "асимметрия_1/2": 0.0,
        "асимметрия_1/10": 0.0,
        "сигма_по_ПШПВ": 0.0
    }
    
    for f in (0.5, 0.2, 0.1):
        e_left = cross_level(E, net, i_top, f * h, -1)
        e_right = cross_level(E, net, i_top, f * h, +1)
        width_f = e_right - e_left
        left_f = E[i_top] - e_left
        right_f = e_right - E[i_top]
        
        if f == 0.5:
            result["ПШПВ"] = width_f
            result["лево_1/2"] = left_f
            result["право_1/2"] = right_f
        elif f == 0.2:
            result["FW1/5"] = width_f   # правка 0: генератор клал сюда ПШПВ, ширина на 1/5 терялась
            result["FW1/5_к_ПШПВ"] = width_f / result["ПШПВ"]
        elif f == 0.1:
            result["FW1/10"] = width_f
            result["лево_1/10"] = left_f
            result["право_1/10"] = right_f
            result["FW1/10_к_ПШПВ"] = width_f / result["ПШПВ"]
    
    result["асимметрия_1/2"] = result["лево_1/2"] / result["право_1/2"]
    result["асимметрия_1/10"] = result["лево_1/10"] / result["право_1/10"]
    result["сигма_по_ПШПВ"] = result["ПШПВ"] / 2.3548
    
    return result

def compare_to_gauss(shape: dict) -> dict:
    return {
        "FW1/5_изб_%": 100 * (shape["FW1/5_к_ПШПВ"] / 1.5239 - 1),
        "FW1/10_изб_%": 100 * (shape["FW1/10_к_ПШПВ"] / 1.8226 - 1)
    }

def format_report(name: str, shape: dict, cmp: dict) -> str:
    return f"""{name}
Энергия вершины: {shape['E_вершины']:.3f} кэВ
ПШПВ: {shape['ПШПВ']:.2f} кэВ ({100 * shape['ПШПВ'] / shape['E_вершины']:.1f}%)
FW1/5: {shape['FW1/5']:.2f} кэВ, избыток над Гауссом: {cmp['FW1/5_изб_%']:.1f}%
FW1/10: {shape['FW1/10']:.2f} кэВ, избыток над Гауссом: {cmp['FW1/10_изб_%']:.1f}%
Асимметрия 1/2: {shape['асимметрия_1/2']:.3f}
Асимметрия 1/10: {shape['асимметрия_1/10']:.3f}"""

if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True)
    parser.add_argument("--dataset", default="sample")
    parser.add_argument("--peak", type=float, default=661.657)
    parser.add_argument("--series", default="net")
    args = parser.parse_args()
    
    with open(args.data, "r", encoding="utf-8") as f:
        text = f.read()
    
    start = text.find("{")
    end = text.rfind("}") + 1
    json_text = text[start:end]
    d = json.loads(json_text)
    
    E = np.array(d[args.dataset]["E_keV"])
    y = np.array(d[args.dataset][args.series])
    
    shape = measure_shape(E, y, args.peak)
    cmp = compare_to_gauss(shape)
    print(format_report(args.dataset, shape, cmp))
