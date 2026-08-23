import numpy as np

def print_window_residuals(e_meas, cps_meas, models, fwhm_fn):
    WINS = [("Bi214 609", 609.3), ("Ac228 911", 911.2),
            ("K40 1460", 1460.8), ("Tl208 2614", 2614.5)]
    
    results = {}
    print("=== NEVYAZKA PO OKNAM LINIJ (gross, znak: +model vyshe izmereniya) ===")
    header = f"{'window':<18} {'izm,cps':>11}"
    for name, _ in models:
        header += f" {name + ',cps':>12} {name + ',%':>9}"
    print(header)
    
    for label, e0 in WINS:
        sigma = fwhm_fn(e0) / 2.35482
        mask = (e_meas >= e0 - 2.5*sigma) & (e_meas <= e0 + 2.5*sigma)
        ym = np.sum(cps_meas[mask])
        row_data = [f"{label:<18} {ym:>11.4f}"]
        results[label] = {}
        for name, model in models:
            rm = np.sum(model[mask])
            if ym <= 0:
                percent = float("nan")
            else:
                percent = 100.0 * (rm - ym) / ym
            row_data.append(f"{rm:>12.4f} {percent:>+9.1f}")
            results[label][name] = percent
        print(" ".join(row_data))
    
    muon_mask = (e_meas >= 2700) & (e_meas < 2828)
    ym = np.sum(cps_meas[muon_mask])
    label = "mu 2700-2828"
    row_data = [f"{label:<18} {ym:>11.4f}"]
    results[label] = {}
    for name, model in models:
        rm = np.sum(model[muon_mask])
        if ym <= 0:
            percent = float("nan")
        else:
            percent = 100.0 * (rm - ym) / ym
        row_data.append(f"{rm:>12.4f} {percent:>+9.1f}")
        results[label][name] = percent
    print(" ".join(row_data))
    
    return results
