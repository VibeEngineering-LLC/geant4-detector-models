import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import fit_coverage as fc
import fit_nuclides as fn
from scipy.optimize import minimize_scalar, minimize
import io, contextlib

def add_tail(model, e_meas, a_frac, beta_kev):
    de = np.median(np.diff(e_meas))
    # Create a matrix of energy differences for broadcasting
    diff = e_meas[:, None] - e_meas[None, :]
    # Only consider contributions from higher energies (j < i)
    mask = diff > 0
    # Compute tail contribution using broadcasting
    tail_contrib = model[:, None] * a_frac * np.exp(-diff / beta_kev) * de / beta_kev
    tail_contrib = np.where(mask, tail_contrib, 0)
    # tail_contrib[i,j] = vklad ot istochnika i (vyshe po energii) v priyomnik j.
    # Nuzhna summa PO ISTOCHNIKAM i DLYA KAZHDOGO priyomnika j -> axis=0, ne axis=1
    # (axis=1 oshibochno summiroval by "iskhodyashchij potok" iz kazhdogo i, chto
    # dobavlyaet kanalu ego zhe sobstvennyj hvost vmesto hvostov ot kanalov vyshe).
    total_tail = np.sum(tail_contrib, axis=0)
    return model + total_tail

def main():
    names, A, e_meas, cps_meas, live = fc.assemble()
    
    with contextlib.redirect_stdout(io.StringIO()):
        amps_dict, a_mu = fn.main()
    
    amp = np.array([amps_dict.get(n, 0.0) for n in names])
    model_base = A @ amp
    
    # Restrict to analysis window
    mask = (e_meas >= 20) & (e_meas < 2830)
    e, y, m0 = e_meas[mask], cps_meas[mask], model_base[mask]
    
    def objective(params):
        a_frac, beta_kev = params
        if a_frac < 0 or beta_kev <= 0:
            return 1e12
        model_tail = add_tail(m0, e, a_frac, beta_kev)
        residuals = (model_tail - y) / np.maximum(y, 1e-6)
        return np.sum(residuals ** 2)

    # Bounds fizicheski osmyslennye: pri beta_kev -> bolshoe znachenie (sravnimoe
    # s diapazonom spektra 2830 keV) exp(-diff/beta) stanovitsya ploskim dlya vsekh
    # diff, i rabotaet tolko otnoshenie a_frac/beta - Nelder-Mead bez granic ubegaet
    # po etomu ploskomu napravleniyu v beskonechnost (naideno mutacionnoj proverkoj,
    # sm. sintetichesky test v chate 22.08). beta_kev <= 300 keV - razumnyj verkhnij
    # predel (neskolko FWHM detektora ~30-100 keV), a_frac <= 20 - shchedryj zapas.
    result = minimize(objective, x0=[0.5, 100.0], method="L-BFGS-B",
                      bounds=[(0.0, 20.0), (1.0, 300.0)])
    
    print(f"{'success':<12} a_frac={result.x[0]:8.4f} beta_kev={result.x[1]:8.1f} chi2={result.fun:10.2f}")
    
    # Metrics
    fc_baseline = fc.fraction_covered(m0, y)
    fr_baseline = fc.form_residual_pct(m0, y)
    
    fc_tail = fc.fraction_covered(add_tail(m0, e, result.x[0], result.x[1]), y)
    fr_tail = fc.form_residual_pct(add_tail(m0, e, result.x[0], result.x[1]), y)
    
    # fraction_covered vozvrashchaet DOLYU [0,1], ne procent - v otlichie ot
    # form_residual_pct (uzhe procent). Ollama-generaciya pechatala fc_baseline/
    # fc_tail bez *100 (0.381 -> "0.4%" vmesto 38.1%) - najdeno vychitkoj 22.08.
    print("%-24s baseline=%6.1f%%  s tail=%6.1f%%" % ("%% zapolneniya formoj:", fc_baseline * 100, fc_tail * 100))
    print("%-24s baseline=%6.1f%%  s tail=%6.1f%%" % ("nevyazka formy:", fr_baseline, fr_tail))
    
    # Band-wise analysis
    bands = ((20, 100), (100, 300), (300, 700), (700, 1500), (1500, 2000),
             (2000, 2400), (2400, 2830))
    
    print("\nBand       Measured   Baseline   With tail   Deficit   Deficit")
    print("           Sum        Sum        Sum         %%        %%")
    
    for e_min, e_max in bands:
        band_mask = (e >= e_min) & (e < e_max)
        if not np.any(band_mask):
            continue
        y_band = y[band_mask]
        m0_band = m0[band_mask]
        m_tail_band = add_tail(m0_band, e[band_mask], result.x[0], result.x[1])
        
        y_sum = np.sum(y_band)
        m0_sum = np.sum(m0_band)
        m_tail_sum = np.sum(m_tail_band)
        
        def_pct_base = 100 * (y_sum - m0_sum) / y_sum
        def_pct_tail = 100 * (y_sum - m_tail_sum) / y_sum
        
        print("%5d-%-6d %10.5f %10.5f %10.5f %8.1f%% %8.1f%%" % (
            e_min, e_max, y_sum, m0_sum, m_tail_sum, def_pct_base, def_pct_tail))
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
