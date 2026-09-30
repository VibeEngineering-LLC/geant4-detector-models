import sys
sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import os

# Constants
ME = 511.0  # keV
Q = 1310.894  # keV
Z = 20
A = 40
fsc = 1 / 137.03599976
alphaZ = fsc * Z
Rnuc = 0.5 * fsc * A ** (1/3)
V0 = 1.13 * fsc**2 * abs(Z)**(4/3)
gamma0 = np.sqrt(1 - alphaZ**2)

gc = [-0.1010678, 0.4245549, -0.6998588, 0.9512363, -0.5748646, 1.0]
imMax = 200.0

def Gamma_func(arg):
    """Port of G4BetaDecayCorrections::Gamma"""
    fac = 1.0
    x = arg - 1.0
    
    loop = 0
    while x > 1.0:
        fac *= x
        x -= 1.0
        loop += 1
        if loop > 1000:
            break
            
    # Polynomial from Abramowitz and Stegun for 0 < arg < 1 (actually x in (-1, 0) usually, but code uses x = arg-1)
    # The C++ code does: sum = gc[0]; for i=1..5: sum = sum*x + gc[i]
    # This is Horner's method for polynomial P(x) = gc[0]*x^5 + ... + gc[5] ? 
    # Let's trace: 
    # i=1: sum = gc[0]*x + gc[1]
    # i=2: sum = (gc[0]*x + gc[1])*x + gc[2] = gc[0]*x^2 + gc[1]*x + gc[2]
    # ...
    # i=5: sum = gc[0]*x^5 + gc[1]*x^4 + gc[2]*x^3 + gc[3]*x^2 + gc[4]*x + gc[5]
    
    sum_val = gc[0]
    for i in range(1, 6):
        sum_val = sum_val * x + gc[i]
        
    return sum_val * fac

def ModSquared(re, im):
    """Port of G4BetaDecayCorrections::ModSquared"""
    im = np.clip(im, -imMax, imMax)
    
    factor1 = ((1 + re)**2 + im**2) ** (re + 0.5)
    factor2 = np.exp(2 * im * np.arctan(im / (1 + re)))
    factor3 = np.exp(2 * (1 + re))
    factor4 = 2 * np.pi
    factor5 = np.exp((1 + re) / ((1 + re)**2 + im**2) / 6)
    factor6 = re**2 + im**2
    
    return factor1 * factor4 * factor5 / factor2 / factor3 / factor6

def FermiFunction(W):
    """Port of G4BetaDecayCorrections::FermiFunction"""
    if Z < 0:
        Wprime = W + V0
    else:
        Wprime = W - V0
        if Wprime <= 1.00001:
            Wprime = 1.00001
            
    p_e = np.sqrt(Wprime**2 - 1)
    eta = alphaZ * Wprime / p_e
    epieta = np.exp(np.pi * eta)
    
    realGamma = Gamma_func(2 * gamma0 + 1)
    mod2Gamma = ModSquared(gamma0, eta)
    
    factor1 = 2 * (1 + gamma0) * mod2Gamma / realGamma / realGamma
    factor2 = epieta * (2 * p_e * Rnuc) ** (2 * (gamma0 - 1))
    factor3 = (Wprime / W) * np.sqrt((Wprime**2 - 1) / (W**2 - 1))
    
    return factor1 * factor2 * factor3

def ShapeFactor(p_e, e_nu):
    """Port of uniqueThirdForbidden shape factor"""
    # eta for shape factor calculation
    W_e = np.sqrt(1 + p_e**2)
    eta = alphaZ * W_e / p_e
    
    gamma1 = np.sqrt(4 - alphaZ**2)
    gamma2 = np.sqrt(9 - alphaZ**2)
    gamma3 = np.sqrt(16 - alphaZ**2)
    
    gamterm0 = Gamma_func(2 * gamma0 + 1)
    gamterm1 = gamterm0 / Gamma_func(2 * gamma1 + 1)
    gamterm2 = gamterm0 / Gamma_func(2 * gamma2 + 1)
    gamterm3 = gamterm0 / Gamma_func(2 * gamma3 + 1)
    
    twoPR = 2 * p_e * Rnuc
    
    # term1
    term1 = (e_nu**6 * (1 + gamma0)) / 1260.0
    
    # term2
    modSq_g1 = ModSquared(gamma1, eta)
    modSq_g0 = ModSquared(gamma0, eta)
    term2 = (2 * (2 + gamma1) * e_nu**4 * p_e**2 * 
             twoPR**(2 * (gamma1 - gamma0 - 1)) * 
             gamterm1**2 * modSq_g1 / modSq_g0) / 5.0
    
    # term3
    modSq_g2 = ModSquared(gamma2, eta)
    term3 = (60 * (3 + gamma2) * p_e**4 * e_nu**2 * 
             twoPR**(2 * (gamma2 - gamma0 - 2)) * 
             gamterm2**2 * modSq_g2 / modSq_g0)
    
    # term4
    modSq_g3 = ModSquared(gamma3, eta)
    term4 = (2240 * p_e**6 * (4 + gamma3) * 
             twoPR**(2 * (gamma3 - gamma0 - 3)) * 
             gamterm3**2 * modSq_g3 / modSq_g0)
    
    return term1 + term2 + term3 + term4

def main():
    # Grid for Geant4 port
    n = 3000
    T_keV = np.linspace(0, Q, n, endpoint=False) + Q / (2 * n) # midpoints
    
    ex = T_keV / ME
    p = np.sqrt(ex * (ex + 2))
    en = Q / ME - ex
    
    # Calculate Fermi Function vectorized
    W = 1 + ex
    F_vals = np.array([FermiFunction(w) for w in W])
    
    # Calculate Shape Factor vectorized
    SF_vals = np.array([ShapeFactor(p[i], en[i]) for i in range(n)])
    
    f = p * (1 + ex) * en**2 * F_vals * SF_vals
    
    # Normalize
    norm = np.trapezoid(f, T_keV)
    Pg = f / norm
    Tk = T_keV
    
    # Import IB spectrum module
    import ib_spectrum_kub as kub
    
    # Get Kub shape
    Tk2, Pk = kub.beta(Q, 20, "u3")
    
    # 1) Mean beta energy and max relative difference
    mean_E_g4 = np.trapezoid(Tk * Pg, Tk)
    mean_E_kub = np.trapezoid(Tk2 * Pk, Tk2)
    
    # Interpolate Kub to G4 grid for comparison
    Pk_interp = np.interp(Tk, Tk2, Pk)
    
    # Max relative difference
    # Guard against zero in denominator
    mask = Pk_interp > 1e-15
    if np.any(mask):
        rel_diff = np.abs((Pg[mask] - Pk_interp[mask]) / Pk_interp[mask])
        max_rel_diff_pct = np.max(rel_diff) * 100
    else:
        max_rel_diff_pct = 0.0
        
    print(f"Средняя энергия бета (G4): {mean_E_g4:.4f} keV")
    print(f"Средняя энергия бета (Kub): {mean_E_kub:.4f} keV")
    print(f"Максимальное относительное различие Pg vs Pk_interp: {max_rel_diff_pct:.2f}%")
    
    # 2) IB Spectrum Calculation
    ks = np.arange(0.5, 1310.0, 1.0) # keV
    
    # dNdk for G4 shape
    # dNdk(k) = integral over T of P(T) * phi(W, k) / ME dT
    # W = 1 + T/ME
    # k in units of ME for phi
    
    dNdk_g4 = np.zeros_like(ks)
    for i, k_keV in enumerate(ks):
        k_ME = k_keV / ME
        # Vectorized over T
        W_vals = 1 + Tk / ME
        # phi(W, k)
        phi_vals = np.array([kub.phi(w, k_ME) for w in W_vals])
        integrand = Pg * phi_vals / ME
        dNdk_g4[i] = np.trapezoid(integrand, Tk)
        
    # dNdk for Kub shape
    dNdk_ours = np.zeros_like(ks)
    for i, k_keV in enumerate(ks):
        k_ME = k_keV / ME
        W_vals = 1 + Tk2 / ME
        phi_vals = np.array([kub.phi(w, k_ME) for w in W_vals])
        integrand = Pk * phi_vals / ME
        dNdk_ours[i] = np.trapezoid(integrand, Tk2)
        
    # Bands
    bands = [
        (25, 150),
        (150, 400),
        (400, 1311),
        (25, 1311)
    ]
    
    print("\nПолосы (keV) | Сумма dNdk G4 | Сумма dNdk Kub | Отношение")
    print("-" * 60)
    
    for low, high in bands:
        mask = (ks >= low) & (ks < high)
        sum_g4 = np.sum(dNdk_g4[mask])
        sum_ours = np.sum(dNdk_ours[mask])
        
        if sum_ours > 1e-15:
            ratio = sum_g4 / sum_ours
        else:
            ratio = 0.0
            
        print(f"{low}-{high} | {sum_g4:.6f} | {sum_ours:.6f} | {ratio:.4f}")
        
    # 3) Save CSV
    script_dir = os.path.dirname(os.path.abspath(__file__))
    audit_dir = os.path.join(script_dir, "..", "audit")
    os.makedirs(audit_dir, exist_ok=True)
    csv_path = os.path.join(audit_dir, "gs45_k40_ib_shapes.csv")
    
    np.savetxt(csv_path, 
               np.column_stack((ks, dNdk_g4, dNdk_ours)), 
               delimiter=",", 
               header="k_keV,dNdk_g4shape,dNdk_ours_u3", 
               comments="")

if __name__ == "__main__":
    main()
