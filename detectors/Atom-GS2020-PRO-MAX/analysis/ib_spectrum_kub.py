# IB spectrum of beta- decay (KUB), photons per decay per keV.
# Sources: Knipp & Uhlenbeck, Physica 3, 425 (1936); Bloch, Phys. Rev. 50, 272 (1936) - phi(); unique-forbidden shape
# factors SH: Konopinski & Uhlenbeck, Phys. Rev. 60, 308 (1941). No Coulomb correction in phi() (Z only in beta()).
# Full list: audit/GS-42-method2-effects.md, section "Pervoistochniki A3-A6".
# Formulas beta/phi/SH copied from verifier ib_calc.py, copied 2026-09-15 (copy forced: scratchpad does not survive session):
#   (scratchpad of the 15.09.2026 session, not preserved)
# Spec: scripts/_spec_ib_template.md, part A.
# Branches/Q from LNHB: K40_tables.txt (beta- 89.56 %, Q 1310.91 keV); Cs137_tables.txt (94.57 % Q 513.97; 5.43 % Q 1175.63; 0.0006 % Q 892.17 omitted per spec).
# CSV grid: k_keV = centre of 1 keV bin [k-0.5, k+0.5], k = kmin+0.5, kmin+1.5, ... while k < Q_max.
# Header Y values: trapezoid on an independent 0.25 keV grid (not the CSV grid).
# Branch file mode (2026-09-15): --branches FILE.json {"shape_label": str, "branches": [{"fraction", "Q_keV", "Z_daughter", "shape"}, ...]}; extra keys ignored.
import sys, argparse, json, numpy as np
sys.stdout.reconfigure(encoding='utf-8')
ME, A = 511.0, 1/137.036
SH = {'a': lambda p, q: 1+0*p, 'u1': lambda p, q: p**2+q**2, 'u2': lambda p, q: p**4+10/3*p**2*q**2+q**4,
      'u3': lambda p, q: p**6+7*p**4*q**2+7*p**2*q**4+q**6}

def beta(E0, Z, s, n=3000):  # Fermi (non-rel) x shape; normalized per keV
    T = (np.arange(n)+.5)*E0/n; W = 1+T/ME; p = np.sqrt(W*W-1); q = (E0-T)/ME; e = 2*np.pi*A*Z*W/p
    P = p*W*q*q*e/(1-np.exp(-e))*SH[s](p, q); return T, P/np.trapezoid(P, T)

def phi(W, k):  # photons per unit k (mc^2) for e- total energy W
    Wp = np.maximum(W-k, 1+1e-12); p, pp = np.sqrt(W*W-1), np.sqrt(Wp*Wp-1)
    return np.where(W-k > 1, A/(np.pi*k)*pp/p*((W*W+Wp*Wp)/(W*pp)*np.log(Wp+pp)-2), 0)

def branches(nuc, shape):
    if nuc == 'K40':
        return [(0.8956, 1310.91, 20, shape)]
    elif nuc == 'Cs137':
        return [(0.9457, 513.97, 56, 'u1'), (0.0543, 1175.63, 56, 'u2')]

def load_branches(path):
    with open(path, encoding='utf-8') as f:
        d = json.load(f)
    br = []
    for b in d['branches']:
        fraction = float(b['fraction'])
        Q_keV = float(b['Q_keV'])
        Z_daughter = int(b['Z_daughter'])
        shape = str(b['shape'])
        if shape not in SH:
            sys.exit(f'Error: shape {shape} not in SH')
        br.append((fraction, Q_keV, Z_daughter, shape))
    if not br:
        sys.exit('Error: empty branches list')
    label = str(d['shape_label'])
    return br, label

def dndk(br, ks):
    result = np.zeros_like(ks)
    for f, Q, Z, s in br:
        T, P = beta(Q, Z, s)
        for i, k in enumerate(ks):
            if k < Q:
                integrand = P * phi(1 + T/ME, k/ME) / ME
                result[i] += f * np.trapezoid(integrand, T)
    return result

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('nuclide')
    parser.add_argument('out_csv')
    parser.add_argument('--branches', default=None)
    parser.add_argument('--shape', choices=['u3', 'a'], default='u3')
    parser.add_argument('--kmin', type=float, default=100.0)
    parser.add_argument('--band2', action='store_true')
    args = parser.parse_args()
    
    if args.branches is not None:
        br, shape_txt = load_branches(args.branches)
    else:
        if args.nuclide not in ('K40', 'Cs137'):
            parser.error('nuclide must be K40 or Cs137 unless --branches is given')
        br = branches(args.nuclide, args.shape)
        shape_txt = args.shape if args.nuclide == 'K40' else 'u1+u2'
    
    qmax = max(b[1] for b in br)
    
    fine = np.arange(args.kmin, qmax + 1e-9, 0.25)
    yf = dndk(br, fine)
    y_tot = np.trapezoid(yf, fine)
    
    m = (fine >= 150) & (fine <= 300)
    y_150 = np.trapezoid(yf[m], fine[m])
    
    y_300 = None
    if args.band2:
        m2 = (fine >= 300) & (fine <= 600)
        y_300 = np.trapezoid(yf[m2], fine[m2])
    
    kc = np.arange(args.kmin + 0.5, qmax, 1.0)
    yc = dndk(br, kc)
    
    header = [
        f'nuclide={args.nuclide}',
        f'shape={shape_txt}',
        f'branches(fraction,Q_keV,Z_daughter,shape)=' + '; '.join(f'{f},{Q},{Z},{s}' for f, Q, Z, s in br),
        f'kmin_keV={args.kmin}',
        f'Y_total_kmin={y_tot:.6e}',
        f'Y_150_300={y_150:.6e}',
        f'Y_300_600={y_300:.6e}' if args.band2 else None,
        'grid: k_keV = centre of 1 keV bin [k-0.5,k+0.5]; header Y by trapezoid on 0.25 keV grid',
        'source: formulas ib_calc.py (verifier 2026-09-15, KUB phi); branches LNHB'
    ]
    
    with open(args.out_csv, 'w', encoding='utf-8', newline='\n') as f:
        for h in header:
            if h is not None:
                f.write('# ' + h + '\n')
        f.write('k_keV,dNdk_per_decay_per_keV\n')
        for k, v in zip(kc, yc):
            f.write(f'{k:.1f},{v:.6e}\n')
    
    for h in header[:7 if args.band2 else 6]:
        print(h)
    print(f'rows={len(kc)} first_k={kc[0]:.1f} last_k={kc[-1]:.1f} -> {args.out_csv}')

if __name__ == '__main__':
    main()
