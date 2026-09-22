import sys, io
sys.stdout.reconfigure(encoding='utf-8')
P = {'neutron':0,'proton':1,'muplus':29,'muminus':30,'electron':31,'positron':32,'photon':33}
g = [l.split() for l in open('dcc/ICRP116.inp').read().splitlines()[2:]]
dcc = {k:[float(r[2+v]) for r in g] for k,v in P.items()}
emid = [float(r[0]) for r in g]; ewid = [float(r[1]) for r in g]
blocks = open('SpecOut/'+sys.argv[1]+'.out').read().split('New condition')[1:]
o = open('OUT/'+sys.argv[1]+'-integrals.csv','w', encoding='utf-8')
o.write('cond,header,particle,flux_total,flux_lt0.5eV,flux_0.5eV-0.1MeV,flux_gt1MeV,flux_gt10MeV,Edep_dose_uSvph\n')
prev = 0.0
for n, b in enumerate(blocks):
    L = [x for x in b.splitlines() if x.strip()]
    hdr = ' '.join(L[1].split()); tot = float(L[2].split()[-1])
    tab = [[float(x) for x in l.split()] for l in L[6:146]]
    for k, v in P.items():
        f = [r[1+v] for r in tab]
        b0 = 0 if k == 'neutron' else 60   # PARMA validity: EM/mu/p from bin 61 (11.3 keV)
        I = lambda lo, hi: sum(f[i]*ewid[i] for i in range(b0, 140) if lo <= emid[i] < hi)
        d = sum(f[i]*ewid[i]*dcc[k][i]*3.6e-3 for i in range(140))
        o.write('%d,"%s",%s,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e\n' % (n, hdr, k, I(0,1e9), I(0,5e-7), I(5e-7,0.1), I(1,1e9), I(10,1e9), d))
    print(n, hdr, 'dose_uSv/h(diff)=%.4f' % (tot-prev)); prev = tot
o.close()
