"""D-018 (11.09.2026): выше последнего репера шкала = шкала файла, сдвинутая к нему;
ниже первого — квадратичная экстраполяция по реперам (форма файла там отвергнута χ²)."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mix_unfold_g1s as g1s
KIT = os.path.join(os.path.dirname(HERE), "reference", "lsrm", "raw_lsrm", "Work", "BG",
                   "Gamma-1S", "Spe - поверки", "Поверка 2016")
spec = g1s.read_lsrm_spe(os.path.join(KIT, "Маринелли", "Смесь_AmTiCsEu_Маринелли.spe"))
e_of, n_refs = g1s.recalibrate_energy(spec, verbose=False)
f = lambda c: float(spec.channel_to_energy(c))
rows = spec.extras["lsrm_peaks_table"]
ch_ref = lambda E: min(rows, key=lambda r: abs(r["energy_keV"] - E))["position_ch"]
lo, hi, bad = ch_ref(g1s.RECAL_REFS[0]), ch_ref(g1s.RECAL_REFS[-1]), 0
p = sorted((ch_ref(E), E) for E in g1s.RECAL_REFS)
q = np.polyfit([x[0] for x in p], [x[1] for x in p], 2)
cases = [(c, hi, f(c) - f(hi)) for c in (hi + 50, hi + 300, spec.n_channels - 1)]
cases += [(c, lo, float(np.polyval(q, c) - np.polyval(q, lo))) for c in (lo - 3, 2)]  # низ: квадратичная, сдвинутая
for c0 in (lo, hi):   # P-029: непрерывность на крайних реперах (пара c0−ε / c0)
    jump = abs(e_of(c0 - 1e-6) - e_of(c0)) + abs(e_of(c0 + 1e-6) - e_of(c0))
    bad += jump > 1e-3
    print(f"непрерывность на кан {c0:8.3f}: скачок {jump:.6f} кэВ {'OK' if jump <= 1e-3 else 'FAIL'}")
for c, c0, want in cases:
    got = e_of(c) - e_of(c0)
    ok = abs(got - want) < 1e-6
    bad += not ok
    print(f"кан {c:8.2f}: сдвиг от репера {got:9.3f}, ожидание {want:9.3f} {'OK' if ok else 'FAIL'}")
print(f"E(кан {hi:.3f}) = {e_of(hi):.3f} (репер {g1s.RECAL_REFS[-1]})")
# #AMT-3: внешние низкие реперы — точное попадание, ниже первого — наклон первого отрезка, монотонность.
XR = [(12.30, 31.69), (14.62, 41.08)]
ex, n2 = g1s.recalibrate_energy(spec, verbose=False, extra_refs=XR)
k0 = (41.08 - 31.69) / (14.62 - 12.30)
for c, want in [(12.30, 31.69), (14.62, 41.08), (10.0, 31.69 - k0 * 2.30), (hi, g1s.RECAL_REFS[-1])]:
    ok = abs(ex(c) - want) < 1e-6
    bad += not ok
    print(f"внешние реперы: кан {c:7.2f} → {ex(c):8.3f}, ожидание {want:8.3f} {'OK' if ok else 'FAIL'}")
mono = all(ex(c + 0.5) > ex(c) for c in np.arange(0, spec.n_channels - 1, 0.5))
bad += not mono
print(f"внешние реперы: монотонность {'OK' if mono else 'FAIL'}, реперов {n2}")
print("RECAL_EXTRAP_FAILURES =", bad)
sys.exit(1 if bad else 0)
