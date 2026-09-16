"""Низ шкалы пробы (<59,5 кэВ) по донорам комплекта: переход каналов донор→проба линеен (кристалл один),
находится по линиям выше 59,5 кэВ через шкалу пробы; низкие пики донора переносятся в каналы пробы."""
import sys, os, numpy as np
A = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, A + "/analysis")
import mix_unfold_g1s as g
K = A + "/reference/lsrm/raw_lsrm/Work/BG/Gamma-1S/Spe - поверки/Поверка 2016/"
smp = g.read_lsrm_spe(K + "Маринелли/Смесь_AmTiCsEu_Маринелли.spe"); e_of, _ = g.recalibrate_energy(smp, verbose=False)
CH = np.arange(0.0, smp.n_channels - 1, 0.01); EE = np.array([e_of(c) for c in CH])
# Только разрешённые в NaI линии: 1085/1112, 411/444, 867 цеплялись за чужие пики (СКО перехода 0,9 кан).
# 964/1408 сняты: у обоих Eu-доноров одинаковые невязки (−1,7…−2,1 / +1,0…+1,1 кан) — в пробе 964
# лежит на комптоновском крае Sc-44 1157 (~944 кэВ). Для переноса вниз важны нижние линии.
EU = [121.78, 244.70, 344.28, 778.90]
BA = [80.90, 302.85, 356.01, 383.85]   # 80,90 — взвешенный дублет 79,61/81,00
DON = [("Точка 5см/Eu-152 #SRC-07_Точечная-5см_5cm.spe", EU), ("Точка 25см/Eu-152 #SRC-07_Точечная-25см_25cm.spe", EU),
       ("Точка 5см/Ba-133 #SRC-07_Точечная-5см_5cm.spe", BA), ("Точка 25см/Ba-133 #SRC-07_Точечная-25см_25cm.spe", BA)]
TR = {}
for f, hi in DON:
    rows = g.read_lsrm_spe(K + f).extras["lsrm_peaks_table"]
    near = lambda E0: min(rows, key=lambda r: abs(r["energy_keV"] - E0))
    pr = [(near(E0)["position_ch"], float(np.interp(E0, EE, CH)), E0) for E0 in hi if abs(near(E0)["energy_keV"] - E0) < max(6, 0.04 * E0)]
    x, y = np.array([p[0] for p in pr]), np.array([p[1] for p in pr]); a = np.polyfit(x, y, 1)
    rk = np.polyval(a, x) - y; TR[f] = a
    print(f"## {f}: линий {len(pr)}, переход кан_пробы = {a[1]:.3f} + {a[0]:.5f}·кан_донора, СКО {np.sqrt(np.mean(rk**2)):.3f} кан (макс {np.abs(rk).max():.3f})")
    print("   невязки, кан: " + "  ".join(f"{p[2]:.1f}:{r:+.3f}" for p, r in zip(pr, rk)))
    for r in sorted((r for r in rows if r["energy_keV"] < 60), key=lambda r: -r["area"]):
        cs = float(np.polyval(a, r["position_ch"]))
        print(f"   пик донора кан {r['position_ch']:6.2f} ±{r.get('d_position_ch') or 0:.2f}, площадь {r['area']:8.0f} → кан пробы {cs:6.2f} → по нынешней шкале пробы {e_of(cs):6.2f} кэВ")

# K-рентген ниже порога поиска пиков прибора (кан < 20) — своя подгонка: гауссиана + линейный фон по всему
# комплексу Kα+Kβ; энергия комплекса — средняя по интенсивностям IAEA (без строки «KB» = KpB1 + KpB2).
from scipy.optimize import curve_fit
KX = {"Eu-152": [(39.522, 20.87), (40.117, 37.8), (45.523, 11.81), (46.575, 3.05)],     # Sm K
      "Ba-133": [(30.625, 32.6), (30.973, 60.2), (35.089, 17.6), (35.818, 4.29)]}      # Cs K
def centroid(y, lo, hi):
    x = np.arange(lo, hi + 1, dtype=float); yy = np.asarray(y, float)[lo:hi + 1]
    f = lambda x, A, c, s, b0, b1: A * np.exp(-0.5 * ((x - c) / s) ** 2) + b0 + b1 * (x - c)
    p, cov = curve_fit(f, x, yy, p0=[yy.max() - yy.min(), x[np.argmax(yy)], 1.5, yy.min(), 0.0],
                       sigma=np.sqrt(np.maximum(yy, 1)), maxfev=20000)
    return p[1], float(np.sqrt(cov[1, 1])), abs(p[2])
print("\n## K-рентген доноров (своя подгонка)")
for f, a in TR.items():
    nk = "Eu-152" if "Eu-152" in f else "Ba-133"; Ek = sum(e * i for e, i in KX[nk]) / sum(i for _, i in KX[nk])
    y = g.read_lsrm_spe(K + f).counts; c0 = 8 + int(np.argmax(np.asarray(y[8:24], float)))
    for w in (3, 4):
        c, dc, s = centroid(y, c0 - w, c0 + w); cs = float(np.polyval(a, c))
        print(f"   {f.split('/')[1][:24]}: макс кан {c0}, окно ±{w}: центр {c:6.2f}±{dc:.2f} (σ {s:.2f}) → кан пробы {cs:6.2f}; "
              f"нынешняя шкала {e_of(cs):6.2f} кэВ против {Ek:.2f} (разница {e_of(cs) - Ek:+.2f})")
