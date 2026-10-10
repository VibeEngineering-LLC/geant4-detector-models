# #GS-68: свет (NPSM) против депозита на монолиниях, обе свёрнуты до ПШПВ карточки; рядом — изм/мод тория/KCl в тех же полосах (в ПШПВ от E).
import sys, json, numpy as np
sys.stdout.reconfigure(encoding="utf-8")
C = json.load(open("<REPOS>/geant4-detector-models/detectors/Atom-GS2020-PRO-MAX/dataset/detector_card.json", encoding="utf-8"))["fwhm"]
P = np.array(C["points_sl"]); fw = lambda e: 1.05 * np.exp(np.interp(np.log(e), np.log(P[:, 0]), np.log(P[:, 1])))
TH = "C:/g4work/gs2020/run_marinelli/out_v5_oisn10/fit_m1_bgw_tail0_fwscale_calshape_calsum.json"; KC = "C:/g4work/gs2020/kcl_1l_v4w85_83/fit_kcl_bgw.json"
D = sys.argv[1] if len(sys.argv) > 1 else "C:/g4work/gs2020/npsm68"; L = [float(v) for v in sys.argv[2:]] or [238.6, 583.2, 911.2, 1460.8, 2614.5]
B = [(-3, -2), (-2, -1.25), (-1.25, -0.5), (-0.5, 0.5), (0.5, 1.25), (1.25, 2), (2, 3)]
def fold(x, y, sig):  # гаусс с σ(x), отсчёты сохраняются
    G = np.exp(-0.5 * ((x[None, :] - x[:, None]) / sig[:, None]) ** 2); G /= G.sum(axis=1, keepdims=True); return G.T @ y
def bands(x, a, b, E):
    f = fw(E); return [a[(x >= E + l * f) & (x < E + h * f)].sum() / max(b[(x >= E + l * f) & (x < E + h * f)].sum(), 1e-9) for l, h in B]
print("полосы в ПШПВ от E:", B)
for E in L:
    t = open(f"{D}/line_E{E}_npsm.csv", encoding="utf-8").read().splitlines(); i = t.index("bin_keV,count_edep,count_light")
    a = np.array([[float(v) for v in r.split(",")] for r in t[i + 1:]]); x, ed, li = a[:, 0], a[:, 1], a[:, 2]
    r = (x > 0.3 * E) & (x < 1.05 * E); p0 = x[r][np.argmax(li[r])]   # вершина света (вес < 1 → шкала сжата)
    w = np.abs(x - p0) < 0.03 * p0; pl = (x[w] * li[w]).sum() / li[w].sum(); k = E / pl   # шкала света по пику
    si = np.sqrt((li[w] * (x[w] - pl) ** 2).sum() / li[w].sum()) * k; print(f'  вершина света {p0:.1f} кэВ-депозита')                      # σ_intr света в пике, кэВ
    lis = np.diff(np.interp(x + 0.5, (x + 0.5) * k, np.cumsum(li)), prepend=0)           # свет в шкале энергии
    sc = fw(x) / 2.3548; fe = fold(x, ed, sc); fl = fold(x, lis, np.sqrt(np.maximum(sc ** 2 - (si * x / E) ** 2, 1.0)))
    d = json.load(open(KC if E == 1460.8 else TH, encoding="utf-8")); me = np.array(d["e"])
    print(f"E {E}: k {k:.4f} σ_intr {si:.2f} кэВ ({100*si/E:.2f} %), ПШПВ карт. {fw(E):.1f}")
    print("  свет/депозит:", " ".join(f"{v:.3f}" for v in bands(x, fl, fe, E)))
    print("  изм/мод     :", " ".join(f"{v:.3f}" for v in bands(me, np.array(d["net"]), np.array(d["model"]), E)), "(KCl)" if E == 1460.8 else "(торий)")
