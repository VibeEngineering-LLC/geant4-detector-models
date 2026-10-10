# #GS-60: чувствительность нетто мягкой зоны к сдвигу шкалы фона (аналог П224 BecqMoni «нуль фона против пробы»).
# Фон сдвигается по энергии на d кэВ (интерполяция плотности отсчётов), печатается изм − фон по полосам.
import json, sys, numpy as np
sys.stdout.reconfigure(encoding="utf-8")
D = r"<WORKDIR>\GEANT4\web\gs2020-th232-page\dist"
BANDS = [(15, 35), (35, 50), (50, 70), (70, 100), (100, 150), (15, 150)]
for fn, pre in (("data.js", "window.G1S="), ("data-k40.js", "window.GS_K40=")):
    s = json.loads(open(D + "\\" + fn, encoding="utf-8").read().strip()[len(pre):].rstrip(";"))["spectrum"]
    E = np.array(s["e_of_ch"], float); c = np.array(s["counts"], float); bg = np.array(s["bg_counts"], float)
    w = np.gradient(E); dens = bg / w
    print(f"== {fn}: изм − фон(E − d), отсчёты")
    print("  d кэВ " + "".join(f"{a:>7}-{b:<5}" for a, b in BANDS))
    for d in (-6, -4, -2, 0, 2, 4, 6):
        b2 = np.interp(E - d, E, dens) * w
        print(f"  {d:+4d}  " + "".join(f"{(c[k] - b2[k]).sum():13.0f}" for k in [((E >= a) & (E < b)) for a, b in BANDS]))
