# #GS-60 (07.10, #CFG-1 «да»): задания внешнего поля для KCl 1 л и Th (вода — готовые #GS-44). Поля: GDML E N OUT SEED HEMI.
E = [30, 40, 50, 60, 80, 100, 120, 150, 200, 250, 300, 400, 500, 609.312, 800, 1000, 1250, 1460.822, 1764.494, 2000, 2614.511]
G = [("kcl1l", "C:/g4work/gs2020/kcl/GS2020_marinelli_kcl_1000ml_v4_w85_83.gdml"), ("th", "C:/g4work/gs2020/gs20_matrix/mx_oisn10.gdml")]
s, L = 79001, []
for tag, g in G:
    for e in E:
        for h in ("up", "down"):
            L.append(f"{g} {e} 50000000 C:/g4work/gs2020/ext/chunks/ext_{tag}_E{e}_{h}.csv {s} {h}"); s += 1
open("C:/g4work/gs2020/ext/jobs_ext_gs60.txt", "w", newline="\n").write("\n".join(L) + "\n")
print(len(L), L[0], L[-1], sep="\n")
