import sys, os, glob, tempfile
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stage2_io

sys.stdout.reconfigure(encoding="utf-8")

def merge_component(items):
    T = sum(item["T"] for item in items)
    total = np.zeros_like(items[0]["total"])
    cats = {}
    first_cats = items[0]["cats"]
    
    for pg, proc, vol in first_cats:
        cats[(pg, proc, vol)] = {
            "names": first_cats[(pg, proc, vol)]["names"],
            "bins": np.zeros_like(first_cats[(pg, proc, vol)]["bins"])
        }
    
    for item in items:
        T_f = item["T"]
        total_f = item["total"]
        cats_f = item["cats"]
        
        # Обновляем суммарные значения
        total[:, 0] += total_f[:, 0] * T_f / T
        total[:, 1] += total_f[:, 1] * T_f**2 / T**2
        total[:, 2] += total_f[:, 2] * T_f / T
        total[:, 3] += total_f[:, 3] * T_f**2 / T**2
        
        # Обновляем категории
        for pg, proc, vol in cats_f:
            if (pg, proc, vol) not in cats:
                cats[(pg, proc, vol)] = {
                    "names": cats_f[(pg, proc, vol)]["names"],
                    "bins": np.zeros_like(cats_f[(pg, proc, vol)]["bins"])
                }
            cats[(pg, proc, vol)]["bins"] += cats_f[(pg, proc, vol)]["bins"] * T_f / T   # скорости категорий — с весом T_f/T
    
    return T, total, cats

def merge_dir(directory, out_prefix, nbin_expected=None):
    meta_files = glob.glob(os.path.join(directory, "*_meta.txt"))
    if not meta_files:
        sys.exit(2)
    
    bases = []
    for mf in meta_files:
        base = os.path.basename(mf).rsplit("_meta.txt", 1)[0]
        bases.append(base)
    
    components = {}
    for base in bases:
        try:
            meta = stage2_io.read_meta(os.path.join(directory, f"{base}_meta.txt"))
            T_sim_s = float(meta["T_sim_s"])
            if T_sim_s <= 0:
                print(f"SKIP {base}: T_sim_s <= 0", file=sys.stderr)
                continue
            total = stage2_io.read_total(os.path.join(directory, f"{base}_total.csv"))
            cats = stage2_io.read_cat(os.path.join(directory, f"{base}_cat.csv"))
        except Exception as e:
            print(f"SKIP {base}: error reading files - {e}", file=sys.stderr)
            continue
        
        if nbin_expected is not None and len(total) != nbin_expected + 1:
            print(f"SKIP {base}: total rows = {len(total)}, expected {nbin_expected + 1}", file=sys.stderr)
            continue
            
        component = base.rsplit("_", 1)[0]
        if component not in components:
            components[component] = []
        components[component].append({
            "T": T_sim_s,
            "total": total,
            "cats": cats
        })
    
    if not components:
        sys.exit(2)
    
    result = {}
    all_total = None
    
    print("component  files  T_sim_s  rate_sum_s^-1  stat_err_s^-1")
    
    for comp, items in components.items():
        T, total, cats = merge_component(items)
        result[comp] = (T, total, cats)
        
        # Записываем файлы
        out_total_path = f"{out_prefix}_total_{comp}.csv"
        stage2_io.write_total(out_total_path, total, f"rates per second, merged T_sim_s={T}")
        
        if all_total is None:
            all_total = np.zeros_like(total)
        all_total += total
        
        # Считаем суммы
        rate_sum = np.sum(total[:-1, 0])
        stat_err = np.sqrt(np.sum(total[:-1, 1]))
        
        print(f"{comp:8}  {len(items):6}  {T:8.2f}  {rate_sum:12.6f}  {stat_err:12.6f}")
    
    # Записываем общий файл
    out_total_all_path = f"{out_prefix}_total_all.csv"
    stage2_io.write_total(out_total_all_path, all_total, "rates per second, sum over components")
    
    # Объединяем категории: сумма скоростей по компонентам (ключи — объединение)
    all_cats = {}
    for comp in components:
        for key, val in result[comp][2].items():
            if key not in all_cats:
                all_cats[key] = {"names": val["names"], "bins": np.zeros_like(val["bins"])}
            all_cats[key]["bins"] += val["bins"]

    out_cat_all_path = f"{out_prefix}_cat_all.csv"
    stage2_io.write_cat(out_cat_all_path, all_cats)
    
    # Добавляем итоговый компонент
    T_total = sum(T for T, _, _ in result.values())
    result["ALL"] = (T_total, all_total, all_cats)
    
    rate_sum_all = np.sum(all_total[:-1, 0])
    stat_err_all = np.sqrt(np.sum(all_total[:-1, 1]))
    print(f"{'ALL':8}  {sum(len(items) for items in components.values()):6}  {T_total:8.2f}  {rate_sum_all:12.6f}  {stat_err_all:12.6f}")
    
    return result

def _w_total(path, rows):
    with open(path, "w") as f:
        f.write("# test\nbin_keV,sumw,sumw2,light_sumw,light_sumw2\n")
        for i in range(6):
            r = rows.get(i, (0, 0, 0, 0))
            f.write("%d,%g,%g,%g,%g\n" % ((i,) + tuple(r)))

def _w_files(d, base, T, rows, cat):
    with open(os.path.join(d, base + "_meta.txt"), "w") as f:
        f.write("T_sim_s=%g\n" % T)
    _w_total(os.path.join(d, base + "_total.csv"), rows)
    with open(os.path.join(d, base + "_cat.csv"), "w") as f:
        f.write("pg,proc,vol,pgname,procname,volname,rate_total,bin_0,bin_1,bin_2,bin_3,bin_4,bin_5\n" + cat + "\n")

def selftest():
    try:
        with tempfile.TemporaryDirectory() as tmp:
            _w_files(tmp, "aa_0", 1, {0: (2, 4, 0, 0)}, "0,1,2,gamma,nCapture,Cargo,2,2,0,0,0,0,0")
            _w_files(tmp, "aa_1", 3, {0: (6, 12, 0, 0)}, "0,1,2,gamma,nCapture,Cargo,6,6,0,0,0,0,0")
            _w_files(tmp, "bb_0", 2, {1: (1, 0.5, 0, 0)}, "1,4,3,neutron,neutronInelastic,Floor,1,0,1,0,0,0,0\n0,1,2,gamma,nCapture,Cargo,1,1,0,0,0,0,0")   # тот же ключ (0,1,2) в другой компоненте: суммируется
            res = merge_dir(tmp, os.path.join(tmp, "out"), nbin_expected=5)
            checks = [
                ("aa rate", res["aa"][1][0, 0], 5.0), ("aa var", res["aa"][1][0, 1], 7.0),
                ("bb rate", res["bb"][1][1, 0], 1.0), ("bb var", res["bb"][1][1, 1], 0.5),
                ("ALL rate0", res["ALL"][1][0, 0], 5.0), ("ALL rate1", res["ALL"][1][1, 0], 1.0),
                ("ALL var0", res["ALL"][1][0, 1], 7.0), ("ALL var1", res["ALL"][1][1, 1], 0.5),
            ]
            tot = stage2_io.read_total(os.path.join(tmp, "out_total_all.csv"))
            cat = stage2_io.read_cat(os.path.join(tmp, "out_cat_all.csv"))
            if set(cat) != {(0, 1, 2), (1, 4, 3)}:
                print("SELFTEST FAIL: category keys", sorted(cat)); sys.exit(1)
            checks += [("file rate0", tot[0, 0], 5.0), ("file rate1", tot[1, 0], 1.0),
                       ("cat shared bin0", cat[(0, 1, 2)]["bins"][0], 5.0 + 1.0 * 2 / 2), ("cat bb bin1", cat[(1, 4, 3)]["bins"][1], 1.0)]
            for name, got, exp in checks:
                if abs(got - exp) > 1e-9:
                    print("SELFTEST FAIL: %s got %r expected %r" % (name, got, exp)); sys.exit(1)
    except SystemExit:
        raise
    except Exception as e:
        print("SELFTEST FAIL: %r" % (e,)); sys.exit(1)
    print("SELFTEST PASS")
    sys.exit(0)

if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--selftest":
        selftest()
    elif len(sys.argv) == 3:
        merge_dir(sys.argv[1], sys.argv[2])
    else:
        print("Usage: stage2_merge.py <dir> <out_prefix>")
        print("       stage2_merge.py --selftest")
        sys.exit(1)
