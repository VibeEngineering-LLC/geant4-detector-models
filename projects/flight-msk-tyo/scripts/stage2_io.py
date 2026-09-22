import numpy as np

def read_meta(path):
    """Читает метаданные из файла, каждая строка key=value."""
    res = {}
    with open(path, "r") as f:
        for line in f:
            if "=" in line:
                k, v = line.split("=", 1)
                res[k.strip()] = v.strip()
    return res

def read_total(path):
    """Читает данные суммирования из CSV."""
    data = np.loadtxt(path, delimiter=",", skiprows=2, ndmin=2)
    return data[:, 1:5]

def write_total(path, arr, note):
    """Записывает данные суммирования в CSV."""
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(f"# {note}\n")
        f.write("bin_keV,sumw,sumw2,light_sumw,light_sumw2\n")
        for i, row in enumerate(arr):
            f.write(f"{i},{','.join('%.10g' % x for x in row)}\n")

def read_cat(path):
    """Читает категории из CSV."""
    res = {}
    with open(path, "r") as f:
        next(f)  # пропустить заголовок
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split(",")
            pg, proc, vol = map(int, parts[:3])
            pgname, procname, volname = parts[3:6]
            rate_total = float(parts[6])
            bins = np.array(parts[7:], dtype=float)
            key = (pg, proc, vol)
            if key in res:
                res[key]["bins"] += bins
            else:
                res[key] = {"names": (pgname, procname, volname), "bins": bins}
    return res

def write_cat(path, cats):
    """Записывает категории в CSV."""
    if not cats:
        return
    n_bins = len(next(iter(cats.values()))["bins"])
    header = ["pg", "proc", "vol", "pgname", "procname", "volname", "rate_total"] + [f"bin_{i}" for i in range(n_bins)]
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(",".join(header) + "\n")
        sorted_cats = sorted(cats.items(), key=lambda x: -sum(x[1]["bins"][:-1]))
        for (pg, proc, vol), data in sorted_cats:
            pgname, procname, volname = data["names"]
            bins = data["bins"]
            rate_total = sum(bins[:-1])
            head = "%d,%d,%d,%s,%s,%s,%.10g" % (pg, proc, vol, pgname, procname, volname, rate_total)
            f.write(head + "," + ",".join("%.10g" % x for x in bins) + "\n")
