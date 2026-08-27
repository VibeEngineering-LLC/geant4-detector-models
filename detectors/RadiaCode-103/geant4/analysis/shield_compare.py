# -*- coding: utf-8 -*-
"""
Сравнение отклика RadiaCode-103 без и с свинцовым домиком.
Считает коэффициенты ослабления по энергетическим полосам, строит график.
Открытый верх — главная причина ограниченной эффективности защиты.
"""
import os, sys
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
import fit_two_criteria as ftc

RUN_DIR = os.path.join(HERE, "..", "run_field", "output")
OUT_DIR = os.path.join(HERE, "..", "verify")
os.makedirs(OUT_DIR, exist_ok=True)

PAIRS = [
    ("K-40", "est_K40_noshield.csv", "est_K40_shield.csv"),
    ("Tl-208", "est_Tl208_noshield.csv", "est_Tl208_shield.csv"),
]

BANDS = [(20,100),(100,300),(300,700),(700,1500),(1500,2400),(2400,2999)]

def load(fname):
    path = os.path.join(RUN_DIR, fname)
    if not os.path.exists(path):
        return None
    meta, cps, counts = ftc.read_template(path)
    return {"meta": meta, "cps": cps, "counts": counts}

def analyze(name, f_no, f_sh):
    d_no = load(f_no)
    d_sh = load(f_sh)
    if d_no is None or d_sh is None:
        print(f"Нет данных для {name}")
        return None

    meta_no = d_no["meta"]
    meta_sh = d_sh["meta"]
    cps_no = d_no["cps"]
    cps_sh = d_sh["cps"]
    counts_no = d_no["counts"]
    counts_sh = d_sh["counts"]

    print(f"\n{name}")
    print(f"Без домика: hits={meta_no['n_hits_in_crystal']}, cps_total={meta_no['cps_total']:.2f}")
    print(f"С домиком:  hits={meta_sh['n_hits_in_crystal']}, cps_total={meta_sh['cps_total']:.2f}")

    atten_total = meta_no["cps_total"] / meta_sh["cps_total"]
    print(f"Полное ослабление: {atten_total:.2f} раз")

    res = {"name": name, "no": cps_no, "sh": cps_sh, "bands": []}
    for i, (e1, e2) in enumerate(BANDS):
        # ОТНОШЕНИЕ СЧИТАЕМ ПО cps, НЕ по сырым counts: прогоны идут на разном
        # числе первичных (1e8 без домика против 3e8 с домиком), и сырые
        # отсчёты занизили бы ослабление ровно втрое. cps уже нормированы
        # на эквивалентное время t_run каждого прогона.
        c_no = np.sum(counts_no[e1:e2])   # только для статистики и надёжности
        c_sh = np.sum(counts_sh[e1:e2])
        s_no = np.sum(cps_no[e1:e2])
        s_sh = np.sum(cps_sh[e1:e2])
        ratio = s_no / s_sh if s_sh > 0 else np.inf
        rel = np.sqrt(1/c_no + 1/c_sh) if c_no > 0 and c_sh > 0 else np.inf
        rel_err = ratio * rel if np.isfinite(rel) else np.inf  # АБСОЛЮТНАЯ ошибка
        unreliable = c_sh < 20
        res["bands"].append({
            "e1": e1, "e2": e2, "c_no": c_no, "c_sh": c_sh,
            "ratio": ratio, "rel_err": rel_err, "unreliable": unreliable
        })
        if unreliable:
            print(f"  {e1}-{e2} кэВ: без={c_no:.0f}, с={c_sh:.0f}, отношение={ratio:.2f} (ненадёжно)")
        else:
            print(f"  {e1}-{e2} кэВ: без={c_no:.0f}, с={c_sh:.0f}, отношение={ratio:.2f} ± {rel_err:.2f}")
    return res

def plot(results):
    n_plots = len(results)
    fig = plt.figure(figsize=(13, 5.5 * n_plots))
    gs = fig.add_gridspec(n_plots, 2, hspace=0.45)

    for i, r in enumerate(results):
        name = r["name"]
        ax1 = fig.add_subplot(gs[i, 0])
        ax2 = fig.add_subplot(gs[i, 1])

        ax1.plot(r["no"], label="без домика", lw=1)
        ax1.plot(r["sh"], label="с домиком", lw=1)
        ax1.set_yscale("log")
        ax1.set_xlim(20, 2999)
        ax1.set_xlabel("Энергия, кэВ")
        ax1.set_ylabel("cps")
        # bbox справа выносил легенду ПОВЕРХ соседней панели со столбиками.
        # В правом верхнем углу спектра данных нет (кривые спадают влево-вниз).
        ax1.legend(loc="upper right", fontsize=8, framealpha=0.95)
        ax1.grid(True, alpha=0.3)

        ratios = [b["ratio"] for b in r["bands"]]
        errors = [b["rel_err"] for b in r["bands"]]
        labels = [f"{b['e1']}-{b['e2']}" for b in r["bands"]]
        x = np.arange(len(ratios))
        ax2.bar(x, ratios, yerr=errors, capsize=3, alpha=0.7)
        ax2.axhline(y=1.0, color="k", linestyle="--")
        ax2.set_xticks(x)
        ax2.set_xticklabels(labels, rotation=20)
        ax2.set_ylabel("во сколько раз ослаблено")
        ax2.grid(True, alpha=0.3)
        ax2.set_title(name)

    fig.suptitle("RadiaCode-103: влияние свинцового домика (Pb 50 мм, полость 150x150x385 мм, верх ОТКРЫТ)")
    path = os.path.join(OUT_DIR, "RC103_shield_effect.png")
    plt.savefig(path, bbox_inches="tight")
    plt.close()
    return path

def main():
    results = []
    for name, f_no, f_sh in PAIRS:
        r = analyze(name, f_no, f_sh)
        if r is not None:
            results.append(r)

    if not results:
        print("Прогоны ещё не готовы")
        raise SystemExit

    path = plot(results)
    print(f"График сохранён в {path}")

    # ТРЕБУЕТ ТОЛКОВАНИЯ
    print("\nТРЕБУЕТ ТОЛКОВАНИЯ:")
    trigger = False
    for r in results:
        atten_total = r["no"].sum() / r["sh"].sum()
        if atten_total < 2:
            print("  (а) полное ослабление меньше 2 — защита почти не работает")
            trigger = True

        # Сравнивать можно только НАДЁЖНЫЕ полосы: у полосы без отсчётов
        # ratio = inf, и она «больше» чего угодно — это давало ложный триггер.
        soft = [b["ratio"] for b in r["bands"] if b["e1"] == 20 and not b["unreliable"]]
        hard = [b["ratio"] for b in r["bands"] if b["e1"] == 1500 and not b["unreliable"]]
        if soft and hard and soft[0] < hard[0]:
            print("  (б) ослабление в полосе 20-100 кэВ меньше, чем в 1500-2400 — физически подозрительно")
            trigger = True

        for b in r["bands"]:
            if b["c_sh"] < 20:
                print(f"  (в) [{r['name']}] в полосе {b['e1']}-{b['e2']} с домиком меньше 20 отсчётов — статистически не обеспечен")
                trigger = True

    if not trigger:
        print("  (пусто — ни один триггер не сработал)")

if __name__ == "__main__":
    main()
