# -*- coding: utf-8 -*-
"""
Рисует результаты двух разложений (A — chi2/ndf, B — форма) для RadiaCode-103.
Данные берутся импортом из fit_two_criteria, повторной реализации подгонки нет.
Легенда вынесена за пределы поля данных — в этой линии уже был дефект,
когда непрозрачная легенда перекрыла часть спектра и график читался как обрыв данных.
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

matplotlib.rcParams["font.size"] = 9
matplotlib.rcParams["figure.dpi"] = 130
matplotlib.rcParams["axes.grid"] = True
matplotlib.rcParams["grid.alpha"] = 0.3

OUT_DIR = os.path.join(HERE, "..", "verify")
os.makedirs(OUT_DIR, exist_ok=True)
# Метка входных шаблонов в имени файла — чтобы прогон на комнатных шаблонах
# не затирал прежний PNG со сферическими (у обоих одна логика сохранения).
_TAG = os.environ.get("G4MODELS_PLOT_TAG", "")

COLORS = {
    "K40": "#1f77b4",
    "Ra226": "#ff7f0e",
    "Pb214": "#2ca02c",
    "Bi214": "#d62728",
    "Pb212": "#9467bd",
    "Ac228": "#8c564b",
    "Bi212": "#e377c2",
    "Tl208": "#17becf",
    "mu": "#7f7f7f",     # мюоны — единицы другие (мюон/с), не Бк/кг
}
UNITS = {"mu": "мюон/с"}

def prepare():
    meas_path = os.path.join(ftc.MEAS_DIR, ftc.MEAS_NAME)
    if not os.path.exists(meas_path):
        raise SystemExit(f"Файл измерения не найден: {meas_path}")
    
    smp = ftc.read_rcxml.read(meas_path)[0]
    cnt = smp.counts[:-1].astype(float)
    ch = np.arange(len(cnt))
    e_meas = np.asarray(sum(c * ch ** i for i, c in enumerate(ftc.CAL_ROOM)))
    live = float(smp.live)

    # cols — СПИСОК массивов (не словарь), сетка каждого шаблона своя: 1 кэВ/бин.
    names, cols, metas, varis = ftc.load_templates()
    A = np.zeros((len(e_meas), len(names)))
    for k in range(len(names)):
        c = cols[k]
        A[:, k] = ftc.fl.rebin_model_to_meas(np.arange(len(c)) + 0.5, c, e_meas)

    sel = (e_meas >= ftc.E_LO) & (e_meas < ftc.E_HI)
    A = A[sel]
    cnt = cnt[sel]
    e_meas = e_meas[sel]

    A_counts = A * live
    weights_a = 1.0 / np.sqrt(np.maximum(cnt, 1.0))
    weights_b = 1.0 / np.maximum(cnt, 1.0)

    # Сигнатура: fit(A_counts, meas_counts, weights, names, title, note)
    amp_a, err_a, pred_a, chi2_ndf_a, form_resid_a = ftc.fit(
        A_counts, cnt, weights_a, names, "A — критерий chi2/ndf", "(для графика)")
    amp_b, err_b, pred_b, chi2_ndf_b, form_resid_b = ftc.fit(
        A_counts, cnt, weights_b, names, "B — критерий невязки формы", "(для графика)")

    # Сырые (НЕсвёрнутые) шаблоны — нужны, чтобы показать линии рентгена до
    # размытия аппаратным разрешением.
    raw = {}
    for name in names:
        if name == "mu":
            continue          # мюонный шаблон лежит отдельно, не в TEMPLATE_DIR
        p = os.path.join(ftc.TEMPLATE_DIR, ftc.TEMPLATE_FMT % name)
        raw[name] = ftc.read_template(p)[1]

    return {
        "names": names,
        "raw": raw,
        "cols": cols,
        "e_sel": e_meas,
        "meas_counts": cnt,
        "A_counts": A_counts,
        "live": live,
        "amp_a": amp_a,
        "err_a": err_a,
        "pred_a": pred_a,
        "chi2_ndf_a": chi2_ndf_a,
        "form_resid_a": form_resid_a,
        "amp_b": amp_b,
        "err_b": err_b,
        "pred_b": pred_b,
        "chi2_ndf_b": chi2_ndf_b,
        "form_resid_b": form_resid_b
    }

def plot_decomposition(data, only_b=False):
    """only_b=True — отдельный файл только с критерием B (просьба оператора)."""
    suffix = "_B" if only_b else ""
    path = os.path.join(OUT_DIR, "RC103_bg_decomposition%s%s.png" % (_TAG, suffix))
    nrow = 2 if only_b else 3
    heights = [3, 2] if only_b else [3, 3, 2]
    fig = plt.figure(figsize=(13, 8.0 if only_b else 11.5))
    # hspace большой намеренно: при 0.25 заголовок нижней панели наезжал на
    # подпись оси X верхней (тот же класс дефекта, что P-003).
    gs = fig.add_gridspec(nrow, hspace=0.62, height_ratios=heights)
    
    meas_cps = data["meas_counts"] / data["live"]
    e = data["e_sel"]
    
    # Обе подгонки решают ОДНО И ТО ЖЕ матричное уравнение A*x = y (x >= 0,
    # неотрицательный МНК), различаясь ТОЛЬКО весовой матрицей W: min ||W(Ax-y)||.
    crit_a = (data["amp_a"], data["pred_a"], data["chi2_ndf_a"], data["form_resid_a"],
              "A: NNLS, min ||W(Ax−y)||,  W = diag(1/√N)  —  минимум χ² (пуассоновские веса)")
    crit_b = (data["amp_b"], data["pred_b"], data["chi2_ndf_b"], data["form_resid_b"],
              "B: NNLS, min ||W(Ax−y)||,  W = diag(1/N)  —  минимум относительной невязки формы")
    for i, (amp, pred, chi2_ndf, form_resid, name) in enumerate(
            [crit_b] if only_b else [crit_a, crit_b]):
        ax = fig.add_subplot(gs[i, 0])
        
        ax.step(e, meas_cps, where="mid", color="0.35", lw=1.0, label="измерение")
        ax.plot(e, pred / data["live"], color="crimson", lw=1.6, label="сумма модели")
        
        zero_names = []
        for k, name_k in enumerate(data["names"]):
            if amp[k] > 0:
                ax.plot(e, amp[k] * data["A_counts"][:, k] / data["live"], 
                        color=COLORS[name_k], lw=0.9, alpha=0.85, 
                        ls=("--" if name_k == "mu" else "-"),
                        label="%s %.1f %s" % (name_k, amp[k],
                                              UNITS.get(name_k, "Бк/кг")))
            else:
                zero_names.append(name_k)
        
        ax.set_yscale("log")
        ax.set_xlim(ftc.E_LO, ftc.E_HI)
        ax.set_ylim(max(meas_cps)*1e-5, max(meas_cps)*2)
        ax.set_xlabel("Энергия, кэВ")
        ax.set_ylabel("Скорость счёта, 1/(с·канал)")
        
        title = (f"Критерий {name}\n"
                 f"невязка формы = {form_resid:.4f},  χ²/ndf = {chi2_ndf:.1f}")
        if zero_names:
            title += f"\nобнулены NNLS: {', '.join(zero_names)}"
        ax.set_title(title, fontsize=10)
        
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8, framealpha=1.0)

    ax3 = fig.add_subplot(gs[nrow - 1, 0])
    
    band_indices = np.arange(len(ftc.BANDS))
    x_labels = [f"{b[0]}-{b[1]}" for b in ftc.BANDS]
    
    ratios_a = []
    ratios_b = []
    # Отношение по полосе = СУММА модели / СУММА измерения. Поканальное среднее
    # отношений — другая величина: при малых counts отдельные каналы дают выбросы.
    for lo, hi in ftc.BANDS:
        mask = (e >= lo) & (e < hi)
        denom = data["meas_counts"][mask].sum()
        if denom > 0:
            ratios_a.append(data["pred_a"][mask].sum() / denom)
            ratios_b.append(data["pred_b"][mask].sum() / denom)
        else:
            ratios_a.append(np.nan)
            ratios_b.append(np.nan)
    
    if not only_b:
        ax3.plot(band_indices, ratios_a, "o-", color="blue", lw=1.4,
                 label="A: W = diag(1/√N), минимум χ²")
    ax3.plot(band_indices, ratios_b, "s-", color="green", lw=1.4,
             label="B: W = diag(1/N), минимум невязки формы")
    
    ax3.axhline(1.0, color="0.5", ls="--", lw=1)
    ax3.set_xticks(band_indices)
    ax3.set_xticklabels(x_labels, rotation=20)
    ax3.set_ylabel("Модель / Измерение")
    ax3.set_ylim(0.4, 1.2)
    ax3.set_xlabel("Полосы энергии, кэВ")
    
    ax3.set_title("Отношение модель/измерение по полосам "
                  "(мюонный столбец включён)")
    # Внизу слева — единственная свободная зона: кривые идут около 1.0 и
    # падают только в последней полосе. В upper right легенда перекрывала точки.
    ax3.legend(loc="lower left", fontsize=8, framealpha=0.95)
    
    fig.suptitle("RadiaCode-103, фон помещения: нуклидное разложение по полным "
                 "спектрам (метод 1, новая GDML-модель)", y=0.995)
    # tight_layout НЕ вызываем: он пересчитывает отступы и съедает заданный
    # hspace, возвращая наложение заголовков.
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path

def plot_templates(data):
    path = os.path.join(OUT_DIR, "RC103_bg_templates%s.png" % _TAG)
    fig = plt.figure(figsize=(13, 7))
    gs = fig.add_gridspec(2, hspace=0.25)
    
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    
    # Шаблоны рисуем на ИХ СОБСТВЕННОЙ сетке (1 кэВ/бин), без ремешка: здесь
    # нужна форма шаблона как есть, а не проекция на каналы прибора.
    for k, name in enumerate(data["names"]):
        c = data["cols"][k]
        e_full = np.arange(len(c)) + 0.5
        norm = c.sum()
        if norm > 0:
            ax1.plot(e_full, c / norm, color=COLORS[name], lw=1.0,
                     alpha=0.85, label=name)

    ax1.set_yscale("log")
    ax1.set_ylim(1e-9, 5e-2)   # иначе автомасштаб уходит в 1e-12 и всё сжимается
    ax1.set_xlim(20, 2830)
    ax1.set_ylabel("Нормированная форма (доля на канал)")
    ax1.set_title("Формы континуума почти совпадают — именно поэтому поканальное разделение звеньев вырождено")
    ax1.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8, framealpha=1.0)
    
    for k, name in enumerate(data["names"]):
        c = data["cols"][k]
        e_full = np.arange(len(c)) + 0.5
        norm = c.sum()
        if norm > 0:
            ax2.plot(e_full, c / norm, color=COLORS[name], lw=1.0,
                     alpha=0.85, label=name)

    # Пунктир — те же шаблоны ДО свёртки с разрешением прибора (только два,
    # иначе каша). Смысл: линии рентгена в расчёте ЕСТЬ, но RC-103 их не
    # разрешает — FWHM на 30 кэВ около 40 %, поэтому это плечо, а не пик.
    for name in ("Pb214", "Bi214"):
        if name not in data["names"]:
            continue
        rw = data["raw"][name]
        if rw.sum() > 0:
            ax2.plot(np.arange(len(rw)) + 0.5, rw / rw.sum(), color=COLORS[name],
                     lw=0.9, ls=":", alpha=0.9, label=name + " до свёртки")

    ax2.set_yscale("log")
    ax2.set_ylim(1e-5, 5e-2)
    ax2.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8)
    ax2.set_xlim(20, 300)
    ax2.set_xlabel("Энергия, кэВ")
    ax2.set_ylabel("Нормированная форма (доля на канал)")
    ax2.set_title("Расчётные шаблоны звеньев (ионный источник + nucleusLimits, метод 1)")
    
    # Подписи особенностей. Текст ставим в координатах ОСЕЙ (0..1 по Y), иначе
    # при логарифмической шкале с плавающими пределами он уезжает за кадр.
    for x, txt, yfrac in ((30.0, "K-рентген Cs/I\n(рождается в кристалле)", 0.93),
                          (80.0, "K-рентген Pb/Bi\n(приходит из стены)", 0.72),
                          (225.0, "обратное рассеяние", 0.93)):
        ax2.axvline(x, color="black", ls="--", lw=0.8, alpha=0.6)
        ax2.annotate(txt, xy=(x, yfrac), xycoords=("data", "axes fraction"),
                     xytext=(6, 0), textcoords="offset points",
                     fontsize=8, ha="left", va="top",
                     bbox=dict(boxstyle="round,pad=0.25", fc="white",
                               ec="0.6", alpha=0.85))
    
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path

def main():
    data = prepare()
    path1 = plot_decomposition(data)
    path1b = plot_decomposition(data, only_b=True)
    path2 = plot_templates(data)
    print(os.path.abspath(path1))
    print(os.path.abspath(path1b))
    print(os.path.abspath(path2))

if __name__ == "__main__":
    main()
