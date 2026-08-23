# -*- coding: utf-8 -*-
"""
RadiaCode-103: функция отклика, разложенная по каналам взаимодействия

Канал ставится в момент события по истории процессов (geometry/main.cc, enum Chan)
и из готового спектра не восстанавливается.

Чтение, свёртка и порядок каналов заимствованы у донора Nano16 импортом,
потому что правило приоритета у обоих приборов одно.

Своё здесь — ширина линии из rcspec, нормировка на квант поля и подписи.

Запуск:
python analysis/draw_channels.py [E1 E2 ...]
"""
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

_HERE = os.path.dirname(os.path.abspath(__file__))
_DONOR = os.path.normpath(os.path.join(
    _HERE, "..", "..", "AtomSpectra-Nano-16-PRO", "analysis"))
sys.path.insert(0, _DONOR)
import draw_channels as dc
sys.path.insert(0, _HERE)
import rcspec

# Подмена зависимости, а не копия формулы — dc.broaden зовёт dc.fwhm,
# и после подмены донорская свёртка работает с разрешением RadiaCode-103;
# коэффициенты живут в rcspec и уточняются там.
dc.fwhm = lambda e: float(rcspec.fwhm(max(e, 1.0)))

# bare: фон 7 дней снимался прибором БЕЗ сосуда (оператор, 21.08)
SRC = os.path.normpath(os.path.join(_HERE, "..", "results", "bare", "response"))
OUT = os.path.normpath(os.path.join(_HERE, "..", "results", "figures",
                                    "rc103_channels.png"))
E_MAX = 3200.0

def step_for(e0):
    """Возвращает шаг канала отображения в кэВ."""
    return max(5.0, round(dc.fwhm(e0) / 4.0 / 5.0) * 5.0)

def main():
    wanted = [float(x) for x in sys.argv[1:]] or [150.0, 662.0, 1461.0, 2614.0]
    
    files = [f for f in os.listdir(SRC) if f.endswith("_chan.csv")]
    if not files:
        raise SystemExit("в %s нет файлов *_chan.csv" % SRC)
    files.sort()
    
    index = {}
    for fname in files:
        path = os.path.join(SRC, fname)
        head, names, rows = dc.read_chan(path)
        e0 = float(head["E_prim_keV"])
        index[e0] = (head, names, rows)
    
    dc.E_MAX = E_MAX
    
    fig, axes = plt.subplots(len(wanted), 1, figsize=(12.6, 4.0 * len(wanted)), sharex=False)
    if len(wanted) == 1:
        axes = [axes]
    
    for ax, e_want in zip(axes, wanted):
        e0 = min(index, key=lambda x: abs(x - e_want))
        head, names, rows = index[e0]
        stamp = head.get("src_sha1", "?")
        
        # Вес на ОДИН квант, вошедший внутрь цилиндра поля;
        # телесного угла здесь нет, источник уже задан поверхностью,
        # охватывающей прибор.
        w = 1.0 / float(head["N_primaries"])
        
        step = step_for(e0)
        nch = int(E_MAX / step)
        cols = [(k + 0.5) * step for k in range(nch)]
        
        grand = sum(sum(v) for _, v in rows)
        
        # Полный отклик
        total = dc.broaden([(e, sum(v)) for e, v in rows], w, nch, step)
        ax.step(cols, total, where="mid", lw=1.6, color="#111111",
                label="полный отклик", zorder=5)
        
        # По каждому каналу
        for name, label, colour in dc.ORDER:
            if name not in names:
                continue
            j = names.index(name)
            raw = sum(v[j] for _, v in rows)
            if raw == 0:
                continue
            cur = dc.broaden([(e, v[j]) for e, v in rows], w, nch, step)
            ax.step(cols, cur, where="mid", lw=1.0, color=colour,
                    label="%s — %.1f %%" % (label, 100.0 * raw / grand))
        
        # Оформление
        ax.set_yscale("log")
        ax.set_xlim(0, min(E_MAX, e0 * 1.15))
        top = max(total)
        ax.set_ylim(top / 2.0e3, top * 40)
        ax.set_ylabel("вероятность на квант поля")
        
        # Локаторы — только ось X: ось Y логарифмическая, линейный
        # MultipleLocator на ней — дефект генерации, убран при вычитке
        ax.xaxis.set_major_locator(MultipleLocator(250))
        ax.xaxis.set_minor_locator(MultipleLocator(50))
        
        # Сетка
        ax.grid(True, which="major", alpha=0.26, lw=0.6)
        ax.grid(True, which="minor", axis="x", alpha=0.11, lw=0.4)
        
        # Легенда
        ax.legend(fontsize=7.8, loc="upper right", ncol=2, framealpha=0.94)
        
        # Заголовок
        ax.set_title("падающий квант %.0f кэВ, канал отображения %.0f кэВ, событий с сигналом %d" % (e0, step, grand), fontsize=10.5, pad=4)
    
    axes[-1].set_xlabel("Энерговыделение, кэВ")
    
    # Общий заголовок
    fig.suptitle("RadiaCode-103: функция отклика, разложенная по каналам взаимодействия\nгеометрия фона БЕЗ домика (поле комнаты, цилиндрическая поверхность, прибор без сосуда)", fontsize=12, y=0.995)
    
    # Подпись внизу
    fig.text(0.5, 0.004, "Свёрнуто с ПШПВ(E) прибора из rcspec.fwhm. Штамп исходников %s. Проценты в легенде — доля канала во ВСЕХ событиях с сигналом на этом узле. Правило приоритета каналов общее с Nano16." % stamp,
             fontsize=8.4, ha="center", color="#555555")
    
    # Отступы
    fig.subplots_adjust(left=0.075, right=0.985, top=0.925, bottom=0.062, hspace=0.16)
    
    # Сохранение
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print("записано: %s" % OUT)
    return 0

if __name__ == "__main__":
    sys.exit(main())
