#!/usr/bin/env python3
"""Шаг укладки и зазор WT-20 по красным кончикам (фото 6971217312.webp).
Метод: цветовая маска красного, связные компоненты, PCA каждого кончика
(большая ось = ось прутка, малая ширина = ⌀3,2 мм — масштаб), центры
проецируются на перпендикуляр к средней оси прутков; диффы = шаг.
Масштаб per-pair: средняя ширина двух соседних кончиков."""
import sys
sys.path.insert(0, '.venv/Lib/site-packages')
import numpy as np
from PIL import Image
from scipy import ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IMG = r'C:\Users\Дмитрий\Мой диск\Дозиметрия\Руководства\Атом Нано 16 про\WT-20\6971217312.webp'
DIA = 3.2  # мм, диаметр прутка = ширина красной метки

rgb = np.asarray(Image.open(IMG).convert('RGB'), dtype=float)
R, G, B = rgb[..., 0], rgb[..., 1], rgb[..., 2]
mask = (R > 110) & (R - G > 45) & (R - B > 45)
print("red pixels:", mask.sum())

lab, nlab = ndimage.label(mask)
sizes = ndimage.sum(mask, lab, range(1, nlab + 1))
keep = [i + 1 for i, s in enumerate(sizes) if s >= 0.3 * np.median(sizes[sizes > 50])]
print("components:", nlab, "-> kept:", len(keep))

tips = []
for li in keep:
    ys, xs = np.nonzero(lab == li)
    pts = np.stack([xs, ys], axis=1).astype(float)
    c = pts.mean(axis=0)
    d = pts - c
    cov = d.T @ d / len(d)
    evals, evecs = np.linalg.eigh(cov)
    ax_major = evecs[:, 1]          # ось прутка
    ax_minor = evecs[:, 0]          # поперёк
    t = d @ ax_major                # вдоль оси
    p = d @ ax_minor                # поперёк
    # ширина: медиана (max-min) поперёк по бинам вдоль оси
    nb = max(4, int((t.max() - t.min()) / 6))
    bins = np.linspace(t.min(), t.max(), nb + 1)
    ws = []
    for b0, b1 in zip(bins[:-1], bins[1:]):
        sel = (t >= b0) & (t < b1)
        if sel.sum() > 3:
            ws.append(p[sel].max() - p[sel].min())
    width = np.median(ws)
    tips.append(dict(c=c, axis=ax_major, width=width, n=len(d)))

# средняя ось прутков (знак нормируем: ось "вверх")
axes = np.array([t['axis'] * (1 if t['axis'][1] < 0 else -1) for t in tips])
mean_axis = axes.mean(axis=0)
mean_axis /= np.linalg.norm(mean_axis)
perp = np.array([-mean_axis[1], mean_axis[0]])
ang = np.degrees(np.arctan2(mean_axis[0], -mean_axis[1]))
print(f"mean rod axis tilt: {ang:.1f} deg from vertical")

# сортировка по проекции на перпендикуляр
proj = np.array([t['c'] @ perp for t in tips])
order = np.argsort(proj)
tips = [tips[i] for i in order]
proj = proj[order]
widths = np.array([t['width'] for t in tips])
print("tip widths px:", np.round(widths, 1))

gaps_px = np.diff(proj)
pair_w = (widths[:-1] + widths[1:]) / 2
pitch_mm = DIA * gaps_px / pair_w          # per-pair масштаб
clear_mm = pitch_mm - DIA                  # зазор между поверхностями
print("pitch px :", np.round(gaps_px, 1))
print("pitch mm :", np.round(pitch_mm, 2))
print("clear mm :", np.round(clear_mm, 2))
print(f"\nPITCH  median {np.median(pitch_mm):.2f}  mean {pitch_mm.mean():.2f}  range {pitch_mm.min():.2f}..{pitch_mm.max():.2f}")
print(f"CLEAR  median {np.median(clear_mm):.2f}  mean {clear_mm.mean():.2f}  range {clear_mm.min():.2f}..{clear_mm.max():.2f}")
print(f"scale global: {DIA/np.median(widths):.4f} mm/px")

fig, ax1 = plt.subplots(figsize=(8, 8))
ax1.imshow(rgb.astype(np.uint8))
for t in tips:
    c = t['c']
    ax1.plot(c[0], c[1], 'y+', ms=12, mew=2)
    a = t['axis'] * 40
    ax1.plot([c[0] - a[0], c[0] + a[0]], [c[1] - a[1], c[1] + a[1]], 'c-', lw=0.8)
ax1.set_title(f"tips={len(tips)}  pitch={np.median(pitch_mm):.2f} mm  clear={np.median(clear_mm):.2f} mm")
plt.tight_layout()
plt.savefig("figures/wt20_red_tips.png", dpi=120)
print("fig: figures/wt20_red_tips.png")
