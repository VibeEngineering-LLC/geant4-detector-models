#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Детерминированный обмер X-ray снимков RC-103 (пиксели → мм).

НИКАКОЙ LLM: только пороги, связные компоненты и bbox по numpy. Результат (JSON)
далее уходит в Ollama-хелпер `scripts/ollama/h2_xray_compare.py` для арифметики
сравнения с координатами Blender-сцены.

Калибровка
----------
Фидуциальной линейки в кадре нет, конус-бим даёт РАЗНЫЙ масштаб по осям, поэтому
масштаб берётся отдельно по каждой оси из габарита корпуса (чертёж оператора:
123.0 x 34.0 x 17.5 мм, см. `scripts/build_case_stl.py:18-20`). Измеряются, по
сути, ОТНОСИТЕЛЬНЫЕ положения внутренностей в долях габарита — этого достаточно
для сверки «капсула в U-вырезе / где LCD / где разъём».

Снимки (`RC103/Фото/`)
---------------------
  photo_2026-08-24_13-16-52 (4).jpg — вид сверху, портрет, фон чёрный  [основной]
  photo_2026-08-24_13-16-52 (3).jpg — вид сверху, пейзаж, фон светлый  [контроль]
  photo_2026-08-24_13-16-52 (5).jpg — вид с ребра, портрет, фон чёрный [Z-стек]

Ориентация: в (4) верх кадра = −X устройства (детектор), низ = +X (USB);
в (3) левый край = −X, правый = +X. Ось Y устройства — поперёк кадра.

Запуск: set PYTHONIOENCODING=utf-8 && python xray_measure.py [--png]
Выход:  JSON в stdout + ../verify/xray_measure.json (+ ../verify/xray_annotated.png)
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from PIL import Image

GEANT4_ROOT = Path(__file__).resolve().parents[1]
RC103_ROOT = GEANT4_ROOT.parent
PHOTO = RC103_ROOT / "Фото"
OUT_JSON = GEANT4_ROOT / "verify" / "xray_measure.json"
OUT_PNG = GEANT4_ROOT / "verify" / "xray_annotated.png"

CASE_L, CASE_W, CASE_H = 123.0, 34.0, 17.5

IMG_TOP_DARK = "photo_2026-08-24_13-16-52 (4).jpg"
IMG_TOP_LIGHT = "photo_2026-08-24_13-16-52 (3).jpg"
IMG_SIDE = "photo_2026-08-24_13-16-52 (5).jpg"


# --------------------------------------------------------------------------- #
# примитивы
# --------------------------------------------------------------------------- #
def gray(name: str) -> np.ndarray:
    return np.asarray(Image.open(PHOTO / name).convert("L"), dtype=np.float64) / 255.0


def bbox(mask: np.ndarray) -> tuple[int, int, int, int] | None:
    ys, xs = np.nonzero(mask)
    if ys.size == 0:
        return None
    return int(xs.min()), int(xs.max()), int(ys.min()), int(ys.max())


def largest_blob(mask: np.ndarray) -> np.ndarray:
    """Крупнейшая 4-связная компонента (без scipy)."""
    seen = np.zeros(mask.shape, dtype=bool)
    best: np.ndarray = np.zeros_like(mask)
    best_size = 0
    ys, xs = np.nonzero(mask)
    for y0, x0 in zip(ys, xs):
        if seen[y0, x0]:
            continue
        stack = [(int(y0), int(x0))]
        cells: list[tuple[int, int]] = []
        while stack:
            y, x = stack.pop()
            if y < 0 or x < 0 or y >= mask.shape[0] or x >= mask.shape[1]:
                continue
            if seen[y, x] or not mask[y, x]:
                continue
            seen[y, x] = True
            cells.append((y, x))
            stack.extend(((y + 1, x), (y - 1, x), (y, x + 1), (y, x - 1)))
        if len(cells) > best_size:
            best_size = len(cells)
            best = np.zeros_like(mask)
            for y, x in cells:
                best[y, x] = True
    return best


class Frame:
    """Кадр вида сверху: связывает пиксели с координатами устройства (мм).

    along  — ось кадра, вдоль которой лежит X устройства ('rows' или 'cols')
    flip   — True, если рост индекса вдоль `along` идёт в −X устройства
    """

    def __init__(self, img: np.ndarray, case: tuple[int, int, int, int],
                 along: str, flip_along: bool, flip_across: bool = False) -> None:
        self.img = img
        self.x0, self.x1, self.y0, self.y1 = case
        self.along = along
        self.flip_along = flip_along
        self.flip_across = flip_across
        n_rows = self.y1 - self.y0 + 1
        n_cols = self.x1 - self.x0 + 1
        if along == "rows":
            self.s_along = CASE_L / n_rows
            self.s_across = CASE_W / n_cols
        else:
            self.s_along = CASE_L / n_cols
            self.s_across = CASE_W / n_rows

    def _map(self, i: int, origin: int, scale: float, full: float, flip: bool) -> float:
        d = (i - origin) * scale
        return (full / 2 - d) if flip else (-full / 2 + d)

    def to_device(self, px: tuple[int, int, int, int]) -> dict:
        """px = (col_min, col_max, row_min, row_max) → мм в системе устройства."""
        cmin, cmax, rmin, rmax = px
        if self.along == "rows":
            a_lo = self._map(rmin, self.y0, self.s_along, CASE_L, self.flip_along)
            a_hi = self._map(rmax + 1, self.y0, self.s_along, CASE_L, self.flip_along)
            b_lo = self._map(cmin, self.x0, self.s_across, CASE_W, self.flip_across)
            b_hi = self._map(cmax + 1, self.x0, self.s_across, CASE_W, self.flip_across)
        else:
            a_lo = self._map(cmin, self.x0, self.s_along, CASE_L, self.flip_along)
            a_hi = self._map(cmax + 1, self.x0, self.s_along, CASE_L, self.flip_along)
            b_lo = self._map(rmin, self.y0, self.s_across, CASE_W, self.flip_across)
            b_hi = self._map(rmax + 1, self.y0, self.s_across, CASE_W, self.flip_across)
        xs = sorted((a_lo, a_hi))
        ys = sorted((b_lo, b_hi))
        return {
            "x_mm": [round(xs[0], 2), round(xs[1], 2)],
            "y_mm": [round(ys[0], 2), round(ys[1], 2)],
            "size_x_mm": round(xs[1] - xs[0], 2),
            "size_y_mm": round(ys[1] - ys[0], 2),
            "center_x_mm": round((xs[0] + xs[1]) / 2, 2),
            "center_y_mm": round((ys[0] + ys[1]) / 2, 2),
            "px": [cmin, cmax, rmin, rmax],
        }

    def roi_mask(self, along_frac: tuple[float, float],
                 across_frac: tuple[float, float] = (0.0, 1.0)) -> np.ndarray:
        """Маска ROI, заданная долями габарита корпуса по каждой оси кадра."""
        m = np.zeros(self.img.shape, dtype=bool)
        if self.along == "rows":
            n_a, o_a = self.y1 - self.y0 + 1, self.y0
            n_b, o_b = self.x1 - self.x0 + 1, self.x0
            a0, a1 = o_a + int(n_a * along_frac[0]), o_a + int(n_a * along_frac[1])
            b0, b1 = o_b + int(n_b * across_frac[0]), o_b + int(n_b * across_frac[1])
            m[a0:a1, b0:b1] = True
        else:
            n_a, o_a = self.x1 - self.x0 + 1, self.x0
            n_b, o_b = self.y1 - self.y0 + 1, self.y0
            a0, a1 = o_a + int(n_a * along_frac[0]), o_a + int(n_a * along_frac[1])
            b0, b1 = o_b + int(n_b * across_frac[0]), o_b + int(n_b * across_frac[1])
            m[b0:b1, a0:a1] = True
        return m

    def feature(self, roi: np.ndarray, *, brighter_than: float | None = None,
                darker_than: float | None = None) -> dict | None:
        if brighter_than is not None:
            mask = (self.img > brighter_than) & roi
        elif darker_than is not None:
            mask = (self.img < darker_than) & roi
        else:
            raise ValueError("нужен один из порогов")
        blob = largest_blob(mask)
        b = bbox(blob)
        if b is None:
            return None
        out = self.to_device(b)
        out["npx"] = int(blob.sum())
        return out

    def meta(self) -> dict:
        return {
            "case_px": {"col": [self.x0, self.x1], "row": [self.y0, self.y1]},
            "scale_mm_per_px": {"along_X": round(self.s_along, 5),
                                "across_Y": round(self.s_across, 5)},
            "anisotropy_pct": round(100 * (self.s_across / self.s_along - 1), 1),
            "orientation": f"X устройства вдоль {self.along}, flip={self.flip_along}",
        }


# --------------------------------------------------------------------------- #
# снимки
# --------------------------------------------------------------------------- #
def measure_top_dark() -> dict:
    """(4): портрет, фон чёрный (mean 0.017 / max 0.024 в углу 20x20) → порог 0.05."""
    g = gray(IMG_TOP_DARK)
    case = bbox(largest_blob(g > 0.05))
    if case is None:
        raise RuntimeError(f"{IMG_TOP_DARK}: силуэт корпуса не найден")
    fr = Frame(g, case, along="rows", flip_along=False)
    res = {"image": IMG_TOP_DARK, "px_size": list(g.shape[::-1]), **fr.meta(), "features": {}}
    # кристалл/капсула: самое плотное пятно на конце −X (первые 22 % длины)
    f = fr.feature(fr.roi_mask((0.00, 0.22)), brighter_than=0.70)
    if f:
        res["features"]["detector_dense"] = f
    # плата: средняя плотность по всей длине (ограничена по Y стенками корпуса)
    f = fr.feature(fr.roi_mask((0.02, 0.98), (0.03, 0.97)), brighter_than=0.20)
    if f:
        res["features"]["pcb_plus_parts"] = f
    # разъём USB: плотное пятно на последних 8 % длины
    f = fr.feature(fr.roi_mask((0.92, 1.00), (0.20, 0.80)), brighter_than=0.55)
    if f:
        res["features"]["usb"] = f
    return res


def measure_top_light() -> dict:
    """(3): пейзаж, фон ≈0.80; корпус даёт отклонение в обе стороны → |g−bg|>0.10."""
    g = gray(IMG_TOP_LIGHT)
    h = g.shape[0]
    band = g.copy()
    band[int(h * 0.80) :, :] = float(np.median(g[:20, :20]))  # гасим подставку
    bg = float(np.median(g[:20, :20]))
    case = bbox(largest_blob(np.abs(band - bg) > 0.10))
    if case is None:
        raise RuntimeError(f"{IMG_TOP_LIGHT}: силуэт корпуса не найден")
    fr = Frame(band, case, along="cols", flip_along=False)
    res = {"image": IMG_TOP_LIGHT, "px_size": list(g.shape[::-1]), "bg_level": round(bg, 4),
           **fr.meta(), "features": {}}
    # капсула: почти чёрный квадрат в первых 25 % длины
    f = fr.feature(fr.roi_mask((0.00, 0.25)), darker_than=0.30)
    if f:
        res["features"]["capsule_dark"] = f
    # окно LCD / рамка дисплея: тёмная полоса в 20…45 % длины, верхняя половина Y
    f = fr.feature(fr.roi_mask((0.18, 0.48), (0.05, 0.60)), darker_than=0.42)
    if f:
        res["features"]["lcd_frame"] = f
    # разъём USB: плотный тёмный блок в последних 12 % длины
    f = fr.feature(fr.roi_mask((0.88, 1.00), (0.15, 0.85)), darker_than=0.55)
    if f:
        res["features"]["usb"] = f
    return res


def measure_side() -> dict:
    """(5): вид с ребра. Прибор в кадре наклонён — абсолютный X ненадёжен, читаем Z."""
    g = gray(IMG_SIDE)
    h = g.shape[0]
    band = g.copy()
    band[int(h * 0.90) :, :] = 0.0
    case = bbox(largest_blob(band > 0.05))
    if case is None:
        raise RuntimeError(f"{IMG_SIDE}: силуэт не найден")
    x0, x1, y0, y1 = case
    s_thick = CASE_H / (x1 - x0 + 1)
    s_along = CASE_L / (y1 - y0 + 1)
    roi = np.zeros(band.shape, dtype=bool)
    roi[y0 : y1 + 1, x0 : x1 + 1] = True
    dense = bbox(largest_blob((band > 0.78) & roi))
    out = {
        "image": IMG_SIDE,
        "px_size": list(g.shape[::-1]),
        "case_px": {"col": [x0, x1], "row": [y0, y1]},
        "scale_mm_per_px": {"thickness_Z": round(s_thick, 5), "along_X": round(s_along, 5)},
        "caveat": "прибор в кадре наклонён; Z читаем, абсолютный X — нет",
        "features": {},
    }
    if dense:
        dx0, dx1, dy0, dy1 = dense
        out["features"]["dense_stack"] = {
            "z_mm": [round((dx0 - x0) * s_thick - CASE_H / 2, 2),
                     round((dx1 + 1 - x0) * s_thick - CASE_H / 2, 2)],
            "thickness_mm": round((dx1 - dx0 + 1) * s_thick, 2),
            "px": [dx0, dx1, dy0, dy1],
        }
    return out


# --------------------------------------------------------------------------- #
def annotate(res: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    fig, axes = plt.subplots(1, 3, figsize=(16, 8))
    for ax, key, img, boost in (
        (axes[0], "top_dark", IMG_TOP_DARK, 3.2),
        (axes[1], "top_light", IMG_TOP_LIGHT, 1.0),
        (axes[2], "side", IMG_SIDE, 3.2),
    ):
        g = gray(img)
        ax.imshow(np.clip(g * boost, 0, 1), cmap="gray")
        d = res[key]
        c = d["case_px"]
        ax.add_patch(Rectangle((c["col"][0], c["row"][0]),
                               c["col"][1] - c["col"][0], c["row"][1] - c["row"][0],
                               fill=False, ec="#e74c3c", lw=1.6))
        for name, f in d.get("features", {}).items():
            p = f["px"]
            ax.add_patch(Rectangle((p[0], p[2]), p[1] - p[0], p[3] - p[2],
                                   fill=False, ec="#f1c40f", lw=1.2))
            lbl = name
            if "center_x_mm" in f:
                lbl += f"  X={f['center_x_mm']}"
            ax.text(p[0], p[2] - 5, lbl, color="#f1c40f", fontsize=6)
        ax.set_title(f"{key}\n{img}", fontsize=7)
        ax.axis("off")
    fig.suptitle("RC-103 X-ray: обмер корпуса и внутренностей (красный — габарит-линейка)",
                 fontsize=11)
    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=140, bbox_inches="tight")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--png", action="store_true")
    args = ap.parse_args()

    missing = [n for n in (IMG_TOP_DARK, IMG_TOP_LIGHT, IMG_SIDE) if not (PHOTO / n).exists()]
    if missing:
        print(json.dumps({"error": "нет файлов", "missing": missing}, ensure_ascii=False), flush=True)
        return 1

    res = {
        "calibration": {
            "case_L_mm": CASE_L, "case_W_mm": CASE_W, "case_H_mm": CASE_H,
            "source": "чертёж оператора; RC103/scripts/build_case_stl.py:18-20",
        },
        "top_dark": measure_top_dark(),
        "top_light": measure_top_light(),
        "side": measure_side(),
    }
    if args.png:
        annotate(res)
        res["annotated_png"] = str(OUT_PNG)
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(res, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(res, ensure_ascii=False, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
