#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Собирает 4 кадра rc110_vis_multiview (спереди/сбоку/разрез/три четверти)
в одну картинку-монтаж 2x2 с подписями под каждой четвертью.

Вход:  build/RadiaCode-110/rc110_view_{front,side,section,iso}.png
       (все 4 - 1600x900, пишет geometry/vis_multiview.cc, см. шапку того
       файла; порядок и подписи здесь СВЯЗАНЫ с тем, что реально пишет
       vis_multiview.cc - при переименовании выходных файлов там править и
       здесь).
Выход: build/RadiaCode-110/rc110_multiview.png

Чистый Pillow, вне сборки Geant4 (task-бриф: "можно и Python вне Geant4-
сборки, не обязательно C++") - запускать обычным python после того, как
rc110_vis_multiview.exe отработал и записал все 4 кадра.

Запуск:
    python assemble_multiview.py
"""
from __future__ import annotations

import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

# ---------------------------------------------------------------------------
# Пути. Скрипт лежит в detectors/RadiaCode-110/geant4/geometry/, кадры и
# результат - в build/RadiaCode-110/ (тот же build-каталог, что и .exe).
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[3]  # geometry -> geant4 -> RadiaCode-110 -> detectors -> repo root
BUILD_DIR = REPO_ROOT / "build" / "RadiaCode-110"

FRAMES = [
    ("rc110_view_front.png", "спереди"),
    ("rc110_view_side.png", "сбоку"),
    ("rc110_view_section.png", "разрез Y=0"),
    ("rc110_view_iso.png", "три четверти"),
]

CAPTION_H = 60          # высота полосы подписи под каждой четвертью, px
GUTTER = 8              # зазор между четвертями и по внешнему краю, px
BG = (10, 10, 10)        # тёмный фон монтажа (совпадает с фоном рендеров)
CAPTION_BG = (24, 24, 24)
CAPTION_FG = (230, 230, 230)


def _load_font(size: int) -> ImageFont.ImageFont:
    """Системный TTF если найдётся, иначе Pillow-дефолт (без Cyrillic-риска
    молчаливой замены на квадратики - стандартные Windows-шрифты Arial/
    Segoe UI кириллицу несут, проверено на этой машине)."""
    candidates = [
        r"C:\Windows\Fonts\segoeui.ttf",
        r"C:\Windows\Fonts\arial.ttf",
    ]
    for c in candidates:
        if Path(c).exists():
            try:
                return ImageFont.truetype(c, size)
            except OSError:
                continue
    return ImageFont.load_default()


def main() -> int:
    frame_paths = [(BUILD_DIR / name, caption) for name, caption in FRAMES]
    missing = [str(p) for p, _ in frame_paths if not p.exists()]
    if missing:
        print("ОШИБКА: не найдены кадры:", *missing, sep="\n  ")
        return 1

    images = [Image.open(p).convert("RGB") for p, _ in frame_paths]
    sizes = {im.size for im in images}
    if len(sizes) != 1:
        print(f"ОШИБКА: кадры разного размера {sizes} - vis_multiview.cc "
              "должен писать все 4 в одном разрешении (см. RenderView в "
              "vis_multiview.cc, /vis/open TOOLSSG_OFFSCREEN 1600x900 один "
              "раз на всю сессию)")
        return 1
    cell_w, cell_h = sizes.pop()

    font = _load_font(28)

    canvas_w = cell_w * 2 + GUTTER * 3
    canvas_h = (cell_h + CAPTION_H) * 2 + GUTTER * 3
    canvas = Image.new("RGB", (canvas_w, canvas_h), BG)
    draw = ImageDraw.Draw(canvas)

    positions = [(0, 0), (1, 0), (0, 1), (1, 1)]  # спереди/сбоку сверху, разрез/изо снизу
    for (img, (_, caption)), (col, row) in zip(zip(images, frame_paths), positions):
        x = GUTTER + col * (cell_w + GUTTER)
        y = GUTTER + row * (cell_h + CAPTION_H + GUTTER)
        canvas.paste(img, (x, y))
        cap_box = (x, y + cell_h, x + cell_w, y + cell_h + CAPTION_H)
        draw.rectangle(cap_box, fill=CAPTION_BG)
        bbox = draw.textbbox((0, 0), caption, font=font)
        tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
        tx = x + (cell_w - tw) // 2
        ty = y + cell_h + (CAPTION_H - th) // 2 - bbox[1]
        draw.text((tx, ty), caption, fill=CAPTION_FG, font=font)

    out_path = BUILD_DIR / "rc110_multiview.png"
    canvas.save(out_path)

    colors = canvas.getcolors(maxcolors=2_000_000)
    n_colors = len(colors) if colors else "many(>2e6)"
    print(f"OK: {out_path} ({canvas.size[0]}x{canvas.size[1]}, "
          f"distinct_colors={n_colors}, {out_path.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())