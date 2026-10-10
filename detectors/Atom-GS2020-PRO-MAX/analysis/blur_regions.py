# -*- coding: utf-8 -*-
# Замыливание областей снимка (серийные номера и пр.) перед публикацией; метаданные не переносятся (новый файл).
# python blur_regions.py <in> <out> x0,y0,x1,y1 [x0,y0,x1,y1 ...]
import sys
from PIL import Image, ImageFilter

src, dst, boxes = sys.argv[1], sys.argv[2], [tuple(int(v) for v in b.split(",")) for b in sys.argv[3:]]
im = Image.open(src).convert("RGB")
for b in boxes:
    im.paste(im.crop(b).filter(ImageFilter.GaussianBlur(6)).resize((max(1, (b[2] - b[0]) // 8), max(1, (b[3] - b[1]) // 8))).resize((b[2] - b[0], b[3] - b[1])), b)
clean = Image.new("RGB", im.size)
clean.paste(im)
clean.save(dst, quality=90)
print("замылено областей", len(boxes), "->", dst, im.size)
