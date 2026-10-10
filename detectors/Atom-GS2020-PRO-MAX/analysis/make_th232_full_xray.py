r"""
Build the "2 %" method-2 library for GS2020 Th-232: the donor's gamma-line library (th232.yaml) unchanged,
plus K/L X-ray lines of the chain (#XR-1), reusing the existing add_xray_lines function from the sibling
module make_th232_lib05.py.

Operator 2026-09-25: "рентген K L учитываем всегда во всех методах".
The donor th232.yaml file itself is NOT modified (belongs to another module's zone).
Output: scripts/configs/th232_gs2020_full_xray.yaml.
Run with: python make_th232_full_xray.py
"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")
import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_th232_lib05 import SRC_CFG, add_xray_lines

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "configs", "th232_gs2020_full_xray.yaml")


def main():
    with open(SRC_CFG, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    lines = cfg["library"]["lines"]
    seen = {round(float(l["e_kev"]), 3) for l in lines}
    n_before = len(lines)
    n_xray = add_xray_lines(lines, seen)
    lines.sort(key=lambda x: float(x["e_kev"]))
    cfg["library"]["lines"] = lines

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg, f, allow_unicode=True, sort_keys=False)

    print(f"Донорских γ-линий: {n_before}, добавлено рентгеновских (#XR-1): {n_xray}, итог: {len(lines)}")
    print(f"Выходной файл: {OUT}")


if __name__ == "__main__":
    main()
