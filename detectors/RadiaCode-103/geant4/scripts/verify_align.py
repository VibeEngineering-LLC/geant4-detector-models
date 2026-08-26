#!/usr/bin/env python3
"""Проверка совмещения STL-корпуса и внутренностей. См. ../docs/GEANT4-MODEL.md."""
from __future__ import annotations

import struct
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

GEANT4_ROOT = Path(__file__).resolve().parents[1]
RC103_ROOT = GEANT4_ROOT.parent
sys.path.insert(0, str(RC103_ROOT / "scripts"))
from rc103_coords import COORDS  # noqa: E402

from stl_to_gdml import resolve_stl  # noqa: E402

STL = resolve_stl()
OUT = GEANT4_ROOT / "verify" / "RC103_align_check.png"

CASE_SIZE = COORDS["Case_STL"]["size"]
CASE_X = (-CASE_SIZE[0] / 2, CASE_SIZE[0] / 2)
CASE_Z = (-CASE_SIZE[2] / 2, CASE_SIZE[2] / 2)
CAVITY = (120.0, 31.0, 14.5)

CRYSTAL_C = COORDS["Crystal_CsI"]["center"]
CRYSTAL_SIZE = COORDS["Crystal_CsI"]["size"]
SIPM_C = COORDS["SiPM"]["center"]
SIPM_SIZE = COORDS["SiPM"]["size"]
CAPSULE_C = COORDS["Capsule_body"]["center"]
CAPSULE_SIZE = COORDS["Capsule_body"]["size"]
PCB_C = COORDS["PCB_FR4"]["center"]
PCB_SIZE = COORDS["PCB_FR4"]["size"]
BATT_C = COORDS["Battery_LiPo"]["center"]
BATT_SIZE = COORDS["Battery_LiPo"]["size"]
DISP_LCD_C = COORDS["Display_LCD"]["center"]
DISP_LCD_SIZE = COORDS["Display_LCD"]["size"]
DISP_WIN_C = COORDS["Display_window"]["center"]
DISP_WIN_SIZE = COORDS["Display_window"]["size"]


def remap(x: float, y: float, z: float) -> tuple[float, float, float]:
    return (x, y, z)


def load_stl_sample(path: Path, step: int = 40) -> list:
    data = path.read_bytes()
    n = struct.unpack_from("<I", data, 80)[0]
    tris = []
    off = 84
    for i in range(n):
        v = struct.unpack_from("<12fH", data, off)
        if i % step == 0:
            tri = [remap(v[3 + j * 3], v[4 + j * 3], v[5 + j * 3]) for j in range(3)]
            tris.append(tri)
        off += 50
    return tris


def box_faces(cx, cy, cz, hx, hy, hz):
    x0, x1 = cx - hx, cx + hx
    y0, y1 = cy - hy, cy + hy
    z0, z1 = cz - hz, cz + hz
    c = [
        (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
        (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
    ]
    faces = [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4], [2, 3, 7, 6], [0, 3, 7, 4], [1, 2, 6, 5]]
    return [[c[i] for i in f] for f in faces]


def box_edges(cx, cy, cz, hx, hy, hz):
    x0, x1 = cx - hx, cx + hx
    y0, y1 = cy - hy, cy + hy
    z0, z1 = cz - hz, cz + hz
    corners = np.array([
        [x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
        [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1],
    ])
    edge_idx = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]
    return [corners[[a, b]] for a, b in edge_idx]


def add_box_wire(ax, cx, cy, cz, hx, hy, hz, *, color, lw=1.5):
    ax.add_collection3d(Line3DCollection(box_edges(cx, cy, cz, hx, hy, hz), colors=color, linewidths=lw))


def add_solid_box(ax, cx, cy, cz, hx, hy, hz, *, fc, ec, alpha=0.7):
    ax.add_collection3d(
        Poly3DCollection(box_faces(cx, cy, cz, hx, hy, hz), alpha=alpha, facecolor=fc, edgecolor=ec, linewidths=0.8)
    )


def draw_scene(ax, tris) -> None:
    ax.add_collection3d(
        Poly3DCollection(tris, alpha=0.45, facecolor="#6fa8dc", edgecolor="#2a6496", linewidths=0.08)
    )
    hx, hy, hz = (s / 2 for s in CASE_SIZE)
    add_box_wire(ax, 0, 0, 0, hx, hy, hz, color="#c0392b", lw=2.0)
    px, py, pz = (s / 2 for s in PCB_SIZE)
    add_solid_box(ax, *PCB_C, px, py, pz, fc="#1b7a3d", ec="#0d4d24", alpha=0.8)
    bx, by, bz = (s / 2 for s in BATT_SIZE)
    add_solid_box(ax, *BATT_C, bx, by, bz, fc="#c0c4c8", ec="#7a8088", alpha=0.9)
    cx, cy, cz = (s / 2 for s in CAPSULE_SIZE)
    add_solid_box(ax, *CAPSULE_C, cx, cy, cz, fc="#888888", ec="#555555", alpha=0.35)
    kx, ky, kz = (s / 2 for s in CRYSTAL_SIZE)
    add_solid_box(ax, *CRYSTAL_C, kx, ky, kz, fc="gold", ec="darkorange", alpha=0.95)
    sx, sy, sz = (s / 2 for s in SIPM_SIZE)
    add_solid_box(ax, *SIPM_C, sx, sy, sz, fc="#00bcd4", ec="#00838f", alpha=0.95)
    lx, ly, lz = (s / 2 for s in DISP_LCD_SIZE)
    add_solid_box(ax, *DISP_LCD_C, lx, ly, lz, fc="#880e4f", ec="#4a148c", alpha=0.9)
    wx, wy, wz = (s / 2 for s in DISP_WIN_SIZE)
    add_solid_box(ax, *DISP_WIN_C, wx, wy, wz, fc="#1565c0", ec="#0d47a1", alpha=0.85)
    ax.set_xlabel("X, USB → +X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z, −Z = лицо")
    ax.set_xlim(-72, 72)
    ax.set_ylim(-20, 20)
    ax.set_zlim(-14, 14)
    ax.view_init(elev=20, azim=-65)
    ax.set_box_aspect((2.2, 0.55, 0.35))


def draw_xz_section(ax) -> None:
    ax.add_patch(
        plt.Rectangle(
            (CASE_X[0], CASE_Z[0]),
            CASE_X[1] - CASE_X[0],
            CASE_Z[1] - CASE_Z[0],
            fill=False,
            ec="#c0392b",
            lw=2,
            label="корпус STL bbox",
        )
    )
    cvx, _, cvz = CAVITY
    ax.add_patch(
        plt.Rectangle(
            (-cvx / 2, -cvz / 2),
            cvx,
            cvz,
            fill=False,
            ec="#27ae60",
            lw=1.2,
            ls="--",
            label="полость",
        )
    )
    ax.add_patch(
        plt.Rectangle(
            (PCB_C[0] - PCB_SIZE[0] / 2, PCB_C[2] - PCB_SIZE[2] / 2),
            PCB_SIZE[0],
            PCB_SIZE[2],
            fc="#1b7a3d",
            ec="#0d4d24",
            lw=1,
            alpha=0.8,
            label="плата",
        )
    )
    ax.add_patch(
        plt.Rectangle(
            (BATT_C[0] - BATT_SIZE[0] / 2, BATT_C[2] - BATT_SIZE[2] / 2),
            BATT_SIZE[0],
            BATT_SIZE[2],
            fc="#c0c4c8",
            ec="#7a8088",
            lw=1,
            alpha=0.9,
            label="LiPo 60 mm",
        )
    )
    cap_x, _, cap_z = CAPSULE_C
    ax.add_patch(
        plt.Rectangle(
            (cap_x - CAPSULE_SIZE[0] / 2, cap_z - CAPSULE_SIZE[2] / 2),
            CAPSULE_SIZE[0],
            CAPSULE_SIZE[2],
            fill=False,
            ec="#7f8c8d",
            lw=1,
        )
    )
    cr_x, _, cr_z = CRYSTAL_C
    ax.add_patch(
        plt.Rectangle(
            (cr_x - CRYSTAL_SIZE[0] / 2, cr_z - CRYSTAL_SIZE[2] / 2),
            CRYSTAL_SIZE[0],
            CRYSTAL_SIZE[2],
            fc="gold",
            ec="darkorange",
            lw=1.2,
            alpha=0.85,
            label="CsI 10 mm",
        )
    )
    sp_x, _, sp_z = SIPM_C
    ax.add_patch(
        plt.Rectangle(
            (sp_x - SIPM_SIZE[0] / 2, sp_z - SIPM_SIZE[2] / 2),
            SIPM_SIZE[0],
            SIPM_SIZE[2],
            fc="#00bcd4",
            ec="#00838f",
            lw=1,
            label="SiPM 6×6",
        )
    )
    dl_x, _, dl_z = DISP_LCD_C
    ax.add_patch(
        plt.Rectangle(
            (dl_x - DISP_LCD_SIZE[0] / 2, dl_z - DISP_LCD_SIZE[2] / 2),
            DISP_LCD_SIZE[0],
            DISP_LCD_SIZE[2],
            fc="#880e4f",
            ec="#4a148c",
            lw=1,
            alpha=0.9,
            label="экран LCD",
        )
    )
    dw_x, _, dw_z = DISP_WIN_C
    ax.add_patch(
        plt.Rectangle(
            (dw_x - DISP_WIN_SIZE[0] / 2, dw_z - DISP_WIN_SIZE[2] / 2),
            DISP_WIN_SIZE[0],
            DISP_WIN_SIZE[2],
            fc="#1565c0",
            ec="#0d47a1",
            lw=1,
            alpha=0.85,
            label="дисплей",
        )
    )
    ax.set_xlim(-72, 72)
    ax.set_ylim(-16, 16)
    ax.set_aspect("equal")
    ax.set_xlabel("X (USB → +X)")
    ax.set_ylabel("Z (−Z = лицо)")
    ax.set_title("разрез Y = 0")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, alpha=0.3)


def main() -> None:
    tris = load_stl_sample(STL, step=35)
    fig = plt.figure(figsize=(13, 5))
    fig.suptitle("RC-103: корпус STL + плата + LiPo + детектор", fontsize=12)
    ax1 = fig.add_subplot(121, projection="3d")
    ax1.set_title("3D (общий вид)")
    draw_scene(ax1, tris)
    ax2 = fig.add_subplot(122)
    draw_xz_section(ax2)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=160, bbox_inches="tight")
    print(f"saved: {OUT}")


if __name__ == "__main__":
    main()
