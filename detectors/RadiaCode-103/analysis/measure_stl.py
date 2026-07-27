# -*- coding: utf-8 -*-
"""Промер сосуда Маринелли по STL: габариты, толщины стенок, форма и глубина
колодца, объём полости в зависимости от уровня заполнения.

Метод: вертикальные лучи по сетке (x,y). Отсортированные точки пересечения с
поверхностью дают чередующиеся интервалы «материал / пустота», отсюда сразу и
толщины, и границы полостей, и объём — без предположений о форме тела.

Запуск:  python measure_stl.py <can.stl> [cap.stl]
"""
import struct
import sys

import numpy as np


def read_stl(path):
    with open(path, "rb") as f:
        n = struct.unpack("<I", f.read(84)[80:84])[0]
        buf = f.read(n * 50)
    a = np.frombuffer(buf, dtype=np.uint8).reshape(n, 50)
    return a[:, 12:48].copy().view("<f4").reshape(n, 3, 3).astype(np.float64)


def volume(tri):
    a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]
    return abs(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0) / 1000.0


class Mesh:
    def __init__(self, tri):
        self.tri = tri
        self.a, self.b, self.c = tri[:, 0], tri[:, 1], tri[:, 2]
        self.v0 = self.b[:, :2] - self.a[:, :2]
        self.v1 = self.c[:, :2] - self.a[:, :2]
        self.den = self.v0[:, 0] * self.v1[:, 1] - self.v0[:, 1] * self.v1[:, 0]
        self.ok = np.abs(self.den) > 1e-12
        self.lo = tri[:, :, :2].min(1)
        self.hi = tri[:, :, :2].max(1)

    def hits(self, x, y):
        m = self.ok & (self.lo[:, 0] <= x) & (x <= self.hi[:, 0]) \
                    & (self.lo[:, 1] <= y) & (y <= self.hi[:, 1])
        if not m.any():
            return np.zeros(0)
        a, v0, v1, den = self.a[m], self.v0[m], self.v1[m], self.den[m]
        px, py = x - a[:, 0], y - a[:, 1]
        u = (px * v1[:, 1] - py * v1[:, 0]) / den
        v = (py * v0[:, 0] - px * v0[:, 1]) / den
        ins = (u >= 0) & (v >= 0) & (u + v <= 1)
        if not ins.any():
            return np.zeros(0)
        u, v = u[ins], v[ins]
        az, bz, cz = a[ins, 2], self.b[m][ins, 2], self.c[m][ins, 2]
        z = np.sort(az + u * (bz - az) + v * (cz - az))
        return z[np.concatenate(([True], np.diff(z) > 1e-6))]


def slice_z(tri, z):
    pts = []
    for k in range(3):
        p, q = tri[:, k], tri[:, (k + 1) % 3]
        d1, d2 = p[:, 2] - z, q[:, 2] - z
        m = (d1 * d2) < 0
        if m.any():
            t = (d1[m] / (d1[m] - d2[m]))[:, None]
            pts.append(p[m, :2] + t * (q[m, :2] - p[m, :2]))
    return np.vstack(pts) if pts else np.zeros((0, 2))


def main():
    if len(sys.argv) < 2:
        raise SystemExit(
            "Укажите файл STL: python measure_stl.py <корпус.stl> "
            "[<кристалл.stl>]. Модели прибора лежат в ../drawings/.")
    can_path = sys.argv[1]
    tri = read_stl(can_path)
    lo, hi = tri.reshape(-1, 3).min(0), tri.reshape(-1, 3).max(0)
    print("=== %s" % can_path.split("\\")[-1])
    print("  треугольников %d, объём пластика %.2f см³" % (len(tri), volume(tri)))
    print("  габарит X %.2f  Y %.2f  Z %.2f мм" % tuple(hi - lo))
    print("  Z от %.2f до %.2f" % (lo[2], hi[2]))
    ztop = hi[2]

    mesh = Mesh(tri)

    # --- наружный и внутренний радиусы, форма колодца по сечениям
    print("\n  сечения: Rmax тела, ближайшая к оси поверхность, габарит колодца")
    for frac in (0.05, 0.2, 0.4, 0.6, 0.8, 0.95):
        z = lo[2] + frac * (hi[2] - lo[2])
        p = slice_z(tri, z)
        if len(p) == 0:
            continue
        r = np.hypot(p[:, 0], p[:, 1])
        near = p[r < 0.6 * r.max()]
        rn = np.hypot(near[:, 0], near[:, 1]) if len(near) else np.array([np.nan])
        bb = (near[:, 0].max() - near[:, 0].min(),
              near[:, 1].max() - near[:, 1].min()) if len(near) else (np.nan, np.nan)
        print("   z=%7.2f  Rmax=%6.2f  Rmin(бл.)=%6.2f  колодец %6.2f x %6.2f"
              % (z, r.max(), rn.min(), bb[0], bb[1]))

    # --- вертикальные разрезы
    print("\n  вертикальные разрезы (интервалы материала по z):")
    rmax = np.hypot(tri.reshape(-1, 3)[:, 0], tri.reshape(-1, 3)[:, 1]).max()
    for x, y, tag in [(0, 0, "ось колодца"), (0, 12, "стенка колодца"),
                      (0, 0.55 * rmax, "проба"), (0.8 * rmax, 0, "проба у стенки"),
                      (0.97 * rmax, 0, "стенка стакана")]:
        z = mesh.hits(x, y)
        if len(z) == 0:
            print("   (%6.1f,%6.1f) %-18s мимо тела" % (x, y, tag))
            continue
        segs = ["%.2f-%.2f" % (z[i], z[i + 1]) for i in range(0, len(z) - 1, 2)]
        print("   (%6.1f,%6.1f) %-18s %s" % (x, y, tag, "  |  ".join(segs)))

    # --- объём полости в зависимости от уровня
    rin = None
    zmid = lo[2] + 0.75 * (hi[2] - lo[2])
    p = slice_z(tri, zmid)
    if len(p):
        r = np.hypot(p[:, 0], p[:, 1])
        rin = np.sort(r)[len(r) // 20]        # внутренняя поверхность стакана
    print("\n  внутренний радиус (оценка по сечению): %.2f мм" % rin)

    step = 0.6
    g = np.arange(-rin + step / 2, rin, step)
    cols = []
    for x in g:
        for y in g:
            if x * x + y * y > rin * rin:
                continue
            z = mesh.hits(x, y)
            if len(z) == 0:
                cols.append([(lo[2], ztop)])
                continue
            zz = z[z <= ztop + 1e-9]
            if len(zz) % 2:
                zz = np.append(zz, ztop)
            gaps = []
            for i in range(1, len(zz) - 1, 2):
                gaps.append((zz[i], zz[i + 1]))
            if len(zz) and zz[-1] < ztop - 1e-6:
                gaps.append((zz[-1], ztop))
            cols.append(gaps)

    def vol_to(level):
        t = 0.0
        for gg in cols:
            for z0, z1 in gg:
                a, b = z0, min(z1, level)
                if b > a:
                    t += b - a
        return t * step * step / 1000.0

    print("\n  уровень z, мм   объём полости, см³")
    for lv in np.linspace(lo[2] + 0.4 * (ztop - lo[2]), ztop, 8):
        print("   %8.2f        %8.1f" % (lv, vol_to(lv)))
    print("   до среза       %8.1f  <-- полный объём пробы" % vol_to(ztop))

    if len(sys.argv) > 2:
        t2 = read_stl(sys.argv[2])
        l2, h2 = t2.reshape(-1, 3).min(0), t2.reshape(-1, 3).max(0)
        print("\n=== %s" % sys.argv[2].split("\\")[-1])
        print("  объём пластика %.2f см³, габарит %.2f x %.2f x %.2f"
              % (volume(t2), *(h2 - l2)))


if __name__ == "__main__":
    main()
