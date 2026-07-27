# -*- coding: utf-8 -*-
"""Чтение спектров RadiaCode из XML приложения.

В одном файле может лежать несколько наборов (проба и её фон). Калибровка
квадратичная: E = c0 + c1*ch + c2*ch².
"""
import xml.etree.ElementTree as ET

import numpy as np


class Spec:
    def __init__(self, name, device, live, coef, counts, weight=None, volume=None,
                 bgname=None):
        self.name = name
        self.device = device
        self.live = live            # с, «живое» время
        self.coef = coef            # (c0, c1, c2)
        self.counts = counts        # отсчёты по каналам
        self.weight = weight        # кг
        self.volume = volume        # л
        self.bgname = bgname

    @property
    def energy(self):
        ch = np.arange(len(self.counts))
        c0, c1, c2 = self.coef
        return c0 + c1 * ch + c2 * ch ** 2

    def channel_of(self, E):
        """Обратная задача к квадратичной калибровке."""
        c0, c1, c2 = self.coef
        if abs(c2) < 1e-12:
            return (E - c0) / c1
        d = c1 ** 2 + 4 * c2 * (E - c0)
        return (-c1 + np.sqrt(d)) / (2 * c2)

    def kev_per_channel(self, E):
        ch = self.channel_of(E)
        return self.coef[1] + 2 * self.coef[2] * ch

    def __repr__(self):
        return ("<%s %s: %d кан., живое %d с, %d отсчётов, %.0f..%.0f кэВ>"
                % (self.device, self.name, len(self.counts), self.live,
                   self.counts.sum(), self.energy[0], self.energy[-1]))


def _txt(node, path, default=None):
    e = node.find(path)
    return e.text if e is not None and e.text is not None else default


def read(path):
    """-> список Spec из файла."""
    root = ET.parse(path).getroot()
    out = []
    for rd in root.iter("ResultData"):
        es = rd.find("EnergySpectrum")
        if es is None:
            continue
        coef = [float(c.text) for c in es.find("EnergyCalibration/Coefficients")]
        while len(coef) < 3:
            coef.append(0.0)
        counts = np.array([float(d.text) for d in es.find("Spectrum")])
        w = _txt(rd, "SampleInfo/Weight")
        v = _txt(rd, "SampleInfo/Volume")
        out.append(Spec(
            name=_txt(rd, "SampleInfo/Name", "?"),
            device=_txt(rd, "DeviceConfigReference/Name", "?"),
            live=float(_txt(es, "MeasurementTime", "0")),
            coef=tuple(coef[:3]),
            counts=counts,
            weight=float(w) if w else None,
            volume=float(v) if v else None,
            bgname=_txt(rd, "BackgroundSpectrumFile"),
        ))
    return out


if __name__ == "__main__":
    import sys
    for p in sys.argv[1:]:
        print("===", p.split("\\")[-1])
        for s in read(p):
            print("   ", s)
            if s.weight:
                print("      масса %.3f кг, объём %.1f л, фон: %s"
                      % (s.weight, s.volume or 0, s.bgname))
