# -*- coding: utf-8 -*-
# Width of the Tl-208 2614.5 keV peak vs assumed centroid and window size.
# Reason: BecqMoni reports centroid 2648.86 keV, not 2614.51 - a shifted
# window inflates any half-height width estimate.
import os
import sys
_H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _H)
import fwhm_curve as fc

counts, energy, live = fc.load("field")
ecal = fc.ecal_from(energy)
print("E0      wf    FWHM      sig")
for E in (2614.51, 2630.0, 2648.86, 2665.0):
    for wf in (1.25, 1.6, 2.0):
        try:
            r = fc.measure_fwhm(counts, energy_keV=E, energy_cal=ecal,
                                window_factor=wf)
            w = fc.val(r, "fwhm_keV")
            s = fc.val(r, "significance_sigma") or fc.val(r, "significance")
            print("%7.2f %5.2f %8s %8s"
                  % (E, wf, ("%.2f" % w) if w else "-", ("%.2f" % s) if s else "-"))
        except Exception as exc:
            print("%7.2f %5.2f refused: %s" % (E, wf, exc))
