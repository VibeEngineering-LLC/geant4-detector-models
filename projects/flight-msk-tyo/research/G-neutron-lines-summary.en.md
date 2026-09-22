# G. Neutron-induced lines (n,g) and (n,n') — summary (22.09.2026, rev. 4)

Reaction source: `research/B-gamma-lines-airframe.en.md`, section 9. CANDIDATE significance is computed strictly at the tabulated level energy (`scripts/lines_report.py`, 3-channel window, 4σ threshold) — not the same as "is there a real peak nearby": the tabulated energy and the centre of the observed peak may not coincide by a few keV. TWO significances are given: in the SUMMED raw spectrum (what is actually plotted on the page) and separately in the component produced by field neutrons (`prod3_total_neutron.csv`), because the gamma continuum of other components can mask the line in the sum even when the reaction itself is statistically robust.
"Found" — the candidate passed the significance threshold in the sum (that is what is shown on the page). "Seen in the neutron component, not resolved in the sum" — the threshold is not passed in the sum but is passed with margin in the neutron component (typical case — a narrow line on a high continuum of other particles). "Significant peak nearby, unidentified" — not significant at the tabulated point in either the sum or the component, but an independent peak search found a significant maximum in the sum within tolerance; the program did not link them automatically. "Candidate" — nothing significant at the point or nearby, in either the sum or the component. "Outside search window" — the adjacent energy is occupied by another line.
The "Yield per 100 captures" column is the probability that a neutron captured by THIS SPECIFIC nuclide emits a gamma ray of exactly this energy (not the probability of capture itself); source — the IAEA PGAA adopted database, `research/data/nglist_a.dat`/`promptgammas.xls` (see `research/B-gamma-lines-airframe.md` line 53). For (n,n′) — a dash: this is inelastic scattering, no capture occurs.

| E, keV | Element | Reaction | Material | Status | Significance in the sum, σ | Significance in neutron component, σ | Nearest peak | Yield per 100 captures |
|---|---|---|---|---|---|---|---|---|
| 27.36 | I | ¹²⁷I(n,γ) | the crystal itself (NaI/CsI) | outside search window | — | — | — | 6.94/100 capt. |
| 29.83 | K | ³⁹K(n,γ) | people, luggage | outside search window | — | — | — | 65.71 |
| 30.64 | Al | ²⁷Al(n,γ) | skin, structural framing | seen in the neutron component, not resolved in the sum | 2.7 | 17.9 | — | 34.55 |
| 57.61 | I | ¹²⁷I(n,n′) (1st level 57.608 keV, 7/2⁺, ENSDF) and ¹²⁷I(n,γ) 58.11 | NaI/CsI crystal — the classic "iodine escape bump" | seen in the neutron component, not resolved in the sum | -2.2 | 18.0 | — | 4.52 (capture) |
| 58.11 | I | ¹²⁷I(n,n′) (1st level 57.608 keV, 7/2⁺, ENSDF) and ¹²⁷I(n,γ) 58.11 | NaI/CsI crystal — the classic "iodine escape bump" | seen in the neutron component, not resolved in the sum | 1.1 | 18.2 | — | 4.52 (capture) |
| 78.08 | P | ³¹P(n,γ) | passengers' bones | seen in the neutron component, not resolved in the sum | -0.8 | 5.1 | — | 34.30 |
| 78.91 | Ag | ¹⁰⁷Ag(n,γ) | contacts, solder | seen in the neutron component, not resolved in the sum | 3.4 | 7.0 | — | 10.37 |
| 90.99 | Na | ²³Na(n,γ) | people, NaI | candidate | -2.2 | 1.0 | — | 44.34 |
| 115.22 | Zn | ⁶⁴Zn(n,γ) | 7075 alloy | seen in the neutron component, not resolved in the sum | 2.1 | 7.3 | — | 15.18 |
| 116.61 | Cs | ¹³³Cs(n,γ) | CsI crystal | seen in the neutron component, not resolved in the sum | 1.4 | 6.7 | — | 4.75 |
| 133.61 | I | ¹²⁷I(n,γ) | NaI/CsI crystal | seen in the neutron component, not resolved in the sum | -3.1 | 16.2 | — | 22.90 |
| 137.85 | I | ¹²⁷I(n,γ)¹²⁸I, isomeric level 4⁻, T½=845 ns (populated via a cascade through the 167.4 keV level, T½=175 ns; nearly all cascade energy is absorbed in the crystal — Yoon et al., NIM A 848 (2017) 162, arXiv:1608.02409) | NaI/CsI crystal | found | 5.0 | 36.5 | — | ≥0.07 (direct population per our emitter database; higher when the 167.4 keV cascade is included) |
| 159.28 | Cu | ⁶³Cu(n,γ) | wiring, 2024 alloy | seen in the neutron component, not resolved in the sum | -0.9 | 11.6 | — | 14.34 |
| 176.40 | Cs | ¹³³Cs(n,γ) | CsI crystal | seen in the neutron component, not resolved in the sum | 0.0 | 9.8 | — | 8.15 |
| 185.96 | Cu | ⁶⁵Cu(n,γ) | wiring | candidate | 2.0 | 0.4 | — | 11.24 |
| 205.62 | Cs | ¹³³Cs(n,γ) | CsI | seen in the neutron component, not resolved in the sum | 0.4 | 13.9 | — | 5.15 |
| 278.25 | Cu | ⁶³Cu(n,γ) | wiring | candidate | 0.1 | -15.8 | — | 19.76 |
| 339.42 | Ni | ⁵⁸Ni(n,γ) | steel, heat-resistant alloys | candidate | -2.4 | -2.9 | — | 3.71 |
| 341.71 | Ti | ⁴⁸Ti(n,γ) | titanium parts | seen in the neutron component, not resolved in the sum | 0.4 | 4.0 | — | 23.35 |
| 352.35 | Fe | ⁵⁶Fe(n,γ) | steel, fasteners | candidate | -4.2 | -0.4 | — | 10.54 |
| 464.98 | Ni | ⁵⁸Ni(n,γ) | steel | seen in the neutron component, not resolved in the sum | 3.6 | 4.3 | — | 18.73 |
| 465.14 | Cu | ⁶⁵Cu(n,γ) | wiring | found | 4.3 | 5.4 | — | 6.22 |
| 512.65 | P | ³¹P(n,γ) | bones | found | 33.3 | 9.9 | — | 45.93 |
| 517.07 | Cl | ³⁵Cl(n,γ) | PVC wiring, body salt | candidate | 1.6 | 2.4 | — | 17.41 |
| 519.66 | Ca | ⁴⁰Ca(n,γ) | passengers' skeleton | candidate | -0.5 | 2.7 | — | 12.27 |
| 558.32 | Cd | ¹¹³Cd(n,γ) | cadmium-plated fasteners | candidate | 1.7 | 0.7 | — | 9.03 (σ_elem 227 b!) |
| 564.05 | Cr | ⁵²Cr(n,γ) | stainless steel | seen in the neutron component, not resolved in the sum | 0.7 | 4.5 | — | 14.87 |
| 585.00 | Mg | ²⁴Mg(n,γ) | Mg in 2024/7075 alloys | candidate | 0.4 | -3.1 | — | 58.58 |
| 608.77 | Cu | ⁶³Cu(n,γ) | wiring | candidate | -1.2 | 1.5 | — | 5.97 |
| 636.66 | P | ³¹P(n,γ) | bones | candidate | 0.6 | -0.5 | — | 18.08 |
| 669.72 | Cu | ⁶³Cu(n,n′) | wiring, busbars | candidate | 1.1 | 0.3 | — | — |
| 770.30 | K | ³⁹K(n,γ) | people, luggage | candidate | 0.6 | -1.1 | — | 43.00 |
| 786.30 | Cl | ³⁵Cl(n,γ) | PVC, body | seen in the neutron component, not resolved in the sum | 1.2 | 6.0 | — | 7.85 / 12.45 |
| 788.43 | Cl | ³⁵Cl(n,γ) | PVC, body | candidate | 1.0 | -0.0 | — | 7.85 / 12.45 |
| 840.99 | S | ³²S(n,γ) | body, rubber | seen in the neutron component, not resolved in the sum | 0.2 | 4.3 | — | 63.32 |
| 843.76 | Al | ²⁷Al(n,n′) σ 100 mb @2 MeV and ²⁷Mg decay (71.8 %) and μ⁻ capture on Al | skin | significant peak nearby, unidentified | -0.9 | -7.4 | 837.5 keV, 4.6σ (Δ=-6.3) | — |
| 846.78 | Fe | ⁵⁶Fe(n,n′) σ 528 mb @2 MeV, peak 1999 mb @1.38 MeV; and ⁵⁶Mn decay (98.85 %) | steel, fasteners, engines | candidate | 1.5 | -1.8 | — | — |
| 855.69 | Zn | ⁶⁴Zn(n,γ) | 7075 alloy | candidate | 0.8 | -0.7 | — | 6.00 |
| 869.21 | Na | ²³Na(n,γ) | body, NaI | candidate | 3.9 | 2.7 | — | 20.38 |
| 870.68 | O | ¹⁶O(n,γ) | air, water, people | candidate | 1.0 | 2.5 | — | 93.16, but σ₀ = 0.19 mb |
| 874.39 | Na | ²³Na(n,γ) | body, NaI | candidate | -0.0 | 3.1 | — | 14.34 |
| 877.98 | Ni | ⁵⁸Ni(n,γ) | steel | candidate | 2.1 | -0.9 | — | 5.24 |
| 983.53 | Ti | ⁴⁸Ti(n,n′) | titanium parts | candidate | 1.8 | -5.9 | — | — |
| 1014.52 | Al | ²⁷Al(n,n′) σ 234 mb @2 MeV and ²⁷Mg decay (28.2 %) | skin | candidate | 1.0 | 1.3 | — | — |
| 1071.22 | P | ³¹P(n,γ) | bones | candidate | -7.7 | -11.5 | — | 14.48 |
| 1087.75 | O | ¹⁶O(n,γ) | water, people | candidate | -0.9 | -2.3 | — | 83.16 (σ₀ tiny) |
| 1158.89 | K | ³⁹K(n,γ) | people | candidate | 2.8 | 2.1 | — | 7.62 |
| 1164.87 | Cl | ³⁵Cl(n,γ) | PVC insulation, body salt | candidate | -0.5 | -1.3 | — | 20.46 |
| 1249.62 | Sn | ¹¹⁸Sn(n,γ) | solders | candidate | 0.5 | 2.0 | — | 2.26 |
| 1261.77 | C | ¹²C(n,γ) | fuel, composite, people | candidate | 2.8 | 1.4 | — | 35.13 |
| 1273.35 | Si | ²⁸Si(n,γ) | electronics, glass | found | 4.2 | -0.9 | — | 16.33 |
| 1326.01 | Cu | ⁶³Cu(n,n′) | wiring | candidate | -5.1 | -5.3 | — | — |
| 1368.63 | Na | ²⁴Na decay (99.99 %) and ²⁴Mg(n,n′) 1368.67 | Al alloys, activation | candidate | -9.1 | -1.3 | — | — |
| 1381.74 | Ti | ⁴⁸Ti(n,γ) | titanium parts, engines | candidate | -5.8 | 1.0 | — | 65.74, σ_elem 3.8 b |
| 1434.06 | V | ⁵²V decay (100 %) and ⁵²Cr(n,n′) σ 1003 mb @2 MeV | stainless steel | candidate | 0.5 | 2.5 | — | — |
| 1434.09 | V | ⁵²V decay (100 %) and ⁵²Cr(n,n′) σ 1003 mb @2 MeV | stainless steel | candidate | 0.5 | 2.5 | — | — |
| 1454.21 | Ni | ⁵⁸Ni(n,n′) σ 578 mb @2 MeV | heat-resistant alloys | significant peak nearby, unidentified | 0.8 | 0.8 | 1457.5 keV, 5.0σ (Δ=3.3) | — |
| 1460.85 | Ar | ⁴⁰Ar(n,n′) — overlaps with ⁴⁰K! | cabin air | found | 8.2 | 0.1 | — | — |
| 1585.94 | Ti | ⁴⁸Ti(n,γ) | titanium | candidate | 0.3 | -0.5 | — | 7.92 |
| 1612.79 | Fe | ⁵⁶Fe(n,γ) | steel | candidate | -2.2 | 2.1 | — | 5.91 |
| 1618.97 | K | ³⁹K(n,γ) | people | candidate | 1.5 | 0.1 | — | 6.19 |
| 1622.88 | Al | ²⁷Al(n,γ) | skin | candidate | 0.6 | 0.4 | — | 4.28 |
| 1635.20 | N | ¹⁴N(n,n′) (3948.1→2312.8), σ(MT52) 100 mb @6.75 MeV | air, people | candidate | -0.1 | -4.6 | — | — |
| 1725.29 | Fe | ⁵⁶Fe(n,γ) | steel | candidate | -0.6 | -8.3 | — | 6.99 |
| 1779.03 | Si | ²⁸Si(n,n′), σ 579 mb @5 MeV, peak 1032 mb @2.86 MeV | glass, electronics, fillers | candidate | 0.3 | 3.9 | — | — |
| 1884.82 | N | ¹⁴N(n,γ) | air, people | candidate | 1.6 | 2.5 | — | 18.42 |
| 1942.67 | Ca | ⁴⁰Ca(n,γ) | passengers' skeleton | candidate | 0.2 | -0.3 | — | 85.85 |
| 1951.14 | Cl | ³⁵Cl(n,γ) | PVC, body | candidate | 2.8 | 1.2 | — | 14.54 |
| 2001.31 | Ca | ⁴⁰Ca(n,γ) | skeleton | candidate | -2.9 | -0.7 | — | 16.07 / 9.98 |
| 2009.84 | Ca | ⁴⁰Ca(n,γ) | skeleton | candidate | -1.5 | -0.5 | — | 16.07 / 9.98 |
| 2073.79 | K | ³⁹K(n,γ) | people | candidate | 0.2 | -1.6 | — | 6.52 |
| 2085.10 | Fe | ⁵⁶Fe(n,n′) 2nd level | steel | candidate | -2.8 | -4.5 | — | — |
| 2092.90 | Si | ²⁸Si(n,γ) | glass, electronics | candidate | -2.8 | -2.5 | — | 18.70 |
| 2184.42 | O | ¹⁶O(n,γ) | water | candidate | 0.9 | 1.6 | — | 86.32 (σ₀ tiny) |
| 2223.25 | H | ¹H(n,γ)²H | people, fuel, luggage, plastic, water — the main capture source | found | 7.6 | 18.6 | — | 100 / 100 captures, σ = 0.3326 b |
| 2312.59 | N | ¹⁴N(n,n′), σ(MT51) 19.7 mb @4.9 MeV | air, people | candidate | -3.9 | 2.9 | — | — |
| 2320.80 | Cr | ⁵²Cr(n,γ) | stainless steel | candidate | 0.3 | 0.8 | — | 17.89 |
| 2379.66 | S | ³²S(n,γ) | body, rubber | candidate | 0.5 | 2.2 | — | 37.96 |
| 2517.81 | Na | ²³Na(n,γ) | body, NaI | candidate | 0.2 | -1.7 | — | 13.19 |
| 2614.53 | Pb | ²⁰⁷Pb(n,γ) (0.06) and ²⁰⁸Tl 2614.51 (99.75 %) | instrument lead shield / Th alloys | candidate | -0.4 | -0.2 | — | — |
| 2741.50 | O | ¹⁶O(n,n′) 8871.9→6129.89 and ¹⁶N decay (0.82 %) | air, water, people | candidate | 1.6 | 1.0 | — | — |
| 2752.27 | Na | ²³Na(n,γ) | body, NaI | candidate | 1.9 | 1.8 | — | 12.34 |
| 2828.17 | Mg | ²⁴Mg(n,γ) | Mg in alloys | candidate | 0.9 | 1.6 | — | 44.78 |
| 2930.67 | S | ³²S(n,γ) | body | candidate | -0.9 | -4.3 | — | 15.18 |
| 3033.90 | Al | ²⁷Al(n,γ) | skin | candidate | -2.5 | -3.0 | — | 7.75 |
| 3220.59 | S | ³²S(n,γ) | body | candidate | -0.2 | 0.2 | — | 21.35 |
| 3272.02 | O | ¹⁶O(n,γ) | water | candidate | 1.6 | 2.1 | — | 18.58 |
| 3465.06 | Al | ²⁷Al(n,γ) | skin | candidate | -2.4 | -2.8 | — | 6.32 |
| 3538.97 | Si | ²⁸Si(n,γ) | glass, electronics | candidate | -0.7 | -1.1 | — | 67.23 |
| 3587.46 | Na | ²³Na(n,γ) | body, NaI | candidate | 0.2 | -3.6 | — | 11.25 |
| 3683.92 | C | ¹²C(n,γ) 34.56 and ¹³C* from ¹⁶O(n,α), σ(MT107) 506 mb @4.6 MeV | fuel, composite, people / oxygen | candidate | -3.7 | -7.0 | — | — |
| 3899.89 | P | ³¹P(n,γ) | bones | candidate | -1.2 | -4.6 | — | 17.09 |
| 3916.84 | Mg | ²⁴Mg(n,γ) | Mg in alloys | candidate | -6.0 | 1.6 | — | 59.70 |
| 3981.45 | Na | ²³Na(n,γ) | body, NaI | candidate | 1.3 | -0.3 | — | 12.77 |
| 4133.41 | Al | ²⁷Al(n,γ) | skin | candidate | 1.6 | -2.0 | — | 6.45 |
| 4259.53 | Al | ²⁷Al(n,γ) | skin | candidate | 1.4 | 1.1 | — | 6.62 |
| 4438.90 | C | ¹²C(n,n′), σ 373 mb @8 MeV, peak 460 mb @8.14 MeV; Doppler-broadened to ~180 keV | fuel, composite, plastic, people | candidate | -1.2 | 1.0 | — | — |
| 4444.03 | B | ¹⁰B(n,γ) | borosilicate glass | candidate | -0.5 | -2.1 | — | 7.80 |
| 4508.73 | N | ¹⁴N(n,γ) | air, people | candidate | 1.5 | 0.0 | — | 16.54 |
| 4690.68 | Al | ²⁷Al(n,γ) | skin | candidate | -1.6 | -5.5 | — | 4.72 / 5.45 |
| 4733.84 | Al | ²⁷Al(n,γ) | skin | candidate | 2.1 | 0.3 | — | 4.72 / 5.45 |
| 4869.61 | S | ³²S(n,γ) | body | candidate | 2.0 | 0.0 | — | 11.86 |
| 4933.89 | Si | ²⁸Si(n,γ) | glass, electronics | candidate | 0.8 | 0.0 | — | 63.28 |
| 4945.30 | C | ¹²C(n,γ) = Sₙ(¹³C) | fuel, composite, people | candidate | -3.6 | 0.0 | — | 73.94 (but σ₀ = 3.53 mb) |
| 5104.89 | N | ¹⁴N(n,n′) 5105.89→g.s. | air, people | candidate | 1.7 | 0.0 | — | — |
| 5269.16 | N | ¹⁴N(n,γ) | air, people | candidate | 1.1 | -3.5 | — | 29.57 / 21.05 |
| 5297.82 | N | ¹⁴N(n,γ) | air, people | candidate | -0.5 | 0.0 | — | 29.57 / 21.05 |
| 5380.02 | K | ³⁹K(n,γ) | people | candidate | -0.0 | 0.0 | — | 6.95 |
| 5420.57 | S | ³²S(n,γ) | body, rubber | candidate | 0.1 | 0.0 | — | 56.20 |
| 5533.40 | N | ¹⁴N(n,γ) | air, people | candidate | 1.6 | 0.8 | — | 19.42 / 10.53 |
| 5562.06 | N | ¹⁴N(n,γ) | air, people | candidate | 0.7 | 0.0 | — | 19.42 / 10.53 |
| 5617.90 | Cr | ⁵²Cr(n,γ) | stainless steel | candidate | 0.1 | 0.0 | — | 17.37 |
| 5920.45 | Fe | ⁵⁶Fe(n,γ) | steel | candidate | -1.4 | 0.0 | — | 8.69 |
| 6018.53 | Fe | ⁵⁶Fe(n,γ) | steel | candidate | 1.1 | 0.0 | — | 8.76 |
| 6110.84 | Cl | ³⁵Cl(n,γ) | PVC, body | candidate | 1.5 | -0.5 | — | 15.13 |
| 6128.63 | O | ¹⁶O(n,n′), σ(MT52) 149 mb @8 MeV, peak 362 mb @8.38 MeV; and ¹⁶N decay (67 %) | air, water, people, fuel | candidate | -5.7 | 0.0 | — | — |
| 6322.43 | N | ¹⁴N(n,γ) | air | candidate | 0.9 | 0.0 | — | 18.17 |
| 6379.80 | Si | ²⁸Si(n,γ) | glass | candidate | -4.6 | 0.0 | — | 11.69 |
| 6395.48 | Na | ²³Na(n,γ) | body, NaI | candidate | 2.2 | 0.0 | — | 18.87 |
| 6418.43 | Ti | ⁴⁸Ti(n,γ) | titanium | candidate | -2.0 | 0.0 | — | 24.87 |
| 6419.59 | Ca | ⁴⁰Ca(n,γ) | skeleton — overlap with Ti! | candidate | 0.6 | 0.0 | — | 42.93 |
| 6555.91 | Ti | ⁴⁸Ti(n,γ) | titanium | candidate | -0.7 | 0.0 | — | 4.24 |
| 6600.63 | Cu | ⁶⁵Cu(n,γ) | wiring | candidate | 1.3 | 0.0 | — | 3.92 / 3.73 |
| 6619.61 | Cl | ³⁵Cl(n,γ) | PVC | candidate | -0.7 | 0.0 | — | 5.81 |
| 6680.00 | Cu | ⁶⁵Cu(n,γ) | wiring | candidate | 2.2 | 0.0 | — | 3.92 / 3.73 |
| 6760.08 | Ti | ⁴⁸Ti(n,γ) | titanium | candidate | -7.8 | 0.0 | — | 37.69 |
| 6785.50 | P | ³¹P(n,γ) | bones | candidate | 1.6 | 0.0 | — | 15.52 |
| 6915.50 | O | ¹⁶O(n,n′) 6917.1→g.s., σ(MT53) 59 mb @8 MeV | air, water | candidate | 0.8 | 0.0 | — | — |
| 7069.20 | Zn | ⁶⁴Zn(n,γ) | 7075 alloy | candidate | 0.2 | 0.0 | — | 1.85 / 1.80 |
| 7111.95 | Zn | ⁶⁴Zn(n,γ) | 7075 alloy | candidate | -2.9 | 0.0 | — | 1.85 / 1.80 |
| 7115.15 | O | ¹⁶O(n,n′) 7116.85→g.s.; and ¹⁶N decay (4.9 %) | air, water | candidate | 2.0 | 0.0 | — | — |
| 7199.20 | Si | ²⁸Si(n,γ) | glass | candidate | 0.2 | 0.0 | — | 7.06 |
| 7306.93 | Cu | ⁶³Cu(n,γ) | wiring | candidate | -0.0 | 0.0 | — | 7.10 |
| 7367.78 | Pb | ²⁰⁷Pb(n,γ) = Sₙ | instrument lead shield | candidate | 1.6 | 0.0 | — | 21.92 |
| 7374.49 | Cr | ⁵²Cr(n,γ) | stainless steel | candidate | -1.4 | 0.0 | — | 10.53 |
| 7413.97 | Cl | ³⁵Cl(n,γ) | PVC, body | candidate | 1.3 | 0.0 | — | 7.55 / 6.11 |
| 7631.14 | Fe | ⁵⁶Fe(n,γ), doublet | steel, fasteners, engines | candidate | -4.7 | 0.0 | — | 25.21 + 21.20 |
| 7637.40 | Cu | ⁶³Cu(n,γ) | wiring | candidate | 1.9 | 0.0 | — | 11.95 |
| 7724.03 | Al | ²⁷Al(n,γ) = Sₙ(²⁸Al) | skin | candidate | 1.5 | 0.0 | — | 21.34 |
| 7768.92 | K | ³⁹K(n,γ) | people | candidate | 1.3 | 0.0 | — | 5.57 |
| 7790.33 | Cl | ³⁵Cl(n,γ) | PVC, body | candidate | 1.2 | 0.0 | — | 7.55 / 6.11 |
| 7819.52 | Ni | ⁶⁰Ni(n,γ) | steel | candidate | 1.0 | 0.0 | — | 11.59 |
| 7863.55 | Zn | ⁶⁴Zn(n,γ) | 7075 alloy | candidate | 1.7 | 0.0 | — | 12.82 |
| 7915.62 | Cu | ⁶³Cu(n,γ) = Sₙ | wiring, busbars | candidate | 1.5 | 0.0 | — | 19.23 |
| 7938.46 | Cr | ⁵²Cr(n,γ) = Sₙ | stainless steel | candidate | 0.2 | 0.0 | — | 55.79 |
| 8120.57 | Ni | ⁵⁸Ni(n,γ) | steel | candidate | 2.1 | 0.0 | — | 2.96 |
| 8533.51 | Ni | ⁵⁸Ni(n,γ) | steel | candidate | 0.8 | 0.0 | — | 16.02 |
| 8869.30 | O | ¹⁶O(n,n′) 8871.9→g.s. (rel. 9.3) | air, water | candidate | -1.5 | 0.0 | — | — |
| 8998.41 | Ni | ⁵⁸Ni(n,γ) = Sₙ(⁵⁹Ni) | heat-resistant alloys, engines | candidate | 1.3 | 0.0 | — | 33.11 |
| 10829.12 | N | ¹⁴N(n,γ) = Sₙ(¹⁵N) | cabin and outside air | outside search window | — | — | — | 14.16 |
