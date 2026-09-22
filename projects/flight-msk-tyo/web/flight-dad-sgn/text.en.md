@@HEAD@@

## Input data: the radiation field at cruise altitude

The cosmic-ray field at 9800 m outside the aircraft is taken from the PARMA/EXPACS 4.10 model (Sato 2015, 2016) for the flight conditions: solar activity from the Oulu neutron monitor on 20.09.2026 (W = 76), geomagnetic cutoff rigidity 17.28 GV (midpoint of the Da Nang–Ho Chi Minh City route; one of the highest cutoff rigidities on Earth — near the equator, Earth's magnetic field deflects most primary cosmic particles, and the flux at the same altitude is noticeably lower than at mid and high latitudes; by how much was not calculated separately for this flight), atmospheric depth above the aircraft 278.7 g/cm². The aircraft itself is not included in the PARMA model (a "clear-sky" field) — its effect is accounted for separately in Geant4.

The table gives the flux density (summed over all directions) and the energy distribution within each component.

@@FIELD@@

What this means for the spectrum:

- **By particle number the field is dominated by γ-rays (89 %)**, about 40 % of them softer than 100 keV. They give the bulk of the instrument's count rate. Their origin: primary cosmic protons and nuclei collide with nitrogen and oxygen nuclei at altitudes of about 15–20 km and produce pions; the neutral pion decays immediately into two high-energy γ-rays. Each such photon starts an electron-photon shower: pair production and bremsstrahlung alternately split the energy among an ever-growing number of ever-softer particles. Below a few MeV the shower "cools": photons scatter repeatedly on air electrons (Compton scattering) and lose energy, piling up in the tens-to-hundreds-of-keV range until they are absorbed by the photoelectric effect on nitrogen and oxygen (noticeable below a few tens of keV). Hence the abundance of soft photons. The 511 keV line is positron annihilation — positrons born in the same showers and in π⁺ → μ⁺ → e⁺ decay, in the air around the aircraft. On top of this there are weak nuclear γ lines from excited nitrogen and oxygen nuclei (e.g. 4.44 and 6.13 MeV) — the tabulated PARMA field does not resolve them as separate lines except for 511 keV. The energy distribution of γ-rays in PARMA is a single analytical fit to PHITS shower calculations (Sato 2015 [1]); the model does not separate photons by production mechanism, so the contribution of atmospheric π⁰ decay is already folded into the photon spectrum. Angular distributions follow Sato 2016 [2].
- **Neutrons — 6 % of the flux**, but the energetically broadest component: from thermal (below 1 eV) to hundreds of MeV. Energetic neutrons (above 10 MeV, 23 % of the neutron flux) produce secondary particles in the structure and in people.
- **Electrons and positrons — 4 %**; positrons are half as numerous as electrons, and when they stop they give the 511 keV line inside the aircraft.
- **Muons and protons — under 1 % of the flux**, but almost all of them are more energetic than 10 MeV and fly straight through the crystal — they are what forms the count rate at 3–10 MeV.
- **Not included in the PARMA field** (and not in the calculation): helium nuclei and heavier, pions, kaons. The literature suggests their contribution to the count below 10 MeV is small, but this has not been checked by our own calculation.

## Calculation model

**The aircraft** is an Airbus A320 or A321 (from livery photos, not established more precisely). Both belong to the same family and have the same fuselage cross-section: narrow-body, outer diameter 3.95 m, cabin width 3.70 m, single aisle, three seats on each side; the fuselage is aluminum (skin and structural framing of 2024/7075-type alloys), not composite. The A321 is 6.9 m longer than the A320, but the model uses only an 8 m section around the instrument's location, so the two are not distinguished for the calculation.

Cosmic particles at 9800 m (the PARMA/EXPACS field: neutrons, protons, muons, electrons, positrons, gamma rays) strike an 8 m section of the A320 fuselage, pass through the skin, trim and the layer of people, and reach the instrument. The instrument sits by the window at the level of a seated passenger. The count is computed in two stages: stage I — particle transport through the cabin up to a small "tube" around the instrument (recording all particles entering it), stage II — replaying the recorded particles on a detailed model of the instrument.

@@SCHEME@@

| Element | Dimensions and position (fuselage cross-section) | Material |
|---|---|---|
| Skin and structural framing | outer radius 197.5 cm; 1.85 g/cm² (estimate ±30 %) | aluminum alloy (2024-type: Al 93.5 %, Cu, Mg, Mn) |
| Thermal-acoustic insulation | 0.5 mm film | PET (Mylar), 1.38 g/cm³ |
| Cabin trim | 5 mm, inner radius 185 cm | polyethylene (CH₂), 0.96 g/cm³ |
| Floor | top 60 cm below the axis, 1.1 cm | aluminum alloy, 2.78 g/cm³ (beam-equivalent) |
| People | solid layer 340 cm wide, from the floor to 60 cm above the axis (120 cm tall), 0.11 g/cm³ (153 passengers at 84 kg, 85 % load factor) | ICRU-44 soft tissue (H, C, N, O, Na, P, S, Cl, K) |
| Cargo hold | trapezoid below the floor, 75–150 cm half-width, 0.078 g/cm³ | cellulose (mass-equivalent) |
| Air | cabin 0.9 mg/cm³ (cabin pressure), outside 0.364 mg/cm³ (9800 m altitude) | air |
| Instrument region | cylinder ("tube") of radius 14 cm, length 300 cm, axis along the cabin, centre 150 cm left of the axis and 15 cm above it, inside the layer of people | cabin air |
| Along-axis segment | 8 m along the axis (±4 m) | — |

@@RENDER@@

The stage-II instrument model is the `ASN16Detector` geometry (donor `geant4-detector-models/detectors/AtomSpectra-Nano-16-PRO`): a CsI(Tl) crystal 18×15×60 mm in an aluminum housing 86×42×25 mm with 1.50–1.95 mm walls, a reflective wrap (PTFE + Al foil) on five faces of the crystal, and an electronics board underneath — reproduced in full, without simplification.

**Why a "tube" in the cabin instead of the instrument itself.** The tube is not a model of the instrument, but an imaginary recording surface. The instrument itself (a detailed АтомНано 16 model: crystal, housing, reflector, board, photomultiplier) is used in full — at the second stage. Splitting into two stages is needed for statistics. A 16 cm³ crystal in a cabin 3.7 m wide is a tiny target: of the particles launched at the fuselage cross-section, a negligible fraction reaches the crystal, and almost all the counting time would go to particles that never touch the instrument. A tube of radius 14 cm and length 3 m has a side area of about 26 400 cm² — roughly 600 times the crystal's surface area (45 cm²). So for the same time it receives hundreds of times more particles. Using the tube's length along the cabin exploits the fact that the cabin is uniform along its length: the instrument can be placed at any point along the tube, and all these positions are equivalent. At the second stage the recorded particles (type, energy, direction, entry point and time) are replayed on the instrument model placed at a random point along the tube; each event is used K times at different instrument positions (K = 10 for γ-rays, 30 for electrons and positrons, 100 for neutrons, protons and muons).

**The two-stage scheme in detail.**

1. *Stage I — the field in the cabin.* PARMA particles are launched from a disc of radius 4.8 m onto the fuselage section, and pass through the structure and the people with all secondary processes. Every particle (primary or secondary) crossing the side surface of the tube inward is recorded: type, energy, direction, angle and coordinate along the axis, time, event number, and the process and material of birth. The tube's end caps are not recorded. The "equivalent flight time" T = N / (Φ·πR²), where N is the number of launched particles, Φ is the flux of the PARMA component, R is the disc radius.
2. *Stage II — the instrument.* For each stage-I event, an instrument position along the tube is chosen at random: one of the event's records is picked and a random offset within a window ℓ = 120 cm is added to its coordinate; all records of that event that fall within the window are replayed on the instrument, with their energies, directions and times. This preserves the simultaneity of particles from a single shower (important for summing within one trigger). The weight of each replay is w = n·ℓ / ((L − ℓ)·m), where n is the number of records for the event, m is the number of records within the window, L = 300 cm is the tube length; it returns the correct count rate per unit length. Each event is replayed K times at different instrument positions (the weight is divided by K; K values — see the section on dips in the plots). Count rate = Σ weights / T.
3. *Triggering.* Energy depositions in the crystal separated by less than 1 μs are summed into a single count. Delayed decays (e.g. ¹²⁸I, 25 min) are weighted by 1 − t/T_flight, where T_flight = 80 min: the fraction of the flight remaining after the decay moment.
4. *Why a 120 cm window and not smaller.* A photon entering the tube at a shallow angle travels tens of centimetres along the tube axis before reaching the crystal. With a 30 cm window such photons were lost and the count was underestimated by 22 %; this was caught by the closure check. At 120 cm the residual discrepancy is 2.6 %.

The validity of this substitution was checked with a closure calculation: an empty cabin, direct photon counting at the instrument versus the "tube → instrument" scheme — a 2.6 % discrepancy. The price of this approximation: inside the tube there is air, not people — particles born in tissue closer than 14 cm to the instrument, and backscatter from the body of the instrument's owner, are not accounted for (see caveats).

The particle source is a disc of radius 4.8 m facing the fuselage section; directions and energies are sampled from PARMA tables for 9800 m altitude, W = 76, cutoff rigidity 17.28 GV. The "equivalent flight time" is obtained by dividing the number of particles by the flux: 100 s for each field component. The instrument is АтомНано 16: an 18 × 15 × 60 mm CsI crystal in an aluminum housing, a PTFE reflector, foil, a board with copper and solder, a silicon photomultiplier; the calculation gives the energy deposited in the crystal, to which the energy resolution is then applied.

Physics: Geant4 11.4.2, the FTFP_BERT_HPT list (neutrons up to 20 MeV via evaluated nuclear data with thermal scattering on hydrogen and aluminum), electromagnetic physics option4 (with fluorescence, Auger and PIXE in the instrument), radioactive decay, capture of stopped muons, photonuclear reactions. ⁴⁰K decay inside people is sampled separately.

## Spectrum and decomposition

The black line is the full spectrum (count rate 71.19 ± 0.20 counts/s, about 256,290 ± 730 counts/h in 20–10,000 keV; a further 2.0 s⁻¹ below 20 keV and about 10,565 counts/h above 10 MeV are not shown on the plot). Components are broken down by the type of the primary field particle (neutrons, protons, muons, electrons, positrons, gamma rays) and by ⁴⁰K in people's bodies; "origin" is the particle type / process / material where the particle that hit the instrument was born. Toggle the rows you need.

@@CHART@@

## Detection mechanism: where the counts come from

**What the instrument registers.** The scintillator records only the energy deposited in the CsI crystal during a single trigger (a 1 μs window in the model). The instrument cannot determine the type of the primary particle: it registers only the total energy delivered to the crystal by that particle and its secondary products. All energy depositions within one trigger are summed into a single count at the corresponding energy.

**Two cuts through the same spectrum.**

- *By field particle* — by the type of cosmic-ray particle that arrived from outside the aircraft and started the whole chain. For example: a neutron was captured in a passenger's body, the hydrogen emitted a 2223 keV γ-ray, it reached the instrument — this is a count in the "neutrons" layer.
- *By origin* — by the particle that entered a small region around the instrument (a cylinder of radius 14 cm), and by which process and where it was born. The same example here falls under "γ-rays: neutron capture · in the passenger layer". If several particles deposited energy within one trigger, the count is assigned to the one that contributed the most energy.

**How to read the origin labels.**

- "arrived from outside air (primary) · air outside the airframe" — a cosmic-field particle that arrived from outside and did not turn into other particles before reaching the instrument's region. Along the way it could have scattered in the skin or in people and lost energy — the label is kept regardless.
- "electron bremsstrahlung · in the passenger layer / in the skin" — a γ-ray emitted by a fast electron or positron decelerating in passengers' bodies or in aluminum.
- "positron annihilation · …" — a pair of 511 keV γ-rays from a positron that stopped in that material.
- "neutron capture · in the passenger layer" — capture γ-rays from a slowed-down neutron absorbed by a nucleus, mostly hydrogen (2223 keV).
- "pair production · …" — an electron or positron born from a high-energy γ-ray.
- "other origins" — all "particle + process + material" combinations that did not make the top dozen layers. Broken down below.

**What is in "other origins".** This is 177 small combinations out of 189, together 5.0 % of the total count (about 12,720 of 256,290 counts/h in 20–10,000 keV). Each one individually contributes less than 0.6 % of the count. The largest single one is bremsstrahlung γ-rays from electrons in the floor, 0.5 %. By the type of particle that reached the instrument:

| Type | counts/h | share of "other" | what it is |
|---|---|---|---|
| γ-rays | 8,830 | 69.4 % | secondary photons born in the structure, in people and in luggage |
| electrons and positrons | 2,545 | 20.0 % | secondary electrons knocked out in the structure |
| protons | 760 | 6.0 % | primary cosmic-ray protons and secondary protons from nuclear reactions |
| neutrons | 480 | 3.8 % | secondary fast neutrons from nuclear reactions |
| other (muons and other particles) | 100 | 0.8 % | muons and rare secondary particles (not split by type in the calculation) |

By the process that produced the particle (for γ-rays — what emitted them, for electrons — what knocked them out):

| Process | counts/h | share of "other" | what it is |
|---|---|---|---|
| positron annihilation (γ) | 2,570 | 20.2 % | pairs of 511 keV γ-rays from positrons that stopped in the floor, luggage, trim, insulation — outside the top dozen by material |
| bremsstrahlung (γ) | 2,430 | 19.1 % | γ-rays from electrons in minor materials (luggage, floor, trim) |
| neutron inelastic scattering (γ) | 1,750 | 13.7 % | γ-rays from excited nuclei after neutron scattering; over half in the passenger layer (carbon, oxygen, nitrogen nuclei), the rest in the skin and floor |
| pair production (electrons and positrons) | 1,390 | 11.0 % | pairs produced by γ-rays above 1.022 MeV; mostly in the skin |
| radioactive decay (γ) | 1,280 | 10.1 % | ⁴⁰K γ-rays in people and decays of nuclei produced by cosmic particles |
| primary protons | 670 | 5.2 % | cosmic-ray protons reaching the instrument without interacting |
| electron ionization (electrons) | 560 | 4.4 % | delta electrons from fast electrons and positrons |
| neutron inelastic scattering (neutrons) | 420 | 3.3 % | neutrons scattered with energy loss by nuclei in the structure |
| Compton scattering (electrons) | 370 | 2.9 % | recoil electrons knocked out by γ-rays; mostly in the passenger layer |
| neutron capture (γ) | 320 | 2.6 % | capture γ-rays from neutrons absorbed by nuclei in luggage and trim (in "other" — not in the top dozen) |
| proton inelastic scattering (γ) | 290 | 2.3 % | γ-rays from excited nuclei after proton collisions |
| muon ionization (electrons) | 160 | 1.2 % | delta electrons from muons |
| the rest | about 500 | 4.0 % | photonuclear reactions, hadron elastic scattering and other processes |

By the material where the particle was born: people — 27 % of "other", floor — 22 %, skin — 17 %, cargo hold — 11 %, trim — 9 %, air outside the airframe — 8 %, insulation — 3 %, the rest — 2 %.

By energy, "other" gives 3 % of counts in the 20–100 keV band, 5 % in 100–500 keV, 9 % in 500–1000 keV, and 14 and 13 % in the 1–3 and 3–10 MeV bands. In other words, the higher the energy, the larger the contribution of numerous minor secondary processes: at low energies primary radiation dominates, at high energies hundreds of different nuclear and electromagnetic cascades do.

**Indirect neutron detection.** A neutron does not ionize or excite the scintillator directly. In the spectrum it shows up only through secondary particles. In this model neutrons give 4.2 % of the count (about 10,800 counts/h in 20–10,000 keV) via three routes:

1. **Neutrons from outside reach the instrument itself (about 6,200 counts/h, more than half of the neutron contribution).** They interact with the crystal's nuclei. Capture of a neutron by iodine ¹²⁷I or caesium ¹³³Cs gives a cascade of γ-rays and conversion electrons directly in the crystal: peaks around 137 keV (levels of ¹²⁸I at 133.6 and 137.9 keV and ¹³⁴Cs at 138.7 keV), 175 keV (¹³⁴Cs at 173.8 and 176.4 keV) and 205 keV. Inelastic scattering on ¹²⁷I excites the 57.6 and 202.9 keV levels — hence the peak near 59 keV. The resulting ¹²⁸I decays with a 25 min half-life, giving a β continuum up to 2.1 MeV; the model accounts for it with a correction for the fraction of decays occurring after landing. Another contribution is recoil nuclei from elastic scattering: in the model their energy is counted in full (light quenching is not modelled), so this contribution is overestimated. Peak positions were identified from product-nucleus levels; a breakdown by specific reaction inside the crystal was not computed separately.
2. **Neutrons slow down and are captured in people and the structure.** Then a capture γ-ray reaches the instrument: 2223 keV from hydrogen, lines from nitrogen, chlorine and aluminum. In the "by origin" cut these counts are attributed to γ-rays, in the "by field particle" cut — to neutrons.
3. **Inelastic scattering of fast neutrons on aluminum, oxygen, carbon** gives γ-rays from excited nuclei (skin and floor — about 170 counts/h, people — about 190 counts/h).

A real instrument cannot distinguish a neutron count from a gamma count: a CsI(Tl) detector without pulse-shape discrimination does not register neutrons separately. So "neutrons" on the plot is a contribution to the spectrum that exists only in the model. In a measured spectrum it could only be noticed through indirect signatures, chiefly the peaks near 58 and 137 keV.

## Contribution of field components and structural materials by energy band

The count rate in each energy band is broken down in two independent ways. The first table gives the breakdown by the type of primary cosmic-ray particle that started the chain of events leading to the count (neutrons, protons, μ⁺, μ⁻, electrons, positrons, γ-rays, and γ-rays from ⁴⁰K decay in people's bodies); values are fractions of the total count rate in the band, given in the "total" column. The second table gives the breakdown by the medium where the track of the registered particle started (for secondary particles — the medium where it was born): "air outside the airframe" — the track started in the atmosphere or in the air outside the fuselage, i.e. the particle is a primary field particle; it may then have scattered repeatedly in the skin and the passenger layer, but is identified as primary, so a significant fraction of such photons already have energy below the primary value; "people" — the soft-tissue layer modelling the passengers; "skin" — the skin and structural framing made of aluminum alloy; "trim" and "floor" — the cabin's inner trim and floor; "insulation" — the polyethylene terephthalate thermal-acoustic insulation; "luggage" — the cargo hold (material mass-equivalent to cellulose); "air near the instrument" — the air in the recording region around the instrument; "other" — the remaining volumes, including cabin air. Table rows are energy bands; the count is assigned to whichever source contributed the largest share of energy deposition in the crystal for that trigger.

@@BANDS@@

## Lines and spectral regions — what they are

1. **20–60 keV, with a maximum around 61 keV.** Not a line, but a continuum. Gamma rays of the PARMA field (primary atmospheric photons, 10.7 cm⁻²·s⁻¹) reach the instrument after passing through the skin, trim and passenger layer; rising absorption toward low energies (photoelectric effect in aluminum and soft tissue) competing with rising flux gives a maximum around 60 keV; below 30 keV the count drops to a third of the maximum, and this is mostly Compton and X-ray background rather than photons of those energies themselves. The response shape below 60 keV is set by the iodine (33.17 keV) and caesium (35.99 keV) absorption edges and the escape of iodine K X-rays (28.6 keV) from the crystal.
2. **100–500 keV.** The Compton continuum of the same photons and bremsstrahlung from secondary electrons and positrons (the eBrem process in the passenger layer and the skin); neutrons contribute via inelastic scattering and capture — up to 9 % in the 300–500 keV band.
3. **511 keV.** Positron annihilation: the primary field's own 511 keV line (about 58 % of the peak area; area about 2900–3000 counts/h, different baseline-fitting methods give 2870–2970), annihilation of positrons stopped in people's bodies (about 26 %) and in the skin. The capture line on phosphorus ³¹P (512.65 keV), with a FWHM of 36.6 keV at 511 keV, cannot be separated from annihilation; the Ga, Br, Ca lines from the airframe table are not in the model — these elements are not present in the geometry.
4. **1461 keV.** ⁴⁰K in people's bodies (54.8 Bq/kg per UNSCEAR, about 196 kBq for the whole layer): 0.38 s⁻¹ across the whole spectrum. In the pure K-40 component the peak sits at 1460.5 keV (area about 50 counts/h); in the full spectrum it sits on the primary-particle continuum (about 7 counts/(h·keV)); the peak area there is about 50 counts/h (48–57 across independent fits; the table below gives 78 ± 14 for a ±3 keV window and a parabolic baseline), significance about 5σ; on the smeared spectrum the peak is barely visible.
5. **2223 keV.** Neutron capture by hydrogen ¹H(n,γ)²H: about 92 % of the peak area comes from neutron capture (the nCapture process): mostly in the passenger layer (68–81 % depending on the window), the rest in trim and the cargo hold; neutrons are slowed down and captured on hydrogen (thermal S(α,β) scattering is included).
6. **3–10 MeV.** Mostly not lines. By primary field particle: 40 % of the count is charged particles (electrons and positrons 26 %, muons 11 %, protons 3.5 %) flying straight through the crystal and leaving ionization losses; 56 % is gamma rays, 4 % is neutrons (see the first table). Muons and protons have a maximum around 8.5–9 MeV — the "minimum ionizing" peak (the Landau tail extends above 10 MeV; about 10,565 counts/h of all particles lie beyond the plotted range): 1.24 MeV·cm²/g × 4.51 g/cm³ × a crystal chord length of 1.5–1.8 cm. Capture lines on Al (7724 keV) and Cu (7637 keV; both elements are in the skin alloy and the board) are not resolved here: a 1.5 cm crystal does not fully absorb γ-rays of these energies. Above 10 MeV — beyond the plotted range.

## Lines found in the spectrum

The search runs on the raw energy-deposition spectrum (1 keV per channel, before convolution): a peak is found in the spectrum, and lines from the airframe line table are matched to it within 4 keV; the baseline under the peak is a parabola fitted to the logarithm of the count on the flanks, ±3 keV window. Candidates from elements and decay chains not present in the model (Ga, Br, Ca, Gd, Fe, Ti, Zn, Ni; the Th/U series) are removed from the table. Significance is computed from the calculation's own statistics (Σw²), threshold 4σ; on the smeared spectrum, with a FWHM of 26–60 keV, narrow lines merge into the continuum, so lines are not searched for there. Several rows for one peak are unresolved candidates sharing one combined area. Three outcomes are possible, not two: (1) the peak is significant AND coincides with a tabulated line — the "Lines found in the spectrum" section; (2) the peak is significant, but no tabulated line lands on it within 4 keV (either there is no line at that energy in the database, or there is one but its centre is shifted further than the tolerance) — the "Peaks outside the list" section below, with the nearest database candidate given for each, if one was found in a wider window; (3) a tabulated line exists, but at its exact energy there is no statistically significant excess — such a line appears in neither section (example — 57.6 keV, discussed in the text above). The Σw² statistical error for peaks outside the list is underestimated (a few histories carry disproportionate weight), so a coincidence with a candidate is a hint, not a confirmed identification. For some significant peaks outside the list there is no candidate nearby at all — for example, the peak near 662 keV (7.3σ, not to be confused with Cs-137: this nuclide is not in the model) has no entry in the reaction database; such peaks remain unidentified.

@@LINES@@

## Neutron-induced lines — full summary

Only lines that passed the significance filter in the summed spectrum are shown further up the page. Below are all 150 (n,γ) and (n,n′) reactions from source database [1] for structural, crystal, human-body and cargo materials, with a calculation status: found and labelled, seen but unidentified, a candidate below the statistical threshold in the summed spectrum, or outside the search window (the adjacent energy is occupied). Method details are at the start of the "Lines found in the spectrum" section above.

@@NEUTRONLINES@@

## Comparison with a measurement from the same flight (preliminary)

A measurement from this same flight exists from another, smaller-volume instrument with a CsI crystal (the instrument's specifications are not published), recorded on 20.09.2026 (a 7395 s spectrogram). A cruise-altitude segment (3700–4400 s of the recording) was compared with the calculation for АтомНано 16 (16.2 cm³). The crystals differ, so the comparison is by count per 1 cm³ of volume. For hard particles (above 1 MeV) this is a reasonable approximation — they pass straight through the crystal. For soft photons (below 300 keV) it is crude: they are absorbed near the surface, and a smaller crystal has more surface area per unit volume. An accurate comparison would require a separate calculation for that instrument's geometry.

| Energy band, keV (instrument scale) | Measured, counts/s per 1 cm³ | Calculated, counts/s per 1 cm³ | Calc./measured |
|---|---|---|---|
| 30–100 | 1.42 | 1.82 | 1.28 |
| 100–300 | 1.77 | 1.59 | 0.90 |
| 300–1000 | 0.34 | 0.44 | 1.29 |
| 1000–3000 | 0.20 | 0.19 | 0.95 |
| above 3300 (the instrument's last channel) | 0.38 | 0.42 | 1.11 |

What this means: in every band the calculation differs from the measurement by no more than 30 %. The "above 3300 keV" band is interesting because the calculation gives 6.8 counts/s for 16.2 cm³, about half of it gamma rays, and such a large hard-gamma contribution raises a question. The measurement confirms it at the level: 1.15 counts/s in the second instrument's last channel. The reason for the large hard-gamma contribution is the field itself: the flux of γ-rays above 3 MeV in this calculation's source is about 1.5 cm⁻²·s⁻¹ (0.67 in 3–10 MeV, 0.74 in 10–100 MeV, 0.11 above 100 MeV), three times the flux of charged particles of the same energies. Source of these numbers: the γ-ray field file `tables/dadsgn/w76_h98_ip33.tab`, built by the PARMA 4.10 program [3] for the flight conditions (W = 76, Rc = 17.28 GV, depth 278.7 g/cm²) using model [1]; the numbers were obtained by integrating this table with the script `scripts/field_table.py` and are not given in this form in the papers themselves [1, 2]. The total γ flux (> 10 keV) of 10.74 cm⁻²·s⁻¹ matched an independent run of the contour's `parma_tables.exe`. A crystal up to 6 cm long interacts with a noticeable fraction of such γ-rays (the photoelectric effect is not the dominant process, but pair production and subsequent absorption give counts in the MeV region and in pile-up). We have no direct calculated or measured separation of γ and charged particles above 3 MeV — agreement at the summed level does not prove the breakdown is correct.

## Contribution of wings, fuel and engines (estimate)

The main calculation does not include the wings, fuel or engines. To estimate what they change, a separate comparative calculation was made: two identical setups with a source disc of radius 9 m (needed to cover the wings), one with the fuselage only, the other with two fuel tanks added (blocks 650 × 400 × 15 cm, kerosene at 0.8 g/cm³, 3.1 t each; centres at 5.25 m from the axis and 1.7 m below it), wing skin of aluminum alloy 7 mm thick, and two engines (simplified as cylinders of radius 90 cm and length 250 cm, made of titanium at density 0.36 g/cm³, 2.3 t each). Three field components were computed — neutrons, γ-rays, electrons and positrons — with 15 s of equivalent flight time per setup; protons and muons were not computed. The table gives the ratio of the count rate "with fuel and wings" to "without them" by band, and its statistical error.

@@FUEL@@

Conclusion: the total count rate in 20–10,000 keV is practically unchanged (0.99 ± 0.01). Only γ-rays show a noticeable change: minus 3 % (1 % error) — the fuel tanks partially shield γ radiation. The neutron component in the table shows plus 9 % (6 % error) and plus 25 % ± 16 % in 20–100 keV — this is within two errors, so the excess cannot be considered significant. In total neutrons contribute only about 4 % of the count, so even a noticeable rise in the neutron part does not shift the overall result. Caveats: a short history (15 s versus 100 s in the main calculation), a simplified shape for the tanks and engines, protons and muons not included, and a different source-disc size in this setup; so the table shows ratios, not absolute values, and does not replace a calculation with the real wing geometry.

## Dose rate for a passenger (estimate)

The instrument computes the spectrum in the crystal, not the dose in the body. To estimate the dose rate, a sphere of soft tissue (ICRU-4) of radius 12 cm was computed at the same point along the tube instead of the instrument — the same transport method, the same stage-I events. Three layers were taken: skin (0.07 mm from the surface, as the individual equivalent Hp(0.07)), a depth of 10 mm (as Hp(10)), and a depth of 30 mm; "whole sphere" is the volume average. The equivalent dose is weighted by the wR factor per ICRP 103 (for neutrons — by the energy of the particle entering the tube, not by its energy in the body). ⁴⁰K in people's bodies is not included in the dose: this is background γ radiation from potassium-40 in passengers' own bodies, independent of flight altitude and not part of the cosmic-radiation dose that is the subject of this section.

@@DOSE@@

The absorbed dose from all components together is about 0.70 μGy/h, the equivalent dose about 1.07 μSv/h (at a depth of 10 mm; the spread between layers is within 5 %). At the instrument's count rate of 71.2 counts/s, about 8300 counts/cm³ accumulate in an hour if converted to the 16.2 cm³ crystal — but that is a count, not a dose: the units differ, there is no direct conversion. Half of the equivalent dose comes from γ-rays and electrons/positrons; neutrons contribute about a third — noticeably more than their share of the instrument's count (6 % of the field flux), because neutrons have a high wR factor (5–20 versus 1 for photons). Caveats: the body is a uniform sphere without organs and without self-shielding; orientation and position match the instrument's, not averaged over posture; the effective dose (organ-weighted) was not computed — only the equivalent dose at a point. The statistical error is the same as for the instrument spectrum (fractions of a percent for γ, a few percent for the other components); the geometric systematics of the body are not included in it.

## How it was calculated

1. Source — PARMA/EXPACS 4.10 (Sato 2015, 2016): neutrons, protons, μ±, e±, γ at 9800 m; W = 76 (Oulu, 20.09.2026), cutoff rigidity 17.28 GV (Da Nang 17.07, Ho Chi Minh City 17.45), depth 278.7 g/cm².
2. Stage I (Geant4 11.4.2, FTFP_BERT_HPT, EM option4, radioactive decay, thermal S(α,β) scattering on the hydrogen of people and plastic and on aluminum): the field across the A320 fuselage cross-section — skin and framing 1.85 g/cm², trim, floor, a people layer of 0.11 g/cm³, cargo hold, cabin air 0.9 mg/cm³; recording of particles entering the instrument's region. 100 s of field history; ⁴⁰K in people — a separate run.
3. Stage II: replaying the recorded particles on the АтомНано 16 model (CsI 18×15×60 mm, housing, board, SiPM), with fluorescence, Auger and PIXE; production threshold 0.05 mm; for each recorded event the instrument position is sampled K times (weight divided by K). The replay window along the tube axis is 120 cm (at 30 cm the count was underestimated by 22 % due to lost shallow-angle photons). Closure check in an empty cabin: the replay is 2.6 % below direct photon counting at the instrument (z = −0.9); the arithmetic matches the fraction of photons at less than 13° to the tube axis (1 − cos 13° = 2.6 %), but this is a hypothesis, not confirmed by a separate run.
4. Neutron-capture gamma cascades — a custom emitter built from a database of 58 isotopes, used instead of the standard one.

## AI-use disclosures

- **Who did this.** The problem setup, the choice of physical assumptions, the code (Geant4 programs, processing scripts), the text and the page were produced by an AI assistant (Claude) on the operator's instructions; part of the code was generated by a local language model and edited by hand. A human physicist did not independently check the calculation.
- **How it was checked.** The checks were also done by AI agents (two "sterile" passes with a clean context, recomputation of numbers by independent scripts, mutation testing of the acceptance tests) — this reduces the risk of errors but does not replace expert review. The first pass found 19 incorrect statements on the page, the second found 4 more; the 22 % count underestimate was found only by the closure check. So the calculation and the text may still contain errors that these checks do not catch.
- **Reference values** (K-edges, energy losses, potassium content, line energies) are partly given from the model's memory rather than from an opened primary source; where a source was not opened, this is stated in the audit reports (`audit/`).
- **The physics model was chosen by the AI:** the Geant4 physics list, the simplified cabin geometry, the people layer, substituting a custom emitter for neutron-capture gamma cascades, the "record → replay" scheme — these are authorial assumptions, not a standard validated against measurement. The closure check covers only photons.
- **How to use this.** The result is an approximate model for understanding what the spectrum is made of, not a measurement and not a basis for dose assessment, radiation-safety conclusions, or certification. Such purposes require independent review by a specialist and comparison with an in-flight instrument measurement.

## Model limitations

- **On the "dips" in the plots.** For minor components (neutrons above 6 MeV, muons and protons in the thousands-of-keV region) the field history was too short: the calculation contained only a handful of events per channel, and after convolution with the instrument response function, dips remained between them, dropping on the logarithmic axis. This was calculation statistics, not a physical dip. Fixed by two means: (1) stage I was recomputed with a longer equivalent flight time for neutrons, protons, muons, electrons and positrons (200–1600 s instead of 100 s, as independent histories with different seeds; the total time is summed, as for the other components); (2) stage II was repeated with more draws of the instrument position per recorded event (K = 100 for the main history of neutrons, protons and muons, K = 25–30 for the additional histories and for electrons and positrons, K = 10 for γ-rays and K-40; the error falls as 1/√K, checked on the neutron file). Result: within a 41 keV window (the convolution width at 662 keV) the relative statistical error is now ≤15 % everywhere for neutrons, protons, muons, electrons and positrons in 20–3000 keV — no statistical dips remain. For γ-rays and K-40 the error in this window is above 15 % in places, but this is no longer statistics (the γ statistics themselves are small, the error on wide bands is ≤3 %): the dips there are real spectral structure — lines and edges — not noise. Over wide bands (20–100, 100–300, …, 3000–10000 keV) the statistical error is now no higher than 4–5 % for any component. The per-bin error at 1 keV before convolution remains larger than in the 41 keV window — expected for a narrow channel.
- **No measurement.** The instrument did not record on this flight; the model has not been compared against an instrument. Comparing PARMA with an "official" EXPACS build for this rigidity is not possible — there is no dose-coefficient file.
- **Not included in the source:** heavy ions (He–Ni), pions and kaons (PARMA does not provide them), primaries above 100 GeV. For the range up to 10 MeV, the literature suggests their contribution is small, but this has not been checked by our own count.
- **The structure is simplified:** only an 8 m fuselage section (±4 m along the axis); wings, fuel and engines are not included in the main calculation (a separate estimate is in the wings-and-fuel section: the total count changes by about 1 %, and for the neutron component the effect is within the error); the cockpit is not included. Seats, overhead bins, panels, windows — mass-equivalent; the framing mass 1.85 g/cm² is an estimate ±30 %, affecting 20–100 keV. People are a solid layer, not individual bodies.
- **The tray table and particle backscatter from the operator's own body are not modelled:** the tray-table flag is enabled at stage II, but the donor geometry only builds a tray table inside the "penal" block (wt20), which is switched off — there is no tray table among the stage-II volumes; around the instrument there is a 14 cm-radius tube of air, no nearby tissue. This contribution has not been estimated.
- **⁴⁰K activity in people** of 54.8 Bq/kg was obtained by calculation from a potassium content of 1.8 g/kg (in the project — with a reference to UNSCEAR 2000, the value from the primary source was not opened); the ICRP reference value (140 g potassium per 70 kg) gives about 63 Bq/kg, i.e. K-40 could be up to 15 % higher. "Whole layer" means the ±4 m section (3.57 t), not the entire cabin.
- **Channels below 20 keV** (2.0 s⁻¹ in the model) are not within the plotted range, but on convolution with the instrument response function they leak into the lower channels (up to +55 % at 20.5 keV); the instrument's real threshold is not set in the model.
- **Line boundaries:** the significance of lines and peaks is from the calculation's own statistics (Σw²) without baseline uncertainty; for the neutron component the error is underestimated; lines from the airframe table whose elements are not in the model are not shown. Areas in the table are measured on the unsmeared spectrum (line positions in channels: 510.5, 1460.5, 2223.5 keV), while the plot is smeared — there a peak of the same area is wider and lower.
- **The uncertainty is statistics only.** "± 0.23 s⁻¹" is the statistics of a 100 s field history. Systematics are not included: framing mass ±30 %, PARMA/EXPACS uncertainty (no stated tolerance, estimated 20–30 %), W and cutoff rigidity, instrument orientation (long axis along the cabin — an assumption; the 3–10 MeV band also depends on it), K-40 ±15 %, closure ±3 %. The overall uncertainty of the absolute level is no better than 20–30 %.
- **The closure check is limited:** photons in an empty cabin, statistical error about 3 %; by band, 20–100 keV −5.1 ± 4.2 %, 1–3 MeV −22 ± 22 % (no significant confirmation). Neutrons and charged particles are not covered by this check. Photons below 10 keV entering the tube are discarded, and the tube's end caps are not recorded — "2.0 s⁻¹ below 20 keV" is a lower bound.
- **Quenching (kB = 0):** for recoil nuclei and α particles the light yield is not suppressed, so the neutron contribution in 100 keV – 1 MeV is overestimated. A real CsI(Tl) with a SiPM is noisy below 20–30 keV; there is no noise in the model.
- **The instrument's own background and K-40 outside people** (crystal, board, SiPM, luggage, seats, glass; cosmogenic activation of the crystal) are not included.
- **The instrument response is idealised:** the light yield is proportional to the energy deposit (no CsI non-linearity), there is no electronics threshold or dead time; FWHM ∝ √E from a single point — 41.6 keV at 661.7 keV.
- **9800 m throughout** is an upper bound: during climb and descent the flux is lower; the flight lasts about 1 h 20 min.
- **Physics approximations:** the thermal kernel for cellulose and PET is that of polyethylene; for Ar-36 and Ar-38 in air, Ar-37 and Ar-39 data are substituted; γ-cascade angular correlations are switched off.

## Sources

Calculation code and data: [github.com/Verter73/GEANT4/tree/main/flight-msk-tyo](https://github.com/Verter73/GEANT4/tree/main/flight-msk-tyo) (private repository).

1. Sato T. Analytical model for estimating terrestrial cosmic ray fluxes nearly anytime and anywhere in the world: extension of PARMA/EXPACS. *PLOS ONE* 10(12): e0144679 (2015). [doi:10.1371/journal.pone.0144679](https://doi.org/10.1371/journal.pone.0144679) ([full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4682948/)) — energy spectra of n, p, ions, μ±, e±, γ; relation of parameters to W and cutoff rigidity.
2. Sato T. Analytical model for estimating the zenith angle dependence of terrestrial cosmic ray fluxes. *PLOS ONE* 11(8): e0160390 (2016). [doi:10.1371/journal.pone.0160390](https://doi.org/10.1371/journal.pone.0160390) ([full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4973932/)) — angular distributions.
3. [PARMA/EXPACS](https://phits.jaea.go.jp/expacs/), program version 4.10 (21.03.2021), JAEA — the program from which the field tables were obtained.
4. [Oulu neutron monitor](https://cosmicrays.oulu.fi/) (Oulu NM) — the 20.09.2026 count, from which the index W = 76 was determined by the formula in [1].
5. In-flight spectra and cabin field models: [Fisicas, radiacode-flight-spectroscopy (GitHub, 2026)](https://github.com/Fisicas/radiacode-flight-spectroscopy); [Kochkin et al., ILDAS, JGR Atmospheres (2017/2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5815387/); [Matthiä et al., Space Weather 12 (2014)](https://doi.org/10.1002/2013SW001022); [Takada et al., ApJ 733:13 (2011)](https://arxiv.org/abs/1103.3436); Dachev Ts. et al., Fundamental Space Research 2009 (full-text link not found).
6. Geant4 11.4.2: [Agostinelli et al., NIM A 506 (2003) 250](https://doi.org/10.1016/S0168-9002(03)01368-8); [Allison et al., NIM A 835 (2016) 186](https://doi.org/10.1016/j.nima.2016.06.125) — particle transport, the FTFP_BERT_HPT list, electromagnetic physics option4.
