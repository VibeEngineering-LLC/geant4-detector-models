# Gamma-1S Reference Spectra in BecqMoni Format

A set of certified sources measured on the **Gamma-1S** spectrometric
system (detector UDS-GC-63×63-USB №SN-01, NaI(Tl) 63×63 mm, Lsrm SpectraLine
software). Each file is self-contained: the sample spectrum and the
background of the same geometry are stored in one XML — the background
is recorded as `<BackgroundEnergySpectrum>`, so there is no need to load
a separate file.

Total records: **40** across 5 geometries.
The source `.spe` files are in the sibling folder `../reference_kits/`;
the conversion is done by `scripts/convert/reference_kits_to_becqmoni.py`
(idempotent).

## How to read the tables

**Activity** is taken from the `COMMENT` field of the LSRM header — this
is the source's passport value as of the reference date, not a
measurement result. Where "Bq/kg" is given, this is the specific
activity of the filler material; where "Bq" is given, this is the
activity of a point source. Decay has not been recalculated to the
measurement date.

A separate note on the РИСН-379 mix: in the `.spe` file, the `COMMENT`
field holds only **one** line, so the table lists only one nuclide out
of four (Am-241 for Marinelli, Cs-137 for Denta and Petri-60). The full
source composition is not recorded in the spectra.

**Units.** BecqMoni stores `Weight` in kilograms and `Volume` in liters,
while LSRM `.spe` uses grams and milliliters. The conversion divides by
1000; the "mass" column below is the original value, in grams.

**Dead time** is calculated as `(real − live) / real`, from the times
in the header.

## Denta_120mL

Denta vial 120 mL, 0 cm. For the РИСН-379 mix, the header shows
`GEOMETRY=Дента-100` — an operator typo; the Denta vessel is always
120 mL. The 100 g in `RAWMASS` is the fill mass, not the volume.

| Nuclide | Source | Passport activity | Mass, g | Measurement date | Live, s | Dead, % | Channels | File |
|---|---|---|---:|---|---:|---:|---:|---|
| Cs-137 | Cs137_420-7-14 | Cs-137 1 760 Bq/kg ±5 % as of 2002-05-24 | 68 | 2024-10-31 | 4820 | 0.04 | 1024 | `sample_Cs137_420-7-14_Дента-120мл_0cm.xml` |
| K-40 | K40_420-7-20 | K-40 2 530 Bq/kg ±6 % as of 2002-05-24 | 79 | 2024-10-31 | 57055 | 0.04 | 1024 | `sample_K40_420-7-20_Дента-120мл_0cm.xml` |
| Mix_AmTiCsEu | №SRC-04_Am-Ti-Eu-Cs | Cs-137 2 210 Bq/kg ±5 % as of 2002-05-31 | 100 | 2016-05-31 | 3600 | 0.07 | 1024 | `sample_№SRC-04_Am-Ti-Eu-Cs_Дента-100.xml` |
| Ra-226 | Ra226_420-7-18 | Ra-226 1 780 Bq/kg ±6 % as of 2002-05-24 | 74 | 2024-10-31 | 6938 | 0.04 | 1024 | `sample_Ra226_420-7-18_Дента-120мл_0cm.xml` |
| Th-232 | Th232_420-7-17 | Th-232 1 940 Bq/kg ±6 % as of 2007-09-17 | 192 | 2024-10-31 | 6309 | 0.06 | 1024 | `sample_Th232_420-7-17_Дента-120мл_0cm.xml` |

Background: `background_bg_2016_empty_shield_point5cm.spe`, live time 54000 s.

## Marinelli_1L

Marinelli vessel 1 L, in direct contact with the detector.

| Nuclide | Source | Passport activity | Mass, g | Measurement date | Live, s | Dead, % | Channels | File |
|---|---|---|---:|---|---:|---:|---:|---|
| Cs-137 | M_cs_легкий_2001-2005 | Cs-137 1 890 Bq/kg ±5 % as of 1997-05-30 | 570 | 1999-08-04 | 2156 | 2 | 1023 | `sample_M_cs_легкий_2001-2005.xml` |
| K-40 | M_k_легкий_2001-2005 | K-40 2 540 Bq/kg ±10 % | 665 | 1999-08-04 | 2930 | 1.94 | 1023 | `sample_M_k_легкий_2001-2005.xml` |
| Mix_AmTiCsEu | Смесь_AmTiCsEu | Am-241 4 200 Bq/kg ±10 % as of 2002-05-31 | 1000 | 2016-05-30 | 3600 | 0.31 | 1024 | `sample_Смесь_AmTiCsEu_Маринелли.xml` |
| Ra-226 | M_ra_легкий_2001-2007 | Ra-226 1 850 Bq/kg ±10 % | 622 | 1999-08-04 | 3814 | 2.08 | 1023 | `sample_M_ra_легкий_2001-2007.xml` |
| Th-232 | Th232_420-7-17 | Th-232 1 940 Bq/kg ±6 % as of 2007-09-17 | 1600 | 2024-10-24 | 11359 | 0.21 | 1024 | `Th232_420-7-17_Маринелли_0cm.xml` |

Background: `background_bg_2016_marinelli_water_marinelli.spe`, live time 54000 s.
Background: `Фон закр кр вода_13.spe`, live time 46800 s.

## Petri_60mL

Petri dish 60 mL, 0 cm.

| Nuclide | Source | Passport activity | Mass, g | Measurement date | Live, s | Dead, % | Channels | File |
|---|---|---|---:|---|---:|---:|---:|---|
| Cs-137 | Cs137_420-7-14 | Cs-137 1 760 Bq/kg ±5 % as of 2002-05-24 | 34 | 2024-10-30 | 9074 | 0.04 | 1024 | `sample_Cs137_420-7-14_Петри-60мл_0cm.xml` |
| K-40 | K40_420-7-20 | K-40 2 530 Bq/kg ±6 % as of 2002-05-24 | 40 | 2024-10-28 | 58043 | 0.04 | 1024 | `sample_K40_420-7-20_Петри-60мл_0cm.xml` |
| Mix_AmTiCsEu | №SRC-04_Am-Ti-Eu-Cs | Cs-137 2 210 Bq/kg ±5 % as of 2002-05-31 | 60 | 2016-05-31 | 3600 | 0.06 | 1024 | `sample_№SRC-04_Am-Ti-Eu-Cs_Петри-60.xml` |
| Ra-226 | Ra226_420-7-18 | Ra-226 1 780 Bq/kg ±6 % as of 2002-05-24 | 37 | 2024-10-29 | 9051 | 0.04 | 1024 | `sample_Ra226_420-7-18_Петри-60мл_0cm.xml` |
| Th-232 | Th232_420-7-17 | Th-232 1 940 Bq/kg ±6 % as of 2007-09-17 | 96 | 2024-10-24 | 61929 | 0.05 | 1024 | `sample_Th232_420-7-17_Петри-60мл_0cm.xml` |

Background: `background_bg_2016_empty_shield_point5cm.spe`, live time 54000 s.

## Point_25cm

Point source, 25 cm from the detector end face.

| Nuclide | Source | Passport activity | Mass, g | Measurement date | Live, s | Dead, % | Channels | File |
|---|---|---|---:|---|---:|---:|---:|---|
| Am-241 | Am-241 42.13 | Am-241 118 000 Bq ±5 % as of 2013-12-03 | — | 2016-05-17 | 600 | 0.07 | 1024 | `sample_Am-241 42.13_Точечная-25см_25cm.xml` |
| Ba-133 | Ba-133 #SRC-07 | Ba-133 44 100 Bq ±2 % as of 2008-10-01 | — | 2016-05-18 | 1800 | 0.09 | 1024 | `sample_Ba-133 #SRC-07_Точечная-25см_25cm.xml` |
| Cd-109 | Cd-109 #SRC-07 | Cd-109 1 033 000 Bq ±2 % as of 2008-10-01 | — | 2016-05-18 | 1800 | 0.05 | 1024 | `sample_Cd-109 #SRC-07_Точечная-25см_25cm.xml` |
| Ce-139 | Ce-139_591 | Ce-139 191 000 Bq ±3 % as of 2013-12-01 | — | 2016-05-18 | 1800 | 0.05 | 1024 | `sample_Ce-139_591_Точечная-25см_25cm.xml` |
| Co-60 | Co-60 #SRC-07 | Co-60 107 800 Bq ±2 % as of 2008-10-01 | — | 2016-05-18 | 1800 | 0.11 | 1024 | `sample_Co-60 #SRC-07_Точечная-25см_25cm.xml` |
| Cs-137 | Cs-137 №SRC-01 | Cs-137 106 000 Bq ±3 % as of 2017-05-19 | — | 2024-10-22 | 3435 | 0.12 | 1024 | `sample_Cs-137 №SRC-01_Точечная-25см_25cm.xml` |
| Eu-152 | Eu-152 #SRC-07 | Eu-152 46 700 Bq ±2 % as of 2008-10-01 | — | 2016-05-18 | 3600 | 0.1 | 1024 | `sample_Eu-152 #SRC-07_Точечная-25см_25cm.xml` |
| Mn-54 | Mn-54_587 | Mn-54 224 000 Bq ±3 % as of 2013-12-01 | — | 2016-05-18 | 1800 | 0.07 | 1024 | `sample_Mn-54_587_Точечная-25см_25cm.xml` |
| Na-22 | Na-22 #SRC-12 | Na-22 229 000 Bq ±5 % as of 2022-11-14 | — | 2024-10-22 | 3187 | 0.44 | 1024 | `sample_Na-22 #SRC-12_Точечная-25см_25cm.xml` |
| Th-228 | Th-228 №SRC-03 | Th-228 100 000 Bq ±3 % as of 2021-04-26 | — | 2024-10-24 | 7634 | 0.1 | 1024 | `sample_Th-228 №SRC-03_Точечная-25см_25cm.xml` |
| Y-88 | Y-88 №SRC-02 | Y-88 350 000 Bq ±3 % as of 2023-10-09 | — | 2024-10-22 | 2905 | 0.1 | 1024 | `sample_Y-88 №SRC-02_Точечная-25см_25cm.xml` |

Background: `background_bg_2016_open_lid_point25cm.spe`, live time 54000 s.

## Point_5cm

Point source, 5 cm from the detector end face.

| Nuclide | Source | Passport activity | Mass, g | Measurement date | Live, s | Dead, % | Channels | File |
|---|---|---|---:|---|---:|---:|---:|---|
| Am-241 | Am-241 42.13 | Am-241 118 000 Bq ±5 % as of 2013-12-03 | — | 2016-05-17 | 300 | 0.44 | 1024 | `sample_Am-241 42.13_Точечная-5см_5cm.xml` |
| Ba-133 | Ba-133 #SRC-07 | Ba-133 44 100 Bq ±2 % as of 2008-10-01 | — | 2016-05-17 | 1800 | 0.56 | 1024 | `sample_Ba-133 #SRC-07_Точечная-5см_5cm.xml` |
| Bi-207 | Bi-207 #SRC-11 04.2017 | Bi-207 97 000 Bq ±3 % as of 2017-05-25 | — | 2024-10-22 | 1800 | 2.02 | 1024 | `sample_Bi-207__176_04_2017_Точечная-5см_5cm.xml` |
| Cd-109 | Cd-109 #SRC-07 | Cd-109 1 033 000 Bq ±2 % as of 2008-10-01 | — | 2016-05-17 | 1800 | 0.06 | 1024 | `sample_Cd-109 #SRC-07_Точечная-5см_5cm.xml` |
| Ce-139 | Ce-139_591 | Ce-139 191 000 Bq ±3 % as of 2013-12-01 | — | 2016-05-17 | 1800 | 0.07 | 1024 | `sample_Ce-139_591_Точечная-5см_5cm.xml` |
| Co-57 | Co-57 #SRC-07 | Co-57 99 500 Bq ±2 % as of 2008-10-01 | — | 2016-05-17 | 1800 | 0.04 | 1024 | `sample_Co-57 #SRC-07_Точечная-5см_5cm.xml` |
| Co-60 | Co-60 #SRC-07 | Co-60 107 800 Bq ±2 % as of 2008-10-01 | — | 2016-05-17 | 1200 | 0.72 | 1024 | `sample_Co-60 #SRC-07_Точечная-5см_5cm.xml` |
| Cs-137 | Cs-137 #SRC-07 | Cs-137 94 200 Bq ±2 % as of 2008-10-01 | — | 2016-05-16 | 1800 | 0.74 | 1024 | `sample_Cs-137 #SRC-07_Точечная-5см_5cm.xml` |
| Eu-152 | Eu-152 #SRC-07 | Eu-152 46 700 Bq ±2 % as of 2008-10-01 | — | 2016-05-17 | 1800 | 0.68 | 1024 | `sample_Eu-152 #SRC-07_Точечная-5см_5cm.xml` |
| Mn-54 | Mn-54_587 | Mn-54 224 000 Bq ±3 % as of 2013-12-01 | — | 2016-05-16 | 600 | 0.34 | 1024 | `sample_Mn-54_587_Точечная-5см_5cm.xml` |
| Na-22 | Na-22_585 | Na-22 133 000 Bq ±3 % as of 2013-12-01 | — | 2016-05-16 | 600 | 2.32 | 1024 | `sample_Na-22_585_Точечная-5см_5cm.xml` |
| Th-228 | Th-228 #SRC-07 | Th-228 37 700 Bq ±2 % as of 2008-10-01 | — | 2016-05-16 | 1800 | 0.09 | 1024 | `sample_Th-228 #SRC-07_Точечная-5см_5cm.xml` |
| Y-88 | Y-88_589 | Y-88 234 000 Bq ±3 % as of 2013-12-01 | — | 2016-05-16 | 1200 | 0.05 | 1024 | `sample_Y-88_589_Точечная-5см_5cm.xml` |
| Zn-65 | Zn-65 #SRC-08.2019 | Zn-65 3 100 Bq ±3 % as of 2023-08-02 | — | 2024-10-22 | 2072 | 0.04 | 1024 | `sample_Zn-65__342_2019_Точечная-5см_5cm.xml` |

Background: `background_bg_2016_empty_shield_point5cm.spe`, live time 54000 s.

## Efficiency curves

Stored separately, in the `efficiency/` subfolder. Lsrm SpectraLine
format: `.efr` — efficiency points from reference sources, `.efa` —
description of the geometry and sample matrix.

**The curves are tied to a specific detector unit** (№SN-01) and cannot
be transferred to another instrument, even of the identical model — PMT
gain, crystal light yield, and the electronics chain all differ.

| Kit geometry | File | Contents | Volume, mL | Density, g/cm³ | Distance, cm | Assembled |
|---|---|---|---|---|---|---|
| Denta_120mL | `УДС-ГЦ-63х63-USB__SN-01_-_Дента.efa` | geometry and matrix description | 120 | 1.6 | 0 | 05-11-2024 |
| Denta_120mL | `УДС-ГЦ-63х63-USB__SN-01_-_Дента.efr` | efficiency points from reference sources | 120 | 1.65833333333333 | 0 | 05-11-2024 |
| Marinelli_1L | `УДС-ГЦ-63х63-USB__SN-01_-_Маринелли.efa` | geometry and matrix description | 1000 | 1.6 | 0 | 05-11-2024 |
| Marinelli_1L | `УДС-ГЦ-63х63-USB__SN-01_-_Маринелли.efr` | efficiency points from reference sources | 1000 | 1.6 | 0 | 05-11-2024 |
| Petri_60mL | `УДС-ГЦ-63х63-USB__SN-01_-_Петри.efa` | geometry and matrix description | 60 | 1.6 | 0 | 05-11-2024 |
| Petri_60mL | `УДС-ГЦ-63х63-USB__SN-01_-_Петри.efr` | efficiency points from reference sources | 60 | 1.6 | 0 | 05-11-2024 |
| Point_25cm | `УДС-ГЦ-63х63-USB__SN-01_-_Точечная-25см.efa` | geometry and matrix description | 0 | not essential | 25 | 06-11-2024 |
| Point_25cm | `УДС-ГЦ-63х63-USB__SN-01_-_Точечная-25см.efr` | efficiency points from reference sources | 0 | not essential | 25 | 06-11-2024 |
| Point_5cm | `УДС-ГЦ-63х63-USB__SN-01_-_Точечная-5см.efr` | efficiency points from reference sources | 0 | 0 | 5 | 06-11-2024 |

## Machine-readable index

`INDEX.json` — the same data plus the energy calibration coefficients,
real time, the number of peaks found by SpectraLine, and the instrument
GUID recorded in the XML.
