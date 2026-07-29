# External sources

Literature and documentation underpinning the decisions on geometry, physics and
methodology. A reference listed here means the source was actually used when
building the model or when choosing a particular number — not merely that it is
related to the topic.

The LSRM methodology (spectrum processing, dead time, intrinsic activity) is
kept separately, in `detectors/Gamma-1S/reference/lsrm/references/`: it is tied
to one specific verification chain rather than to modelling in general.

## Geant4

- **Application Developers Guide** —
  https://geant4.web.cern.ch/documentation/dev/bfad_html/ForApplicationDevelopers/index.html
  The main guide: construction of geometry and materials, physics lists,
  primary generation (GPS), user actions, UI macro commands. The `dev` version;
  where it disagrees with the behaviour of the 11.2.1 build, trust the build.

- **Source code** — https://github.com/Geant4/geant4
  The final authority when the documentation is ambiguous or silent. What
  actually gets checked against the code rather than the guide: NIST material
  compositions (`source/materials/src/G4NistMaterialBuilder.cc` — the origin of
  `G4_MAGNESIUM_OXIDE` and the rest), the set of processes in the physics lists
  (`G4EmStandardPhysics_option4`), and the behaviour of GPS commands under
  combinations of `/gps/ang/*`. Check the branch of the build tag — `v11.2.1`
  here, not `master`: a version mismatch is the usual reason why "the
  documentation says otherwise".

## Photon attenuation coefficients

The reference against which the toolkit's own cross sections are checked:
without it, the analysis of the low-energy end silently rests on the untested
assumption that Geant4 computes the attenuation in the entrance-face layers
correctly. Cross-check — `analysis/compare_xcom.py`, result —
`results/mu_geant4_vs_xcom.csv`.

- **NIST XCOM** — https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html
  Computation of μ/ρ by component (coherent, incoherent, photoelectric, pair
  production) at an arbitrary energy. It gives two totals — with and without
  coherent scattering; the two must not be confused, at 59.5 keV they differ by
  12%. The main attenuation calculation in `mucalc.cc` runs WITHOUT coherent
  scattering, and the "without coherent" column is the one that corresponds
  to it.
- **Hubbell J.H., Seltzer S.M.** Tables of X-Ray Mass Attenuation Coefficients
  and Mass Energy-Absorption Coefficients. NISTIR 5632.
  https://physics.nist.gov/PhysRefData/XrayMassCoef/cover.html
  Tables on a fixed grid: Table 3 — elements, Table 4 — compounds and mixtures.
  The quantity tabulated there is μ/ρ **with** coherent scattering (Section 2
  defines the total cross section as a sum in which σ_coh enters explicitly;
  verified on Al: the node 2.778·10⁻¹ at 60 keV matches precisely that column).
  MgO, NaI and rubbers are absent from Table 4 — for them XCOM is the only
  primary source.
- **NIST STAR, material compositions** —
  https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=243
  Natural rubber (ICRU-37): H 0.118371 / C 0.881629, ρ = 0.92 g/cm³. NIST has
  no ready-made μ/ρ table for rubber, so it is computed in XCOM in mixture mode
  from this composition.

## MgO reflector density

The reflector is a packed, pressed powder rather than solid oxide, and this is
what matters: its mass thickness governs the soft (low-energy) end of the
efficiency curve. A summary of the discussion is in the comment in
`detectors/Gamma-1S/geometry/G1SDetector.hh` next to the `mgoDensity` parameter.

- **Mendes B.M. et al.** Monitoring internal contamination from Occupationally
  Exposed Workers of an ¹⁸F-FDG production plant. *Braz. J. Rad. Sci.* 07-03A
  (2019) 01-12. Table 1 — composition and densities of the components of an
  MCNP model of NaI(Tl): ρ(MgO reflector) = 2.0 g/cm³ (after Mouhti 2017 and
  Salgado 2012).
- **Appl. Radiat. Isot.** — Detection efficiency evaluation for low energy of a
  NaI(Tl) scintillation detector.
  https://www.sciencedirect.com/science/article/abs/pii/S0969806X22003681
- **Appl. Radiat. Isot.** — A computational modelling of low-energy gamma ray
  detection efficiency of a cylindrical NaI(Tl) detector.
  https://www.sciencedirect.com/science/article/abs/pii/S0969806X21002310
- **Appl. Radiat. Isot.** — Optimization of the Monte Carlo simulation model of
  NaI(Tl) detector by Geant4 code.
  https://www.sciencedirect.com/science/article/abs/pii/S0969804317307479

  The three works above describe, among other things, the technique of
  CALIBRATING the reflector density from the ratio of the efficiencies of two
  soft lines. The technique is well known, but in this repository it is
  **deliberately not used**: a density fitted in this way absorbs everything
  else that is wrong at the entrance face (the Al thickness, air gaps, whatever
  is not in the technical drawing) and ceases to be a density, remaining an
  effective mass thickness under someone else's name. After such a fit,
  agreement on those same lines stops being a check. A physical value is taken,
  and the result is verified.

- **Bell S. et al.** MgO reflectance data for Monte Carlo simulation of
  LaBr₃:Ce scintillation crystals. *NIM A* 701 (2013) 44–53.
  https://www.sciencedirect.com/science/article/abs/pii/S0168900212012065
  On OPTICAL reflectance as a function of layer thickness; there is no density
  there. Useful in one respect: reflectance saturates with thickness — which
  means the 6 mm at the entrance face are not there for the sake of light
  collection.
