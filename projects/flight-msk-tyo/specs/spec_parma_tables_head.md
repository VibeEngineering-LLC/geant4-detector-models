You are a C++ engineer. Write ONE complete C++17 source file `parma_tables.cpp` and output ONLY the code (no markdown fences, no prose).

## Goal
A command-line tool that dumps the tabulated PARMA/EXPACS cosmic-ray source distribution for ONE particle and ONE explicit condition
(solar modulation W, vertical cut-off rigidity Rc, atmospheric depth d, geometry g) as a plain text file. The tool is derived from the
official `main-generator.cpp` (reproduced in full at the bottom of this message). It must NOT generate particles. It must reuse the
same mesh, the same PARMA functions and the same formulas as `main-generator.cpp`, changed only where this spec says so.

## Command line (exactly this)
    parma_tables <ip> <W> <Rc_GV> <depth_gcm2> <g> <out_file>
* `ip` integer particle id: 0 neutron, 29 mu+, 30 mu-, 31 e-, 32 e+, 33 photon, 1 proton (ids 1-28 are H..Ni; only `IangPart[ip] != 0` is allowed, otherwise print an error to stderr and exit with code 1).
* `W`, `Rc_GV`, `depth_gcm2`, `g` are doubles and are passed DIRECTLY as the arguments `s, r, d, g` of the PARMA functions. Do NOT call `getHPcpp`, `getrcpp` or `getdcpp`.
* Wrong argument count -> usage line to stderr, exit code 2.

## Mesh
Identical to `main-generator.cpp`: `nebin = 1000` energy bins, logarithmic, `nabin = 100` cosine bins, linear from -1 to +1.
`emin = 1.0e-2` MeV, `emax = 1.0e5` MeV for every particle except the neutron, for which `emin = 1.0e-8` MeV. Use the same `ehigh[]`, `emid[]`, `ahigh[]`, `amid[]` construction.
Convention: cosine +1 means particle coming from above (moving downward), exactly as in the original.

## What to write to `out_file` (text, UTF-8 without BOM, all floating numbers printed with `%.9e`)
Line 1: `# ip=<ip> W=<W> Rc=<Rc> d=<d> g=<g>`
Line 2: `nebin nabin ie511 flux511_cont flux511_line total_flux`  (six fields, single spaces) where
* `total_flux` = TotalFlux computed EXACTLY like the original (`etable[nebin]` before normalisation, including the `annihratio` factor in the 511-keV bin for photons), in /cm2/s;
* `ie511` = the index (1-based) of the energy bin that contains the electron mass `0.51099895` MeV for photons, else `0`;
* `flux511_cont` = `getSpecCpp(ip,s,r,d,emid[ie511],g)*(ehigh[ie511]-ehigh[ie511-1])` for photons, else `0`;
* `flux511_line` = `get511fluxCpp(s,r,d)` for photons, else `0`.
Line 3: the `nebin+1` values `ehigh[0..nebin]`, separated by single spaces.
Then `nebin` lines, line k (k = 1..nebin) has `nabin` values: for ia = 1..nabin
`D[k][ia] = getSpecCpp(ip,s,r,d,emid[k],g) * getSpecAngFinalCpp(IangPart[ip],s,r,d,emid[k],g,amid[ia]) * (2*pi) * (ahigh[ia]-ahigh[ia-1])`
(units /cm2/s/MeV per angular bin). Do NOT apply `annihratio` to `D`; the 511-keV line is carried separately by `flux511_line`.
Finally print to stderr one line: `total_flux_continuum=<sum over k of sum over ia of D[k][ia]*(ehigh[k]-ehigh[k-1])>` (continuum only, without the line) so that the caller can check it.

## Requirements
* Declare the PARMA functions exactly as the original does (they are defined in `subroutines.cpp`, which is linked separately). Do not define them.
* `static` arrays for the big tables; no dynamic allocation needed.
* Buffer output with one `FILE*` (`fopen`/`fprintf`), check that the file opened, return 0 on success.
* No random numbers, no `getGenerationCpp`, no particle loop, no `GeneOut` files.
* Comments in Russian, short.
* The code must compile with `cl /EHsc /O2 parma_tables.cpp subroutines.cpp` and with g++ -std=c++17.

## Original `main-generator.cpp` (reference, do not copy its particle generation)
