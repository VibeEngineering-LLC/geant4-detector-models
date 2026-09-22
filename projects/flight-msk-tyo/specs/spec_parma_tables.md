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
#include <iostream>
#include <cmath>
#include <random>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

double getHPcpp(int,int,int);
double getrcpp(double,double);
double getdcpp(double,double);
double getSpecCpp(int,double,double,double,double,double);
double getSpecAngFinalCpp(int,double,double,double,double,double,double);
double getGenerationCpp(double,double,int*,int,double*,double*);
double get511fluxCpp(double, double, double);

int main()
{
const int npart = 33;
static int IangPart[npart+1] = {1,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,5,5,6};
const int nebin=1000; // number of energy mesh (divided by log)
const int nabin=100;  // number of angle mesh (divided by linear)
static double ehigh[nebin+1],emid[nebin+1]; // higher and middle point of energy bin
static double ahigh[nabin+1],amid[nabin+1]; // higher and middle point of angular bin
static double etable[nebin+1] = {}; // probability table (0.0 for 0, 1.0 for nebin)
static double atable[nabin+1][nebin+1] = {}; // probability table (0.0 for 0, 1.0 for nabin) 
double atable2[nabin+1]; // temporary used dimension for anguluar probability
static double Flux[nabin+1][nebin+1] = {}; // Monte Carlo generated flux

int ia,ie,i,ia2,ie511;
double e,phi,u,v,w,sx,cx,xd,yd,zd,x,y,z,annihratio,flux511bin;

double emass = 0.51099895e0; // mass of electron in MeV

 // Set condition 
int nevent=1000; // number of particles to be generated
int ip = 1; // Particle ID (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)
int iyear = 2019; // Year
int imonth = 2;   // Month
int iday = 1;    // Day
double glat = 30.5; // Latitude (deg), -90 =< glat =< 90
double glong = -76.2; // Longitude (deg), -180 =< glong =< 180
double alti = 0.0; // Altitude (km)
double g = 0.15; // Local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin
double radi = 100.0; // radius of the target area in cm (put your target inside this area)

// Set energy and angle ranges for generation
double emin = 1.0e0;  // Minimum energy of particle
double emax = 1.0e5;  // Maximum energy of particle
double amin = -1.0;   // Minimum cosine of particle
double amax =  1.0;   // Maximum cosine of particle

// calculate parameters
double s = getHPcpp(iyear,imonth,iday); // solar modulation potential
double r = getrcpp(glat,glong);  // Vertical cut-off rigidity (GV)
double d = getdcpp(alti,glat);   // Atmospheric depth (g/cm2), set glat = 100 for use US Standard Atmosphere 1976.

if(IangPart[ip] == 0) {
 cout << "Angular distribution is not available for the particle";
 exit(1);
}

if(ip==0 && emin<1.0e-8) {emin=1.0e-8;} // Minimum energy for neutron is 10 meV
if(ip!=0 && emin<1.0e-2) {emin=1.0e-2;} // Minimum energy for other particle is 10 keV

// Make energy and angle mesh
double elog=log10(emin);
double estep=(log10(emax)-log10(emin))/nebin;
for(ie=0;ie<=nebin;ie++) {
 ehigh[ie]=pow(10,elog);
 if(ie!=0) emid[ie]=sqrt(ehigh[ie]*ehigh[ie-1]);
 elog=elog+estep;
}

// setup for 511 keV annihilation gamma
if (ip == npart) {
	for (ie = 1; ie <= nebin; ie++) {
		if (ehigh[ie - 1] < emass && emass <= ehigh[ie]) {
			ie511 = ie;
			flux511bin = getSpecCpp(ip, s, r, d, emid[ie], g) * (ehigh[ie] - ehigh[ie - 1]);
			annihratio = (flux511bin + get511fluxCpp(s, r, d)) / flux511bin;
			break;
		}
	}
} else {
	ie511 = 0;
}

double astep=(amax-amin)/nabin;
for(ia=0;ia<=nabin;ia++) {
 ahigh[ia]=amin+astep*ia;
 if(ia!=0) amid[ia]=(ahigh[ia]+ahigh[ia-1])*0.5;
}

// Make probability table (absolute value)
for(ie=1;ie<=nebin;ie++) {
 for(ia=1;ia<=nabin;ia++) {
  atable[ia][ie]=atable[ia-1][ie]+getSpecCpp(ip,s,r,d,emid[ie],g)*getSpecAngFinalCpp(IangPart[ip],s,r,d,emid[ie],g,amid[ia])*(2.0*acos(-1.0))*(ahigh[ia]-ahigh[ia-1]); // angular integrated value
 }
}
for (ie = 1; ie <= nebin; ie++) {
	if (ip == npart && ie == ie511) {
		etable[ie] = etable[ie - 1] + atable[nabin][ie] * (ehigh[ie] - ehigh[ie - 1]) * annihratio; // energy integrated value
	}
	else {
		etable[ie] = etable[ie - 1] + atable[nabin][ie] * (ehigh[ie] - ehigh[ie - 1]); // energy integrated value
	}
}
double TotalFlux=etable[nebin]; // Total Flux (/cm2/s), used for normalization

// Make probability table (normalized to 1)
for(ie=1;ie<=nebin;ie++) {
 etable[ie]=etable[ie]/etable[nebin];
 for(ia=1;ia<=nabin;ia++) {
  atable[ia][ie]=atable[ia][ie]/atable[nabin][ie];
 }
}

// Particle Generation
mt19937 engine; // Mersenne Twister
uniform_real_distribution<double> rand01(0.0,1.0); // random number between 0 to 1

ofstream sf("GeneOut/generation.out",ios::out);
sf << "ip= " << ip << " ,W-index= " << s << " ,Rc(GV)= " << r << " ,depth(g/cm2)= " << d << " ,g= " << g << "\n";
sf << "Total Flux (/cm2/s)= " << TotalFlux << "\n";
sf << "  Energy(MeV/n)              u              v              w              x              y              z\n"; 
for(i=1;i<=nevent;i++) {
 e=getGenerationCpp(rand01(engine),rand01(engine),&ie,nebin,ehigh,etable); // energy
 if (ip == npart && ie == ie511) { // may be annihilation gamma
	 if (rand01(engine) > 1.0 / annihratio) { e = emass; }
 }
 phi=2.0*acos(-1.0)*(rand01(engine)-0.5); // azimuth angle (rad)
 for(ia2=0;ia2<=nabin;ia2++) {atable2[ia2]=atable[ia2][ie];} 
 cx=getGenerationCpp(rand01(engine),rand01(engine),&ia,nabin,ahigh,atable2); // z direction, -1.0:upward, 0.0:horizontal, 1.0:downward
 do {
  xd = (rand01(engine)-0.5)*2.0*radi;
  yd = (rand01(engine)-0.5)*2.0*radi;
 } while (sqrt(xd*xd+yd*yd) > radi);
 zd = radi;

 sx=sqrt(1-cx*cx); // sin(theta)

 x=xd*cx*cos(phi)-yd*sin(phi)+zd*sx*cos(phi);
 y=xd*cx*sin(phi)+yd*cos(phi)+zd*sx*sin(phi);
 z=-xd*sx+zd*cx;
 u = -sx*cos(phi);
 v = -sx*sin(phi);
 w = -cx;
 sf << scientific << setw(15) << e << setw(15) << u << setw(15) << v << setw(15) << w << setw(15) << x << setw(15) << y << setw(15) << z <<"\n";
 Flux[ia][ie]=Flux[ia][ie]+TotalFlux/(ehigh[ie]-ehigh[ie-1])/((ahigh[ia]-ahigh[ia-1])*2.0*acos(-1.0)); // /cm2/s/sr/MeV
 Flux[0][ie]=Flux[0][ie]+TotalFlux/(ehigh[ie]-ehigh[ie-1]); // /cm2/s/MeV
}

ofstream of("GeneOut/flux.out",ios::out);
of << "ip= " << ip << " ,W-index= " << s << " ,Rc(GV)= " << r << " ,depth(g/cm2)= " << d << " ,g= " << g << "\n";
of << "Total Flux (/cm2/s)= " << TotalFlux << "\n";
of << "Angular and energy differential fluxes in /cm2/s/(MeV/n)/sr\n"; 
of << "   E_low(MeV/n) /cm2/s/(MeV/n)";
for(ia=1;ia<=nabin;ia++) {of << setw(15) << ahigh[ia-1];}
of << "\n";
for(ie=1;ie<=nebin;ie++) {
 of << scientific << setw(15) << ehigh[ie-1];
 for(ia=0;ia<=nabin;ia++) {of << setw(15) << Flux[ia][ie]/nevent;}
 of << "\n";
} 

} // End of main

// ***********************************************************
double getGenerationCpp(double rand1, double rand2, int* ibin, int nbin, double* high, double* table)
// ***********************************************************
{
double getGeneration = 0.0;

int i;

for(i=1;i<=nbin-1;i++) {
 if(rand1<=table[i]) {break;}
}
*ibin=i; // bin ID

getGeneration=high[*ibin-1]*rand2+high[*ibin]*(1.0-rand2);

return getGeneration;
}
