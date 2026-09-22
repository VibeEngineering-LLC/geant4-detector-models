#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;

// PARMA функции (определены в subroutines.cpp)
double getSpecCpp(int ip, double s, double r, double d, double e, double g);
double getSpecAngFinalCpp(int iang, double s, double r, double d, double e, double g, double a);
double get511fluxCpp(double s, double r, double d);

// Таблицы
const int nebin = 1000;
const int nabin = 100;
static double ehigh[nebin+1], emid[nebin+1];
static double ahigh[nabin+1], amid[nabin+1];
static double D[nebin+1][nabin+1];

// Массивы для частиц
const int npart = 33;
static int IangPart[npart+1] = {1,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,5,5,6};

int main(int argc, char* argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s <ip> <W> <Rc_GV> <depth_gcm2> <g> <out_file>\n", argv[0]);
        return 2;
    }

    int ip = atoi(argv[1]);
    double s = atof(argv[2]);
    double r = atof(argv[3]);
    double d = atof(argv[4]);
    double g = atof(argv[5]);
    const char* out_file = argv[6];

    if (ip < 0 || ip > npart || IangPart[ip] == 0) {
        fprintf(stderr, "Error: Invalid particle ID or no angular distribution available\n");
        return 1;
    }

    // Установка энергетической сетки
    double emin = (ip == 0) ? 1.0e-8 : 1.0e-2;
    double emax = 1.0e5;
    double elog = log10(emin);
    double estep = (log10(emax) - log10(emin)) / nebin;

    for (int ie = 0; ie <= nebin; ie++) {
        ehigh[ie] = pow(10, elog);
        if (ie != 0) emid[ie] = sqrt(ehigh[ie] * ehigh[ie-1]);
        elog += estep;
    }

    // Установка угловой сетки
    double amin = -1.0;
    double amax = 1.0;
    double astep = (amax - amin) / nabin;
    for (int ia = 0; ia <= nabin; ia++) {
        ahigh[ia] = amin + astep * ia;
        if (ia != 0) amid[ia] = (ahigh[ia] + ahigh[ia-1]) * 0.5;
    }

    // Вычисление таблицы D
    double total_flux_continuum = 0.0;
    int ie511 = 0;
    double flux511_cont = 0.0;
    double flux511_line = 0.0;

    if (ip == npart) { // Фотон
        double emass = 0.51099895;
        for (int ie = 1; ie <= nebin; ie++) {
            if (ehigh[ie-1] < emass && emass <= ehigh[ie]) {
                ie511 = ie;
                flux511_cont = getSpecCpp(ip, s, r, d, emid[ie], g) * (ehigh[ie] - ehigh[ie-1]);
                flux511_line = get511fluxCpp(s, r, d);
                break;
            }
        }
    }

    double total_flux = 0.0;

    for (int ie = 1; ie <= nebin; ie++) {
        double energy_bin_width = ehigh[ie] - ehigh[ie-1];
        double spec_val = getSpecCpp(ip, s, r, d, emid[ie], g);
        double sum_angular = 0.0;

        for (int ia = 1; ia <= nabin; ia++) {
            double angular_bin_width = ahigh[ia] - ahigh[ia-1];
            double ang_val = getSpecAngFinalCpp(IangPart[ip], s, r, d, emid[ie], g, amid[ia]);
            D[ie][ia] = spec_val * ang_val * 2.0 * M_PI * angular_bin_width;
            sum_angular += D[ie][ia];
        }

        if (ip == npart && ie == ie511) {
            total_flux += sum_angular * energy_bin_width * (1.0 + flux511_line / flux511_cont);
        } else {
            total_flux += sum_angular * energy_bin_width;
        }
        total_flux_continuum += sum_angular * energy_bin_width;
    }

    // PARMA читает input/... от текущего каталога и при неудаче молча даёт NaN:
    // нечисловой или нулевой поток = отказ, а не файл с мусором.
    if (!std::isfinite(total_flux) || !std::isfinite(total_flux_continuum) || total_flux <= 0.0) {
        fprintf(stderr, "Error: total_flux=%g is not a positive finite number "
                        "(run from the directory that contains input/)\n", total_flux);
        return 3;
    }

    FILE* fp = fopen(out_file, "w");
    if (!fp) {
        perror("Error opening output file");
        return 1;
    }

    fprintf(fp, "# ip=%d W=%.9e Rc=%.9e d=%.9e g=%.9e\n", ip, s, r, d, g);
    fprintf(fp, "%d %d %d %.9e %.9e %.9e\n", nebin, nabin, ie511, flux511_cont, flux511_line, total_flux);

    // Запись энергетических границ
    for (int ie = 0; ie <= nebin; ie++) {
        fprintf(fp, "%.9e%c", ehigh[ie], (ie < nebin) ? ' ' : '\n');
    }

    // Запись значений D
    for (int ie = 1; ie <= nebin; ie++) {
        for (int ia = 1; ia <= nabin; ia++) {
            fprintf(fp, "%.9e%c", D[ie][ia], (ia < nabin) ? ' ' : '\n');
        }
    }

    fclose(fp);

    fprintf(stderr, "total_flux_continuum=%.9e\n", total_flux_continuum);

    return 0;
}
