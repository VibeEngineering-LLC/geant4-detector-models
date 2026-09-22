// Rc (ГВ) и атмосферная глубина (г/см2) по PARMA для списка точек: rc_probe <lat> <lon> <alt_km> [<lat> <lon> <alt_km> ...]
// Запуск строго из каталога parma/ (читает input/CORdata.inp).
#include <cstdio>
#include <cstdlib>
double getrcpp(double, double);
double getdcpp(double, double);
int main(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "usage: rc_probe lat lon alt_km ...\n"); return 2; }
    for (int i = 1; i + 2 < argc; i += 3) {
        double la = std::atof(argv[i]), lo = std::atof(argv[i + 1]), h = std::atof(argv[i + 2]);
        std::printf("lat=%.3f lon=%.3f alt_km=%.2f Rc_GV=%.4f depth_gcm2=%.4f\n", la, lo, h, getrcpp(la, lo), getdcpp(h, la));
    }
    return 0;
}
