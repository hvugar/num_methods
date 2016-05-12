#include "borderparabolic.h"

void BorderParabolic::main(int argc, char **argv)
{
    BorderParabolic bp;
    DoubleVector u;
    bp.calculateU(u, bp.hx, bp.ht, bp.N, bp.M);

    FILE* file = fopen("20160512.txt", "w");
    IPrinter::printVector(u, NULL, bp.N, 0, 0, file);
    IPrinter::printVector(u, NULL, 10, 0, 0, stdout);
    fclose(file);
}

BorderParabolic::BorderParabolic()
{
    x0 = 0.0;
    x1 = 1.0;
    t0 = 0.0;
    t1 = 1.0;
    hx = 0.001;
    ht = 0.001;
    N  = 1000;
    M  = 1000;
}

BorderParabolic::~BorderParabolic() {}

double BorderParabolic::initial(unsigned int i) const
{
    return 1.0;
}

double BorderParabolic::boundary(Boundary type, unsigned int j) const
{
    return 0.0;
}

double BorderParabolic::f(unsigned int i, unsigned int j) const
{
    return 0.0;
}
