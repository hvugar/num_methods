#include "borderparabolic.h"

void BorderParabolic::main(int argc, char **argv)
{
    C_UNUSED(argc);
    C_UNUSED(argv);
    BorderParabolic bp;
    DoubleMatrix u;
    bp.calculateU(u, bp.hx, bp.ht, bp.N, bp.M, 0.6);

    IPrinter::printMatrix(u, 10, 10, NULL, stdout);

    FILE* file = fopen("20160516.txt", "w");
    //    IPrinter::printVector(u, NULL, bp.N, 0, 0, file);
    //    IPrinter::printVector(u, NULL, 10, 0, 0, stdout);
    IPrinter::printMatrix(u, bp.M, bp.N, NULL, file);
    fclose(file);
}

BorderParabolic::BorderParabolic()
{
    x0 = 0.0;
    x1 = 1.0;
    t0 = 0.0;
    t1 = 3.0;
    hx = 0.001;
    ht = 0.001;
    N  = 1000;
    M  = 3000;
}

BorderParabolic::~BorderParabolic() {}

double BorderParabolic::initial(unsigned int i) const
{
    C_UNUSED(i);
    return 0.0;
}

double BorderParabolic::boundary(Boundary type, unsigned int j) const
{
    C_UNUSED(type);
    C_UNUSED(j);
    return 5.0;
}

double BorderParabolic::f(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}
