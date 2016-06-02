#include "sampleborderhyperbolic.h"

void SampleBorderHyperBolic::main(int argc, char **argv)
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    SampleBorderHyperBolic bp;
    DoubleMatrix m;
    bp.calculateU(m, bp.hx, bp.ht, bp.M, bp.N);

    FILE *file = fopen("matrix.txt", "w");
    IPrinter::printMatrix(m, bp.M, bp.N, NULL, file);
    fclose(file);
}

SampleBorderHyperBolic::SampleBorderHyperBolic()
{
    x0 = t0 = 0.0;
    x1 = t1 = 1.0;
    a = 1.0;
    hx = 0.001;
    ht = 0.001;
    N = 1000;
    M = 1000;
}

double SampleBorderHyperBolic::initial1(unsigned int i) const
{
    double x = i*hx;
    return -4.0*(x-0.5)*(x-0.5)+1.0;
}

double SampleBorderHyperBolic::initial2(unsigned int i) const
{
    return 1.0;
}

double SampleBorderHyperBolic::boundary(Boundary type, unsigned int j) const
{
    return 0.0;
}

double SampleBorderHyperBolic::f(unsigned int i, unsigned int j) const
{
    return 0.0;
}
