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
    x1 = 1.0;
    t1 = 2.0;
    a = 1.0;
    hx = 0.001;
    ht = 0.001;
    N = 1000;
    M = 2000;
}

double SampleBorderHyperBolic::initial1(unsigned int i UNUSED_PARAM) const
{
    double x = i*hx;
    //return -4.0*(x-0.5)*(x-0.5)+1.0;
    return sin(2.0*M_PI*x);
}

double SampleBorderHyperBolic::initial2(unsigned int i UNUSED_PARAM) const
{
    //double x = i*hx;
    return -1.0;
}

double SampleBorderHyperBolic::boundary(Boundary type, unsigned int j) const
{
    C_UNUSED(type);
    C_UNUSED(j);
    return 0.0;
}

double SampleBorderHyperBolic::f(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}
