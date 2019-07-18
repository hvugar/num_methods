#include "conjugategradinettest.h"

void ConjugateGradinetTest::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    ConjugateGradinetTest f;

    ConjugateGradient g;
    g.setFunction(&f);
    g.setGradient(&f);
    g.setPrinter(&f);
    g.setOptimalityTolerance(0.00001);
    g.setFunctionTolerance(0.00001);
    g.setStepTolerance(0.00001);
    g.setR1MinimizeEpsilon(0.1, 0.00000001);
    g.setMaxIterationCount(200);
    g.setNormalize(false);
    g.showExitMessage(true);

    DoubleVector x(10);
    for (unsigned int i=0; i<x.length(); i++) x[i] = (rand() / 1000) * 0.001;

    g.calculate(x);
}

double ConjugateGradinetTest::fx(const DoubleVector &x) const
{
    double ret = 0.0;
    for (unsigned int i=0; i<x.length(); i++) ret += ((i+1)*0.2*x[i]-(i+1)*0.1)*((i+1)*0.2*x[i]-(i+1)*0.1);
    return ret;
}

void ConjugateGradinetTest::gradient(const DoubleVector &x, DoubleVector &g) const
{
    g.resize(x.length());
    for (unsigned int i=0; i<x.length(); i++)
        g[i] = 2.0*(i+1)*0.2*((i+1)*0.2*x[i]-(i+1)*0.1);
}

void ConjugateGradinetTest::print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const
{
    printf("I[%4d] %14.10f %10.6f: ", iteration, f, alpha);
    IPrinter::print(x, x.length());
}

