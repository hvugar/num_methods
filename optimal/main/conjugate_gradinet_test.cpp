#include "conjugate_gradinet_test.h"
#include <float.h>

void ConjugateGradinetTest::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    ConjugateGradinetTest f;

    ConjugateGradient g;
    g.setFunction(&f);
    g.setGradient(&f);
    g.setPrinter(&f);
    g.setOptimalityTolerance(0.00000001);
    g.setFunctionTolerance(0.00000001);
    g.setStepTolerance(0.00000001);
    g.setR1MinimizeEpsilon(0.01, 10.0*DBL_EPSILON);
    g.setMaxIterationCount(200);
    g.setNormalize(false);
    g.setResetIteration(false);
    //g.setAlgorithm(ConjugateGradient::Algorithm::POLAK_RIBIERE);
    g.showExitMessage(true);

    DoubleVector x(20);
    for (unsigned int i=0; i<x.length(); i++) x[i] = 1.0;//(rand() / 1000) * 0.001;

    g.calculate(x);
}

double ConjugateGradinetTest::fx(const DoubleVector &x) const
{
    double ret = 0.0;
    //for (unsigned int i=0; i<x.length(); i++) ret += ( (0.2*(i+1)*x[i]-0.1*(i+5)) * (0.2*(i+1)*x[i]-0.1*(i+5)) );
    for (unsigned int i=0; i<x.length(); i++) ret += ( ((i+1)*x[i]-(i+5)) * ((i+1)*x[i]-(i+5)) );
    //for (unsigned int i=0; i<x.length(); i++) ret += x[i]*x[i];
    //for (unsigned int i=0; i<x.length(); i++) ret += (2.0*x[i]-5.0)*(2.0*x[i]-5.0);
    return ret;
}

void ConjugateGradinetTest::gradient(const DoubleVector &x, DoubleVector &g) const
{
    g.resize(x.length());
    //for (unsigned int i=0; i<x.length(); i++) g[i] = 2.0 * ( 0.2*(i+1) ) * ( 0.2*(i+1)*x[i]-0.1*(i+5) );
    for (unsigned int i=0; i<x.length(); i++) g[i] = 2.0 * (i+1) * ((i+1)*x[i]-(i+5));
    //for (unsigned int i=0; i<x.length(); i++) g[i] = 2.0*x[i];
    //for (unsigned int i=0; i<x.length(); i++) g[i] = 2.0*2.0*(2.0*x[i]-5.0);
}

void ConjugateGradinetTest::print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const
{
    //printf("I[%4d] %14.10f %10.6f: \n", iteration, f, alpha);
    //IPrinter::print(x, x.length());

    printf("%4d %20.14f %20.14f: \n", iteration, fx(x), alpha); //IPrinter::print(x, x.length());

}

