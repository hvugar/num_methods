#include <stdio.h>
#include <function.h>
#include <r1minimize.h>
#include <gradient.h>
#include <fpgradient.h>
#include <cjtgradient.h>
#include <methods.h>
#include "cfunction.h"
#include "samplecontrol.h"

struct Rosenbrock : public RnFunction
{
    virtual double fx(const std::vector<double>& x);
};

double Rosenbrock::fx(const std::vector<double>& x)
{
    double x1 = x[0];
    double x2 = x[1];
    return ((1 - x1) * (1 - x1)) + 100 * (x2 - x1 * x1) * (x2 - x1 * x1);
}

#define first
int main()
{
#ifdef first
    /* Function */
    Rosenbrock r;

    /* initial point */
    std::vector<double> x;
    x.push_back(-1.0);
    x.push_back(+1.2);

    /* Minimization */
    FastProximalGradient fg;
    fg.setF(&r);
    fg.setEpsilon(0.000001);
    fg.setGradientStep(0.000001);
    fg.setR1MinimizeEpsilon(0.1, 0.000001);
    fg.setX(x);
    fg.calculate();

#else
    CFunction* f = new CFunction(0.0, 1.0, 0.001);
    std::vector<double> u(f->N);
    for (int i=0; i<f->N; i++) u[i] = 0.01;

    SampleControl sc;
    sc.setF(f);
    sc.setX(u);
    sc.setEpsilon(0.0000001);
    sc.setGradientStep(0.000001);
    sc.setR1MinimizeEpsilon(0.01, 0.000001);
    sc.calculate();
#endif
    return 0;
}

