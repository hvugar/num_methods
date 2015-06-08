#include <stdio.h>
#include <function.h>
#include <r1minimize.h>
#include <gradient.h>
#include <fpgradient.h>
#include <cjtgradient.h>
#include <prjgradient.h>
#include <methods.h>

#include "cfunction1.h"
#include "cfunction2.h"
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

#define first1
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
    ProjectionGradient fg;
    fg.a = 0.4;
    fg.b = 0.7;
    fg.setF(&r);
    fg.setEpsilon(0.000001);
    fg.setGradientStep(0.000001);
    fg.setR1MinimizeEpsilon(0.1, 0.000001);
    fg.setX(x);
    fg.calculate();

#else
    CFunction2* f = new CFunction2(0.0, 1.0, 0.01);
    std::vector<double> u(f->n);
    for (int i=0; i<f->n; i++) u[i] = 0.00001;//3*f->t[i];//*f->t[i];

    SampleControl sc;
    sc.setF(f);
    sc.setX(u);
    sc.setEpsilon(0.01);
    sc.setGradientStep(0.000001);
    sc.setR1MinimizeEpsilon(0.01, 0.000001);
    sc.calculate();
#endif
    return 0;
}

