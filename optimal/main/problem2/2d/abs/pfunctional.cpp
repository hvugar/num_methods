#include "pfunctional.h"

PFunctional::PFunctional()
{
    rfactor = 10.0;
}

void PFunctional::calculate(DoubleVector &x, double r)
{
    while (true)
    {
        jfunc->setPenaltyCoefficient(r);

        ConjugateGradient g;
        g.setFunction(jfunc);
        g.setGradient(jfunc);
        g.setPrinter(jfunc);
        g.setProjection(jfunc);
        g.setEpsilon1(0.0001);
        g.setEpsilon2(0.0000001);
        g.setEpsilon3(0.0000001);
        g.setR1MinimizeEpsilon(1.0, 0.001);
        //g.setNormalize(true);
        g.showEndMessage(true);
        //g.setResetIteration(false);

        g.calculate(x);

        r *= rfactor;
    }
}

void PFunctional::calculate(DoubleVector &x)
{
    double r[] = {0.01, 0.10, 1.0, 10.0, 20.0, 100.0, 200.0, 500.0};
    //double e[] = {1.00, 1.00, 1.0, 1.00, 0.10, 0.100, 0.010, 0.010};
    double e[] = {0.00, 0.00, 0.0, 0.00, 0.00, 0.000, 0.000, 0.000};
    for (unsigned int i=0; i<6; i++)
    {
        jfunc->setPenaltyCoefficient(r[i]);
        jfunc->setEpsilon(e[i]);

        ConjugateGradient g;
        g.setFunction(jfunc);
        g.setGradient(jfunc);
        g.setPrinter(jfunc);
        g.setProjection(jfunc);
        g.setEpsilon1(0.0000001);
        g.setEpsilon2(0.0000001);
        g.setEpsilon3(0.0000001);
        g.setR1MinimizeEpsilon(1.0, 0.001);
        g.setNormalize(true);
        g.showEndMessage(true);
        //g.setResetIteration(false);

        g.calculate(x);
    }
}
