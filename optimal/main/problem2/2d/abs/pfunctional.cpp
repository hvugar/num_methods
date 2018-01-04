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
        g.setEpsilon1(0.000);
        g.setEpsilon2(0.000);
        g.setEpsilon3(0.000);
        g.setR1MinimizeEpsilon(1.0, 0.001);
        g.setNormalize(true);
        g.showEndMessage(true);
        //g.setResetIteration(false);

        g.calculate(x);

        r *= rfactor;
    }
}
