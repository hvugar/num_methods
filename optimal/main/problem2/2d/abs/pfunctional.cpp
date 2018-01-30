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
        g.setEpsilon2(0.0001);
        g.setEpsilon3(0.0001);
        g.setR1MinimizeEpsilon(0.1, 0.01);
        //g.setNormalize(true);
        g.showEndMessage(true);
        //g.setResetIteration(false);

        g.calculate(x);

        r *= rfactor;
    }
}

void PFunctional::calculate(DoubleVector &x, const DoubleVector &r, const DoubleVector &e)
{
    for (unsigned int i=0; i<r.length(); i++)
    {
        jfunc->setPenaltyCoefficient(r[i]);
        jfunc->setEpsilon(e[i]);

        ConjugateGradient g;
        //SteepestDescentGradient g;
        g.setFunction(jfunc);
        g.setGradient(jfunc);
        g.setPrinter(jfunc);
        g.setProjection(jfunc);
        g.setEpsilon1(0.0000001);
        g.setEpsilon2(0.0000001);
        g.setEpsilon3(0.0000001);
        g.setR1MinimizeEpsilon(2.0, 0.0001);
        g.setNormalize(true);
        g.showEndMessage(true);
        //g.setResetIteration(false);
        g.calculate(x);
    }
}
