#include "pfunctional.h"

PFunctional::PFunctional()
{
    rfactor = 10.0;
}

void PFunctional::calculate(DoubleVector &x, double r)
{
    while (true)
    {
        func->setPenaltyCoefficient(r);

        ConjugateGradient g;
        g.setFunction(func);
        g.setGradient(func);
        g.setPrinter(func);
        g.setProjection(func);
        g.setOptimalityTolerance(0.0001);
        g.setStepTolerance(0.0001);
        g.setFunctionTolerance(0.0001);
        g.setR1MinimizeEpsilon(0.1, 0.01);
        //g.setNormalize(true);
        g.showExitMessage(true);
        //g.setResetIteration(false);

        g.calculate(x);

        r *= rfactor;
    }
}

void PFunctional::calculate(DoubleVector &x, const DoubleVector &r, const DoubleVector &e, const DoubleVector &s)
{
    grad->setFunction(func);
    grad->setGradient(func);
    grad->setPrinter(func);
    grad->setProjection(func);

    for (unsigned int i=0; i<r.length(); i++)
    {
        func->setPenaltyCoefficient(r[i]);
        func->setRegEpsilon(e[i]);
        grad->setR1MinimizeEpsilon(s[i], 0.00001);
        grad->calculate(x);
        IPrinter::printSeperatorLine();
    }
}
