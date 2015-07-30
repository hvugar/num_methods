#include "rosenbrock.h"
#include <gradient_cjt.h>
#include <gradient_sd.h>

double Rosenbrock::fx(const DoubleVector& x) const
{
    double x1 = x[0];
    double x2 = x[1];
    return ((1 - x1) * (1 - x1)) + 100 * (x2 - x1 * x1) * (x2 - x1 * x1);
}

void Rosenbrock::gradient(DoubleVector &g) const
{
    double h = grad_step;

//    for (unsigned i=0; i<g.size(); i++)
//    {
//        m_x[i] = m_x[i] - h;
//        double f1 = m_fn->fx(m_x);
//        m_x[i] = m_x[i] + 2*h;
//        double f2 = m_fn->fx(m_x);
//        m_x[i] = m_x[i] - h;

//        g[i] = (f2 - f1) / (2 * h);
//    }
}

void Rosenbrock::Main()
{
    /* Function */
    Rosenbrock r;

    /* initial point */
    DoubleVector x;
    x.push_back(-1.2);
    x.push_back(+1.0);

    /* Minimization */
    ConjugateGradient g;
    g.setFunction(&r);
    g.setEpsilon(0.000001);
    g.setGradientStep(0.000001);
    g.setR1MinimizeEpsilon(0.1, 0.000001);
    g.setX(x);
    g.calculate();
}

