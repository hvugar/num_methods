#include "rosenbrock.h"
#include <gradient_cjt.h>
#include <gradient_sd.h>

double Rosenbrock::fx(const DoubleVector& x)
{
    double x1 = x[0];
    double x2 = x[1];
    return ((1 - x1) * (1 - x1)) + 100 * (x2 - x1 * x1) * (x2 - x1 * x1);
}

void Rosenbrock::gradient(double gradient_step, const DoubleVector& x, DoubleVector &g)
{
    RnFunction::Gradient(this, gradient_step, x, g);
}

void RosenbrockPrinter::print(unsigned int iterationCount, const DoubleVector& m_x, const DoubleVector &s, double m_alpha, RnFunction* f) const
{
    if (iterationCount == 1)
    {
        printf("No\t|x1      \t|x2      \t|f(x)      \t|s1      \t|s2      \t|grad_norm  \t|alpha  \t");
        printf("\n--------+---------------+---------------+---------------+---------------+---------------+---------------+-------------\n");
    }

    double y = f->fx(m_x);
    double nr = m_x.EuclideanNorm();

    printf("%d\t", iterationCount);
    m_x[0]>=0.0 ? printf("|+%.10f\t", fabs(m_x[0])) : printf("|%.10f\t", m_x[0]);
    m_x[1]>=0.0 ? printf("|+%.10f\t", fabs(m_x[1])) : printf("|%.10f\t", m_x[1]);
    y>=0.0 ? printf("|%+.10f\t", y) : printf("|%.10f\t", y);
    s[0]>=0.0 ? printf("|%+.10f\t", s[0]) : printf("|%0.10f\t", s[0]);
    s[1]>=0.0 ? printf("|%+.10f\t", s[1]) : printf("|%0.10f\t", s[1]);
    nr>=0.0 ? printf("|%+.10f\t", nr) : printf("|%.10f\t", nr);
    m_alpha>=0.0 ? printf("|%+.10f\t", m_alpha) : printf("|%.10f\t", m_alpha);
    printf("\n");
}

void Rosenbrock::Main()
{
    /* Function */
    Rosenbrock r;

    /* initial point */
    DoubleVector x0;
    x0.push_back(-1.2);
    x0.push_back(+1.0);

    /* Minimization */
    ConjugateGradient g;
    g.setFunction(&r);
    g.setEpsilon(0.000001);
    g.setGradientStep(0.000001);
    g.setR1MinimizeEpsilon(0.1, 0.000001);
    g.setX(x0);
    g.setPrinter(new RosenbrockPrinter);
    g.calculate();
}

