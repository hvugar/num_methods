#include "rosenbrock.h"
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <gradient_cs.h>

double Rosenbrock::fx(const DoubleVector& x)
{
    Rosenbrock* r = const_cast<Rosenbrock*>(this);
    r->count++;

    double x1 = x[0];
    double x2 = x[1];
    return ((1 - x1) * (1 - x1)) + 100 * (x2 - x1 * x1) * (x2 - x1 * x1);
}

void Rosenbrock::gradient(const DoubleVector& x, DoubleVector &g, double gradient_step)
{
    RnFunction::Gradient(this, gradient_step, x, g);
}

void Rosenbrock::gradient(const DoubleVector &x, DoubleVector &g)
{
    RnFunction::Gradient(this, grad_step, x, g);
}

void Rosenbrock::print(unsigned int iterationCount, const DoubleVector& m_x, const DoubleVector &s, double m_alpha, RnFunction* f) const
{
    if (iterationCount == 1)
    {
        printf("No\t|x1      \t|x2      \t|f(x)      \t|s1      \t|s2      \t|grad_norm  \t|alpha  \t");
        printf("\n--------+---------------+---------------+---------------+---------------+---------------+---------------+-------------\n");
    }

    double y = f->fx(m_x);
    double nr = s.EuclideanNorm();

    printf("%d\t", iterationCount);
    m_x[0]>=0.0 ? printf("|+%.10f\t", fabs(m_x[0])) : printf("|%.10f\t", m_x[0]);
    m_x[1]>=0.0 ? printf("|+%.10f\t", fabs(m_x[1])) : printf("|%.10f\t", m_x[1]);
    y>=0.0 ? printf("|%+.10f\t", y) : printf("|%.10f\t", y);
    s[0]>=0.0 ? printf("|+%.6f\t", s[0]) : printf("|%.6f\t", s[0]);
    s[1]>=0.0 ? printf("|+%.6f\t", s[1]) : printf("|%.6f\t", s[1]);
    nr>=0.0 ? printf("|+%.6f\t", nr) : printf("|%.6f\t", nr);
    m_alpha>=0.0 ? printf("|%+.10f\t", m_alpha) : printf("|%.10f\t", m_alpha);
    printf("\n");
}

void Rosenbrock::main()
{
    /* Function */
    Rosenbrock r;

    /* initial point */
    DoubleVector x0(2);
    x0[0] = -1.2;
    x0[1] = +1.0;

    /* Minimization */
//    SteepestDescentGradient g1;
//    g1.setFunction(&r);
//    g1.setEpsilon1(0.000001);
//    g1.setEpsilon2(0.000001);
//    g1.setGradientStep(0.000001);
//    g1.setR1MinimizeEpsilon(0.1, 0.000001);
//    g1.setPrinter(&r);
//    g1.calculate(x0);

    puts("-----------------------------------------------------------------");
    x0[0] = -1.2;
    x0[1] = +1.0;
    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&r);
    g2.setEpsilon1(0.000001);
    g2.setEpsilon2(0.000001);
    g2.setGradientStep(0.000001);
    g2.setR1MinimizeEpsilon(0.1, 0.000001);
    g2.setPrinter(&r);
    g2.setNormalize(true);
    g2.calculate(x0);

    printf("Function call count: %u\n", r.count);

//    puts("-----------------------------------------------------------------");
//    x0[0] = -1.2;
//    x0[1] = +1.0;
//    /* Minimization */
//    ConstStepGradient g3;
//    g3.setFunction(&r);
//    g3.setEpsilon1(0.000001);
//    g3.setEpsilon2(0.000001);
//    g3.setGradientStep(0.000001);
//    g3.setR1MinimizeEpsilon(0.1, 0.000001);
//    g3.setPrinter(&r);
//    g3.setNormalize(false);
//    g3.calculate(x0);
}

