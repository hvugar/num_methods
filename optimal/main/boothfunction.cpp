#include "boothfunction.h"
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <gradient_cs.h>

double BoothFunction::fx(const DoubleVector& x)
{
    double x1 = x[0];
    double x2 = x[1];
    return (x1 + 2.0*x2 - 7.0)*(x1 + 2.0*x2 - 7.0) + (2.0*x1 + x2 - 5.0)*(2.0*x1 + x2 - 5.0);
}

void BoothFunction::gradient(const DoubleVector& x, DoubleVector &g, double gradient_step)
{
    RnFunction::Gradient(this, gradient_step, x, g);
}

void BoothFunction::print(unsigned int iterationCount, const DoubleVector& m_x, const DoubleVector &s, double m_alpha, RnFunction* f) const
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
    y>=0.0 ? printf("|+%.6f\t", y) : printf("|%.6f\t", y);
    s[0]>=0.0 ? printf("|+%.6f\t", s[0]) : printf("|%.6f\t", s[0]);
    s[1]>=0.0 ? printf("|+%.6f\t", s[1]) : printf("|%.6f\t", s[1]);
    nr>=0.0 ? printf("|+%.6f\t", nr) : printf("|%.6f\t", nr);
    m_alpha>=0.0 ? printf("|%+.10f\t", m_alpha) : printf("|%.10f\t", m_alpha);
    printf("\n");
}

void BoothFunction::project(DoubleVector &x, int index)
{
}

void BoothFunction::main()
{
    /* Function */
    BoothFunction func;

    func.a = -10.0;
    func.b = +10.0;

    /* initial point */
    DoubleVector x0(2);
    x0[0] = -1.2;
    x0[1] = +1.0;

    /* Minimization */
    SteepestDescentGradient g1;
    g1.setFunction(&func);
    g1.setEpsilon1(0.000001);
    g1.setEpsilon2(0.000001);
    g1.setGradientStep(0.000001);
    g1.setR1MinimizeEpsilon(0.1, 0.000001);
    g1.setPrinter(&func);
//    g1.calculate(x0);

    puts("-----------------------------------------------------------------");
    x0[0] = +15.2;
    x0[1] = -14.0;
    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&func);
    g2.setEpsilon1(0.000001);
    g2.setEpsilon2(0.000001);
    g2.setGradientStep(0.000001);
    g2.setR1MinimizeEpsilon(0.1, 0.000001);
    g2.setPrinter(&func);
    g2.setProjection(&func);
    g2.setNormalize(false);
    g2.calculate(x0);

    puts("-----------------------------------------------------------------");
    x0[0] = -1.2;
    x0[1] = +1.0;
    /* Minimization */
    ConstStepGradient g3;
    g3.setFunction(&func);
    g3.setEpsilon1(0.000001);
    g3.setEpsilon2(0.000001);
    g3.setGradientStep(0.000001);
    g3.setR1MinimizeEpsilon(0.1, 0.000001);
    g3.setPrinter(&func);
    g3.setNormalize(false);
//    g3.calculate(x0);
}

