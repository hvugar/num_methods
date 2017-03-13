#include "boothfunction.h"
#include <math.h>

void BoothFunction::main(int argc, char ** argv)
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    /* Function */
    BoothFunction func;
    func.grad_step = 0.000001;

    func.a = -10.0;
    func.b = +10.0;

    /* initial point */
    DoubleVector x0(2);
    x0[0] = -1.2;
    x0[1] = +1.0;

    /* Minimization */
    SteepestDescentGradient g1;
    g1.setGradient(&func);
    g1.setFunction(&func);
    g1.setEpsilon1(0.000001);
    g1.setEpsilon2(0.000001);
    g1.setR1MinimizeEpsilon(0.1, 0.000001);
    g1.setPrinter(&func);
//    g1.calculate(x0);

    puts("-----------------------------------------------------------------");
    x0[0] = +15.2;
    x0[1] = -14.0;
    /* Minimization */
    ConjugateGradient g2;
    g2.setGradient(&func);
    g2.setFunction(&func);
    g2.setEpsilon1(0.000001);
    g2.setEpsilon2(0.000001);
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
    g3.setR1MinimizeEpsilon(0.1, 0.000001);
    g3.setPrinter(&func);
    g3.setNormalize(false);
//    g3.calculate(x0);
}

double BoothFunction::fx(const DoubleVector& x) const
{
    double x1 = x[0];
    double x2 = x[1];
    return (x1 + 2.0*x2 - 7.0)*(x1 + 2.0*x2 - 7.0) + (2.0*x1 + x2 - 5.0)*(2.0*x1 + x2 - 5.0);
}

void BoothFunction::gradient(const DoubleVector& x, DoubleVector &g)
{
    IGradient::Gradient(this, grad_step, x, g);
}

void BoothFunction::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double, GradientMethod::MethodResult) const
{
    if (i == 0)
    {
        printf("No\t|x1      \t|x2      \t|f(x)      \t|s1      \t|s2      \t|grad_norm  \t|alpha  \t");
        printf("\n--------+---------------+---------------+---------------+---------------+---------------+---------------+-------------\n");
    }

    BoothFunction *pm = const_cast<BoothFunction*>(this);
    double y = pm->fx(x);
    double nr = g.EuclideanNorm();

    printf("%d\t", i);
    x[0]>=0.0 ? printf("|+%.10f\t", fabs(x[0])) : printf("|%.10f\t", x[0]);
    x[1]>=0.0 ? printf("|+%.10f\t", fabs(x[1])) : printf("|%.10f\t", x[1]);
    y>=0.0 ? printf("|+%.6f\t", y) : printf("|%.6f\t", y);
    g[0]>=0.0 ? printf("|+%.6f\t", g[0]) : printf("|%.6f\t", g[0]);
    g[1]>=0.0 ? printf("|+%.6f\t", g[1]) : printf("|%.6f\t", g[1]);
    nr>=0.0 ? printf("|+%.6f\t", nr) : printf("|%.6f\t", nr);
    printf("\n");
}

void BoothFunction::project(DoubleVector &x, int index)
{
    C_UNUSED(x);
    C_UNUSED(index);
}
