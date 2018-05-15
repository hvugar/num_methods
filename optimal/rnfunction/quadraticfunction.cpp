#include "quadraticfunction.h"

void QuadraticFunction::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    /* Function */
    QuadraticFunction qf;
    qf.grad_step = 0.000001;

    /* initial point */
    puts("-----------------------------------------------------------------");
    DoubleVector x0(2);
    x0[0] = +0.0;
    x0[1] = +0.0;

    /* Minimization */
    ConjugateGradient g2;
    g2.setGradient(&qf);
    g2.setFunction(&qf);
    g2.setEpsilon1(0.00001);
    g2.setEpsilon2(0.00001);
    g2.setEpsilon3(0.00001);
    g2.setR1MinimizeEpsilon(0.001, 0.0000001);
    g2.setPrinter(&qf);
    g2.setNormalize(true);
    g2.setAlgorithm(ConjugateGradient::FLETCHER_REEVES);
    g2.setResetIteration(true);
    g2.calculate(x0);

    printf("Function call count: %u\n", qf.count);
}

QuadraticFunction::QuadraticFunction()
{
    count = 0;
    grad_step = 0.000001;
}

double QuadraticFunction::fx(const DoubleVector& x) const
{
    QuadraticFunction* qf = const_cast<QuadraticFunction*>(this);
    qf->count++;

    double x1 = x[0];
    double x2 = x[1];
    return (x1-5.0)*(x1-5.0)*(x1-5.0)*(x1-5.0) + (x2+4.0)*(x2+4.0);
    //return 4.0*x1*x1 + 3.0*x2*x2 - 4.0*x1*x2 + x1;
}

void QuadraticFunction::gradient(const DoubleVector& x, DoubleVector &g) const
{
    IGradient::Gradient(this, grad_step, x, g);
}

void QuadraticFunction::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double, GradientMethod::MethodResult) const
{
    if (i == 0)
    {
        printf("No\t|x1      \t|x2      \t|f(x)      \t|s1      \t|s2      \t|grad_norm  \t|alpha  \t");
        printf("\n--------+---------------+---------------+---------------+---------------+---------------+---------------+-------------\n");
    }

    double y = const_cast<QuadraticFunction*>(this)->fx(x);
    double nr = g.EuclideanNorm();

    printf("%d\t", i);
    x[0]>=0.0 ? printf("|+%.10f\t", fabs(x[0])) : printf("|%.10f\t", x[0]);
    x[1]>=0.0 ? printf("|+%.10f\t", fabs(x[1])) : printf("|%.10f\t", x[1]);
    y>=0.0 ? printf("|%+.10f\t", y) : printf("|%.10f\t", y);
    g[0]>=0.0 ? printf("|+%.10f\t", g[0]) : printf("|%.10f\t", g[0]);
    g[1]>=0.0 ? printf("|+%.10f\t", g[1]) : printf("|%.10f\t", g[1]);
    nr>=0.0 ? printf("|+%.10f\t", nr) : printf("|%.10f\t", nr);
    printf("\n");
}
