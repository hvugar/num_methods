#include "rosenbrock.h"
#include <math.h>

void Rosenbrock::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    /* Function */
    Rosenbrock r;
    r.grad_step = 0.000001;

    /* initial point */
    DoubleVector x0(2);
    x0[0] = -1.2;
    x0[1] = +1.0;

    puts("-----------------------------------------------------------------");
    x0[0] = -1.2;
    x0[1] = +1.0;
    /* Minimization */
    ConjugateGradient g2;
    g2.setGradient(&r);
    g2.setFunction(&r);
    g2.setEpsilon1(0.000001);
    g2.setEpsilon2(0.000001);
    g2.setEpsilon3(0.000001);
    g2.setR1MinimizeEpsilon(0.1, 0.000001);
    g2.setPrinter(&r);
    g2.setNormalize(true);
    g2.calculate(x0);

    printf("Function call count: %u\n", r.count);

//    puts("-----------------------------------------------------------------");
//    x0[0] = -1.2;
//    x0[1] = +1.0;
//    /* Minimization */
//    SteepestDescentGradient g3;
//    g3.setGradient(&r);
//    g3.setFunction(&r);
//    g3.setEpsilon1(0.000001);
//    g3.setEpsilon2(0.000001);
//    g3.setEpsilon3(0.000001);
//    g3.setR1MinimizeEpsilon(0.1, 0.000001);
//    g3.setPrinter(&r);
//    g3.setNormalize(true);
//    g3.calculate(x0);
}

Rosenbrock::Rosenbrock()
{
    count = 0;
    grad_step = 0.000001;
}

double Rosenbrock::fx(const DoubleVector& x) const
{
    Rosenbrock* r = const_cast<Rosenbrock*>(this);
    r->count++;

    double x1 = x[0];
    double x2 = x[1];
    return ((x1-1.0) * (x1-1.0)) + 100 * (x2 - x1 * x1) * (x2 - x1 * x1);
}

void Rosenbrock::gradient(const DoubleVector& x, DoubleVector &g)
{
    IGradient::Gradient(this, grad_step, x, g);
}

void Rosenbrock::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double, GradientMethod::MethodResult) const
{
    if (i == 0)
    {
        printf("No\t|x1      \t|x2      \t|f(x)      \t|s1      \t|s2      \t|grad_norm  \t|alpha  \t");
        printf("\n--------+---------------+---------------+---------------+---------------+---------------+---------------+-------------\n");
    }

    double y = const_cast<Rosenbrock*>(this)->fx(x);
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



