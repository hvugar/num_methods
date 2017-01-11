#include "borderparabolicd.h"
#include <math.h>

void BorderParabolicD::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    BorderParabolicD bp;

    // Real solution
    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;
        DoubleMatrix ru(bp.M+1,bp.N+1);
        IPrinter::printSeperatorLine("Real solution");
        for (unsigned int i=0; i<=bp.M; i++)
        {
            for (unsigned int j=0; j<=bp.N; j++)
                ru.at(i,j) = bp.U(j,i);
        }
        IPrinter::printMatrix(14, 10, ru, 10, 10, NULL);
        ru.clear();
    }

    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;

        DoubleMatrix u;
        IPrinter::printSeperatorLine("calculateU");
        bp.calculateU(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
        u.clear();
    }

    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;
        DoubleMatrix u;
        IPrinter::printSeperatorLine("calculateN4L2RM");
        bp.calculateN4L2RD(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
    }

    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;
        DoubleMatrix u;
        IPrinter::printSeperatorLine("calculateN4R2LM");
        bp.calculateN4R2LD(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
    }

    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;
        DoubleMatrix u;
        IPrinter::printSeperatorLine("calculateN6L2RM");
        bp.calculateN6L2RD(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
    }

    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;
        DoubleMatrix u;
        IPrinter::printSeperatorLine("calculateN6R2LM");
        bp.calculateN6R2LD(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
    }
}

double BorderParabolicD::initial(unsigned int i UNUSED_PARAM) const
{
    return U(i,0);
}

double BorderParabolicD::boundary(Boundary type UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    if (type == Left)  return U(0,j);
    if (type == Right) return U(N,j);
    return 0.0;
}

double BorderParabolicD::f(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    double t UNUSED_PARAM = j*ht;
    double x UNUSED_PARAM = i*hx;
#ifdef SAMPLE_1
    return 2.0*t - 2.0*a*a;
#endif
#ifdef SAMPLE_2
    return x - a*a*(40.0*cos(20.0*x) - 400.0*x*sin(20.0*x));
#endif
#ifdef SAMPLE_3
    return (sin(10.0*x) + 2.0*x*exp(2.0*x*t) - a*a*(4.0*t*t*exp(2.0*x*t) - 100.0*t*sin(10.0*x)));
#endif
#ifdef SAMPLE_4
    return 1.0 - 12.0*a*a*x*x;
#endif
#ifdef SAMPLE_5
    return 1.0 + 100.0*sin(10.0*x)*a*a;
#endif
#ifdef SAMPLE_6
    unsigned int m = 3;
    unsigned int n = 3;
    return pow(x,m) - a*a*(m*(m-1)*pow(x,m-2)*t + n*(n-1)*pow(x,n-2));
#endif
}

double BorderParabolicD::U(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    double t = j*ht;
#ifdef SAMPLE_1
    return x*x + t*t;
#endif
#ifdef SAMPLE_2
    return x*sin(20.0*x) + t*x;
#endif
#ifdef SAMPLE_3
    return t*sin(10.0*x) + exp(2.0*x*t);
#endif
#ifdef SAMPLE_4
    return x*x*x*x + t;
#endif
#ifdef SAMPLE_5
    return sin(10.0*x) + t;
#endif
#ifdef SAMPLE_6
    unsigned int m = 3;
    unsigned int n = 3;
    return pow(x,m)*t + pow(x,n);
#endif
}
