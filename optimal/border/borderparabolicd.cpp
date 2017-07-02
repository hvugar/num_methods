#include "borderparabolicd.h"
#include <math.h>

void BorderParabolicD::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    BorderParabolicD bp;
    bp.a = 1.0;

    // Real solution
    {
        bp.hx = 0.1;
        bp.ht = 0.1;
        bp.N = 10;
        bp.M = 10;
        DoubleMatrix ru(bp.M+1,bp.N+1);
        IPrinter::printSeperatorLine("Real solution");
        for (unsigned int i=0; i<=bp.M; i++)
        {
            for (unsigned int j=0; j<=bp.N; j++)
                ru.at(i,j) = bp.U(j,i);
        }
        //IPrinter::printMatrix(14, 10, ru, 10, 10, NULL);
        IPrinter::printVector(14,10,ru.row(1));
        ru.clear();
    }

    bp.hx = 0.01;
    bp.ht = 0.0001;
    bp.N = 100;
    bp.M = 10000;

    {
        //bp.hx = 0.01;
        //bp.ht = 0.01;
        //bp.N = 100;
        //bp.M = 100;
        DoubleMatrix u;
        IPrinter::printSeperatorLine("calculateU");
        bp.calculateU(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        IPrinter::printVector(14,10,u.row(1));
        //IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
        u.clear();
    }

    {
        //bp.hx = 0.01;
        //bp.ht = 0.01;
        //bp.N = 100;
        //bp.M = 100;
        DoubleMatrix u;
        IPrinter::printSeperatorLine("calculateN2L2RD");
        bp.calculateN2L2RD(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        //IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
        IPrinter::printVector(14,10,u.row(1));
    }

    {
        //bp.hx = 0.01;
        //bp.ht = 0.01;
        //bp.N = 100;
        //bp.M = 100;
        DoubleMatrix u;
        IPrinter::printSeperatorLine("calculateN4L2RD");
        bp.calculateN4L2RD(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        //IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
        IPrinter::printVector(14,10,u.row(1));
    }

    return;

    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;
        DoubleMatrix u;
        IPrinter::printSeperatorLine("calculateN4L2RD");
        bp.calculateN4L2RD(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
    }

    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;
        DoubleMatrix u;
        IPrinter::printSeperatorLine("calculateN4R2LD");
        bp.calculateN4R2LD(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
    }

    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;
        DoubleMatrix u;
        IPrinter::printSeperatorLine("calculateN6L2RD");
        bp.calculateN6L2RD(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
    }

    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;
        DoubleMatrix u;
        IPrinter::printSeperatorLine("calculateN6R2LD");
        bp.calculateN6R2LD(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
    }
}

double BorderParabolicD::initial(unsigned int i UNUSED_PARAM) const
{
#ifndef SAMPLE_0
    return U(i,0);
#endif
#ifndef SAMPLE_8
    return U(i,0);
#endif
#ifdef SAMPLE_8
    return 1.0;
#endif
}

double BorderParabolicD::boundary(Boundary type UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
#ifndef SAMPLE_0
    if (type == Left)  return U(0,j);
    if (type == Right) return U(N,j);
#endif
#ifndef SAMPLE_8
    if (type == Left)  return U(0,j);
    if (type == Right) return U(N,j);
#endif
#ifdef SAMPLE_8
    if (type == Left)  return 4.0;
    if (type == Right) return 3.0;
#endif
    return 0.0;
}

double BorderParabolicD::f(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    double t UNUSED_PARAM = j*ht;
    double x UNUSED_PARAM = i*hx;
#ifdef SAMPLE_0
    return x*x - 2.0*a*a*t;
#endif
#ifdef SAMPLE_1
    return 2.0*t - 2.0*a*a;
#endif
#ifdef SAMPLE_11
    return 1.0 - 2.0*a*a;
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
#ifdef SAMPLE_7
    return x*x - 2.0*t;
#endif
#ifdef SAMPLE_8
    return 0.0;
#endif
#ifdef SAMPLE_9
    double k = 20.0;
    return exp(k*x)*(1.0 - a*a*k*k*t);
#endif
}

double BorderParabolicD::U(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    double t = j*ht;
#ifdef SAMPLE_0
    return x*x*t;
#endif
#ifdef SAMPLE_1
    return x*x + t*t;
#endif
#ifdef SAMPLE_11
    return x*x + t;
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
#ifdef SAMPLE_7
    return x*x*t;
#endif
#ifdef SAMPLE_8
    return 0.0;
#endif
#ifdef SAMPLE_9
    double k = 20.0;
    return exp(k*x)*t;
#endif
}
