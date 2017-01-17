#include "borderparabolicd.h"
#include <math.h>

void BorderParabolicD::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    BorderParabolicD bp;
    bp.a = 1.0;

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
        IPrinter::printSeperatorLine("calculateN4L2RD");
        bp.calculateN2L2RD(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
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
#ifndef SAMPLE_8
    return U(i,0);
#endif
#ifdef SAMPLE_8
    return 1.0;
#endif
}

double BorderParabolicD::boundary(Boundary type UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
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
}

double BorderParabolicD::U(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    double t = j*ht;
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
}

void BorderParabolicD::calculateN2L2RD(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a)
{
    unsigned int k = 2;
    u.clear();
    u.resize(M+1, N+1);

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector x(k, 0.0);
    DoubleMatrix ems(N-k, k);

    /* initial condition */
    for (unsigned int i=0; i<=N; i++) u.at(0,i) = initial(i);

    /* border conditions */
    for (unsigned int j=1; j<=M; j++)
    {
        u.at(j,0) = boundary(Left, j);
        u.at(j,N) = boundary(Right, j);
    }

    //double alpha = -(a*a*ht)/(hx*hx);
    //double beta  = 1.0 + (2.0*a*a*ht)/(hx*hx);

    double alpha = (a*a*ht)/(hx*hx);

    for (unsigned int m=1; m<=M; m++)
    {
        A[0][0] = -2.0*alpha - 1.0;
        A[0][1] = +1.0*alpha;
        b[0]    = -u.at(m-1,1) - alpha*u.at(m,0) - ht*f(1,m);

        A[0][1] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;

        ems.at(0,0) = A[0][1];
        ems.at(0,1) = b[0];

        for (unsigned int n=2; n<=N-2; n++)
        {
            double g1 = +1.0*alpha;
            double g2 = -2.0*alpha - 1.0;
            double g3 = +1.0*alpha;
            double fi = -u.at(m-1,n) - ht*f(n,m);

            g2 /= -g1;
            g3 /= -g1;
            fi /= +g1;
            g1 = 1.0;

            A[0][0] = A[0][1] + g2;
            A[0][1] = g3;
            b[0]    = b[0] - fi;
            \
            A[0][1] /= A[0][0];
            b[0]    /= A[0][0];
            A[0][0] = 1.0;

            ems.at(n-1,0) = A[0][1];
            ems.at(n-1,1) = b[0];
        }

        A[1][0] = +alpha;
        A[1][1] = -2.0*alpha - 1.0;
        b[1]    = -u.at(m-1,N-1) - alpha*u.at(m,N) - ht*f(N-1,m);

        GaussianElimination(A, b, x);

        //printf("%14.10f %14.10f\n", x[0], x[1]);
        //break;

        u.at(m, N-1) = x.at(1);
        u.at(m, N-2) = x.at(0);
        for (unsigned int i=N-3; i>=1; i--)
        {
            u.at(m,i) = -ems.at(i-1,0)*u.at(m,i+1) + ems.at(i-1,1);
            u.at(m,i) = -ems.at(i-1,0)*u.at(m,i+1) + ems.at(i-1,1);
        }

//        for (unsigned int i=N-2; i>=1; i--)
//        {
//            u.at(m,i-1) = (1.0+2.0*alpha)*u.at(m,i) - alpha*u.at(m,i+1) - u.at(m-1,i);
//            u.at(m,i-1) /= alpha;
//            if (i == N-6) break;
//        }
        printf("%14.10f %14.10f %14.10f %14.10f %14.10f\n", u.at(m,N-5), u.at(m,N-4), u.at(m,N-3), u.at(m,N-2), u.at(m,N-1));
        printf("%14.10f %14.10f %14.10f %14.10f %14.10f\n", U(N-5,m), U(N-4,m), U(N-3,m), U(N-2,m), U(N-1,m));

        IPrinter::printVector(u.row(m));
        break;
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}
