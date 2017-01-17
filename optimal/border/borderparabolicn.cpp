#include "borderparabolicn.h"
#include <math.h>

void BorderParabolicN::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    BorderParabolicN bp;
    bp.a = 1.0;

    // Real solution
    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;
        DoubleMatrix ru(bp.M+1,bp.N+1);
        for (unsigned int i=0; i<=bp.M; i++)
        {
            for (unsigned int j=0; j<=bp.N; j++)
                ru.at(i,j) = bp.U(j,i);
        }
        IPrinter::printMatrix(14, 10, ru, 10, 10, NULL);
        ru.clear();
        IPrinter::printSeperatorLine();
    }

    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;

        DoubleMatrix u;
        bp.calculateN2(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
        u.clear();
        IPrinter::printSeperatorLine();
    }

    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;

        DoubleMatrix u;
        bp.calculateN3(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
        u.clear();
        IPrinter::printSeperatorLine();
    }

    return;

    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;
        DoubleMatrix u;
        bp.calculateN4L2RN(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
        u.clear();
        IPrinter::printSeperatorLine();
    }

    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;
        DoubleMatrix u;
        bp.calculateN4R2LN(u, bp.hx, bp.ht, bp.N, bp.M, bp.a);
        IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
        IPrinter::printSeperatorLine();
    }
}

double BorderParabolicN::initial(unsigned int i UNUSED_PARAM) const
{
    return U(i,0);
}

double BorderParabolicN::boundary(Boundary type UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    double t UNUSED_PARAM = j*ht;

#ifdef SAMPLE_1
    if (type == Left)  return 0.0;
    if (type == Right) return 2.0;
#endif
#ifdef SAMPLE_2
    if (type == Left)  return 0.0;
    if (type == Right) return 2.0;
#endif
#ifdef SAMPLE_3
    if (type == Left)  return 1.0;
    if (type == Right) return 1.0;
#endif
#ifdef SAMPLE_4
    if (type == Left)  return t;
    if (type == Right) return sin(20.0) + 20.0*cos(20.0) + t;
#endif
#ifdef SAMPLE_5
    if (type == Left)  return 12.0*t;
    if (type == Right) return 10.0*t*cos(10.0) + 2.0*t*exp(2.0*t);
#endif
    return 0.0;
}

double BorderParabolicN::f(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    double t UNUSED_PARAM = j*ht;
    double x UNUSED_PARAM = i*hx;
#ifdef SAMPLE_1
    return 2.0*t - 2.0*a*a;
#endif
#ifdef SAMPLE_2
    return 1.0 - 2.0*a*a;
#endif
#ifdef SAMPLE_3
    return 1.0;
#endif
#ifdef SAMPLE_4
    return x - a*a*(40.0*cos(20.0*x) - 400.0*x*sin(20.0*x));
#endif
#ifdef SAMPLE_5
    return (sin(10.0*x) + 2.0*x*exp(2.0*x*t) - a*a*(4.0*t*t*exp(2.0*x*t) - 100.0*t*sin(10.0*x)));
#endif
}

double BorderParabolicN::U(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    double t = j*ht;
#ifdef SAMPLE_1
    return x*x + t*t;
#endif
#ifdef SAMPLE_2
    return x*x + t;
#endif
#ifdef SAMPLE_3
    return x + t;
#endif
#ifdef SAMPLE_4
    return x*sin(20.0*x) + t*x;
#endif
#ifdef SAMPLE_5
    return t*sin(10.0*x) + exp(2.0*x*t);
#endif
}

void BorderParabolicN::calculateN2(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a)
{
    u.clear();
    u.resize(M+1, N+1);

    for (unsigned int i=0; i<=N; i++) u[0][i] = initial(i);

    DoubleMatrix A(2,2,0.0);
    DoubleVector b(2,0.0);
    DoubleVector x(2,0.0);

    DoubleMatrix ems(N+1, 3);

    for (unsigned int m=1; m<=M; m++)
    {
        A[0][0] = 1.0 + (a*a*ht)/(2.0*hx*hx);
        A[0][1] = -(a*a*ht)/(2.0*hx*hx);
        b[0]    = u[m-1][0] + ht * f(0, m) - (a*a*ht)/(hx) * boundary(Left, m);

        A[0][1] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;

        A[1][0] = -(a*a*ht)/(2.0*hx*hx);
        A[1][1] = 1.0 + (a*a*ht)/(2.0*hx*hx);
        b[1]    = u.at(m-1,N) + ht * f(N, m) + (a*a*ht)/(hx) * boundary(Right, m);

        A[1][0] /= A[1][1];
        b[1]    /= A[1][1];
        A[1][1] = 1.0;

        ems.at(0, 0) = A[0][0];
        ems.at(0, 1) = A[0][1];
        ems.at(0, 2) = b[0];

        ems.at(N, 0) = A[1][0];
        ems.at(N, 1) = A[1][1];
        ems.at(N, 2) = b[1];

        for (unsigned int n=1; n<=N-1; n++)
        {
            double g1 = (1.0 + (2.0*a*a*ht)/(hx*hx)) / ((a*a*ht)/(hx*hx));
            double g2 = -1.0;
            double g0 = (u[m-1][n] + ht * f(n, m)) / (-(a*a*ht)/(hx*hx));

            if (n==1)
            {
                A[0][0] = g1;
                A[0][1] = A[0][1] + g2;
                b[0] = b[0] - g0;

                A[0][1] /= A[0][0];
                b[0]    /= A[0][0];
                A[0][0] = 1.0;

                ems.at(n, 0) = A[0][0];
                ems.at(n, 1) = A[0][1];
                ems.at(n, 2) = b[0];
            }
            else if (n == N-1)
            {
                g1 = -g1/g2;
                g0 = -g0/g2;
                g2 = +1.0/g2;

                A[1][0] = A[1][0] + g2;
                A[1][1] = g1;
                b[1] = b[1] - g0;

                A[1][0] /= A[1][1];
                b[1]    /= A[1][1];
                A[1][1] = 1.0;

                ems.at(n, 0) = A[1][0];
                ems.at(n, 1) = A[1][1];
                ems.at(n, 2) = b[1];
            }
            else
            {
                A[0][0] = A[0][1] + g1;
                A[0][1] = g2;
                b[0] = b[0] - g0;

                A[0][1] /= A[0][0];
                b[0]    /= A[0][0];
                A[0][0] = 1.0;

                ems.at(n, 0) = A[0][0];
                ems.at(n, 1) = A[0][1];
                ems.at(n, 2) = b[0];
            }
        }

        GaussianElimination(A, b, x);

        u.at(m, N-1) = x[1];
        u.at(m, N-2) = x[0];

        u.at(m, N) = -ems.at(N,0)*u.at(m,N-2) + ems.at(N,2);
        for (unsigned int n=N-3; n>=1; n--)
        {
            u.at(m, n) = -ems.at(n,1)*u.at(m,n+1) + ems.at(n,2);
            u.at(m, n) /= ems.at(n, 0);
        }
        u.at(m, 0) = -ems.at(0,1)*u.at(m,2) + ems.at(0,2);
    }
}

void BorderParabolicN::calculateN3(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a)
{
    u.clear();
    u.resize(M+1, N+1);

    for (unsigned int i=0; i<=N; i++) u[0][i] = initial(i);

    DoubleMatrix A(2,2,0.0);
    DoubleVector b(2,0.0);
    DoubleVector x(2,0.0);

    DoubleMatrix ems(N+1, 3);

    for (unsigned int m=1; m<=M; m++)
    {
        A[0][0] = 1.0 + (a*a*ht)/(hx*hx);
        A[0][1] = -(a*a*ht)/(hx*hx);
        b[0]    = u[m-1][0] + ht * f(0, m) - (a*a*ht)/(hx) * boundary(Left, m);

        A[0][1] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;

        ems.at(0, 0) = A[0][0];
        ems.at(0, 1) = A[0][1];
        ems.at(0, 2) = b[0];

        //printf("%4u %14.10f %14.10f\n", 0, A[0][0]*U(0,m)+A[0][1]*U(1,m), b[0]);

        for (unsigned int n=1; n<=N-1; n++)
        {
            double g1 = (1.0 + (2.0*a*a*ht)/(hx*hx)) / ((a*a*ht)/(hx*hx));
            double g2 = -1.0;
            double g0 = (u[m-1][n] + ht * f(n, m)) / (-(a*a*ht)/(hx*hx));

            A[0][0] = A[0][1] + g1;
            A[0][1] = g2;
            b[0] = b[0] - g0;

            A[0][1] /= A[0][0];
            b[0]    /= A[0][0];
            A[0][0] = 1.0;

            ems.at(n, 0) = A[0][0];
            ems.at(n, 1) = A[0][1];
            ems.at(n, 2) = b[0];
        }

        A[1][0] = -(a*a*ht)/(hx*hx);
        A[1][1] = 1.0 + (a*a*ht)/(hx*hx);
        b[1]    = u.at(m-1,N) + ht * f(N, m) + (a*a*ht)/(hx) * boundary(Right, m);

        A[1][0] /= A[1][1];
        b[1]    /= A[1][1];
        A[1][1] = 1.0;

        ems.at(N, 0) = A[1][0];
        ems.at(N, 1) = A[1][1];
        ems.at(N, 2) = b[1];

        //printf("%4u %14.10f %14.10f\n", N, A[1][0]*U(N-1,m)+A[1][1]*U(N,m), b[1]);

        GaussianElimination(A, b, x);

        u.at(m, N-0) = x[1];
        u.at(m, N-1) = x[0];

        //printf("%14.10f %14.10f\n", x[0], x[1]);
        //printf("%14.10f %14.10f\n", U(N-1,m), U(N-0,m));

        for (unsigned int n=N-2; n!=UINT_MAX; n--)
        {
            u.at(m, n) = -ems.at(n,1)*u.at(m,n+1) + ems.at(n,2);
            u.at(m, n) /= ems.at(n, 0);
        }
        //u.at(m, 0) = -ems.at(0,1)*u.at(m,2) + ems.at(0,2);
        //IPrinter::printVector(u.row(m));
        //break;
    }
}

