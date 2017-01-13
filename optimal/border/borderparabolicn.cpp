#include "borderparabolicn.h"
#include <math.h>

void BorderParabolicN::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    BorderParabolicN bp;

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

void BorderParabolicN::calculateN2(DoubleMatrix &u, double hx, double ht, double N, double M, double a)
{
    u.clear();
    u.resize(M+1, N+1);

    DoubleVector da(N+1);
    DoubleVector db(N+1);
    DoubleVector dc(N+1);
    DoubleVector dd(N+1);
    DoubleVector rx(N+1);

    for (unsigned int i=0; i<=N; i++) u[0][i] = initial(i);

    DoubleMatrix A(2,2,0.0);
    DoubleVector b(2,0.0);
    DoubleVector x(2,0.0);

    for (unsigned int m=1; m<=M; m++)
    {
        A[0][0] = 1.0 + (a*a*ht)/(4.0*hx*hx);
        A[0][1] = -(a*a*ht)/(4.0*hx*hx);
        b[0]    = u[m-1][0] + ht * f(0, m) - (a*a*ht)/(2.0*hx) * boundary(Left, m);

        printf("%14.10f %14.10f\n", A[0][0]*U(0,m) + A[0][1]*U(1,m), b[0]);

        A[0][1] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;

        for (unsigned int n=1; n<=N-2; n++)
        {
            double g1 = (1.0 + (2.0*a*a*ht)/(hx*hx)) / ((a*a*ht)/(hx*hx));
            double g2 = -1.0;
            double g0 = (u[m-1][n] + ht * f(n, m)) / ((a*a*ht)/(hx*hx));

            A[0][0] = A[0][1] + g1;
            A[0][1] = g2;
            b[0] = b[0] - g0;

            A[0][1] /= A[0][0];
            b[0]    /= A[0][0];
            A[0][0] = 1.0;

            //if (n<10)
            //printf("%14.10f %14.10f %14.10f %14.10f\n", u[m-1][0], boundary(Left, m), A[0][0]*U(0,m) + A[0][1]*U(1,m), b[0]);

        }

        A[1][0] = -(a*a*ht)/(hx*hx);
        A[1][1] = 1.0 + (a*a*ht)/(hx*hx);
        b[1]    = u.at(m-1,N) + ht * f(N, m) + (a*a*ht)/(hx) * boundary(Right, m);
        printf("%14.10f %14.10f\n", A[1][0]*U(N-1,m) + A[1][1]*U(N,m), b[1]);

        GaussianElimination(A, b, x);

        printf("%14.10f %14.10f\n", x[0], x[1]);

        break;

//        for (unsigned int i=0; i<=N; i++)
//        {
//            u[j][i] = rx[i];
//        }
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();
}

void BorderParabolicN::calculateN4L2RN(DoubleMatrix &u, double hx, double ht, double N, double M, double a)
{
    unsigned int k = 4;
    double alpha0 = (a*a*ht)/hx;
    double alpha1 = (a*a*ht)/(12.0*hx*hx);
    double alpha2 = (a*a*ht)/(24.0*hx*hx);
    u.resize(M+1, N+1);

    double D1[5][5] =
    {
        {-25.0, +48.0, -36.0, +16.0, -3.0},
        {-3.0,  -10.0, +18.0, -6.0,  +1.0},
        {+1.0,  -8.0,  +0.0,  +8.0,  -1.0},
        {-1.0,  +6.0, -18.0,  +10.0, +3.0},
        {+3.0,  -16.0, +36.0, -48.0, +25.0}
    };

    double D2[k+1][k+1] =
    {
        {+70.0, -208.0, +228.0, -112.0, +22.0},
        {+22.0, -40.0,  +12.0,  +8.0,   -2.0},
        {-2.0,  +32.0,  -60.0,  +32.0,  -2.0},
        {-2.0,  +8.0,   +12.0,  -40.0,  +22.0},
        {+22.0, -112.0, +228.0, -208.0, +70.0}
    };

    DoubleMatrix A(k+1, k+1, 0.0);
    DoubleVector b(k+1, 0.0);
    DoubleVector x(k+1, 0.0);
    DoubleMatrix ems(N-(k-1), k+1);

    /* initial condition */
    for (unsigned int i=0; i<=N; i++) u.at(0,i) = initial(i);

    for (unsigned int m=1; m<=M; m++)
    {
        //0
        A[0][0] = -3.0*alpha1 - 1.0;
        A[0][1] = -10.0*alpha1;
        A[0][2] = +18.0*alpha1;
        A[0][3] = -6.0*alpha1;
        A[0][4] = +alpha1;
        b[0]    = -u.at(m-1,0) - ht*f(0,m) + alpha0*boundary(Left,m);

        A[0][1] /= A[0][0];
        A[0][2] /= A[0][0];
        A[0][3] /= A[0][0];
        A[0][4] /= A[0][0];
        b.at(0) /= A[0][0];
        A[0][0] = 1.0;

        ems.at(0,0) = A[0][1];
        ems.at(0,1) = A[0][2];
        ems.at(0,2) = A[0][3];
        ems.at(0,3) = A[0][4];
        ems.at(0,4) = b[0];

        for (unsigned int n=0; n<=N-5; n++)
        {
            double g0 = +70.0*alpha2 - 1.0;
            double g1 = -208.0*alpha2;
            double g2 = +228.0*alpha2;
            double g3 = -112.0*alpha2;
            double g4 = +22.0*alpha2;
            double fi = -u.at(m-1,n) - ht*f(n,m);

            g1 /= -g0;
            g2 /= -g0;
            g3 /= -g0;
            g4 /= -g0;
            fi /= +g0;
            g0  = 1.0;

            double a00 = A[0][0];
            A[0][0] = A[0][1] + a00 * g1;
            A[0][1] = A[0][2] + a00 * g2;
            A[0][2] = A[0][3] + a00 * g3;
            A[0][3] = A[0][4] + a00 * g4;
            A[0][4] = 0.0;
            b.at(0)  = b.at(0) - a00*fi;

            A[0][1] /= A[0][0];
            A[0][2] /= A[0][0];
            A[0][3] /= A[0][0];
            A[0][4] /= A[0][0];
            b.at(0) /= A[0][0];
            A[0][0] = 1.0;

            ems.at(n+1,0) = A[0][1];
            ems.at(n+1,1) = A[0][2];
            ems.at(n+1,2) = A[0][3];
            ems.at(n+1,3) = A[0][4];
            ems.at(n+1,4) = b.at(0);
        }

        //N-3
        A[1][0] = +22.0*alpha2;
        A[1][1] = -40.0*alpha2 - 1.0;
        A[1][2] = +12.0*alpha2;
        A[1][3] = +8.0*alpha2;
        A[1][4] = -2.0*alpha2;
        b[1]    = -u.at(m-1,N-3) - ht*f(N-3,m);

        //N-2
        A[2][0] = -2.0*alpha2;
        A[2][1] = +32.0*alpha2;
        A[2][2] = -60.0*alpha2 - 1.0;
        A[2][3] = +32.0*alpha2;
        A[2][4] = -2.0*alpha2;
        b[2]    = -u.at(m-1,N-2) - ht*f(N-2,m);

        //N-1
        A[3][0] = -2.0*alpha2;
        A[3][1] = +8.0*alpha2;
        A[3][2] = +12.0*alpha2;
        A[3][3] = -40.0*alpha2 - 1.0;
        A[3][4] = +22.0*alpha2;
        b[3]    = -u.at(m-1,N-1) - ht*f(N-1,m);

        //N
        A[4][0] = -alpha1;
        A[4][1] = +6.0*alpha1;
        A[4][2] = -18.0*alpha1;
        A[4][3] = +10.0*alpha1;
        A[4][4] = +3.0*alpha1 + 1.0;
        b[4]    = +u.at(m-1,N) + ht*f(N,m) + alpha0*boundary(Right,m);

        GaussianElimination(A, b, x);

        u.at(m, N-4) = x.at(0);
        u.at(m, N-3) = x.at(1);
        u.at(m, N-2) = x.at(2);
        u.at(m, N-1) = x.at(3);
        u.at(m, N-0) = x.at(4);
        //printf("%4u %14.10f %14.10f %14.10f %14.10f %14.10f\n", m, U(N-4,m), U(N-3,m), U(N-2,m), U(N-1,m), U(N,m));
        //printf("%4u %14.10f %14.10f %14.10f %14.10f %14.10f\n", m, u.at(m,N-4), u.at(m,N-3), u.at(m,N-2), u.at(m,N-1), u.at(m,N));
        for (unsigned int n=N-5; n!=UINT32_MAX; n--)
        {
            u.at(m,n) = -ems.at(n,0)*u.at(m,n+1)
                    -ems.at(n,1)*u.at(m,n+2)
                    -ems.at(n,2)*u.at(m,n+3)
                    -ems.at(n,3)*u.at(m,n+4)
                    +ems.at(n,4);
            //printf("%14.10f\n", u1.at(m,n));
        }

        //IPrinter::printVector(u1.row(m));
        //if (m>10) break;
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void BorderParabolicN::calculateN4R2LN(DoubleMatrix &u, double hx, double ht, double N, double M, double a)
{
    unsigned int k = 4;
    double alpha0 = (a*a*ht)/hx;
    double alpha1 = (a*a*ht)/(12.0*hx*hx);
    double alpha2 = (a*a*ht)/(24.0*hx*hx);
    u.resize(M+1, N+1);

    double D2[k+1][k+1] =
    {
        {+70.0, -208.0, +228.0, -112.0, +22.0},
        {+22.0, -40.0,  +12.0,  +8.0,   -2.0},
        {-2.0,  +32.0,  -60.0,  +32.0,  -2.0},
        {-2.0,  +8.0,   +12.0,  -40.0,  +22.0},
        {+22.0, -112.0, +228.0, -208.0, +70.0}
    };

    double D1[5][5] =
    {
        {-25.0, +48.0, -36.0, +16.0, -3.0},
        {-3.0,  -10.0, +18.0, -6.0,  +1.0},
        {+1.0,  -8.0,  +0.0,  +8.0,  -1.0},
        {-1.0,  +6.0, -18.0,  +10.0, +3.0},
        {+3.0,  -16.0, +36.0, -48.0, +25.0}
    };

    DoubleMatrix A(k+1, k+1, 0.0);
    DoubleVector b(k+1, 0.0);
    DoubleVector x(k+1, 0.0);
    DoubleMatrix ems(N-(k-1), k+1);

    /* initial condition */
    for (unsigned int i=0; i<=N; i++) u.at(0,i) = initial(i);

    for (unsigned int m=1; m<=M; m++)
    {
        //0
        A[0][0] = -3.0*alpha1 - 1.0;
        A[0][1] = -10.0*alpha1;
        A[0][2] = +18.0*alpha1;
        A[0][3] = -6.0*alpha1;
        A[0][4] = +alpha1;
        b[0]    = -u.at(m-1,0) - ht*f(0,m) + alpha0*boundary(Left,m);

        //1
        A[1][0] = +22.0*alpha2;
        A[1][1] = -40.0*alpha2 - 1.0;
        A[1][2] = +12.0*alpha2;
        A[1][3] = +8.0*alpha2;
        A[1][4] = -2.0*alpha2;
        b[1]    = -u.at(m-1,1) - ht*f(1,m);

        //2
        A[2][0] = -2.0*alpha2;
        A[2][1] = +32.0*alpha2;
        A[2][2] = -60.0*alpha2 - 1.0;
        A[2][3] = +32.0*alpha2;
        A[2][4] = -2.0*alpha2;
        b[2]    = -u.at(m-1,2) - ht*f(2,m);

        //3
        A[3][0] = -2.0*alpha2;
        A[3][1] = +8.0*alpha2;
        A[3][2] = +12.0*alpha2;
        A[3][3] = -40.0*alpha2 - 1.0;
        A[3][4] = +22.0*alpha2;
        b[3]    = -u.at(m-1,3) - ht*f(3,m);

        //N
        A[4][0] = -alpha1;
        A[4][1] = +6.0*alpha1;
        A[4][2] = -18.0*alpha1;
        A[4][3] = +10.0*alpha1;
        A[4][4] = +3.0*alpha1 + 1.0;
        b[4]    = +u.at(m-1,N) + ht*f(N,m) + alpha0*boundary(Right,m);

        A[4][0] /= A[4][4];
        A[4][1] /= A[4][4];
        A[4][2] /= A[4][4];
        A[4][3] /= A[4][4];
        b[4]    /= A[4][4];
        A[4][4] = 1.0;

        ems.at(N-4,0) = A[4][0];
        ems.at(N-4,1) = A[4][1];
        ems.at(N-4,2) = A[4][2];
        ems.at(N-4,3) = A[4][3];
        ems.at(N-4,4) = b[4];

        for (unsigned int n=N; n>=5; n--)
        {
            double g1 = +22.0*alpha2;
            double g2 = -112.0*alpha2;
            double g3 = +228.0*alpha2;
            double g4 = -208.0*alpha2;
            double g5 = +70.0*alpha2 - 1.0;
            double fi = -u.at(m-1,n) - ht*f(n,m);

            g4 /= -g5;
            g3 /= -g5;
            g2 /= -g5;
            g1 /= -g5;
            fi /= +g5;
            g5  = 1.0;

            double a00 = A[4][4];
            A[4][4] = A[4][3] + a00 * g4;
            A[4][3] = A[4][2] + a00 * g3;
            A[4][2] = A[4][1] + a00 * g2;
            A[4][1] = A[4][0] + a00 * g1;
            b[4]    = b[4] - a00*fi;
            A[4][0] = 0.0;

            A[4][3] /= A[4][4];
            A[4][2] /= A[4][4];
            A[4][1] /= A[4][4];
            A[4][0] /= A[4][4];
            b[4]    /= A[4][4];
            A[4][4] = 1.0;

            ems.at(n-5,0) = A[4][0];
            ems.at(n-5,1) = A[4][1];
            ems.at(n-5,2) = A[4][2];
            ems.at(n-5,3) = A[4][3];
            ems.at(n-5,4) = b[4];
        }

        GaussianElimination(A, b, x);

        u.at(m, 0) = x.at(0);
        u.at(m, 1) = x.at(1);
        u.at(m, 2) = x.at(2);
        u.at(m, 3) = x.at(3);
        u.at(m, 4) = x.at(4);
        for (unsigned int n=5; n<=N; n++)
        {
            u.at(m,n) = -ems.at(n-4,3)*u.at(m,n-1) - ems.at(n-4,2)*u.at(m,n-2) - ems.at(n-4,1)*u.at(m,n-3) - ems.at(n-4,0)*u.at(m,n-4) + ems.at(n-4,4);
        }
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}
