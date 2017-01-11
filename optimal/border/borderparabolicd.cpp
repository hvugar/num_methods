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
        for (unsigned int i=0; i<=bp.M; i++)
        {
            for (unsigned int j=0; j<=bp.N; j++)
                ru.at(i,j) = bp.U(j,i);
        }
        IPrinter::printSeperatorLine("Real solution");
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
        bp.calculateN4L2RM(u);
        IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
    }

    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;
        DoubleMatrix u;
        IPrinter::printSeperatorLine("calculateN4R2LM");
        bp.calculateN4R2LM(u);
        IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
    }

    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;
        DoubleMatrix u;
        IPrinter::printSeperatorLine("calculateN6L2RM");
        bp.calculateN6L2RM(u);
        IPrinter::printMatrix(14, 10, u, 10, 10, NULL);
    }

    {
        bp.hx = 0.01;
        bp.ht = 0.01;
        bp.N = 100;
        bp.M = 100;
        DoubleMatrix u;
        IPrinter::printSeperatorLine("calculateN6R2LM");
        bp.calculateN6R2LM(u);
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

void BorderParabolicD::calculateN4L2RM(DoubleMatrix &u)
{
    unsigned int k = 4;
    double alpha = (ht*a*a)/(24.0*hx*hx);
    u.resize(M+1, N+1);

    double D[k+1][k+1] =
    {
        {+70.0, -208.0, +228.0, -112.0, +22.0},
        {+22.0, -40.0,  +12.0,  +8.0,   -2.0},
        {-2.0,  +32.0,  -60.0,  +32.0,  -2.0},
        {-2.0,  +8.0,   +12.0,  -40.0,  +22.0},
        {+22.0, -112.0, +228.0, -208.0, +70.0}
    };

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

    for (unsigned int m=1; m<=M; m++)
    {
        A.at(0,0) = D[1][1]*alpha - 1.0;
        A.at(0,1) = D[1][2]*alpha;
        A.at(0,2) = D[1][3]*alpha;
        A.at(0,3) = D[1][4]*alpha;
        b.at(0)   = -u.at(m-1,1) - (D[1][0]*alpha)*u.at(m,0) - ht*f(1,m);

        A.at(0,1) /= A.at(0,0);
        A.at(0,2) /= A.at(0,0);
        A.at(0,3) /= A.at(0,0);
        b.at(0)   /= A.at(0,0);
        A.at(0,0) = 1.0;

        ems.at(0,0) = A.at(0,1);
        ems.at(0,1) = A.at(0,2);
        ems.at(0,2) = A.at(0,3);
        ems.at(0,3) = b.at(0);

        // + * * * *
        for (unsigned int n=1; n<=N-(k+1); n++)
        {
            double g1 = D[0][0]*alpha-1.0;
            double g2 = D[0][1]*alpha;
            double g3 = D[0][2]*alpha;
            double g4 = D[0][3]*alpha;
            double g5 = D[0][4]*alpha;
            double g0 = u.at(m-1,n) + ht*f(n,m);

            g2 /= -g1;
            g3 /= -g1;
            g4 /= -g1;
            g5 /= -g1;
            g0 /= -g1;
            g1 = 1.0;

            A.at(0,0) = A.at(0,1) + g2;
            A.at(0,1) = A.at(0,2) + g3;
            A.at(0,2) = A.at(0,3) + g4;
            A.at(0,3) = g5;
            b.at(0)   = b.at(0) - g0;
            \
            A.at(0,1) /= A.at(0,0);
            A.at(0,2) /= A.at(0,0);
            A.at(0,3) /= A.at(0,0);
            b.at(0)   /= A.at(0,0);
            A.at(0,0) = 1.0;

            ems.at(n,0) = A.at(0,1);
            ems.at(n,1) = A.at(0,2);
            ems.at(n,2) = A.at(0,3);
            ems.at(n,3) = b.at(0);
        }

        A.at(1,0) = D[1][0]*alpha;
        A.at(1,1) = D[1][1]*alpha - 1.0;
        A.at(1,2) = D[1][2]*alpha;
        A.at(1,3) = D[1][3]*alpha;
        b.at(1)   = -u.at(m-1,N-3) - (D[1][4]*alpha)*u.at(m,N) - ht*f(N-3,m);

        A.at(2,0) = D[2][0]*alpha;
        A.at(2,1) = D[2][1]*alpha;
        A.at(2,2) = D[2][2]*alpha - 1.0;
        A.at(2,3) = D[2][3]*alpha;
        b.at(2)   = -u.at(m-1,N-2) - (D[2][4]*alpha)*u.at(m,N) - ht*f(N-2,m);

        A.at(3,0) = D[3][0]*alpha;
        A.at(3,1) = D[3][1]*alpha;
        A.at(3,2) = D[3][2]*alpha;
        A.at(3,3) = D[3][3]*alpha - 1.0;
        b.at(3)   = -u.at(m-1,N-1) - (D[3][4]*alpha)*u.at(m,N) - ht*f(N-1,m);

        GaussianElimination(A, b, x);

        u.at(m, N-1) = x.at(3);
        u.at(m, N-2) = x.at(2);
        u.at(m, N-3) = x.at(1);
        u.at(m, N-4) = x.at(0);
        for (unsigned int i=N-(k+1); i>=1; i--)
        {
            u.at(m,i) = -ems.at(i-1,0)*u.at(m,i+1) - ems.at(i-1,1)*u.at(m,i+2) - ems.at(i-1,2)*u.at(m,i+3) + ems.at(i-1,3);
            //u.at(m,i) = -208.0*alpha*u.at(m,i+1) + 228.0*alpha*u.at(m,i+2) - 112.0*alpha*u.at(m,i+3) + 22.0*alpha*u.at(m,i+4) + (u.at(m-1,i)+ht*f(i,m));
            //u.at(m,i) /= -(70.0*alpha-1.0);
        }

        if (m==0)
        {
            printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, U(N-4,m), U(N-3,m), U(N-2,m), U(N-1,m));
            printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, x[0], x[1], x[2], x[3]);
            //printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, u1.at(m,N-4), u1.at(m,N-3), u1.at(m,N-2), u1.at(m,N-1));
            IPrinter::printVector(14,10, u.row(m));
            DoubleVector v(N+1); for (unsigned int n=0; n<=N; n++) v.at(n) = U(n,m);
            IPrinter::printVector(14,10, v);
            v.clear();
            IPrinter::printSeperatorLine();
            break;
        }
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void BorderParabolicD::calculateN4R2LM(DoubleMatrix &u)
{
    unsigned int k = 4;
    double alpha = (a*a*ht)/(24.0*hx*hx);
    u.resize(M+1, N+1);

    double D[k+1][k+1] =
    {
        {+70.0, -208.0, +228.0, -112.0, +22.0},
        {+22.0, -40.0,  +12.0,  +8.0,   -2.0},
        {-2.0,  +32.0,  -60.0,  +32.0,  -2.0},
        {-2.0,  +8.0,   +12.0,  -40.0,  +22.0},
        {+22.0, -112.0, +228.0, -208.0, +70.0}
    };

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector x(k, 0.0);
    DoubleMatrix ems(N-k, k);

    /* initial condition */
    for (unsigned int i=0; i<=N; i++) u.at(0,i) = initial(i);

    /* border conditions */
    for (unsigned int j=1; j<=M; j++)
    {
        u.at(j,0) = boundary(Left,j);
        u.at(j,N) = boundary(Right, j);
    }

    for (unsigned int m=1; m<=M; m++)
    {
        A.at(0,0) = D[1][1]*alpha - 1.0;
        A.at(0,1) = D[1][2]*alpha;
        A.at(0,2) = D[1][3]*alpha;
        A.at(0,3) = D[1][4]*alpha;
        b.at(0)   = -u.at(m-1,1) - (D[1][0]*alpha)*u.at(m,0) - ht*f(1,m);

        A.at(1,0) = D[2][1]*alpha;
        A.at(1,1) = D[2][2]*alpha - 1.0;
        A.at(1,2) = D[2][3]*alpha;
        A.at(1,3) = D[2][4]*alpha;
        b.at(1)   = -u.at(m-1,2) - (D[2][0]*alpha)*u.at(m,0) - ht*f(2,m);

        A.at(2,0) = D[3][1]*alpha;
        A.at(2,1) = D[3][2]*alpha;
        A.at(2,2) = D[3][3]*alpha - 1.0;
        A.at(2,3) = D[3][4]*alpha;
        b.at(2)   = -u.at(m-1,3) - (D[3][0]*alpha)*u.at(m,0) - ht*f(3,m);

        A.at(3,0) = D[3][0]*alpha;
        A.at(3,1) = D[3][1]*alpha;
        A.at(3,2) = D[3][2]*alpha;
        A.at(3,3) = D[3][3]*alpha - 1.0;
        b.at(3)   = -u.at(m-1,N-1) - (D[3][4]*alpha)*u.at(m,N) - ht*f(N-1,m);

        A.at(3,0) /= A.at(3,3);
        A.at(3,1) /= A.at(3,3);
        A.at(3,2) /= A.at(3,3);
        b.at(3)   /= A.at(3,3);
        A.at(3,3) = 1.0;

        ems.at(N-5,0) = A.at(3,0);
        ems.at(N-5,1) = A.at(3,1);
        ems.at(N-5,2) = A.at(3,2);
        ems.at(N-5,3) = b.at(3);

        for (unsigned int n=N-1; n>=k+1; n--)
        {
            double g1 = D[4][0]*alpha;
            double g2 = D[4][1]*alpha;
            double g3 = D[4][2]*alpha;
            double g4 = D[4][3]*alpha;;
            double g5 = D[4][4]*alpha - 1.0;
            double g0  = u.at(m-1,n) + ht*f(n,m);

            g4 /= -g5;
            g3 /= -g5;
            g2 /= -g5;
            g1 /= -g5;
            g0 /= -g5;
            g5 = 1.0;

            A.at(3,3) = A.at(3,2) + g4;
            A.at(3,2) = A.at(3,1) + g3;
            A.at(3,1) = A.at(3,0) + g2;
            A.at(3,0) = g1;
            b.at(3)   = b.at(3) - g0;

            A.at(3,2) /= A.at(3,3);
            A.at(3,1) /= A.at(3,3);
            A.at(3,0) /= A.at(3,3);
            b.at(3)   /= A.at(3,3);
            A.at(3,3) = 1.0;

            ems.at(n-5,0) = A.at(3,0);
            ems.at(n-5,1) = A.at(3,1);
            ems.at(n-5,2) = A.at(3,2);
            ems.at(n-5,3) = b.at(3);
        }

        GaussianElimination(A, b, x);

        u.at(m, 1) = x.at(0);
        u.at(m, 2) = x.at(1);
        u.at(m, 3) = x.at(2);
        u.at(m, 4) = x.at(3);
        for (unsigned int i=k+1; i<=N-1; i++)
        {
            u.at(m,i) = -ems.at(i-4,2)*u.at(m,i-1) - ems.at(i-4,1)*u.at(m,i-2) - ems.at(i-4,0)*u.at(m,i-3) + ems.at(i-4,3);
            //u1.at(m,i) = -208.0*alpha*u1.at(m,i-1) + 228.0*alpha*u1.at(m,i-2) - 112.0*alpha*u1.at(m,i-3) + 22.0*alpha*u1.at(m,i-4) + (u1.at(m-1,i)+ht*f(i,m));
            //u1.at(m,i) /= -(70.0*alpha-1.0);
        }

        if (m==0)
        {
            printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, U(1,m), U(2,m), U(3,m), U(4,m));
            printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, x[0], x[1], x[2], x[3]);
            //printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, u1.at(m,N-4), u1.at(m,N-3), u1.at(m,N-2), u1.at(m,N-1));
            IPrinter::printVector(14,10, u.row(m));
            DoubleVector v(N+1); for (unsigned int n=0; n<=N; n++) v.at(n) = U(n,m);
            IPrinter::printVector(14,10, v);
            v.clear();
            IPrinter::printSeperatorLine();
            //break;
        }
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void BorderParabolicD::calculateN6L2RM(DoubleMatrix &u)
{
    unsigned int k = 6;
    double alpha = (ht*a*a)/(180.0*hx*hx);
    u.resize(M+1, N+1);

    double D[k+1][k+1] =
    {
        {+812.0, -3132.0, +5265.0, -5080.0, +2970.0, -972.0,  +137.0},
        {+137.0, -147.0,  -255.0,  +470.0,  -285.0,  +93.0,   -13.0},
        {-13.0,  +228.0,  -420.0,  +200.0,  +15.0,   -12.0,   +2.0},
        {+2.0,   -27.0,   +270.0,  -490.0,  +270.0,  -27.0,   +2.0},
        {+2.0,   -12.0,   +15.0,   +200.0,  -420.0,  +228.0,  -13.0},
        {-13.0,  +93.0,   -285.0,  +470.0,  -255.0,  -147.0,  +137.0},
        {+137.0, -972.0,  +2970.0, -5080.0, +5265.0, -3132.0, +812.0},
    };

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector x(k);
    DoubleMatrix ems(N-k, k);

    /* initial condition */
    for (unsigned int i=0; i<=N; i++) u.at(0,i) = initial(i);

    /* border conditions */
    for (unsigned int j=1; j<=M; j++)
    {
        u.at(j,0) = boundary(Left, j);
        u.at(j,N) = boundary(Right, j);
    }

    for (unsigned int m=1; m<=M; m++)
    {
        A.at(0,0) = D[1][1]*alpha - 1.0;
        A.at(0,1) = D[1][2]*alpha;
        A.at(0,2) = D[1][3]*alpha;
        A.at(0,3) = D[1][4]*alpha;
        A.at(0,4) = D[1][5]*alpha;
        A.at(0,5) = D[1][6]*alpha;
        b.at(0)   = -u.at(m-1,1) - (D[1][0]*alpha)*u.at(m,0) - ht*f(1,m);

        A.at(0,1) /= A.at(0,0);
        A.at(0,2) /= A.at(0,0);
        A.at(0,3) /= A.at(0,0);
        A.at(0,4) /= A.at(0,0);
        A.at(0,5) /= A.at(0,0);
        b.at(0)   /= A.at(0,0);
        A.at(0,0) = 1.0;

        ems.at(0,0) = A.at(0,1);
        ems.at(0,1) = A.at(0,2);
        ems.at(0,2) = A.at(0,3);
        ems.at(0,3) = A.at(0,4);
        ems.at(0,4) = A.at(0,5);
        ems.at(0,5) = b.at(0);

        // + * * * *
        for (unsigned int n=1; n<=N-(k+1); n++)
        {
            double g1 = D[0][0]*alpha - 1.0;
            double g2 = D[0][1]*alpha;
            double g3 = D[0][2]*alpha;
            double g4 = D[0][3]*alpha;
            double g5 = D[0][4]*alpha;
            double g6 = D[0][5]*alpha;
            double g7 = D[0][6]*alpha;
            double g0 = u.at(m-1,n) + ht*f(n,m);

            g2 /= -g1;
            g3 /= -g1;
            g4 /= -g1;
            g5 /= -g1;
            g6 /= -g1;
            g7 /= -g1;
            g0 /= -g1;
            g1 = 1.0;

            A.at(0,0) = A.at(0,1) + g2;
            A.at(0,1) = A.at(0,2) + g3;
            A.at(0,2) = A.at(0,3) + g4;
            A.at(0,3) = A.at(0,4) + g5;
            A.at(0,4) = A.at(0,5) + g6;
            A.at(0,5) = g7;
            b.at(0)   = b.at(0) - g0;

            A.at(0,1) /= A.at(0,0);
            A.at(0,2) /= A.at(0,0);
            A.at(0,3) /= A.at(0,0);
            A.at(0,4) /= A.at(0,0);
            A.at(0,5) /= A.at(0,0);
            b.at(0)   /= A.at(0,0);
            A.at(0,0) = 1.0;

            ems.at(n,0) = A.at(0,1);
            ems.at(n,1) = A.at(0,2);
            ems.at(n,2) = A.at(0,3);
            ems.at(n,3) = A.at(0,4);
            ems.at(n,4) = A.at(0,5);
            ems.at(n,5) = b.at(0);
        }

        A[1][0] = D[1][0]*alpha;
        A[1][1] = D[1][1]*alpha - 1.0;
        A[1][2] = D[1][2]*alpha;
        A[1][3] = D[1][3]*alpha;
        A[1][4] = D[1][4]*alpha;
        A[1][5] = D[1][5]*alpha;
        b[1]    = -u.at(m-1,N-5) - (D[1][6]*alpha)*u.at(m,N) - ht*f(N-5,m);

        A[2][0] = D[2][0]*alpha;
        A[2][1] = D[2][1]*alpha;
        A[2][2] = D[2][2]*alpha - 1.0;
        A[2][3] = D[2][3]*alpha;
        A[2][4] = D[2][4]*alpha;
        A[2][5] = D[2][5]*alpha;
        b[2]    = -u.at(m-1,N-4) - (D[2][6]*alpha)*u.at(m,N) - ht*f(N-4,m);

        A[3][0] = D[3][0]*alpha;
        A[3][1] = D[3][1]*alpha;
        A[3][2] = D[3][2]*alpha;
        A[3][3] = D[3][3]*alpha - 1.0;
        A[3][4] = D[3][4]*alpha;
        A[3][5] = D[3][5]*alpha;
        b[3]    = -u.at(m-1,N-3) - (D[3][6]*alpha)*u.at(m,N) - ht*f(N-3,m);

        A[4][0] = D[4][0]*alpha;
        A[4][1] = D[4][1]*alpha;
        A[4][2] = D[4][2]*alpha;
        A[4][3] = D[4][3]*alpha;
        A[4][4] = D[4][4]*alpha - 1.0;
        A[4][5] = D[4][5]*alpha;
        b[4]    = -u.at(m-1,N-2) - (D[4][6]*alpha)*u.at(m,N) - ht*f(N-2,m);

        A[5][0] = D[5][0]*alpha;
        A[5][1] = D[5][1]*alpha;
        A[5][2] = D[5][2]*alpha;
        A[5][3] = D[5][3]*alpha;
        A[5][4] = D[5][4]*alpha;
        A[5][5] = D[5][5]*alpha - 1.0;
        b[5]    = -u.at(m-1,N-1) - (D[5][6]*alpha)*u.at(m,N) - ht*f(N-1,m);

        GaussianElimination(A, b, x);

        u.at(m, N-1) = x.at(5);
        u.at(m, N-2) = x.at(4);
        u.at(m, N-3) = x.at(3);
        u.at(m, N-4) = x.at(2);
        u.at(m, N-5) = x.at(1);
        u.at(m, N-6) = x.at(0);
        for (unsigned int i=N-(k+1); i>=1; i--)
        {
            u.at(m,i) = -ems.at(i-1,0)*u.at(m,i+1) - ems.at(i-1,1)*u.at(m,i+2) - ems.at(i-1,2)*u.at(m,i+3) - ems.at(i-1,3)*u.at(m,i+4) - ems.at(i-1,4)*u.at(m,i+5) + ems.at(i-1,5);
            //u.at(m,i) = -208.0*alpha*u.at(m,i+1) + 228.0*alpha*u.at(m,i+2) - 112.0*alpha*u.at(m,i+3) + 22.0*alpha*u.at(m,i+4) + (u.at(m-1,i)+ht*f(i,m));
            //u.at(m,i) /= -(70.0*alpha-1.0);
        }

        if (m==0)
        {
            printf("%d %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n", m, U(N-6,m), U(N-5,m), U(N-4,m), U(N-3,m), U(N-2,m), U(N-1,m));
            printf("%d %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n", m, x[0], x[1], x[2], x[3], x[4], x[5]);
            //printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, u1.at(m,N-4), u1.at(m,N-3), u1.at(m,N-2), u1.at(m,N-1));
            IPrinter::printVector(14,10, u.row(m));
            DoubleVector v(N+1); for (unsigned int n=0; n<=N; n++) v.at(n) = U(n,m);
            IPrinter::printVector(14,10, v);
            v.clear();
            IPrinter::printSeperatorLine();
            break;
        }
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void BorderParabolicD::calculateN6R2LM(DoubleMatrix &u)
{
    unsigned int k = 6;
    double alpha = (ht*a*a)/(180.0*hx*hx);
    u.resize(M+1, N+1);

    double D[k+1][k+1] =
    {
        {+812.0, -3132.0, +5265.0, -5080.0, +2970.0, -972.0,  +137.0},
        {+137.0, -147.0,  -255.0,  +470.0,  -285.0,  +93.0,   -13.0},
        {-13.0,  +228.0,  -420.0,  +200.0,  +15.0,   -12.0,   +2.0},
        {+2.0,   -27.0,   +270.0,  -490.0,  +270.0,  -27.0,   +2.0},
        {+2.0,   -12.0,   +15.0,   +200.0,  -420.0,  +228.0,  -13.0},
        {-13.0,  +93.0,   -285.0,  +470.0,  -255.0,  -147.0,  +137.0},
        {+137.0, -972.0,  +2970.0, -5080.0, +5265.0, -3132.0, +812.0},
    };

    DoubleMatrix A(k, k, 0.0);
    DoubleVector b(k, 0.0);
    DoubleVector x(k);
    DoubleMatrix ems(N-k, k);

    /* initial condition */
    for (unsigned int i=0; i<=N; i++) u.at(0,i) = initial(i);

    /* border conditions */
    for (unsigned int j=1; j<=M; j++)
    {
        u.at(j,0) = boundary(Left, j);
        u.at(j,N) = boundary(Right, j);
    }

    for (unsigned int m=1; m<=M; m++)
    {
        A.at(0,0) = D[1][1]*alpha - 1.0;
        A.at(0,1) = D[1][2]*alpha;
        A.at(0,2) = D[1][3]*alpha;
        A.at(0,3) = D[1][4]*alpha;
        A.at(0,4) = D[1][5]*alpha;
        A.at(0,5) = D[1][6]*alpha;
        b.at(0)   = -u.at(m-1,1) - (D[1][0]*alpha)*u.at(m,0) - ht*f(1,m);

        A.at(1,0) = D[2][1]*alpha;
        A.at(1,1) = D[2][2]*alpha - 1.0;
        A.at(1,2) = D[2][3]*alpha;
        A.at(1,3) = D[2][4]*alpha;
        A.at(1,4) = D[2][5]*alpha;
        A.at(1,5) = D[2][6]*alpha;
        b.at(1)   = -u.at(m-1,2) - (D[2][0]*alpha)*u.at(m,0) - ht*f(2,m);

        A.at(2,0) = D[3][1]*alpha;
        A.at(2,1) = D[3][2]*alpha;
        A.at(2,2) = D[3][3]*alpha - 1.0;
        A.at(2,3) = D[3][4]*alpha;
        A.at(2,4) = D[3][5]*alpha;
        A.at(2,5) = D[3][6]*alpha;
        b.at(2)   = -u.at(m-1,3) - (D[3][0]*alpha)*u.at(m,0) - ht*f(3,m);

        A.at(3,0) = D[4][1]*alpha;
        A.at(3,1) = D[4][2]*alpha;
        A.at(3,2) = D[4][3]*alpha;
        A.at(3,3) = D[4][4]*alpha - 1.0;
        A.at(3,4) = D[4][5]*alpha;
        A.at(3,5) = D[4][6]*alpha;
        b.at(3)   = -u.at(m-1,4) - (D[4][0]*alpha)*u.at(m,0) - ht*f(4,m);

        A.at(4,0) = D[5][1]*alpha;
        A.at(4,1) = D[5][2]*alpha;
        A.at(4,2) = D[5][3]*alpha;
        A.at(4,3) = D[5][4]*alpha;
        A.at(4,4) = D[5][5]*alpha - 1.0;
        A.at(4,5) = D[5][6]*alpha;
        b.at(4)   = -u.at(m-1,5) - (D[5][0]*alpha)*u.at(m,0) - ht*f(5,m);

        A.at(5,0) = D[5][0]*alpha;
        A.at(5,1) = D[5][1]*alpha;
        A.at(5,2) = D[5][2]*alpha;
        A.at(5,3) = D[5][3]*alpha;
        A.at(5,4) = D[5][4]*alpha;
        A.at(5,5) = D[5][5]*alpha - 1.0;
        b.at(5)   = -u.at(m-1,N-1) - (D[5][6]*alpha)*u.at(m,N) - ht*f(N-1,m);

        A.at(5,0) /= A.at(5,5);
        A.at(5,1) /= A.at(5,5);
        A.at(5,2) /= A.at(5,5);
        A.at(5,3) /= A.at(5,5);
        A.at(5,4) /= A.at(5,5);
        b.at(5)   /= A.at(5,5);
        A.at(5,5) = 1.0;

        ems.at(N-7,0) = A.at(5,0);
        ems.at(N-7,1) = A.at(5,1);
        ems.at(N-7,2) = A.at(5,2);
        ems.at(N-7,3) = A.at(5,3);
        ems.at(N-7,4) = A.at(5,4);
        ems.at(N-7,5) = b.at(5);

        for (unsigned int n=N-1; n>=k+1; n--)
        {
            double g1 = D[6][0]*alpha;
            double g2 = D[6][1]*alpha;
            double g3 = D[6][2]*alpha;
            double g4 = D[6][3]*alpha;
            double g5 = D[6][4]*alpha;
            double g6 = D[6][5]*alpha;
            double g7 = D[6][6]*alpha - 1.0;
            double g0 = +u.at(m-1,n) + ht*f(n,m);

            g6 /= -g7;
            g5 /= -g7;
            g4 /= -g7;
            g3 /= -g7;
            g2 /= -g7;
            g1 /= -g7;
            g0 /= -g7;
            g7 = 1.0;

            A.at(5,5) = A.at(5,4) + g6;
            A.at(5,4) = A.at(5,3) + g5;
            A.at(5,3) = A.at(5,2) + g4;
            A.at(5,2) = A.at(5,1) + g3;
            A.at(5,1) = A.at(5,0) + g2;
            A.at(5,0) = g1;
            b.at(5)   = b.at(5) - g0;

            A.at(5,4) /= A.at(5,5);
            A.at(5,3) /= A.at(5,5);
            A.at(5,2) /= A.at(5,5);
            A.at(5,1) /= A.at(5,5);
            A.at(5,0) /= A.at(5,5);
            b.at(5)   /= A.at(5,5);
            A.at(5,5) = 1.0;

            ems.at(n-7,0) = A.at(5,0);
            ems.at(n-7,1) = A.at(5,1);
            ems.at(n-7,2) = A.at(5,2);
            ems.at(n-7,3) = A.at(5,3);
            ems.at(n-7,4) = A.at(5,4);
            ems.at(n-7,5) = b.at(5);
        }

        GaussianElimination(A, b, x);

        u.at(m, 1) = x.at(0);
        u.at(m, 2) = x.at(1);
        u.at(m, 3) = x.at(2);
        u.at(m, 4) = x.at(3);
        u.at(m, 5) = x.at(4);
        u.at(m, 6) = x.at(5);
        for (unsigned int i=k+1; i<=N-1; i++)
        {
            u.at(m,i) = - ems.at(i-6,4)*u.at(m,i-1) - ems.at(i-6,3)*u.at(m,i-2) - ems.at(i-6,2)*u.at(m,i-3) - ems.at(i-6,1)*u.at(m,i-4) - ems.at(i-6,0)*u.at(m,i-5) + ems.at(i-6,5);
            //u1.at(m,i) = -208.0*alpha*u1.at(m,i-1) + 228.0*alpha*u1.at(m,i-2) - 112.0*alpha*u1.at(m,i-3) + 22.0*alpha*u1.at(m,i-4) + (u1.at(m-1,i)+ht*f(i,m));
            //u1.at(m,i) /= -(70.0*alpha-1.0);
        }

        if (m==0)
        {
            printf("%d %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n", m, U(1,m), U(2,m), U(3,m), U(4,m), U(5,m), U(6,m));
            printf("%d %18.10f %18.10f %18.10f %18.10f %18.10f %18.10f\n", m, x[0], x[1], x[2], x[3], x[4], x[5]);
            //printf("%d %18.10f %18.10f %18.10f %18.10f\n", m, u1.at(m,N-4), u1.at(m,N-3), u1.at(m,N-2), u1.at(m,N-1));
            IPrinter::printVector(14,10, u.row(m));
            DoubleVector v(N+1); for (unsigned int n=0; n<=N; n++) v.at(n) = U(n,m);
            IPrinter::printVector(14,10, v);
            v.clear();
            IPrinter::printSeperatorLine();
            //break;
        }
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}
