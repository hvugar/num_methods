#include "parabolicequation.h"
#include "cmethods.h"
#include <rungekutta.h>

double MIN1 = +10000.0;
double MAX1 = -10000.0;

void saveData1(const DoubleMatrix& m, int i, unsigned int N2, unsigned int N1)
{
    double min = m.min();
    double max = m.max();
    if (MIN1 > min) MIN1 = min;
    if (MAX1 < max) MAX1 = max;

    char buffer[20];
    int n = 0;
    if (i<10) n = sprintf(buffer, "data/0000000%d.txt", i);
    if (i<100 && i>=10) n = sprintf(buffer, "data/000000%d.txt", i);
    if (i<1000 && i>=100) n = sprintf(buffer, "data/00000%d.txt", i);
    if (i<10000 && i>=1000) n = sprintf(buffer, "data/0000%d.txt", i);
    if (i<100000 && i>=10000) n = sprintf(buffer, "data/000%d.txt", i);
    buffer[n] = '\0';
    FILE *file = fopen(buffer, "w");
    IPrinter::printMatrix(m, N2, N1, NULL, file);
    fclose(file);

    printf("File: %s min: %.10f max: %.10f min: %.10f max: %.10f\n", buffer, MIN1, MAX1, min, max);
}

void IParabolicEquation::calculateU(DoubleVector &u, double hx, double ht, unsigned int N, unsigned int M, double a) const
{
    u.clear();
    u.resize(N+1);

    DoubleVector da(N-1);
    DoubleVector db(N-1);
    DoubleVector dc(N-1);
    DoubleVector dd(N-1);
    DoubleVector rx(N-1);

    double alpha = -(a*a*ht)/(hx*hx);
    double beta  = 1.0 + (2.0*a*a*ht)/(hx*hx);

    for (unsigned int j=0; j<=M; j++)
    {
        if (j == 0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u[i] = initial(i);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = alpha;
                db[i-1] = beta;
                dc[i-1] = alpha;
                dd[i-1] = u[i] + ht * f(i, j);
            }

            da[0]   = 0.0;
            dc[N-2] = 0.0;

            u[0] = boundary(Left, j);//m1(j);
            u[N] = boundary(Right, j);//m2(j);

            dd[0]   -= alpha * u[0];
            dd[N-2] -= alpha * u[N];

            tomasAlgorithm(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.size());

            for (unsigned int i=1; i<=N-1; i++)
            {
                u[i] = rx[i-1];
            }
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();
}

void IParabolicEquation::calculateU(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a) const
{
    u.clear();
    u.resize(M+1, N+1);

    DoubleVector da(N-1);
    DoubleVector db(N-1);
    DoubleVector dc(N-1);
    DoubleVector dd(N-1);
    DoubleVector rx(N-1);

    double alpha = -(a*a*ht)/(hx*hx);
    double beta  = 1.0 + (2.0*a*a*ht)/(hx*hx);

    for (unsigned int j=0; j<=M; j++)
    {
        if (j == 0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u[j][i] = initial(i);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = alpha;
                db[i-1] = beta;
                dc[i-1] = alpha;
                dd[i-1] = u[j-1][i] + ht * f(i, j);
            }

            da[0]   = 0.0;
            dc[N-2] = 0.0;

            u[j][0] = boundary(Left, j);//m1(j);
            u[j][N] = boundary(Right, j);//m2(j);

            dd[0]   -= alpha * u[j][0];
            dd[N-2] -= alpha * u[j][N];

            tomasAlgorithm(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.size());

            for (unsigned int i=1; i<=N-1; i++)
            {
                u[j][i] = rx[i-1];
            }
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();
}

void IParabolicEquation::calculateU1(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a) const
{
    u.clear();
    u.resize(M+1, N+1);

    DoubleVector da(N-1);
    DoubleVector db(N-1);
    DoubleVector dc(N-1);
    DoubleVector dd(N-1);
    DoubleVector rx(N-1);

    double alpha = 1.0;
    double beta  = (1.0 + (2.0*a*a*ht)/(hx*hx))/(-(a*a*ht)/(hx*hx));

    for (unsigned int j=0; j<=M; j++)
    {
        if (j == 0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u[j][i] = initial(i);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = alpha;
                db[i-1] = beta;
                dc[i-1] = alpha;
                dd[i-1] = (u.at(j-1,i) + ht * f(i, j))/(-(a*a*ht)/(hx*hx));
            }

            da[0]   = 0.0;
            dc[N-2] = 0.0;

            u[j][0] = boundary(Left, j);
            u[j][N] = boundary(Right, j);

            dd[0]   -= alpha * u.at(j,0);
            dd[N-2] -= alpha * u.at(j,N);

            tomasAlgorithm(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.size());

            for (unsigned int i=1; i<=N-1; i++)
            {
                u.at(j,i) = rx.at(i-1);
            }
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();
}

void IParabolicEquation::calculateN(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a) const
{
    u.clear();
    u.resize(M+1, N+1);

    DoubleVector da(N-1);
    DoubleVector db(N-1);
    DoubleVector dc(N-1);
    DoubleVector dd(N-1);
    DoubleVector rx(N-1);

    double alpha = -(a*a*ht)/(hx*hx);
    double beta  = 1.0 + (2.0*a*a*ht)/(hx*hx);

    for (unsigned int j=0; j<=M; j++)
    {
        if (j == 0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u[j][i] = initial(i);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = alpha;
                db[i-1] = beta;
                dc[i-1] = alpha;
                dd[i-1] = u[j-1][i] + ht * f(i, j);
            }

            da[0]   = 0.0;
            dc[N-2] = 0.0;

            db[0]   = alpha+beta;
            db[N-2] = alpha+beta;

            dd[0]   += alpha * hx * boundary(Left, j);
            dd[N-2] -= alpha * hx * boundary(Right, j);

            tomasAlgorithm(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.size());

            for (unsigned int i=1; i<=N-1; i++)
            {
                u[j][i] = rx[i-1];
            }

            u[j][0] = u[j][1]   - hx * boundary(Left, j);
            u[j][N] = u[j][N-1] + hx * boundary(Right, j);
        }
    }
    IPrinter::printSeperatorLine();

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();
}

void IParabolicEquation::calculateN4L2RD(DoubleMatrix &u, double hx, double ht, double N, double M, double a)
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
        A[0][0] = D[1][1]*alpha - 1.0;
        A[0][1] = D[1][2]*alpha;
        A[0][2] = D[1][3]*alpha;
        A[0][3] = D[1][4]*alpha;
        b[0]    = -u.at(m-1,1) - (D[1][0]*alpha)*u.at(m,0) - ht*f(1,m);

        A[0][1] /= A[0][0];
        A[0][2] /= A[0][0];
        A[0][3] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;

        ems.at(0,0) = A[0][1];
        ems.at(0,1) = A[0][2];
        ems.at(0,2) = A[0][3];
        ems.at(0,3) = b[0];

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

            A[0][0] = A[0][1] + g2;
            A[0][1] = A[0][2] + g3;
            A[0][2] = A[0][3] + g4;
            A[0][3] = g5;
            b[0]    = b[0] - g0;
            \
            A[0][1] /= A[0][0];
            A[0][2] /= A[0][0];
            A[0][3] /= A[0][0];
            b[0]    /= A[0][0];
            A[0][0] = 1.0;

            ems.at(n,0) = A[0][1];
            ems.at(n,1) = A[0][2];
            ems.at(n,2) = A[0][3];
            ems.at(n,3) = b[0];
        }

        A[1][0] = D[1][0]*alpha;
        A[1][1] = D[1][1]*alpha - 1.0;
        A[1][2] = D[1][2]*alpha;
        A[1][3] = D[1][3]*alpha;
        b[1]    = -u.at(m-1,N-3) - (D[1][4]*alpha)*u.at(m,N) - ht*f(N-3,m);

        A[2][0] = D[2][0]*alpha;
        A[2][1] = D[2][1]*alpha;
        A[2][2] = D[2][2]*alpha - 1.0;
        A[2][3] = D[2][3]*alpha;
        b[2]    = -u.at(m-1,N-2) - (D[2][4]*alpha)*u.at(m,N) - ht*f(N-2,m);

        A[3][0] = D[3][0]*alpha;
        A[3][1] = D[3][1]*alpha;
        A[3][2] = D[3][2]*alpha;
        A[3][3] = D[3][3]*alpha - 1.0;
        b[3]    = -u.at(m-1,N-1) - (D[3][4]*alpha)*u.at(m,N) - ht*f(N-1,m);

        GaussianElimination(A, b, x);

        u.at(m, N-1) = x.at(3);
        u.at(m, N-2) = x.at(2);
        u.at(m, N-3) = x.at(1);
        u.at(m, N-4) = x.at(0);
        for (unsigned int i=N-(k+1); i>=1; i--)
        {
            u.at(m,i) = -ems.at(i-1,0)*u.at(m,i+1) - ems.at(i-1,1)*u.at(m,i+2) - ems.at(i-1,2)*u.at(m,i+3) + ems.at(i-1,3);
        }
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void IParabolicEquation::calculateN4R2LD(DoubleMatrix &u, double hx, double ht, double N, double M, double a)
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
        A[0][0] = D[1][1]*alpha - 1.0;
        A[0][1] = D[1][2]*alpha;
        A[0][2] = D[1][3]*alpha;
        A[0][3] = D[1][4]*alpha;
        b[0]    = -u.at(m-1,1) - (D[1][0]*alpha)*u.at(m,0) - ht*f(1,m);

        A[1][0] = D[2][1]*alpha;
        A[1][1] = D[2][2]*alpha - 1.0;
        A[1][2] = D[2][3]*alpha;
        A[1][3] = D[2][4]*alpha;
        b[1]    = -u.at(m-1,2) - (D[2][0]*alpha)*u.at(m,0) - ht*f(2,m);

        A[2][0] = D[3][1]*alpha;
        A[2][1] = D[3][2]*alpha;
        A[2][2] = D[3][3]*alpha - 1.0;
        A[2][3] = D[3][4]*alpha;
        b[2]    = -u.at(m-1,3) - (D[3][0]*alpha)*u.at(m,0) - ht*f(3,m);

        A[3][0] = D[3][0]*alpha;
        A[3][1] = D[3][1]*alpha;
        A[3][2] = D[3][2]*alpha;
        A[3][3] = D[3][3]*alpha - 1.0;
        b[3]    = -u.at(m-1,N-1) - (D[3][4]*alpha)*u.at(m,N) - ht*f(N-1,m);

        A[3][0] /= A[3][3];
        A[3][1] /= A[3][3];
        A[3][2] /= A[3][3];
        b[3]    /= A[3][3];
        A[3][3] = 1.0;

        ems.at(N-5,0) = A[3][0];
        ems.at(N-5,1) = A[3][1];
        ems.at(N-5,2) = A[3][2];
        ems.at(N-5,3) = b[3];

        for (unsigned int n=N-1; n>=k+1; n--)
        {
            double g1 = D[4][0]*alpha;
            double g2 = D[4][1]*alpha;
            double g3 = D[4][2]*alpha;
            double g4 = D[4][3]*alpha;
            double g5 = D[4][4]*alpha - 1.0;
            double g0 = u.at(m-1,n) + ht*f(n,m);

            g4 /= -g5;
            g3 /= -g5;
            g2 /= -g5;
            g1 /= -g5;
            g0 /= -g5;
            g5 = 1.0;

            A[3][3] = A[3][2] + g4;
            A[3][2] = A[3][1] + g3;
            A[3][1] = A[3][0] + g2;
            A[3][0] = g1;
            b[3]    = b[3] - g0;

            A[3][2] /= A[3][3];
            A[3][1] /= A[3][3];
            A[3][0] /= A[3][3];
            b[3]    /= A[3][3];
            A[3][3] = 1.0;

            ems.at(n-5,0) = A[3][0];
            ems.at(n-5,1) = A[3][1];
            ems.at(n-5,2) = A[3][2];
            ems.at(n-5,3) = b[3];
        }

        GaussianElimination(A, b, x);

        u.at(m, 1) = x.at(0);
        u.at(m, 2) = x.at(1);
        u.at(m, 3) = x.at(2);
        u.at(m, 4) = x.at(3);
        for (unsigned int i=k+1; i<=N-1; i++)
        {
            u.at(m,i) = -ems.at(i-4,2)*u.at(m,i-1) - ems.at(i-4,1)*u.at(m,i-2) - ems.at(i-4,0)*u.at(m,i-3) + ems.at(i-4,3);
        }
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void IParabolicEquation::calculateN6L2RD(DoubleMatrix &u, double hx, double ht, double N, double M, double a)
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
        A[0][0] = D[1][1]*alpha - 1.0;
        A[0][1] = D[1][2]*alpha;
        A[0][2] = D[1][3]*alpha;
        A[0][3] = D[1][4]*alpha;
        A[0][4] = D[1][5]*alpha;
        A[0][5] = D[1][6]*alpha;
        b[0]    = -u.at(m-1,1) - (D[1][0]*alpha)*u.at(m,0) - ht*f(1,m);

        A[0][1] /= A[0][0];
        A[0][2] /= A[0][0];
        A[0][3] /= A[0][0];
        A[0][4] /= A[0][0];
        A[0][5] /= A[0][0];
        b[0]    /= A[0][0];
        A[0][0] = 1.0;

        ems.at(0,0) = A[0][1];
        ems.at(0,1) = A[0][2];
        ems.at(0,2) = A[0][3];
        ems.at(0,3) = A[0][4];
        ems.at(0,4) = A[0][5];
        ems.at(0,5) = b[0];

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

            A[0][0] = A[0][1] + g2;
            A[0][1] = A[0][2] + g3;
            A[0][2] = A[0][3] + g4;
            A[0][3] = A[0][4] + g5;
            A[0][4] = A[0][5] + g6;
            A[0][5] = g7;
            b[0]   = b[0] - g0;

            A[0][1] /= A[0][0];
            A[0][2] /= A[0][0];
            A[0][3] /= A[0][0];
            A[0][4] /= A[0][0];
            A[0][5] /= A[0][0];
            b[0]    /= A[0][0];
            A[0][0] = 1.0;

            ems.at(n,0) = A[0][1];
            ems.at(n,1) = A[0][2];
            ems.at(n,2) = A[0][3];
            ems.at(n,3) = A[0][4];
            ems.at(n,4) = A[0][5];
            ems.at(n,5) = b[0];
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
        }
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void IParabolicEquation::calculateN6R2LD(DoubleMatrix &u, double hx, double ht, double N, double M, double a)
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
        A[0][0] = D[1][1]*alpha - 1.0;
        A[0][1] = D[1][2]*alpha;
        A[0][2] = D[1][3]*alpha;
        A[0][3] = D[1][4]*alpha;
        A[0][4] = D[1][5]*alpha;
        A[0][5] = D[1][6]*alpha;
        b[0]    = -u.at(m-1,1) - (D[1][0]*alpha)*u.at(m,0) - ht*f(1,m);

        A[1][0] = D[2][1]*alpha;
        A[1][1] = D[2][2]*alpha - 1.0;
        A[1][2] = D[2][3]*alpha;
        A[1][3] = D[2][4]*alpha;
        A[1][4] = D[2][5]*alpha;
        A[1][5] = D[2][6]*alpha;
        b[1]    = -u.at(m-1,2) - (D[2][0]*alpha)*u.at(m,0) - ht*f(2,m);

        A[2][0] = D[3][1]*alpha;
        A[2][1] = D[3][2]*alpha;
        A[2][2] = D[3][3]*alpha - 1.0;
        A[2][3] = D[3][4]*alpha;
        A[2][4] = D[3][5]*alpha;
        A[2][5] = D[3][6]*alpha;
        b[2]    = -u.at(m-1,3) - (D[3][0]*alpha)*u.at(m,0) - ht*f(3,m);

        A[3][0] = D[4][1]*alpha;
        A[3][1] = D[4][2]*alpha;
        A[3][2] = D[4][3]*alpha;
        A[3][3] = D[4][4]*alpha - 1.0;
        A[3][4] = D[4][5]*alpha;
        A[3][5] = D[4][6]*alpha;
        b[3]    = -u.at(m-1,4) - (D[4][0]*alpha)*u.at(m,0) - ht*f(4,m);

        A[4][0] = D[5][1]*alpha;
        A[4][1] = D[5][2]*alpha;
        A[4][2] = D[5][3]*alpha;
        A[4][3] = D[5][4]*alpha;
        A[4][4] = D[5][5]*alpha - 1.0;
        A[4][5] = D[5][6]*alpha;
        b[4]    = -u.at(m-1,5) - (D[5][0]*alpha)*u.at(m,0) - ht*f(5,m);

        A[5][0] = D[5][0]*alpha;
        A[5][1] = D[5][1]*alpha;
        A[5][2] = D[5][2]*alpha;
        A[5][3] = D[5][3]*alpha;
        A[5][4] = D[5][4]*alpha;
        A[5][5] = D[5][5]*alpha - 1.0;
        b[5]    = -u.at(m-1,N-1) - (D[5][6]*alpha)*u.at(m,N) - ht*f(N-1,m);

        A[5][0] /= A[5][5];
        A[5][1] /= A[5][5];
        A[5][2] /= A[5][5];
        A[5][3] /= A[5][5];
        A[5][4] /= A[5][5];
        b[5]    /= A[5][5];
        A[5][5] = 1.0;

        ems.at(N-7,0) = A[5][0];
        ems.at(N-7,1) = A[5][1];
        ems.at(N-7,2) = A[5][2];
        ems.at(N-7,3) = A[5][3];
        ems.at(N-7,4) = A[5][4];
        ems.at(N-7,5) = b[5];

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

            A[5][5] = A[5][4] + g6;
            A[5][4] = A[5][3] + g5;
            A[5][3] = A[5][2] + g4;
            A[5][2] = A[5][1] + g3;
            A[5][1] = A[5][0] + g2;
            A[5][0] = g1;
            b[5]    = b[5] - g0;

            A[5][4] /= A[5][5];
            A[5][3] /= A[5][5];
            A[5][2] /= A[5][5];
            A[5][1] /= A[5][5];
            A[5][0] /= A[5][5];
            b[5]    /= A[5][5];
            A[5][5] = 1.0;

            ems.at(n-7,0) = A[5][0];
            ems.at(n-7,1) = A[5][1];
            ems.at(n-7,2) = A[5][2];
            ems.at(n-7,3) = A[5][3];
            ems.at(n-7,4) = A[5][4];
            ems.at(n-7,5) = b[5];
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
        }
    }

    ems.clear();
    x.clear();
    b.clear();
    A.clear();
}

void IBackwardParabolicEquation::calculateU(DoubleMatrix &psi, double hx, double ht, unsigned int N, unsigned int M, double a) const
{
    psi.clear();
    psi.resize(M+1, N+1);

    DoubleVector da(N-1);
    DoubleVector db(N-1);
    DoubleVector dc(N-1);
    DoubleVector dd(N-1);
    DoubleVector rx(N-1);

    double alpha = -(a*ht)/(hx*hx);
    double beta  = 1.0 + (2.0*a*ht)/(hx*hx);

    for (unsigned int k=0; k<=M; k++)
    {
        unsigned int j = M-k;

        if (j == M)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                psi[j][i] = binitial(i);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = alpha;
                db[i-1] = beta;
                dc[i-1] = alpha;
                dd[i-1] = psi[j+1][i] - ht * bf(i, j);
            }

            da[0]   = 0.0;
            dc[N-2] = 0.0;

            psi[j][0] = bboundary(Left, j);//bm1(j);
            psi[j][N] = bboundary(Right, j);//bm2(j);

            dd[0]   -= alpha * psi[j][0];
            dd[N-2] -= alpha * psi[j][N];

            tomasAlgorithm(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.size());

            for (unsigned int i=1; i<=N-1; i++)
            {
                psi[j][i] = rx[i-1];
            }
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();
}

void IParabolicEquation2D::calculateMVD(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    //cleaning matrix
    u.resize(N2+1, N1+1);

    DoubleMatrix uh(N2+1, N1+1);

    DoubleVector da1(N1-1);
    DoubleVector db1(N1-1);
    DoubleVector dc1(N1-1);
    DoubleVector dd1(N1-1);
    DoubleVector rx1(N1-1);

    DoubleVector da2(N2-1);
    DoubleVector db2(N2-1);
    DoubleVector dc2(N2-1);
    DoubleVector dd2(N2-1);
    DoubleVector rx2(N2-1);

    double x1_a = -(a1*a1*ht)/(2.0*h1*h1);
    double x1_b  = 1.0 + (a1*a1*ht)/(h1*h1);
    double x1_c = (a2*a2*ht)/(2.0*h2*h2);

    double x2_a = -(a2*a2*ht)/(2.0*h2*h2);
    double x2_b  = 1.0 + (a2*a2*ht)/(h2*h2);
    double x2_c = (a1*a1*ht)/(2.0*h1*h1);

    //    double lamda1 = ((a1*a1)*ht)/(h1*h1);
    //    double lamda2 = ((a2*a2)*ht)/(h2*h2);

    //    double x1_a = -0.5 * lamda1;
    //    double x1_b = +1.0 + lamda1;
    //    double x1_c = +0.5 * lamda2;

    //    double x2_a = -0.5 * lamda2;
    //    double x2_b = +1.0 + lamda2;
    //    double x2_c = +0.5 * lamda1;

    for (unsigned int k=0; k<=M; k++)
    {
        if (k==0)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u[j][i] = initial(i, j);
                }
            }
        }
        else
        {
            // Approximation to x1 direction
            for (unsigned int j=1; j<N2; j++)
            {
                for (unsigned int i=1; i<N1; i++)
                {
                    da1[i-1] = x1_a;
                    db1[i-1] = x1_b;
                    dc1[i-1] = x1_a;
                    dd1[i-1] = x1_c*(u[j-1][i] - 2.0*u[j][i] + u[j+1][i]) + u[j][i] + (ht/2.0) * f(i, j, 2*k-1);
                    //dd1[i-1] = x1_c*u[j-1][i] + x1_d*u[j][i] + x1_c*u[j+1][i] + (ht/2.0) * f(i, j, 2*k-1);
                }

                da1[0]     = 0.0;
                dc1[N1-2]  = 0.0;

                uh[j][0]  = boundary(0, j, 2*k-1);
                uh[j][N1] = boundary(N1, j, 2*k-1);

                dd1[0]    -= x1_a * uh[j][0];
                dd1[N1-2] -= x1_a * uh[j][N1];

                tomasAlgorithm(da1.data(), db1.data(), dc1.data(), dd1.data(), rx1.data(), rx1.size());

                for (unsigned int i=1; i<N1; i++)
                {
                    uh[j][i] = rx1[i-1];
                }
            }

            for (unsigned int i=0; i<=N1; i++)
            {
                uh[0][i]  = boundary(i, 0, 2*k-1);
                uh[N2][i] = boundary(i, N2, 2*k-1);
            }

            // Approximation to x2 direction
            for (unsigned int i=1; i<N1; i++)
            {
                for (unsigned int j=1; j<N2; j++)
                {
                    da2[j-1] = x2_a;
                    db2[j-1] = x2_b;
                    dc2[j-1] = x2_a;
                    dd2[j-1] = x2_c*(uh[j][i-1] - 2.0*uh[j][i] + uh[j][i+1]) + uh[j][i] + (ht/2.0) * f(i, j, 2*k);
                    //dd2[j-1] = x2_c*uh[j][i-1] + x2_d*uh[j][i] + x2_c*uh[j][i+1] + (ht/2.0) * f(i, j, 2*k);
                }
                da2[0]     = 0.0;
                dc2[N2-2]  = 0.0;

                u[0][i]  = boundary(i, 0, 2*k);
                u[N2][i] = boundary(i, N2, 2*k);

                dd2[0]    -= x2_a * u[0][i];
                dd2[N2-2] -= x2_a * u[N2][i];

                tomasAlgorithm(da2.data(), db2.data(), dc2.data(), dd2.data(), rx2.data(), rx2.size());

                for (unsigned int j=1; j<N2; j++)
                {
                    u[j][i] = rx2[j-1];
                }
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                u[j][0]  = boundary(0, j, 2*k);
                u[j][N1] = boundary(N1, j, 2*k);
            }
        }

        //if (k%100==0)
        //saveData1(u, k, N2, N1);
    }

    da1.clear();
    db1.clear();
    dc1.clear();
    dd1.clear();
    rx1.clear();

    da2.clear();
    db2.clear();
    dc2.clear();
    dd2.clear();
    rx2.clear();
}

void IParabolicEquation2D::calculateMVD1(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    u.clear();
    u.resize(N2+1, N1+1);

    DoubleMatrix uh(N2+1, N1+1);

    DoubleVector da1(N1-1);
    DoubleVector db1(N1-1);
    DoubleVector dc1(N1-1);
    DoubleVector dd1(N1-1);
    DoubleVector rx1(N1-1);

    DoubleVector da2(N2-1);
    DoubleVector db2(N2-1);
    DoubleVector dc2(N2-1);
    DoubleVector dd2(N2-1);
    DoubleVector rx2(N2-1);

    double x1_a = -(a1*a1*ht)/(h1*h1);
    double x1_b  = 1.0 + 2.0*(a1*a1*ht)/(h1*h1);
    double x1_c = (a2*a2*ht)/(h2*h2);
    double x2_a = -(a2*a2*ht)/(h2*h2);
    double x2_b  = 1.0 + 2.0*(a2*a2*ht)/(h2*h2);
    double x2_c = (a1*a1*ht)/(h1*h1);

    for (unsigned int k=0; k<=M; k++)
    {
        if (k==0)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u[j][i] = initial(i, j);
                }
            }
        }
        else
        {
            // Approximation to x1 direction
            if (k%2==1)
            {
                for (unsigned int j=1; j<N2; j++)
                {
                    for (unsigned int i=1; i<N1; i++)
                    {
                        da1[i-1] = x1_a;
                        db1[i-1] = x1_b;
                        dc1[i-1] = x1_a;
                        dd1[i-1] = x1_c*(u[j-1][i] - 2.0*u[j][i] + u[j+1][i]) + u[j][i] + ht * f(i, j, k);
                    }

                    da1[0]     = 0.0;
                    dc1[N1-2]  = 0.0;

                    uh[j][0]  = boundary(0, j, k);
                    uh[j][N1] = boundary(N1, j, k);

                    dd1[0]    -= x1_a * uh[j][0];
                    dd1[N1-2] -= x1_a * uh[j][N1];

                    tomasAlgorithm(da1.data(), db1.data(), dc1.data(), dd1.data(), rx1.data(), rx1.size());

                    for (unsigned int i=1; i<N1; i++)
                    {
                        uh[j][i] = rx1[i-1];
                    }
                }

                for (unsigned int i=0; i<=N1; i++)
                {
                    uh[0][i]  = boundary(i, 0, k);
                    uh[N2][i] = boundary(i, N2, k);
                }
            }
            else
            {
                // Approximation to x2 direction
                for (unsigned int i=1; i<N1; i++)
                {
                    for (unsigned int j=1; j<N2; j++)
                    {
                        da2[j-1] = x2_a;
                        db2[j-1] = x2_b;
                        dc2[j-1] = x2_a;
                        dd2[j-1] = x2_c*(uh[j][i-1] - 2.0*uh[j][i] + uh[j][i+1]) + uh[j][i] + ht * f(i, j, k);
                    }
                    da2[0]     = 0.0;
                    dc2[N2-2]  = 0.0;

                    u[0][i]  = boundary(i, 0, k);
                    u[N2][i] = boundary(i, N2, k);

                    dd2[0]    -= x2_a * u[0][i];
                    dd2[N2-2] -= x2_a * u[N2][i];

                    tomasAlgorithm(da2.data(), db2.data(), dc2.data(), dd2.data(), rx2.data(), rx2.size());

                    for (unsigned int j=1; j<N2; j++)
                    {
                        u[j][i] = rx2[j-1];
                    }
                }

                for (unsigned int j=0; j<=N2; j++)
                {
                    u[j][0]  = boundary(0, j, k);
                    u[j][N1] = boundary(N1, j, k);
                }
            }
        }
    }

    da1.clear();
    db1.clear();
    dc1.clear();
    dd1.clear();
    rx1.clear();

    da2.clear();
    db2.clear();
    dc2.clear();
    dd2.clear();
    rx2.clear();
}

void IParabolicEquation2D::calculateMFS(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    u.clear();
    u.resize(N2+1, N1+1);

    DoubleMatrix uh(N2+1, N1+1);

    DoubleVector da1(N1-1);
    DoubleVector db1(N1-1);
    DoubleVector dc1(N1-1);
    DoubleVector dd1(N1-1);
    DoubleVector rx1(N1-1);

    DoubleVector da2(N2-1);
    DoubleVector db2(N2-1);
    DoubleVector dc2(N2-1);
    DoubleVector dd2(N2-1);
    DoubleVector rx2(N2-1);

    double lamda1 = (a1*a1)*(ht/(h1*h1));
    double x1_a = -lamda1/2.0;
    double x1_b  = 1.0 + lamda1;

    double lamda2 = (a2*a2)*(ht/(h2*h2));
    double x2_a = -lamda2/2.0;
    double x2_b  = 1.0 + lamda2;

    for (unsigned int k=0; k<=M; k++)
    {
        if (k==0)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u[j][i] = initial(i, j);
                }
            }
        }
        else
        {
            // Approximation to x1 direction
            for (unsigned int j=1; j<N2; j++)
            {
                for (unsigned int i=1; i<N1; i++)
                {
                    da1[i-1] = x1_a;
                    db1[i-1] = x1_b;
                    dc1[i-1] = x1_a;
                    dd1[i-1] = u[j][i] + (ht/2.0) * f(i, j, 2*k-1);
                }

                da1[0]     = 0.0;
                dc1[N1-2]  = 0.0;

                uh[j][0]  = boundary(0, j, 2*k-1);
                uh[j][N1] = boundary(N1, j, 2*k-1);

                dd1[0]    -= x1_a * uh[j][0];
                dd1[N1-2] -= x1_a * uh[j][N1];

                tomasAlgorithm(da1.data(), db1.data(), dc1.data(), dd1.data(), rx1.data(), rx1.size());

                for (unsigned int i=1; i<N1; i++)
                {
                    uh[j][i] = rx1[i-1];
                }
            }

            for (unsigned int i=0; i<=N1; i++)
            {
                uh[0][i]  = boundary(i, 0, 2*k-1);
                uh[N2][i] = boundary(i, N2, 2*k-1);
            }

            // Approximation to x2 direction
            for (unsigned int i=1; i<N1; i++)
            {
                for (unsigned int j=1; j<N2; j++)
                {
                    da2[j-1] = x2_a;
                    db2[j-1] = x2_b;
                    dc2[j-1] = x2_a;
                    dd2[j-1] = uh[j][i] + (ht/2.0) * f(i, j, 2*k);
                }
                da2[0]     = 0.0;
                dc2[N2-2]  = 0.0;

                u[0][i]  = boundary(i, 0, 2*k);
                u[N2][i] = boundary(i, N2, 2*k);

                dd2[0]    -= x2_a * u[0][i];
                dd2[N2-2] -= x2_a * u[N2][i];

                tomasAlgorithm(da2.data(), db2.data(), dc2.data(), dd2.data(), rx2.data(), rx2.size());

                for (unsigned int j=1; j<N2; j++)
                {
                    u[j][i] = rx2[j-1];
                }
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                u[j][0]  = boundary(0, j, 2*k);
                u[j][N1] = boundary(N1, j, 2*k);
            }
        }

        saveData1(u, k, N2, N1);
    }

    da1.clear();
    db1.clear();
    dc1.clear();
    dd1.clear();
    rx1.clear();

    da2.clear();
    db2.clear();
    dc2.clear();
    dd2.clear();
    rx2.clear();
}

void IParabolicEquation2D::calculateMVD(DoubleCube &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    u.clear();
    u.resize(M+1, N2+1, N1+1);

    DoubleMatrix u0(N2+1, N1+1);
    DoubleMatrix u1(N2+1, N1+1);

    DoubleVector da1(N1-1);
    DoubleVector db1(N1-1);
    DoubleVector dc1(N1-1);
    DoubleVector dd1(N1-1);
    DoubleVector rx1(N1-1);

    DoubleVector da2(N2-1);
    DoubleVector db2(N2-1);
    DoubleVector dc2(N2-1);
    DoubleVector dd2(N2-1);
    DoubleVector rx2(N2-1);

    double x1_a = -(a1*a1*ht)/(2.0*h1*h1);
    double x1_b  = 1.0 + (a1*a1*ht)/(h1*h1);
    double x1_c = (a2*a2*ht)/(2.0*h2*h2);
    //double x1_d = 1.0 - (a2*a2*ht)/(h2*h2);

    double x2_a = -(a2*a2*ht)/(2.0*h2*h2);
    double x2_b  = 1.0 + (a2*a2*ht)/(h2*h2);
    double x2_c = (a1*a1*ht)/(2.0*h1*h1);
    //double x2_d = 1.0 - (a1*a1*ht)/(h1*h1);

    for (unsigned int k=0; k<=M; k++)
    {
        if (k==0)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u0[j][i] = initial(i, j);
                }
                u[k] = u0;
            }
        }
        else
        {
            // Approximation to x1 direction
            for (unsigned int j=1; j<N2; j++)
            {
                for (unsigned int i=1; i<N1; i++)
                {
                    da1[i-1] = x1_a;
                    db1[i-1] = x1_b;
                    dc1[i-1] = x1_a;
                    dd1[i-1] = x1_c*(u0[j-1][i] - 2.0*u0[j][i] + u0[j+1][i]) + u0[j][i] + (ht/2.0) * f(i, j, 2*k-1);
                    //dd1[i-1] = x1_c*u0[j-1][i] + x1_d*u0[j][i] + x1_c*u0[j+1][i] + (ht/2.0) * f(i, j, 2*k-1);
                }

                da1[0]     = 0.0;
                dc1[N1-2]  = 0.0;

                u1[j][0]  = boundary(0, j, 2*k-1);
                u1[j][N1] = boundary(N1, j, 2*k-1);

                dd1[0]    -= x1_a * u1[j][0];
                dd1[N1-2] -= x1_a * u1[j][N1];

                tomasAlgorithm(da1.data(), db1.data(), dc1.data(), dd1.data(), rx1.data(), rx1.size());

                for (unsigned int i=1; i<N1; i++)
                {
                    u1[j][i] = rx1[i-1];
                }
            }

            for (unsigned int i=0; i<=N1; i++)
            {
                u1[0][i]  = boundary(i, 0, 2*k-1);
                u1[N2][i] = boundary(i, N2, 2*k-1);
            }

            // Approximation to x2 direction
            for (unsigned int i=1; i<N1; i++)
            {
                for (unsigned int j=1; j<N2; j++)
                {
                    da2[j-1] = x2_a;
                    db2[j-1] = x2_b;
                    dc2[j-1] = x2_a;
                    dd2[j-1] = x2_c*(u1[j][i-1] - 2.0*u1[j][i] + u1[j][i+1]) + u1[j][i] + (ht/2.0) * f(i, j, 2*k);
                    //dd2[j-1] = x2_c*u1[j][i-1] + x2_d*u1[j][i] + x2_c*u1[j][i+1] + (ht/2.0) * f(i, j, 2*k);
                }
                da2[0]     = 0.0;
                dc2[N2-2]  = 0.0;

                u0[0][i]  = boundary(i, 0, 2*k);//m3(i, 2*k);
                u0[N2][i] = boundary(i, N2, 2*k);//m4(i, 2*k);

                dd2[0]    -= x2_a * u0[0][i];
                dd2[N2-2] -= x2_a * u0[N2][i];

                tomasAlgorithm(da2.data(), db2.data(), dc2.data(), dd2.data(), rx2.data(), rx2.size());

                for (unsigned int j=1; j<N2; j++)
                {
                    u0[j][i] = rx2[j-1];
                }
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                u0[j][0]  = boundary(0, j, 2*k);//m1(j, 2*k);
                u0[j][N1] = boundary(N1, j, 2*k);//m2(j, 2*k);
            }

            u[k] = u0;
        }
    }

    da1.clear();
    db1.clear();
    dc1.clear();
    dd1.clear();
    rx1.clear();

    da2.clear();
    db2.clear();
    dc2.clear();
    dd2.clear();
    rx2.clear();
}

void IBackwardParabolicEquation2D::calculateMVD(DoubleCube &psi, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    //cleaning cube
    psi.resize(M+1, N2+1, N1+1);

    DoubleMatrix psi0(N2+1, N1+1);
    DoubleMatrix psi1(N2+1, N1+1);
    \
    DoubleVector da1(N1-1);
    DoubleVector db1(N1-1);
    DoubleVector dc1(N1-1);
    DoubleVector dd1(N1-1);
    DoubleVector rx1(N1-1);

    DoubleVector da2(N2-1);
    DoubleVector db2(N2-1);
    DoubleVector dc2(N2-1);
    DoubleVector dd2(N2-1);
    DoubleVector rx2(N2-1);

    double x1_a = -(a1*a1*ht)/(2.0*h1*h1);
    double x1_b  = 1.0 + (a1*a1*ht)/(h1*h1);
    double x1_c = (a2*a2*ht)/(2.0*h2*h2);

    double x2_a = -(a2*a2*ht)/(2.0*h2*h2);
    double x2_b  = 1.0 + (a2*a2*ht)/(h2*h2);
    double x2_c = (a1*a1*ht)/(2.0*h1*h1);

    for (unsigned int k1=0; k1<=M; k1++)
    {
        unsigned int k = M-k1;

        if (k==M)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    psi0[j][i] = binitial(i, j);
                }
            }
            psi[k] = psi0;
        }
        else
        {
            for (unsigned int j=1; j<N2; j++)
            {
                for (unsigned int i=1; i<N1; i++)
                {
                    da1[i-1] = x1_a;
                    db1[i-1] = x1_b;
                    dc1[i-1] = x1_a;
                    dd1[i-1] = x1_c*(psi0[j-1][i] - 2.0*psi0[j][i] + psi0[j+1][i]) + psi0[j][i] - (ht/2.0) * bf(i, j, 2*k+1);
                    //dd1[i-1] = x1_c*psi0[j-1][i] + x1_d*psi0[j][i] + x1_c*psi0[j+1][i] - (ht/2.0) * bf(i, j, 2*k+1);
                }

                da1[0]     = 0.0;
                dc1[N1-2]  = 0.0;

                psi1[j][0]  = bboundary(0, j, 2*k+1);
                psi1[j][N1] = bboundary(N1, j, 2*k+1);

                dd1[0]    -= x1_a * psi1[j][0];
                dd1[N1-2] -= x1_a * psi1[j][N1];

                tomasAlgorithm(da1.data(), db1.data(), dc1.data(), dd1.data(), rx1.data(), rx1.size());

                for (unsigned int i=1; i<N1; i++)
                {
                    psi1[j][i] = rx1[i-1];
                }
            }

            for (unsigned int i=0; i<=N1; i++)
            {
                psi1[0][i]  = bboundary(i, 0, 2*k+1);
                psi1[N2][i] = bboundary(i, N2, 2*k+1);
            }

            for (unsigned int i=1; i<N1; i++)
            {
                for (unsigned int j=1; j<N2; j++)
                {
                    da2[j-1] = x2_a;
                    db2[j-1] = x2_b;
                    dc2[j-1] = x2_a;
                    dd2[j-1] = x2_c*(psi1[j][i-1] - 2.0*psi1[j][i] + psi1[j][i+1]) + psi1[j][i] - (ht/2.0) * bf(i, j, 2*k);
                    //dd2[j-1] = x2_c*psi1[j][i-1] + x2_d*psi1[j][i] + x2_c*psi1[j][i+1] - (ht/2.0) * bf(i, j, 2*k);
                }
                da2[0]     = 0.0;
                dc2[N2-2]  = 0.0;

                psi0[0][i]  = bboundary(i, 0, 2*k);
                psi0[N2][i] = bboundary(i, N2, 2*k);

                dd2[0]    -= x2_a * psi0[0][i];
                dd2[N2-2] -= x2_a *psi0[N2][i];

                tomasAlgorithm(da2.data(), db2.data(), dc2.data(), dd2.data(), rx2.data(), rx2.size());

                for (unsigned int j=1; j<N2; j++)
                {
                    psi0[j][i] = rx2[j-1];
                }
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                psi0[j][0]  = bboundary(0, j, 2*k);
                psi0[j][N1] = bboundary(N1, j, 2*k);
            }

            psi[k] = psi0;
        }
    }

    da1.clear();
    db1.clear();
    dc1.clear();
    dd1.clear();
    rx1.clear();

    da2.clear();
    db2.clear();
    dc2.clear();
    dd2.clear();
    rx2.clear();
}

void IParabolicEquation::calculateL(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a) const
{
    u.clear();

    u.resize(M+1, N+1);
    for (unsigned int j=0; j<=M; j++)
    {
        //u[j].resize(N+1);
        u[j][0] = boundary(Left, j);
        u[j][N] = boundary(Right, j);
    }
    for (unsigned int i=0; i<=N; i++) u[0][i] = initial(i);

    double betta = (a*a)/(hx*hx);
    double alpha = -2.0*betta;

    class A : public CauchyProblem
    {
    public:
        virtual double f(double t, const DoubleVector &y) const
        {
            unsigned int j=(unsigned int)(t/ht);
            if (i==0)          { return bt * (p->boundary(Left, j) - 2.0*y[i] + y[i+1])        + p->f(i+1, j); }
            else if (i==(N-2)) { return bt * (y[i-1]       - 2.0*y[i] + p->boundary(Right, j)) + p->f(i+1, j); }
            else               { return bt * (y[i-1]       - 2.0*y[i] + y[i+1])                + p->f(i+1, j); }
        }
        unsigned int i;
        double al;
        double bt;
        double ht;
        double hx;
        double a;
        unsigned int N;
        const IParabolicEquation *p;
    };

    std::vector<CauchyProblem*> cps;
    cps.resize(N-1);
    for (unsigned int i=1; i<=N-1; i++)
    {
        A *cp = new A;
        cp->al = alpha;
        cp->bt = betta;
        cp->N = N;
        cp->i = i-1;
        cp->x0 = 0.0;
        cp->y0 = u[0][i];
        cp->p = this;
        //cp->u = &u;
        cp->ht = ht;
        cp->hx = hx;
        cp->a = a;
        cps[i-1] = cp;
    }

    DoubleMatrix m;
    CauchyProblem::rungeKutta(cps, 0.0, ht, M, m);
    for (unsigned int j=1; j<=M; j++)
    {
        for (unsigned int i=1; i<=N-1; i++)
        {
            u[j][i] = m[i-1][j-1];
        }
    }
}

void IParabolicEquation::calculateN1(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a) const
{
    u.clear();
    u.resize(M+1, N+1);

    DoubleVector da(N+1);
    DoubleVector db(N+1);
    DoubleVector dc(N+1);
    DoubleVector dd(N+1);
    DoubleVector rx(N+1);

    double alpha = -(a*a*ht)/(hx*hx);
    double beta  = 1.0 + (2.0*a*a*ht)/(hx*hx);

    for (unsigned int j=0; j<=M; j++)
    {
        if (j == 0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u[j][i] = initial(i);
            }
        }
        else
        {
            da[0] = 0.0;
            db[0] = 1.0 + (a*a*ht)/(hx*hx);
            dc[0] = -(a*a*ht)/(hx*hx);
            dd[0] = u[j-1][0] + ht * f(0, j) - (a*a*ht)/(hx) * boundary(Left, j);

            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i] = alpha;
                db[i] = beta;
                dc[i] = alpha;
                dd[i] = u[j-1][i] + ht * f(i, j);
            }

            da[N] = -(a*a*ht)/(hx*hx);
            db[N] = 1.0 + (a*a*ht)/(hx*hx);
            dc[N] = 0.0;
            dd[N] = u[j-1][N] + ht * f(N, j) + (a*a*ht)/(hx) * boundary(Right, j);

            tomasAlgorithm(da.data(), db.data(), dc.data(), dd.data(), rx.data(), rx.size());

            for (unsigned int i=0; i<=N; i++)
            {
                u[j][i] = rx[i];
            }
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();
}
