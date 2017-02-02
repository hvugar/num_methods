#include "hyperbolicequation.h"

void saveVectorFiles(const DoubleVector& v, int i, unsigned int N)
{
    char buffer[20];
    int n = 0;
    if (i<10) n = sprintf(buffer, "data/0000000%d.txt", i);
    if (i<100 && i>=10) n = sprintf(buffer, "data/000000%d.txt", i);
    if (i<1000 && i>=100) n = sprintf(buffer, "data/00000%d.txt", i);
    if (i<10000 && i>=1000) n = sprintf(buffer, "data/0000%d.txt", i);
    if (i<100000 && i>=10000) n = sprintf(buffer, "data/000%d.txt", i);
    buffer[n] = '\0';
    FILE *file = fopen(buffer, "w");
    IPrinter::printVector(v, NULL, N, 0, 0, file);
    fclose(file);
}

double MIN = +100.0;
double MAX = -100.0;

void saveData(const DoubleMatrix& m, int i, unsigned int N2, unsigned int N1)
{
    double min = m.min();
    double max = m.max();
    if (MIN > min) MIN = min;
    if (MAX < max) MAX = max;

    char buffer[20];
    int n = 0;
    if (i<10) n = sprintf(buffer, "data/0000000%d.txt", i);
    if (i<100 && i>=10) n = sprintf(buffer, "data/000000%d.txt", i);
    if (i<1000 && i>=100) n = sprintf(buffer, "data/00000%d.txt", i);
    if (i<10000 && i>=1000) n = sprintf(buffer, "data/0000%d.txt", i);
    buffer[n] = '\0';
    FILE *file = fopen(buffer, "w");
    IPrinter::printMatrix(m, N2, N1, NULL, file);
    fclose(file);

    printf("File: %s min: %.16f max: %.16f\n", buffer, MIN, MAX);
}

void IHyperbolicEquation::calculateU(DoubleVector &u, double hx, double ht, unsigned int M, unsigned int N, double a, double lamda) const
{
    u.resize(N+1);

    DoubleVector u0(N+1);
    DoubleVector u1(N+1);

    DoubleVector da(N-1);
    DoubleVector db(N-1);
    DoubleVector dc(N-1);
    DoubleVector rd(N-1);
    DoubleVector rx(N-1);

    double alpha1 = -(lamda*a*a)*((ht*ht)/(hx*hx));
    double beta1  = 1.0 + (2.0*lamda*a*a*(ht*ht))/(hx*hx);
    double alpha2 = (1.0-2.0*lamda)*a*a*((ht*ht)/(hx*hx));
    double alpha3 = +(lamda*a*a)*((ht*ht)/(hx*hx));

    for (unsigned int j=0; j<=M-1; j++)
    {
        if (j==0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u0[i] = initial1(i);
                u1[i] = u0[i] + ht*initial2(i);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = alpha1;
                db[i-1] = beta1;
                dc[i-1] = alpha1;
                rd[i-1] = (alpha2*(u1[i+1] - 2.0*u1[i] + u1[i-1]) + 2.0*u1[i])
                        + (alpha3*(u0[i+1] - 2.0*u0[i] + u0[i-1]) - u0[i])
                        + (ht*ht)*f(i, j);
            }

            da[0]   = 0.0;
            dc[N-2] = 0.0;

            u[0] = boundary(Left, j+1); //m1(j+1);
            u[N] = boundary(Right, j+1); //m2(j+1);

            rd[0]   -= alpha1 * u[0];
            rd[N-2] -= alpha1 * u[N];

            tomasAlgorithm(da.data(), db.data(), dc.data(), rd.data(), rx.data(), rx.size());

            for (unsigned int i=1; i<=N-1; i++)
            {
                u[i] = rx[i-1];
            }

            for (unsigned int i=0; i<=N; i++)
            {
                u0[i] = u1[i];
                u1[i] = u[i];
            }
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    rd.clear();
    rx.clear();
}

void IHyperbolicEquation::calculateU(DoubleMatrix &u, double hx, double ht, unsigned int M, unsigned int N, double a, double lamda) const
{
    u.resize(M+1, N+1);

    DoubleVector da(N-1);
    DoubleVector db(N-1);
    DoubleVector dc(N-1);
    DoubleVector rd(N-1);
    DoubleVector rx(N-1);

    double alpha1 = -(lamda*a*a)*((ht*ht)/(hx*hx));
    double beta1  = 1.0 + (2.0*lamda*a*a*(ht*ht))/(hx*hx);
    double alpha2 = (1.0-2.0*lamda)*a*a*((ht*ht)/(hx*hx));
    double alpha3 = +(lamda*a*a)*((ht*ht)/(hx*hx));

    for (unsigned int j=0; j<=M-1; j++)
    {
        if (j==0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u[0][i] = initial1(i);
                u[1][i] = u[0][i] + ht*initial2(i);
            }
        }
        else
        {
            u[j+1][0] = boundary(Left,j+1);
            u[j+1][N] = boundary(Right,j+1);

            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = alpha1;
                db[i-1] = beta1;
                dc[i-1] = alpha1;
                rd[i-1] = (alpha2*(u[j][i+1]   - 2.0*u[j][i]   + u[j][i-1])   + 2.0*u[j][i])
                        + (alpha3*(u[j-1][i+1] - 2.0*u[j-1][i] + u[j-1][i-1]) - u[j-1][i])
                        + (ht*ht)*f(i, j);
            }

            da[0]   = 0.0;
            dc[N-2] = 0.0;
            rd[0]   -= alpha1 * u[j+1][0];
            rd[N-2] -= alpha1 * u[j+1][N];

            tomasAlgorithm(da.data(), db.data(), dc.data(), rd.data(), rx.data(), rx.size());

            for (unsigned int i=1; i<=N-1; i++)
            {
                u[j+1][i] = rx[i-1];
            }
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    rd.clear();
    rx.clear();
}

void IBackwardHyperbolicEquation::calculateU(DoubleVector &p, double hx, double ht, unsigned int M, unsigned int N, double a, double lamda) const
{
    p.resize(N+1);

    DoubleVector p0(N+1);
    DoubleVector p1(N+1);

    DoubleVector da(N-1);
    DoubleVector db(N-1);
    DoubleVector dc(N-1);
    DoubleVector rd(N-1);
    DoubleVector rx(N-1);

    double alpha1 = -(lamda*a*a)*((ht*ht)/(hx*hx));
    double beta1  = 1.0 + (2.0*lamda*a*a*(ht*ht))/(hx*hx);
    double alpha2 = (1.0-2.0*lamda)*a*a*((ht*ht)/(hx*hx));
    double alpha3 = +(lamda*a*a)*((ht*ht)/(hx*hx));

    //const unsigned int k = (unsigned)0-1;
    for (unsigned int j1=0; j1<=M-1; j1++)
    {
        unsigned int j = M-j1;

        if (j==M)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                p0[i] = binitial1(i);
                p1[i] = p0[i] - ht*binitial2(i);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = alpha1;
                db[i-1] = beta1;
                dc[i-1] = alpha1;
                rd[i-1] = (alpha2*(p1[i+1] - 2.0*p1[i] + p1[i-1]) + 2.0*p1[i])
                        + (alpha3*(p0[i+1] - 2.0*p0[i] + p0[i-1]) - p0[i])
                        + (ht*ht)*bf(i,j);
            }

            da[0]   = 0.0;
            dc[N-2] = 0.0;

            p[0] = bboundary(Left, j-1);
            p[N] = bboundary(Right, j-1);

            rd[0]   -= alpha1 * p[0];
            rd[N-2] -= alpha1 * p[N];

            tomasAlgorithm(da.data(), db.data(), dc.data(), rd.data(), rx.data(), rx.size());

            for (unsigned int i=1; i<=N-1; i++)
            {
                p[i] = rx[i-1];
            }

            for (unsigned int i=0; i<=N; i++)
            {
                p0[i] = p1[i];
                p1[i] = p[i];
            }
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    rd.clear();
    rx.clear();
}

void IBackwardHyperbolicEquation::calculateU(DoubleMatrix &p, double hx, double ht, unsigned int M, unsigned int N, double a, double lamda) const
{
    p.resize(M+1, N+1);

    DoubleVector da(N-1);
    DoubleVector db(N-1);
    DoubleVector dc(N-1);
    DoubleVector rd(N-1);
    DoubleVector rx(N-1);

    double alpha1 = -(lamda*a*a)*((ht*ht)/(hx*hx));
    double beta1  = 1.0 + (2.0*lamda*a*a*(ht*ht))/(hx*hx);
    double alpha2 = (1.0-2.0*lamda)*a*a*((ht*ht)/(hx*hx));
    double alpha3 = +(lamda*a*a)*((ht*ht)/(hx*hx));

    //const unsigned int k = (unsigned)0-1;
    for (unsigned int j1=0; j1<=M-1; j1++)
    {
        unsigned int j = M-j1;

        if (j==M)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                p[M][i] = binitial1(i);
                p[M-1][i] = p[M][i] - ht*binitial2(i);
            }
        }
        else
        {
            p[j-1][0] = bboundary(Left, j-1);
            p[j-1][N] = bboundary(Right, j-1);

            for (unsigned int i=1; i<=N-1; i++)
            {
                da[i-1] = alpha1;
                db[i-1] = beta1;
                dc[i-1] = alpha1;
                rd[i-1] = (alpha2*(p[j][i+1]   - 2.0*p[j][i]   + p[j][i-1])   + 2.0*p[j][i])
                        + (alpha3*(p[j+1][i+1] - 2.0*p[j+1][i] + p[j+1][i-1]) - p[j+1][i])
                        + (ht*ht)*bf(i,j);
            }

            da[0]   = 0.0;
            dc[N-2] = 0.0;
            rd[0]   -= alpha1 * p[j-1][0];
            rd[N-2] -= alpha1 * p[j-1][N];

            tomasAlgorithm(da.data(), db.data(), dc.data(), rd.data(), rx.data(), rx.size());

            for (unsigned int i=1; i<=N-1; i++)
            {
                p[j-1][i] = rx[i-1];
            }
        }
    }

    da.clear();
    db.clear();
    dc.clear();
    rd.clear();
    rx.clear();
}

void IHyperbolicEquation2D::calculateMVD(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    u.resize(N2+1, N1+1);

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

    double x1_a = -(a1*a1*ht*ht)/(h1*h1);
    double x1_b  = 1.0 + (2.0*a1*a1*ht*ht)/(h1*h1);
    double x1_c = (a2*a2*ht*ht)/(h2*h2);

    double x2_a = -(a2*a2*ht*ht)/(h2*h2);
    double x2_b  = 1.0 + (2.0*a2*a2*ht*ht)/(h2*h2);
    double x2_c = (a1*a1*ht*ht)/(h1*h1);

    for (unsigned int k=1; k<=M; k++)
    {
        if (k==1)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u0[j][i] = initial1(i, j);
                    u1[j][i] = u0[j][i] + ht*initial2(i, j);
                }
            }
        }
        else
        {
            if ((k % 2) == 0)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u[0][i]  = boundary(i, 0, k);
                    u[N2][i] = boundary(i, N2, k);
                }

                // Approximation to x1 direction
                for (unsigned int j=1; j<N2; j++)
                {
                    for (unsigned int i=1; i<N1; i++)
                    {
                        da1[i-1] = x1_a;
                        db1[i-1] = x1_b;
                        dc1[i-1] = x1_a;
                        dd1[i-1] = x1_c*(u1[j-1][i] - 2.0*u1[j][i] + u1[j+1][i]) + 2.0*u1[j][i] - u0[j][i] + (ht*ht) * f(i, j, k);
                    }

                    da1[0]     = 0.0;
                    dc1[N1-2]  = 0.0;

                    u[j][0]  = boundary(0, j, k);
                    u[j][N1] = boundary(N1, j, k);

                    dd1[0]    -= x1_a * u[j][0];
                    dd1[N1-2] -= x1_a * u[j][N1];

                    tomasAlgorithm(da1.data(), db1.data(), dc1.data(), dd1.data(), rx1.data(), rx1.size());

                    for (unsigned int i=1; i<N1; i++)
                    {
                        u[j][i] = rx1[i-1];
                    }
                }

            }
            else
            {
                for (unsigned int j=0; j<=N2; j++)
                {
                    u[j][0]  = boundary(0, j, k);
                    u[j][N1] = boundary(N1, j, k);
                }

                // Approximation to x2 direction
                for (unsigned int i=1; i<N1; i++)
                {
                    for (unsigned int j=1; j<N2; j++)
                    {
                        da2[j-1] = x2_a;
                        db2[j-1] = x2_b;
                        dc2[j-1] = x2_a;
                        dd2[j-1] = x2_c*(u1[j][i-1] - 2.0*u1[j][i] + u1[j][i+1]) + 2.0*u1[j][i] - u0[j][i] + (ht*ht) * f(i, j, k);
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
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u0[j][i] = u1[j][i];
                    u1[j][i] = u[j][i];
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

void IHyperbolicEquation2D::calculateMVD(DoubleCube &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    //cleaning cube
    u.clear();
    u.resize(M+1, N2+1, N1+1);

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

    double x1_a = -(a1*a1*ht*ht)/(h1*h1);
    double x1_b  = 1.0 + (2.0*a1*a1*ht*ht)/(h1*h1);
    double x1_c = (a2*a2*ht*ht)/(h2*h2);

    double x2_a = -(a2*a2*ht*ht)/(h2*h2);
    double x2_b  = 1.0 + (2.0*a2*a2*ht*ht)/(h2*h2);
    double x2_c = (a1*a1*ht*ht)/(h1*h1);

    for (unsigned int k=1; k<=M; k++)
    {
        if (k==1)
        {
            //u[0].resize(N2+1, N1+1);
            //u[1].resize(N2+1, N1+1);

            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u[0][j][i] = initial1(i, j);
                    u[1][j][i] = u[0][j][i] + ht*initial2(i, j);
                }
            }
        }
        else
        {
            //u[k].resize(N2+1, N1+1);

            if ((k % 2) == 0)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u[k][0][i]  = boundary(i, 0, k);
                    u[k][N2][i] = boundary(i, N2, k);
                }

                // Approximation to x1 direction
                for (unsigned int j=1; j<N2; j++)
                {
                    for (unsigned int i=1; i<N1; i++)
                    {
                        da1[i-1] = x1_a;
                        db1[i-1] = x1_b;
                        dc1[i-1] = x1_a;
                        dd1[i-1] = x1_c*(u[k-1][j-1][i] - 2.0*u[k-1][j][i] + u[k-1][j+1][i]) + 2.0*u[k-1][j][i] - u[k-2][j][i] + (ht*ht) * f(i, j, k);
                    }

                    da1[0]     = 0.0;
                    dc1[N1-2]  = 0.0;

                    u[k][j][0]  = boundary(0, j, k);
                    u[k][j][N1] = boundary(N1, j, k);

                    dd1[0]    -= x1_a * u[k][j][0];
                    dd1[N1-2] -= x1_a * u[k][j][N1];

                    tomasAlgorithm(da1.data(), db1.data(), dc1.data(), dd1.data(), rx1.data(), rx1.size());

                    for (unsigned int i=1; i<N1; i++)
                    {
                        u[k][j][i] = rx1[i-1];
                    }
                }

            }
            else
            {
                for (unsigned int j=0; j<=N2; j++)
                {
                    u[k][j][0]  = boundary(0, j, k);
                    u[k][j][N1] = boundary(N1, j, k);
                }

                // Approximation to x2 direction
                for (unsigned int i=1; i<N1; i++)
                {
                    for (unsigned int j=1; j<N2; j++)
                    {
                        da2[j-1] = x2_a;
                        db2[j-1] = x2_b;
                        dc2[j-1] = x2_a;
                        dd2[j-1] = x2_c*(u[k-1][j][i-1] - 2.0*u[k-1][j][i] + u[k-1][j][i+1]) + 2.0*u[k-1][j][i] - u[k-2][j][i] + (ht*ht) * f(i, j, k);
                    }
                    da2[0]     = 0.0;
                    dc2[N2-2]  = 0.0;

                    u[k][0][i]  = boundary(i, 0, k);
                    u[k][N2][i] = boundary(i, N2, k);

                    dd2[0]    -= x2_a * u[k][0][i];
                    dd2[N2-2] -= x2_a * u[k][N2][i];

                    tomasAlgorithm(da2.data(), db2.data(), dc2.data(), dd2.data(), rx2.data(), rx2.size());

                    for (unsigned int j=1; j<N2; j++)
                    {
                        u[k][j][i] = rx2[j-1];
                    }
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

void IHyperbolicEquation2D::calculateU1(DoubleCube &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2, double qamma) const
{
    //cleaning cube
    u.clear();
    u.resize(M+1, N2+1, N1+1);

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

    double x1_a = -(a1*a1*ht*ht)/(h1*h1);
    double x1_b  = 1.0 + (2.0*a1*a1*ht*ht)/(h1*h1) + qamma*ht;
    double x1_c = (a2*a2*ht*ht)/(h2*h2);

    double x2_a = -(a2*a2*ht*ht)/(h2*h2);
    double x2_b  = 1.0 + (2.0*a2*a2*ht*ht)/(h2*h2) + qamma*ht;
    double x2_c = (a1*a1*ht*ht)/(h1*h1);

    for (unsigned int m=1; m<=M; m++)
    {
        if (m==1)
        {
            for (unsigned int n2=0; n2<=N2; n2++)
            {
                for (unsigned int n1=0; n1<=N1; n1++)
                {
                    u.at(0,n2,n1) = initial1(n1, n2);
                    u.at(1,n2,n1) = u(0,n2,n1) + ht*initial2(n1, n2);
                }
            }
        }
        else
        {
            if ((m % 2) == 0)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u.at(m,0,i)  = boundary(i, 0, m);
                    u.at(m,N2,i) = boundary(i, N2, m);
                }

                // Approximation to x1 direction
                for (unsigned int j=1; j<N2; j++)
                {
                    for (unsigned int i=1; i<N1; i++)
                    {
                        da1[i-1] = x1_a;
                        db1[i-1] = x1_b;
                        dc1[i-1] = x1_a;
                        dd1[i-1] = x1_c*(u.at(m-1,j-1,i) - 2.0*u.at(m-1,j,i) + u.at(m-1,j+1,i)) + 2.0*u.at(m-1,j,i) - u.at(m-2,j,i) + qamma*ht*u.at(m-1,j,i) + (ht*ht) * f(i, j, m);
                    }

                    da1[0]     = 0.0;
                    dc1[N1-2]  = 0.0;

                    u.at(m,j,0)  = boundary(0, j, m);
                    u.at(m,j,N1) = boundary(N1, j, m);

                    dd1[0]    -= x1_a * u.at(m,j,0);
                    dd1[N1-2] -= x1_a * u.at(m,j,N1);

                    tomasAlgorithm(da1.data(), db1.data(), dc1.data(), dd1.data(), rx1.data(), rx1.size());

                    for (unsigned int i=1; i<N1; i++)
                    {
                        u.at(m,j,i) = rx1[i-1];
                    }
                }

            }
            else
            {
                for (unsigned int j=0; j<=N2; j++)
                {
                    u.at(m,j,0)  = boundary(0, j, m);
                    u.at(m,j,N1) = boundary(N1, j, m);
                }

                // Approximation to x2 direction
                for (unsigned int i=1; i<N1; i++)
                {
                    for (unsigned int j=1; j<N2; j++)
                    {
                        da2[j-1] = x2_a;
                        db2[j-1] = x2_b;
                        dc2[j-1] = x2_a;
                        dd2[j-1] = x2_c*(u.at(m-1,j,i-1) - 2.0*u.at(m-1,j,i) + u(m-1,j,i+1)) + 2.0*u.at(m-1,j,i) - u(m-2,j,i) + qamma*ht*u.at(m-1,j,i) + (ht*ht) * f(i, j, m);
                    }
                    da2[0]     = 0.0;
                    dc2[N2-2]  = 0.0;

                    u.at(m,0,i)  = boundary(i, 0, m);
                    u.at(m,N2,i) = boundary(i, N2, m);

                    dd2[0]    -= x2_a * u.at(m,0,i);
                    dd2[N2-2] -= x2_a * u.at(m,N2,i);

                    tomasAlgorithm(da2.data(), db2.data(), dc2.data(), dd2.data(), rx2.data(), rx2.size());

                    for (unsigned int j=1; j<N2; j++)
                    {
                        u.at(m,j,i) = rx2[j-1];
                    }
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

void IHyperbolicEquation2D::calculateU1(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2, double qamma) const
{
    //cleaning matrix
    u.resize(N2+1, N1+1);

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

    double x1_a = -(a1*a1*ht*ht)/(h1*h1);
    double x1_b  = 1.0 + (2.0*a1*a1*ht*ht)/(h1*h1) + qamma*ht;
    double x1_c = (a1*a1*ht*ht)/(h1*h1);

    double x2_a = -(a2*a2*ht*ht)/(h2*h2);
    double x2_b  = 1.0 + (2.0*a2*a2*ht*ht)/(h2*h2) + qamma*ht;
    double x2_c = (a2*a2*ht*ht)/(h2*h2);

    for (unsigned int k=1; k<=M; k++)
    {
        if (k==1)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u0[j][i] = initial1(i, j);
                    u1[j][i] = u0[j][i] + ht*initial2(i, j);
                }
            }

            saveData(u0, 0, N2, N1);
            saveData(u1, 1, N2, N1);
        }
        else
        {
            if ((k % 2) == 0)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u[0][i]  = boundary(i, 0, k);
                    u[N2][i] = boundary(i, N2, k);
                }

                // Approximation to x1 direction
                for (unsigned int j=1; j<N2; j++)
                {
                    for (unsigned int i=1; i<N1; i++)
                    {
                        da1[i-1] = x1_a;
                        db1[i-1] = x1_b;
                        dc1[i-1] = x1_a;
                        dd1[i-1] = x1_c*(u1[j-1][i] - 2.0*u1[j][i] + u1[j+1][i]) + 2.0*u1[j][i] - u0[j][i] + qamma*ht*u1[j][i] + (ht*ht) * f(i, j, k);
                    }

                    da1[0]     = 0.0;
                    dc1[N1-2]  = 0.0;

                    u[j][0]  = boundary(0, j, k);
                    u[j][N1] = boundary(N1, j, k);

                    dd1[0]    -= x1_a * u[j][0];
                    dd1[N1-2] -= x1_a * u[j][N1];

                    tomasAlgorithm(da1.data(), db1.data(), dc1.data(), dd1.data(), rx1.data(), rx1.size());

                    for (unsigned int i=1; i<N1; i++)
                    {
                        u[j][i] = rx1[i-1];
                    }
                }

            }
            else
            {
                for (unsigned int j=0; j<=N2; j++)
                {
                    u[j][0]  = boundary(0, j, k);
                    u[j][N1] = boundary(N1, j, k);
                }

                // Approximation to x2 direction
                for (unsigned int i=1; i<N1; i++)
                {
                    for (unsigned int j=1; j<N2; j++)
                    {
                        da2[j-1] = x2_a;
                        db2[j-1] = x2_b;
                        dc2[j-1] = x2_a;
                        dd2[j-1] = x2_c*(u1[j][i-1] - 2.0*u1[j][i] + u1[j][i+1]) + 2.0*u1[j][i] - u0[j][i] + qamma*ht*u1[j][i] + (ht*ht) * f(i, j, k);
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
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u0[j][i] = u1[j][i];
                    u1[j][i] = u[j][i];
                }
            }

            saveData(u, k, N2, N1);
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

void IBackwardHyperbolicEquation2D::calculateU(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    //cleaning matrix
    u.resize(N2+1, N1+1);
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

    double x1_a = -(a1*a1*ht*ht)/(h1*h1);
    double x1_b  = 1.0 + (2.0*a1*a1*ht*ht)/(h1*h1);
    double x1_c = (a1*a1*ht*ht)/(h1*h1);

    double x2_a = -(a2*a2*ht*ht)/(h2*h2);
    double x2_b  = 1.0 + (2.0*a2*a2*ht*ht)/(h2*h2);
    double x2_c = (a2*a2*ht*ht)/(h2*h2);

    for (unsigned int k=1; k<=M; k++)
    {
        if (k==1)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u0[j][i] = binitial1(i, j);
                    u1[j][i] = u0[j][i] + ht*binitial2(i, j);
                }
            }
        }
        else
        {
            if ((k % 2) == 0)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u[0][i]  = bboundary(i, 0, k);
                    u[N2][i] = bboundary(i, N2, k);
                }

                // Approximation to x1 direction
                for (unsigned int j=1; j<N2; j++)
                {
                    for (unsigned int i=1; i<N1; i++)
                    {
                        da1[i-1] = x1_a;
                        db1[i-1] = x1_b;
                        dc1[i-1] = x1_a;
                        dd1[i-1] = x1_c*(u1[j-1][i] - 2.0*u1[j][i] + u1[j+1][i]) + 2.0*u1[j][i] - u0[j][i] + (ht*ht) * bf(i, j, k);
                    }

                    da1[0]     = 0.0;
                    dc1[N1-2]  = 0.0;

                    u[j][0]  = bboundary(0, j, k);
                    u[j][N1] = bboundary(N1, j, k);

                    dd1[0]    -= x1_a * u[j][0];
                    dd1[N1-2] -= x1_a * u[j][N1];

                    tomasAlgorithm(da1.data(), db1.data(), dc1.data(), dd1.data(), rx1.data(), rx1.size());

                    for (unsigned int i=1; i<N1; i++)
                    {
                        u[j][i] = rx1[i-1];
                    }
                }

            }
            else
            {
                for (unsigned int j=0; j<=N2; j++)
                {
                    u[j][0]  = bboundary(0, j, k);
                    u[j][N1] = bboundary(N1, j, k);
                }

                // Approximation to x2 direction
                for (unsigned int i=1; i<N1; i++)
                {
                    for (unsigned int j=1; j<N2; j++)
                    {
                        da2[j-1] = x2_a;
                        db2[j-1] = x2_b;
                        dc2[j-1] = x2_a;
                        dd2[j-1] = x2_c*(u1[j][i-1] - 2.0*u1[j][i] + u1[j][i+1]) + 2.0*u1[j][i] - u0[j][i] + (ht*ht) * bf(i, j, k);
                    }
                    da2[0]     = 0.0;
                    dc2[N2-2]  = 0.0;

                    u[0][i]  = bboundary(i, 0, k);
                    u[N2][i] = bboundary(i, N2, k);

                    dd2[0]    -= x2_a * u[0][i];
                    dd2[N2-2] -= x2_a * u[N2][i];

                    tomasAlgorithm(da2.data(), db2.data(), dc2.data(), dd2.data(), rx2.data(), rx2.size());

                    for (unsigned int j=1; j<N2; j++)
                    {
                        u[j][i] = rx2[j-1];
                    }
                }
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u0[j][i] = u1[j][i];
                    u1[j][i] = u[j][i];
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

void IBackwardHyperbolicEquation2D::calculateU(DoubleCube &p, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    //cleaning cube
    p.resize(M+1, N2+1, N1+1);

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

    double x1_a = -(a1*a1*ht*ht)/(h1*h1);
    double x1_b  = 1.0 + (2.0*a1*a1*ht*ht)/(h1*h1);
    double x1_c = (a2*a2*ht*ht)/(h2*h2);

    double x2_a = -(a2*a2*ht*ht)/(h2*h2);
    double x2_b  = 1.0 + (2.0*a2*a2*ht*ht)/(h2*h2);
    double x2_c = (a1*a1*ht*ht)/(h1*h1);

    for (unsigned int k1=1; k1<=M; k1++)
    {
        unsigned int k = M - k1;

        if (k==(M-1))
        {
            //p[M].resize(N2+1, N1+1);   //for (unsigned int j=0; j<=N2; j++) p[M][j].resize(N1+1);
            //p[M-1].resize(N2+1, N1+1); //for (unsigned int j=0; j<=N2; j++) p[M-1][j].resize(N1+1);

            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    p[M][j][i] = binitial1(i, j);
                    p[M-1][j][i] = p[M][j][i] - ht*binitial2(i, j);
                }
            }
        }
        else
        {
            //p[k].resize(N2+1, N1+1);

            if ((k % 2) == 0)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    p[k][0][i]  = bboundary(i, 0, k);
                    p[k][N2][i] = bboundary(i, N2, k);
                }

                // Approximation to x1 direction
                for (unsigned int j=1; j<N2; j++)
                {
                    for (unsigned int i=1; i<N1; i++)
                    {
                        da1[i-1] = x1_a;
                        db1[i-1] = x1_b;
                        dc1[i-1] = x1_a;
                        dd1[i-1] = x1_c*(p[k+1][j-1][i] - 2.0*p[k+1][j][i] + p[k+1][j+1][i]) + 2.0*p[k+1][j][i] - p[k+2][j][i] + (ht*ht) * bf(i, j, k);
                    }

                    da1[0]     = 0.0;
                    dc1[N1-2]  = 0.0;

                    p[k][j][0]  = bboundary(0, j, k);
                    p[k][j][N1] = bboundary(N1, j, k);

                    dd1[0]    -= x1_a * p[k][j][0];
                    dd1[N1-2] -= x1_a * p[k][j][N1];

                    tomasAlgorithm(da1.data(), db1.data(), dc1.data(), dd1.data(), rx1.data(), rx1.size());

                    for (unsigned int i=1; i<N1; i++)
                    {
                        p[k][j][i] = rx1[i-1];
                    }
                }

            }
            else
            {
                for (unsigned int j=0; j<=N2; j++)
                {
                    p[k][j][0]  = bboundary(0, j, k);
                    p[k][j][N1] = bboundary(N1, j, k);
                }

                // Approximation to x2 direction
                for (unsigned int i=1; i<N1; i++)
                {
                    for (unsigned int j=1; j<N2; j++)
                    {
                        da2[j-1] = x2_a;
                        db2[j-1] = x2_b;
                        dc2[j-1] = x2_a;
                        dd2[j-1] = x2_c*(p[k+1][j][i-1] - 2.0*p[k+1][j][i] + p[k+1][j][i+1]) + 2.0*p[k+1][j][i] - p[k+2][j][i] + (ht*ht) * bf(i, j, k);
                    }
                    da2[0]     = 0.0;
                    dc2[N2-2]  = 0.0;

                    p[k][0][i]  = bboundary(i, 0, k);
                    p[k][N2][i] = bboundary(i, N2, k);

                    dd2[0]    -= x2_a * p[k][0][i];
                    dd2[N2-2] -= x2_a * p[k][N2][i];

                    tomasAlgorithm(da2.data(), db2.data(), dc2.data(), dd2.data(), rx2.data(), rx2.size());

                    for (unsigned int j=1; j<N2; j++)
                    {
                        p[k][j][i] = rx2[j-1];
                    }
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

void IBackwardHyperbolicEquation2D::calculateU1(DoubleCube &p, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2, double qamma) const
{
    //cleaning cube
    p.resize(M+1, N2+1, N1+1);

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

    double x1_a = -(a1*a1*ht*ht)/(h1*h1);
    double x1_b  = 1.0 + (2.0*a1*a1*ht*ht)/(h1*h1) - qamma*ht;
    double x1_c = (a2*a2*ht*ht)/(h2*h2);

    double x2_a = -(a2*a2*ht*ht)/(h2*h2);
    double x2_b  = 1.0 + (2.0*a2*a2*ht*ht)/(h2*h2) + qamma*ht;
    double x2_c = (a1*a1*ht*ht)/(h1*h1);

    for (unsigned int k1=1; k1<=M; k1++)
    {
        unsigned int k = M - k1;

        if (k==(M-1))
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    p.at(M,j,i) = binitial1(i, j);
                    p.at(M-1,j,i) = p.at(M,j,i) - ht*binitial2(i, j);
                }
            }
        }
        else
        {
            if ((k % 2) == 0)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    p.at(k,0,i)  = bboundary(i, 0, k);
                    p.at(k,N2,i) = bboundary(i, N2, k);
                }

                // Approximation to x1 direction
                for (unsigned int j=1; j<N2; j++)
                {
                    for (unsigned int i=1; i<N1; i++)
                    {
                        da1[i-1] = x1_a;
                        db1[i-1] = x1_b;
                        dc1[i-1] = x1_a;
                        dd1[i-1] = x1_c*(p.at(k+1,j-1,i) - 2.0*p.at(k+1,j,i) + p.at(k+1,j+1,i)) + 2.0*p.at(k+1,j,i) - p.at(k+2,j,i) - qamma*ht*p.at(k+1,j,i) + (ht*ht) * bf(i, j, k);
                    }

                    da1[0]     = 0.0;
                    dc1[N1-2]  = 0.0;

                    p.at(k,j,0)  = bboundary(0, j, k);
                    p.at(k,j,N1) = bboundary(N1, j, k);

                    dd1[0]    -= x1_a * p.at(k,j,0);
                    dd1[N1-2] -= x1_a * p.at(k,j,N1);

                    tomasAlgorithm(da1.data(), db1.data(), dc1.data(), dd1.data(), rx1.data(), rx1.size());

                    for (unsigned int i=1; i<N1; i++)
                    {
                        p.at(k,j,i) = rx1[i-1];
                    }
                }

            }
            else
            {
                for (unsigned int j=0; j<=N2; j++)
                {
                    p.at(k,j,0)  = bboundary(0, j, k);
                    p.at(k,j,N1) = bboundary(N1, j, k);
                }

                // Approximation to x2 direction
                for (unsigned int i=1; i<N1; i++)
                {
                    for (unsigned int j=1; j<N2; j++)
                    {
                        da2[j-1] = x2_a;
                        db2[j-1] = x2_b;
                        dc2[j-1] = x2_a;
                        dd2[j-1] = x2_c*(p.at(k+1,j,i-1) - 2.0*p.at(k+1,j,i) + p.at(k+1,j,i+1)) + 2.0*p.at(k+1,j,i) - p.at(k+2,j,i) + qamma*ht*p.at(k+1,j,i) + (ht*ht) * bf(i, j, k);
                    }
                    da2[0]     = 0.0;
                    dc2[N2-2]  = 0.0;

                    p.at(k,0,i)  = bboundary(i, 0, k);
                    p.at(k,N2,i) = bboundary(i, N2, k);

                    dd2[0]    -= x2_a * p.at(k,0,i);
                    dd2[N2-2] -= x2_a * p.at(k,N2,i);

                    tomasAlgorithm(da2.data(), db2.data(), dc2.data(), dd2.data(), rx2.data(), rx2.size());

                    for (unsigned int j=1; j<N2; j++)
                    {
                        p.at(k,j,i) = rx2[j-1];
                    }
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
