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
    for (unsigned int i=0; i<u.size(); i++) u[i].clear();
    u.clear();

    u.resize(M+1);
    for (unsigned int i=0; i<u.size(); i++) u[i].resize(N+1);

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

void IParabolicEquation::calculateN(DoubleMatrix &u, double hx, double ht, unsigned int N, unsigned int M, double a) const
{
    for (unsigned int i=0; i<u.size(); i++) u[i].clear();
    u.clear();

    u.resize(M+1);
    for (unsigned int i=0; i<u.size(); i++) u[i].resize(N+1);

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

    da.clear();
    db.clear();
    dc.clear();
    dd.clear();
    rx.clear();
}

void IBackwardParabolicEquation::calculateU(DoubleMatrix &psi, double hx, double ht, unsigned int N, unsigned int M, double a) const
{
    for (unsigned int i=0; i<psi.size(); i++) psi[i].clear();
    psi.clear();

    psi.resize(M+1);
    for (unsigned int j=0; j<psi.size(); j++) psi[j].resize(N+1);

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

void IParabolicEquation2D::caluclateMVD(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    //cleaning matrix
    for (unsigned int j=0; j<u.size(); j++) u[j].clear();
    u.clear();
    u.resize(N2+1); for (unsigned int j=0; j<=N2; j++) u[j].resize(N1+1);

    DoubleMatrix uh;
    uh.resize(N2+1); for (unsigned int j=0; j<=N2; j++) uh[j].resize(N1+1);

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

                uh[j][0]  = boundary(0, j, 2*k-1);//m1(j, 2*k-1);
                uh[j][N1] = boundary(N1, j, 2*k-1);//m2(j, 2*k-1);

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
                uh[0][i]  = boundary(i, 0, 2*k-1);//m3(i, 2*k-1);
                uh[N2][i] = boundary(i, N2, 2*k-1);//m4(i, 2*k-1);
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

                u[0][i]  = boundary(i, 0, 2*k);//m3(i, 2*k);
                u[N2][i] = boundary(i, N2, 2*k);//m4(i, 2*k);

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
                u[j][0]  = boundary(0, j, 2*k);//m1(j, 2*k);
                u[j][N1] = boundary(N1, j, 2*k);//m2(j, 2*k);
            }
        }

        //if (k%100==0)
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

void IParabolicEquation2D::caluclateMVD1(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    //cleaning matrix
    for (unsigned int j=0; j<u.size(); j++) u[j].clear();
    u.clear();
    u.resize(N2+1); for (unsigned int j=0; j<=N2; j++) u[j].resize(N1+1);

    DoubleMatrix uh;
    uh.resize(N2+1); for (unsigned int j=0; j<=N2; j++) uh[j].resize(N1+1);

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

void IParabolicEquation2D::caluclateMFS(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    //cleaning matrix
    for (unsigned int j=0; j<u.size(); j++) u[j].clear();
    u.clear();
    u.resize(N2+1); for (unsigned int j=0; j<=N2; j++) u[j].resize(N1+1);

    DoubleMatrix uh;
    uh.resize(N2+1); for (unsigned int j=0; j<=N2; j++) uh[j].resize(N1+1);

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

void IParabolicEquation2D::caluclateMVD(DoubleCube &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    //cleaning cube
    for (unsigned int k=0; k<u.size(); k++)
    {
        unsigned int uk_size = u[k].size();
        for (unsigned int j=0; j<uk_size; j++) u[k][j].clear();
        u[k].clear();
    }
    u.clear();
    u.resize(M+1);

    DoubleMatrix u0(N2+1);
    for (unsigned int j=0; j<=N2; j++) u0[j].resize(N1+1);
    DoubleMatrix u1(N2+1);
    for (unsigned int j=0; j<=N2; j++) u1[j].resize(N1+1);

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

                u1[j][0]  = boundary(0, j, 2*k-1);//m1(j, 2*k-1);
                u1[j][N1] = boundary(N1, j, 2*k-1);//m2(j, 2*k-1);

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
                u1[0][i]  = boundary(i, 0, 2*k-1);//m3(i, 2*k-1);
                u1[N2][i] = boundary(i, N2, 2*k-1);//m4(i, 2*k-1);
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

void IBackwardParabolicEquation2D::caluclateMVD(DoubleCube &psi, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    //cleaning cube
    for (unsigned int k=0; k<psi.size(); k++)
    {
        unsigned int uk_size = psi[k].size();
        for (unsigned int j=0; j<uk_size; j++) psi[k][j].clear();
        psi[k].clear();
    }
    psi.clear();
    psi.resize(M+1);

    DoubleMatrix psi0(N2+1);
    for (unsigned int j=0; j<=N2; j++) psi0[j].resize(N1+1);
    DoubleMatrix psi1(N2+1);
    for (unsigned int j=0; j<=N2; j++) psi1[j].resize(N1+1);

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
                    //                    dd1[i-1] = x1_c*psi0[j-1][i] + x1_d*psi0[j][i] + x1_c*psi0[j+1][i] - (ht/2.0) * bf(i, j, 2*k+1);
                }

                da1[0]     = 0.0;
                dc1[N1-2]  = 0.0;

                psi1[j][0]  = bboundary(0, j, 2*k+1);//bm1(j, 2*k+1);
                psi1[j][N1] = bboundary(N1, j, 2*k+1);//bm2(j, 2*k+1);

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
                psi1[0][i]  = bboundary(i, 0, 2*k+1);//bm3(i, 2*k+1);
                psi1[N2][i] = bboundary(i, N2, 2*k+1);//bm4(i, 2*k+1);
            }

            for (unsigned int i=1; i<N1; i++)
            {
                for (unsigned int j=1; j<N2; j++)
                {
                    da2[j-1] = x2_a;
                    db2[j-1] = x2_b;
                    dc2[j-1] = x2_a;
                    dd2[j-1] = x2_c*(psi1[j][i-1] - 2.0*psi1[j][i] + psi1[j][i+1]) + psi1[j][i] - (ht/2.0) * bf(i, j, 2*k);
                    //                    dd2[j-1] = x2_c*psi1[j][i-1] + x2_d*psi1[j][i] + x2_c*psi1[j][i+1] - (ht/2.0) * bf(i, j, 2*k);
                }
                da2[0]     = 0.0;
                dc2[N2-2]  = 0.0;

                psi0[0][i]  = bboundary(i, 0, 2*k);//bm3(i, 2*k);
                psi0[N2][i] = bboundary(i, N2, 2*k);//bm4(i, 2*k);

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
                psi0[j][0]  = bboundary(0, j, 2*k);//bm1(j, 2*k);
                psi0[j][N1] = bboundary(N1, j, 2*k);//bm2(j, 2*k);
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

    u.resize(M+1);
    for (unsigned int j=0; j<=M; j++)
    {
        u[j].resize(N+1);
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
            unsigned int j=(unsigned int)(t/0.001);
            if (i==1)          { return bt * (p->boundary(Left, j) - 2.0*y[i] + y[i+1])       + p->f(i, j); }
            else if (i==(N-1)) { return bt * (y[i-1]       - 2.0*y[i] + p->boundary(Right, j)) + p->f(i, j); }
            else               { return bt * (y[i-1]       - 2.0*y[i] + y[i+1])       + p->f(i, j); }
        }
        unsigned int i;
        double al;
        double bt;
        DoubleMatrix *u;
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
        cp->i = i;
        cp->x0 = 0.0;
        cp->y0 = u[0][i];
        cp->p = this;
        cp->u = &u;
        cps[i-1] = cp;
    }

    DoubleMatrix m;
    CauchyProblem::euler1(cps, 0.0, ht, M, m);
    for (unsigned int j=1; j<=M; j++)
    {
        for (unsigned int i=1; i<=N-1; i++)
        {
            u[j][i] = m[j-1][i-1];
        }
    }
}
