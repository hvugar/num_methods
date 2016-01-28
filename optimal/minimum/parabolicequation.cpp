#include "parabolicequation.h"
#include "cmethods.h"

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
                u[i] = fi(i);
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

            u[0] = m1(j);
            u[N] = m2(j);

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
                u[j][i] = fi(i);
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

            u[j][0] = m1(j);
            u[j][N] = m2(j);

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
            for (unsigned i=0; i<=N; i++)
            {
                psi[j][i] = bfi(i);
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

            psi[j][0] = bm1(j);
            psi[j][N] = bm2(j);

            dd[0]   -= alpha * bm1(j);
            dd[N-2] -= alpha * bm2(j);

            TomasAlgorithm(da, db, dc, dd, rx);

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

void IParabolicEquation2D::calculateU(DoubleMatrix &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    for (unsigned int j=0; j<u.size(); j++) u[j].clear();
    u.clear();

    DoubleMatrix u0;
    DoubleMatrix u1;

    u0.resize(N2+1); for (unsigned int j=0; j<=N2; j++) u0[j].resize(N1+1);
    u1.resize(N2+1); for (unsigned int j=0; j<=N2; j++) u1[j].resize(N1+1);

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

    double x1_alpha1 = -(a2*ht)/(2.0*h2*h2);
    double x1_beta1  = 1.0 + (a2*ht)/(h2*h2);
    double x1_alpha2 = (a1*ht)/(2.0*h1*h1);
    double x1_beta2  = 1.0 - (a1*ht)/(h1*h1);

    double x2_alpha1 = -(a1*ht)/(2.0*h1*h1);
    double x2_beta1  = 1.0 + (a1*ht)/(h1*h1);
    double x2_alpha2 = (a2*ht)/(2.0*h2*h2);
    double x2_beta2  = 1.0 - (a2*ht)/(h2*h2);

    for (unsigned int k=0; k<=M; k++)
    {
        if (k==0)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u0[j][i] = fi(i, j);
                }
            }
        }
        else
        {
            // Approximation to x1 direction
            for (unsigned int i=1; i<N1; i++)
            {
                for (unsigned int j=1; j<N2; j++)
                {
                    da1[j-1] = x1_alpha1;
                    db1[j-1] = x1_beta1;
                    dc1[j-1] = x1_alpha1;
                    dd1[j-1] = x1_alpha2*u0[j][i-1] + x1_beta2*u0[j][i] + x1_alpha2*u0[j][i+1] + (ht/2.0) * f(i, j, 2*k-1);
                }

                da1[0]     = 0.0;
                dc1[N2-2]  = 0.0;
                dd1[0]    -= x1_alpha1 * m3(i, 2*k-1);
                dd1[N2-2] -= x1_alpha1 * m4(i, 2*k-1);

                TomasAlgorithm(da1, db1, dc1, dd1, rx1);

                u1[0][i]  = m3(i, 2*k-1);
                for (unsigned int j=1; j<N2; j++)
                {
                    u1[j][i] = rx1[j-1];
                }
                u1[N2][i] = m4(i, 2*k-1);
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                u1[j][0]  = m1(j, 2*k-1);
                u1[j][N1] = m2(j, 2*k-1);
            }

            // Approximation to x2 direction
            for (unsigned int j=1; j<N2; j++)
            {
                for (unsigned int i=1; i<N1; i++)
                {
                    da2[i-1] = x2_alpha1;
                    db2[i-1] = x2_beta1;
                    dc2[i-1] = x2_alpha1;
                    dd2[i-1] = x2_alpha2*u1[j-1][i] + x2_beta2*u1[j][i] + x2_alpha2*u1[j+1][i] + (ht/2.0) * f(i, j, 2*k);
                }
                da2[0]     = 0.0;
                dc2[N1-2]  = 0.0;
                dd2[0]    -= x2_alpha1 * m1(j, 2*k);
                dd2[N1-2] -= x2_alpha1 * m2(j, 2*k);

                TomasAlgorithm(da2, db2, dc2, dd2, rx2);

                u0[j][0]  = m1(j, 2*k);
                for (unsigned int i=1; i<N1; i++)
                {
                    u0[j][i] = rx2[i-1];
                }
                u0[j][N1] = m2(j, 2*k);
            }

            for (unsigned int i=0; i<=N1; i++)
            {
                u0[0][i]  = m3(i, 2*k);
                u0[N2][i] = m4(i, 2*k);
            }
        }
        char buffer[20];
        int n = sprintf(buffer, "data/1000/heat_%d.txt", k);
        buffer[n] = 0;
        FILE *file = fopen(buffer, "w");
        IPrinter::printMatrix(u0, N2, N1, NULL, file);
        fclose(file);
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

    u = u0;
}

void IParabolicEquation2D::calculateU(DoubleCube &u, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    u.clear();
    u.resize(M+1);

    DoubleMatrix u0(N2+1);
    DoubleMatrix u1(N2+1);

    for (unsigned int j=0; j<=N2; j++) u0[j].resize(N1+1);
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

    double x1_alpha1 = -(a2*ht)/(2.0*h2*h2);
    double x1_beta1  = 1.0 + (a2*ht)/(h2*h2);
    double x1_alpha2 = (a1*ht)/(2.0*h1*h1);
    double x1_beta2  = 1.0 - (a1*ht)/(h1*h1);

    double x2_alpha1 = -(a1*ht)/(2.0*h1*h1);
    double x2_beta1  = 1.0 + (a1*ht)/(h1*h1);
    double x2_alpha2 = (a2*ht)/(2.0*h2*h2);
    double x2_beta2  = 1.0 - (a2*ht)/(h2*h2);

    for (unsigned int k=0; k<=M; k++)
    {
        if (k==0)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    u0[j][i] = fi(i, j);
                }
                u[k] = u0;
            }
        }
        else
        {
            // Approximation to x1 direction
            for (unsigned int i=1; i<N1; i++)
            {
                for (unsigned int j=1; j<N2; j++)
                {
                    da1[j-1] = x1_alpha1;
                    db1[j-1] = x1_beta1;
                    dc1[j-1] = x1_alpha1;
                    dd1[j-1] = x1_alpha2*u0[j][i-1] + x1_beta2*u0[j][i] + x1_alpha2*u0[j][i+1] + (ht/2.0) * f(i, j, 2*k-1);
                }

                da1[0]     = 0.0;
                dc1[N2-2]  = 0.0;
                dd1[0]    -= x1_alpha1 * m3(i, 2*k-1);
                dd1[N2-2] -= x1_alpha1 * m4(i, 2*k-1);

                TomasAlgorithm(da1, db1, dc1, dd1, rx1);

                u1[0][i]  = m3(i, 2*k-1);
                for (unsigned int j=1; j<N2; j++)
                {
                    u1[j][i] = rx1[j-1];
                }
                u1[N2][i] = m4(i, 2*k-1);
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                u1[j][0]  = m1(j, 2*k-1);
                u1[j][N1] = m2(j, 2*k-1);
            }

            // Approximation to x2 direction
            for (unsigned int j=1; j<N2; j++)
            {
                for (unsigned int i=1; i<N1; i++)
                {
                    da2[i-1] = x2_alpha1;
                    db2[i-1] = x2_beta1;
                    dc2[i-1] = x2_alpha1;
                    dd2[i-1] = x2_alpha2*u1[j-1][i] + x2_beta2*u1[j][i] + x2_alpha2*u1[j+1][i] + (ht/2.0) * f(i, j, 2*k);
                }
                da2[0]     = 0.0;
                dc2[N1-2]  = 0.0;
                dd2[0]    -= x2_alpha1 * m1(j, 2*k);
                dd2[N1-2] -= x2_alpha1 * m2(j, 2*k);

                TomasAlgorithm(da2, db2, dc2, dd2, rx2);

                u0[j][0]  = m1(j, 2*k);
                for (unsigned int i=1; i<N1; i++)
                {
                    u0[j][i] = rx2[i-1];
                }
                u0[j][N1] = m2(j, 2*k);
            }

            for (unsigned int i=0; i<=N1; i++)
            {
                u0[0][i]  = m3(i, 2*k);
                u0[N2][i] = m4(i, 2*k);
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

void IBackwardParabolicEquation2D::calculateU(std::vector<DoubleMatrix> &psi, double h1, double h2, double ht, unsigned int N1, unsigned int N2, unsigned int M, double a1, double a2) const
{
    //for (unsigned int j=0; j<psi.size(); j++) psi[j].clear();
    psi.clear();
    psi.resize(M+1);

    DoubleMatrix psi0;
    DoubleMatrix psi1;

    psi0.resize(N2+1); for (unsigned int j=0; j<=N2; j++) psi0[j].resize(N1+1);
    psi1.resize(N2+1); for (unsigned int j=0; j<=N2; j++) psi1[j].resize(N1+1);

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

    double x1_alpha1 = -(a1*ht)/(2.0*h1*h1);
    double x1_beta1  = 1.0 + (a1*ht)/(h1*h1);
    double x1_alpha2 = +(a2*ht)/(2.0*h2*h2);
    double x1_beta2  = 1.0 - (a2*ht)/(h2*h2);

    double x2_alpha1 = -(a2*ht)/(2.0*h2*h2);
    double x2_beta1  = 1.0 + (a2*ht)/(h2*h2);
    double x2_alpha2 = +(a1*ht)/(2.0*h1*h1);
    double x2_beta2  = 1.0 - (a1*ht)/(h1*h1);

    for (unsigned int k1=0; k1<=M; k1++)
    {
        unsigned int k = M-k1;

        if (k==M)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned i=0; i<=N1; i++)
                {
                    psi0[j][i] = bfi(i, j);
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
                    da1[i-1] = x1_alpha1;
                    db1[i-1] = x1_beta1;
                    dc1[i-1] = x1_alpha1;
                    dd1[i-1] = x1_alpha2*psi0[j-1][i] + x1_beta2*psi0[j][i] + x1_alpha2*psi0[j+1][i] - (ht/2.0) * bf(i, j, 2*k+1);
                }

                da1[0]     = 0.0;
                dc1[N1-2]  = 0.0;
                dd1[0]    -= x1_alpha1 * bm1(j, 2*k+1);
                dd1[N1-2] -= x1_alpha1 * bm2(j, 2*k+1);

                TomasAlgorithm(da1, db1, dc1, dd1, rx1);

                psi1[j][0]  = bm1(j, 2*k+1);
                for (unsigned int i=1; i<N1; i++)
                {
                    psi1[j][i] = rx1[i-1];
                }
                psi1[j][N1] = bm2(j, 2*k+1);
            }

            for (unsigned int i=0; i<=N1; i++)
            {
                psi1[0][i]  = bm3(i, 2*k+1);
                psi1[N2][i] = bm4(i, 2*k+1);
            }

            for (unsigned int i=1; i<N1; i++)
            {
                for (unsigned int j=1; j<N2; j++)
                {
                    da2[j-1] = x2_alpha1;
                    db2[j-1] = x2_beta1;
                    dc2[j-1] = x2_alpha1;
                    dd2[j-1] = x2_alpha2*psi1[j][i-1] + x2_beta2*psi1[j][i] + x2_alpha2*psi1[j][i+1] - (ht/2.0) * bf(i, j, 2*k);
                }
                da2[0]     = 0.0;
                dc2[N2-2]  = 0.0;
                dd2[0]    -= x2_alpha1 * bm3(i, 2*k);
                dd2[N2-2] -= x2_alpha1 * bm4(i, 2*k);

                TomasAlgorithm(da2, db2, dc2, dd2, rx2);

                psi0[0][i]  = bm3(i, 2*k);
                for (unsigned int j=1; j<N2; j++)
                {
                    psi0[j][i] = rx2[j-1];
                }
                psi0[N2][i] = bm4(i, 2*k);
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                psi0[j][0]  = bm1(j, 2*k);
                psi0[j][N1] = bm2(j, 2*k);
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

    //psi = psi0;

    //    for (unsigned int k=0; k<=M; k++)
    //    {
    //        if (k==0)
    //        {
    //            for (unsigned int j=0; j<=N2; j++)
    //            {
    //                for (unsigned int i=0; i<=N1; i++)
    //                {
    //                    u0[j][i] = fi(i, j);
    //                }
    //            }
    //        }
    //        else
    //        {
    //            // Approximation to x1 direction
    //            for (unsigned int i=1; i<N1; i++)
    //            {
    //                for (unsigned int j=1; j<N2; j++)
    //                {
    //                    da1[j-1] = x1_alpha1;
    //                    db1[j-1] = x1_beta1;
    //                    dc1[j-1] = x1_alpha1;
    //                    dd1[j-1] = x1_alpha2*u0[j][i-1] + x1_beta2*u0[j][i] + x1_alpha2*u0[j][i+1] + (ht/2.0) * f(i, j, 2*k-1);
    //                }

    //                da1[0]     = 0.0;
    //                dc1[N2-2]  = 0.0;
    //                dd1[0]    -= x1_alpha1 * m3(i, 2*k-1);
    //                dd1[N2-2] -= x1_alpha1 * m4(i, 2*k-1);

    //                TomasAlgorithm(da1, db1, dc1, dd1, rx1);

    //                u1[0][i]  = m3(i, 2*k-1);
    //                for (unsigned int j=1; j<N2; j++)
    //                {
    //                    u1[j][i] = rx1[j-1];
    //                }
    //                u1[N2][i] = m4(i, 2*k-1);
    //            }

    //            for (unsigned int j=0; j<=N2; j++)
    //            {
    //                u1[j][0]  = m1(j, 2*k-1);
    //                u1[j][N1] = m2(j, 2*k-1);
    //            }

    //            // Approximation to x2 direction
    //            for (unsigned int j=1; j<N2; j++)
    //            {
    //                for (unsigned int i=1; i<N1; i++)
    //                {
    //                    da2[i-1] = x2_alpha1;
    //                    db2[i-1] = x2_beta1;
    //                    dc2[i-1] = x2_alpha1;
    //                    dd2[i-1] = x2_alpha2*u1[j-1][i] + x2_beta2*u1[j][i] + x2_alpha2*u1[j+1][i] + (ht/2.0) * f(i, j, 2*k);
    //                }
    //                da2[0]     = 0.0;
    //                dc2[N1-2]  = 0.0;
    //                dd2[0]    -= x2_alpha1 * m1(j, 2*k);
    //                dd2[N1-2] -= x2_alpha1 * m2(j, 2*k);

    //                TomasAlgorithm(da2, db2, dc2, dd2, rx2);

    //                u0[j][0]  = m1(j, 2*k);
    //                for (unsigned int i=1; i<N1; i++)
    //                {
    //                    u0[j][i] = rx2[i-1];
    //                }
    //                u0[j][N1] = m2(j, 2*k);
    //            }

    //            for (unsigned int i=0; i<=N1; i++)
    //            {
    //                u0[0][i]  = m3(i, 2*k);
    //                u0[N2][i] = m4(i, 2*k);
    //            }
    //        }
    //    }

    //    da1.clear();
    //    db1.clear();
    //    dc1.clear();
    //    dd1.clear();
    //    rx1.clear();

    //    da2.clear();
    //    db2.clear();
    //    dc2.clear();
    //    dd2.clear();
    //    rx2.clear();

    //    u = u0;
}
