#include "hyperbolicequation.h"

void IHyperbolicEquation::calculateU(DoubleVector &u, double hx, double ht, unsigned int M, unsigned int N, double a, double lamda) const
{
    u.clear();
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
                u0[i] = fi1(i);
                u1[i] = u0[i] + ht*fi2(i);
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

            u[0] = m1(j+1);
            u[N] = m2(j+1);

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
    for (unsigned int i=0; i<u.size(); i++) u[i].clear();
    u.clear();

    u.resize(M+1);
    for (unsigned int i=0; i<u.size(); i++) u[i].resize(N+1);

    DoubleVector da;
    DoubleVector db;
    DoubleVector dc;
    DoubleVector rd;
    DoubleVector rx;

    da.resize(N-1);
    db.resize(N-1);
    dc.resize(N-1);
    rd.resize(N-1);
    rx.resize(N-1);

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
                u[0][i] = fi1(i);
                u[1][i] = u[0][i] + ht*fi2(i);
            }
        }
        else
        {
            u[j+1][0] = m1(j+1);
            u[j+1][N] = m2(j+1);

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
    p.clear();
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
                p0[i] = bfi1(i);
                p1[i] = p0[i] - ht*bfi2(i);
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

            p[0] = bm1(j-1);
            p[N] = bm2(j-1);

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
    for (unsigned int i=0; i<p.size(); i++) p[i].clear();
    p.clear();

    p.resize(M+1);
    for (unsigned int i=0; i<p.size(); i++) p[i].resize(N+1);

    DoubleVector da;
    DoubleVector db;
    DoubleVector dc;
    DoubleVector rd;
    DoubleVector rx;

    da.resize(N-1);
    db.resize(N-1);
    dc.resize(N-1);
    rd.resize(N-1);
    rx.resize(N-1);

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
                p[M][i] = bfi1(i);
                p[M-1][i] = p[M][i] - ht*bfi2(i);
            }
        }
        else
        {
            p[j-1][0] = bm1(j-1);
            p[j-1][N] = bm2(j-1);

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
