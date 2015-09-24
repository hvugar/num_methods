#include "heat2d.h"

void printLayer(const DoubleMatrix& x)
{
    unsigned int m = x.size()/10;
    for (unsigned int j=0; j<x.size(); j++)
    {
        if (j%m==0)
        {
            for (unsigned int i=0; i<x[j].size(); i++)
            {
                if (i%m==0) printf("%14.10f ", x[j][i]);
            }
            puts("");
        }
    }
}

Heat2DControl::Heat2DControl()
{
    N1 = 1000;
    N2 = 1000;
    M = 1000;

    t0 = 0.0;
    t1 = 1.0;

    x11 = 0.0;
    x12 = 1.0;

    x21 = 0.0;
    x22 = 1.0;

    h1 = (x12 - x11)/N1;
    h2 = (x22 - x21)/N2;
    dt = (t1 - t0)/M;

    a1 = 1.0;
    a2 = 1.0;

    mu.resize(N2+1);
    u0.resize(N2+1);
    mf.resize(N2+1);
    mp.resize(N2+1);
    mx1.resize(N2+1);
    mx2.resize(N2+1);
    unsigned int i=0;
    for (i=0; i<N2+1; i++)
    {
        mu[i].resize(N1+1);
        u0[i].resize(N1+1);
        mf[i].resize(N1+1);
        mp[i].resize(N1+1);
        mx1[i].resize(N1+1);
        mx2[i].resize(N1+1);
    }
}

Heat2DControl::~Heat2DControl()
{

}

double Heat2DControl::fx(const DoubleVector &x)
{
    return 0.0;
}

void Heat2DControl::gradient(double step, const DoubleVector &x, DoubleVector &g)
{

}

double Heat2DControl::u(double x1, double x2, double t) { return x1*x1 + x2*x2 + t*t; }
double Heat2DControl::U(double x1, double x2) { return x1*x1 + x2*x2 + 1.0; }

double Heat2DControl::f(double x1, double x2, double t) { return 2.0*t - 4.0; }
double Heat2DControl::fi(double x1, double x2) { return x1*x1 + x2*x2; }
double Heat2DControl::m1(double x2, double t) { return x2*x2 + t*t; }
double Heat2DControl::m2(double x2, double t) { return x2*x2 + t*t + 1.0; }
double Heat2DControl::m3(double x1, double t) { return x1*x1 + t*t; }
double Heat2DControl::m4(double x1, double t) { return x1*x1 + t*t + 1.0; }

double Heat2DControl::psi_fi(double x1, double x2) { return 2*u(x1, x2, t1) - U(x1, x2); }
double Heat2DControl::psi_m1(double x2, double t) { return 0.0; }
double Heat2DControl::psi_m2(double x2, double t) { return 0.0; }
double Heat2DControl::psi_m3(double x1, double t) { return 0.0; }
double Heat2DControl::psi_m4(double x1, double t) { return 0.0; }

void Heat2DControl::calculateX(const DoubleMatrix &u, const DoubleMatrix &f, DoubleMatrix &x1, DoubleMatrix &x2)
{
    unsigned int i=0;

    DoubleMatrix u0;
    DoubleMatrix u1;
    DoubleMatrix u2;

    u0.resize(N2+1);
    u1.resize(N2+1);
    u2.resize(N2+1);

    for (i=0; i<N2+1; i++)
    {
        u0[i].resize(N1+1);
        u1[i].resize(N1+1);
        u2[i].resize(N1+1);
    }

    for (unsigned int j=0; j<N2+1; j++)
    {
        for (unsigned int i=0; i<N1+1; i++)
        {
            u0[j][i] = fi((h1*i), (h2*j));
        }
    }

    printLayer(u0);

    DoubleVector a;
    DoubleVector b;
    DoubleVector c;
    DoubleVector d;
    DoubleVector x;

    double alpha1 = 0.0;
    double beta1  = 0.0;
    double alpha2 = 0.0;
    double beta2  = 0.0;

    unsigned int K = 2*M;

    for (unsigned int k=0; k<K; k++)
    {
        // Approximation on x direction
        alpha1 = -(a1*dt)/(2.0*h1*h1);
        beta1  = 1.0 + (a1*dt)/(h1*h1);
        alpha2 = (a2*dt)/(2.0*h2*h2);
        beta2  = 1.0 - (a2*dt)/(h2*h2);
        a.resize(N1-1);
        b.resize(N1-1);
        c.resize(N1-1);
        d.resize(N1-1);
        x.resize(N1-1);

        for (unsigned int j=1; j<N2; j++)
        {
            for (unsigned int i=1; i<N1; i++)
            {
                a[i-1] = alpha1;
                b[i-1] = beta1;
                c[i-1] = alpha1;
                d[i-1] = alpha2*u0[j-1][i] + beta2*u0[j][i] + alpha2*u0[j+1][i] + (dt/2.0)*f(i, j);
            }
            a[0]   = 0.0;
            c[N1-2] = 0.0;
            d[0]    -= alpha1 * m1(h2*j, dt*(k+0.5));
            d[N1-2] -= alpha1 * m2(h2*j, dt*(k+0.5));
            //TomasAlgorithm(a, b, c, d, x);
            for (unsigned int i=1; i<N1; i++)
            {
                u1[j][i] = x[i-1];
            }
            u1[j][0]  = m1(h2*j, dt*(k+0.5));
            u1[j][N1] = m2(h2*j, dt*(k+0.5));
        }

        for (unsigned int i=0; i<=N1; i++)
        {
            u1[0][i]  = m3(h1*i, dt*(k+0.5));
            u1[N2][i] = m4(h1*i, dt*(k+0.5));
        }
        a.clear();
        b.clear();
        c.clear();
        d.clear();
        x.clear();

        printf("\nk = %d\n", 2*k+1);
        printLayer(u1);

        ///////////////////////////////////////////////////////////////////////
        // Approximation on y direktion
        alpha1 = -(a2*dt)/(2.0*h2*h2);
        beta1  = 1.0 + (a2*dt)/(h2*h2);
        alpha2 = (a1*dt)/(2.0*h1*h1);
        beta2  = 1.0 - (a1*dt)/(h1*h1);

        a.resize(N2-1);
        b.resize(N2-1);
        c.resize(N2-1);
        d.resize(N2-1);
        x.resize(N2-1);
        for (unsigned int i=1; i<N1; i++)
        {
            for (unsigned int j=1; j<N2; j++)
            {
                a[j-1] = alpha1;
                b[j-1] = beta1;
                c[j-1] = alpha1;
                d[j-1] = alpha2*u1[j][i-1] + beta2*u1[j][i] + alpha2*u1[j][i+1] + (dt/2.0);//*_f(i*h1, j*h2, (k+0.5)*dt);
            }
            a[0]     = 0.0;
            c[N2-2]  = 0.0;
            d[0]    -= alpha1  * m3(h1*i, dt*(k+1.0));
            d[N2-2] -= alpha1 * m4(h1*i, dt*(k+1.0));
            //TomasAlgorithm(a, b, c, d, x);
            for (unsigned int j=1; j<N2; j++)
            {
                u1[j][i] = x[j-1];
            }
            u2[0][i]  = m3(h1*i, dt*(k+1));
            u2[N2][i] = m4(h1*i, dt*(k+1));
        }
        for (unsigned int j=0; j<=N2; j++)
        {
            u2[j][0]  = m1(h2*j, dt*(k+1));
            u2[j][N1] = m2(h2*j, dt*(k+1));
        }
        a.clear();
        b.clear();
        c.clear();
        d.clear();
        x.clear();

        printf("k = %d\n", 2*k+2);
        printLayer(u2);
    }
    printf("end\n");
}

