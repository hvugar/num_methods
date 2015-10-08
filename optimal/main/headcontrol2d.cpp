#include "headcontrol2d.h"

void TomasAlgorithm1(const DoubleVector &a, const DoubleVector &b, const DoubleVector &c, const DoubleVector &d, DoubleVector &x);
void printLayer1(const DoubleMatrix& x);

HeadControl2D::HeadControl2D()
{
    t0 = 0.0;
    t1 = 1.0;
    x10 = x20 = 0.0;
    x11 = x21 = 1.0;

    L = 3;
    e1.resize(L);
    e2.resize(L);

    e1[0] = 0.2;
    e1[1] = 0.5;
    e1[2] = 0.8;

    e2[0] = 0.4;
    e2[1] = 0.2;
    e2[2] = 0.6;

    N1 = 100;
    N2 = 100;
    M  = 10000;

    h1 = (x11-x10) / N1;
    h2 = (x21-x20) / N2;
    ht  = (t1 - t0) / M;

    a1 = a2 = 1.0;

    U.resize(N2+1);
    for (unsigned j=0; j<=N2; j++) U[j].resize(N1+1);
}

double HeadControl2D::fx(const DoubleVector& x)
{
    return 0.0;
}

void HeadControl2D::gradient(double step, const DoubleVector& x, DoubleVector& g)
{
}

double HeadControl2D::u(double x1, double x2, double t)
{
    return x1*x1 + x1*x2 + x2*x2*x2 + x1*x2*t + t*t;
}

double HeadControl2D::fi(double x1, double x2)
{
    return u(x1, x2, t0);
}

double HeadControl2D::m1(double x1, double t)
{
    return u(x1, x20, t);
}

double HeadControl2D::m2(double x1, double t)
{
    return u(x1, x21, t);
}

double HeadControl2D::m3(double x2, double t)
{
    return u(x10, x2, t);
}

double HeadControl2D::m4(double x2, double t)
{
    return u(x11, x2, t);
}

double HeadControl2D::f(double x1, double x2, double t)
{
    double sum = 0.0;

    if (fabs(x1 - e1[0]) < 0.000001 && fabs(x2 - e2[0]) < 0.000001)
    {
        sum += t*t;
    }
    if (fabs(x1 - e1[1]) < 0.000001 && fabs(x2 - e2[1]) < 0.000001)
    {
        sum += t;
    }
    if (fabs(x1 - e1[2]) < 0.000001 && fabs(x2 - e2[2]) < 0.000001)
    {
        sum += t*t*t;
    }

    return sum;
}

void HeadControl2D::calculateU()
{
    DoubleMatrix u0;
    DoubleMatrix u1;

    u0.resize(N2+1);
    u1.resize(N2+1);

    for (unsigned int j=0; j<=N2; j++) // x2-->
    {
        u0[j].resize(N1+1);
        u1[j].resize(N1+1);
        for (unsigned int i=0; i<=N1; i++) // x1-->
        {
            u0[j][i] = fi((h1*i), (h2*j));
        }
    }

    DoubleVector a;
    DoubleVector b;
    DoubleVector c;
    DoubleVector d;
    DoubleVector x;

    double alpha1 = 0.0;
    double beta1  = 0.0;
    double alpha2 = 0.0;
    double beta2  = 0.0;

    for (unsigned int k=0; k<M; k++)
    {
        // Approximation on x direction
        alpha1 = -(a1*ht)/(2.0*h1*h1);
        beta1  = 1.0 + (a1*ht)/(h1*h1);
        alpha2 = (a2*ht)/(2.0*h2*h2);
        beta2  = 1.0 - (a2*ht)/(h2*h2);

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
                d[i-1] = alpha2*u0[j-1][i] + beta2*u0[j][i] + alpha2*u0[j+1][i] + (ht/2.0) * f(i*h1, j*h2, k*ht);
            }

            a[0]     = 0.0;
            c[N1-2]  = 0.0;
            d[0]    -= alpha1 * m1(h2*j, ht*(k+0.5));
            d[N1-2] -= alpha1 * m2(h2*j, ht*(k+0.5));
            TomasAlgorithm1(a, b, c, d, x);
            for (unsigned int i=1; i<N1; i++)
            {
                u1[j][i] = x[i-1];
            }
            u1[j][0]  = m1(h2*j, ht*(k+0.5));
            u1[j][N1] = m2(h2*j, ht*(k+0.5));
        }

        for (unsigned int i=0; i<=N1; i++)
        {
            u1[0][i]  = m3(h1*i, ht*(k+0.5));
            u1[N2][i] = m4(h1*i, ht*(k+0.5));
        }

        a.clear();
        b.clear();
        c.clear();
        d.clear();
        x.clear();

        //        printf("\nk = %f\n", k+0.5);
        //        printLayer(u1);

        ///////////////////////////////////////////////////////////////////////
        // Approximation on y direktion
        alpha1 = -(a2*ht)/(2.0*h2*h2);
        beta1  = 1.0 + (a2*ht)/(h2*h2);
        alpha2 = (a1*ht)/(2.0*h1*h1);
        beta2  = 1.0 - (a1*ht)/(h1*h1);

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
                d[j-1] = alpha2*u1[j][i-1] + beta2*u1[j][i] + alpha2*u1[j][i+1] + (ht/2.0) * f(i*h1, j*h2, k*ht);
            }
            a[0]     = 0.0;
            c[N2-2]  = 0.0;
            d[0]    -= alpha1 * m3(h1*i, ht*(k+1.0));
            d[N2-2] -= alpha1 * m4(h1*i, ht*(k+1.0));
            TomasAlgorithm1(a, b, c, d, x);
            for (unsigned int j=1; j<N2; j++)
            {
                u0[j][i] = x[j-1];
            }
            u0[0][i]  = m3(h1*i, ht*(k+1));
            u0[N2][i] = m4(h1*i, ht*(k+1));
        }

        for (unsigned int j=0; j<=N2; j++)
        {
            u0[j][0]  = m1(h2*j, ht*(k+1));
            u0[j][N1] = m2(h2*j, ht*(k+1));
        }

        a.clear();
        b.clear();
        c.clear();
        d.clear();
        x.clear();

        //        printf("k = %d\n", k);
        //        printLayer(u0);
    }

    U = u0;
    printLayer1(U);

    for (unsigned int j=0; j<=N2; j++)
    {
        u1[j].clear();
        u0[j].clear();
    }
    u1.clear();
    u0.clear();
}

void TomasAlgorithm1(const DoubleVector &a, const DoubleVector &b, const DoubleVector &c, const DoubleVector &d, DoubleVector &x)
{
    if (x.size() != a.size() || x.size() != b.size() || x.size() != c.size() || x.size() != d.size())
        return;

    int n = x.size();
    DoubleVector p(n);
    DoubleVector q(n);

    for (int i=0; i<n; i++)
    {
        if (i==0)
        {
            p[0] = d[0]/b[0];
            q[0] = -c[0]/b[0];
        } else
            if(i==n-1)
            {
                p[n-1] = (d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
                q[n-1] = 0.0;
            }
            else
            {
                p[i] = (d[i]-a[i]*p[i-1])/(b[i]+a[i]*q[i-1]);
                q[i] = -c[i]/(b[i]+a[i]*q[i-1]);
            }
    }

    for (int i=n-1; i>=0; i--)
    {
        if (i==n-1)
        {
            x[i] = p[i];
        }
        else
        {
            x[i] = p[i] + q[i]*x[i+1];
        }
    }
}

void printLayer1(const DoubleMatrix& x)
{
    int m = x.size()/10;
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
