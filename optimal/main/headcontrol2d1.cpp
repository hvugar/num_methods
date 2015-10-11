#include "headcontrol2d1.h"
#include <printer.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <tomasmethod.h>

HeadControl2D1::HeadControl2D1()
{
    t0 = 0.0;
    t1 = 1.0;
    x10 = x20 = 0.0;
    x11 = x21 = 1.0;

    N1 = 100;
    N2 = 100;
    M  = 1000;

    C = (M+1)*(N1+1)*(N2+1);

    h1 = (x11-x10) / N1;
    h2 = (x21-x20) / N2;
    ht = (t1 - t0) / M;

    a1 = a2 = 1.0;

    U.resize(N2+1);
    for (unsigned j=0; j<=N2; j++) U[j].resize(N1+1);
}

void HeadControl2D1::calculateU(const DoubleVector& E, DoubleMatrix& u)
{
    //    puts("*********************************");
    //    puts("Calculating u...");

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
            TomasAlgorithm(a, b, c, d, x);
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
            TomasAlgorithm(a, b, c, d, x);
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
    }

    u = u0;

    for (unsigned int j=0; j<=N2; j++)
    {
        u1[j].clear();
        u0[j].clear();
    }
    u1.clear();
    u0.clear();

    //    puts("Calculated u.");
    //    puts("*********************************");
}

void HeadControl2D1::calculateP(const DoubleVector& E, DoubleVector& g)
{
    //    puts("*********************************");
    //    puts("Calculating psi...");
    DoubleMatrix psi0;
    DoubleMatrix psi1;

    psi0.resize(N2+1);
    psi1.resize(N2+1);

    for (unsigned int j=0; j<N2+1; j++)
    {
        psi0[j].resize(N1+1);
        psi1[j].resize(N1+1);
        for (unsigned int i=0; i<N1+1; i++)
        {
            psi0[j][i] = -2.0*(mu[j][i] - U[j][i]);

            g[M*(N2+1)*(N1+1)+j*(N1+1)+i] = -psi0[j][i] + 2.0*f(i*h1, j*h2, t1);
        }
    }

    ht = -ht;

    DoubleVector a;
    DoubleVector b;
    DoubleVector c;
    DoubleVector d;
    DoubleVector x;

    double alpha1 = 0.0;
    double beta1  = 0.0;
    double alpha2 = 0.0;
    double beta2  = 0.0;

    for (unsigned int p=0; p<M; p++)
    {
        unsigned int k = M-p;

        // Approximation on x direction
        alpha1 = (a1*ht)/(2.0*h1*h1);
        beta1  = 1.0 - (a1*ht)/(h1*h1);
        alpha2 = -(a2*ht)/(2.0*h2*h2);
        beta2  = 1.0 + (a2*ht)/(h2*h2);

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
                d[i-1] = alpha2*psi0[j-1][i] + beta2*psi0[j][i] + alpha2*psi0[j+1][i];
            }
            a[0]     = 0.0;
            c[N1-2]  = 0.0;
            d[0]    -= alpha1 * psi_m1(h2*j, ht*(k-0.5));
            d[N1-2] -= alpha1 * psi_m2(h2*j, ht*(k-0.5));
            TomasAlgorithm(a, b, c, d, x);
            for (unsigned int i=1; i<N1; i++)
            {
                psi1[j][i] = x[i-1];
            }
            psi1[j][0]  = psi_m1(h2*j, ht*(k-0.5));
            psi1[j][N1] = psi_m2(h2*j, ht*(k-0.5));
        }

        for (unsigned int i=0; i<=N1; i++)
        {
            psi1[0][i]  = psi_m3(h1*i, ht*(k-0.5));
            psi1[N2][i] = psi_m4(h1*i, ht*(k-0.5));
        }
        a.clear();
        b.clear();
        c.clear();
        d.clear();
        x.clear();

        ///////////////////////////////////////////////////////////////////////
        // Approximation on y direktion
        alpha1 = (a2*ht)/(2.0*h2*h2);
        beta1  = 1.0 - (a2*ht)/(h2*h2);
        alpha2 = -(a1*ht)/(2.0*h1*h1);
        beta2  = 1.0 + (a1*ht)/(h1*h1);

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
                d[j-1] = alpha2*psi1[j][i-1] + beta2*psi1[j][i] + alpha2*psi1[j][i+1];
            }
            a[0]     = 0.0;
            c[N2-2]  = 0.0;
            d[0]    -= alpha1 * psi_m3(h1*i, ht*(k-1.0));
            d[N2-2] -= alpha1 * psi_m4(h1*i, ht*(k-1.0));
            TomasAlgorithm(a, b, c, d, x);
            for (unsigned int j=1; j<N2; j++)
            {
                psi0[j][i] = x[j-1];
            }
            psi0[0][i]  = psi_m3(h1*i, ht*(k-1.0));
            psi0[N2][i] = psi_m4(h1*i, ht*(k-1.0));
        }
        for (unsigned int j=0; j<=N2; j++)
        {
            psi0[j][0]  = psi_m1(h2*j, ht*(k-1.0));
            psi0[j][N1] = psi_m2(h2*j, ht*(k-1.0));
        }

        a.clear();
        b.clear();
        c.clear();
        d.clear();
        x.clear();

        for (unsigned int j=0; j<N2+1; j++)
        {
            for (unsigned int i=0; i<N1+1; i++)
            {
                g[(k-1)*(N2+1)*(N1+1)+j*(N1+1)+i] = -psi0[j][i] + 2.0*f(i*h1, j*h2, (k-1)*ht);
            }
        }
    }

    ht = -ht;

    for (unsigned int j=0; j<=N2; j++)
    {
        psi1[j].clear();
        psi0[j].clear();
    }

    psi1.clear();
    psi0.clear();

    //    puts("Calculated psi.");
    //    puts("*********************************");
}

double HeadControl2D1::fx(const DoubleVector& E)
{
    //    puts("*********************************");
    //    puts("Calculating integral...");

    //    printf("E 1: %.10f %.10f\n", E[0], E[1]);
    //    printf("E 2: %.10f %.10f\n", E[2], E[3]);
    //    printf("E 3: %.10f %.10f\n", E[4], E[5]);

    //    calculateU(E, mu);
    //    puts("---");
    //    printMatrix1(mu);

    double sum = 0.0;

    for (unsigned int j=0; j<N2; j++)
    {
        for (unsigned int i=0; i<N1; i++)
        {
            int j1 = j;
            int j2 = j+1;
            int i1 = i;
            int i2 = i+1;

            double f1 = mu[j1][i1] - u(i1*h1, j1*h2, 1.0);//U[j1][i1];
            double f2 = mu[j1][i2] - u(i2*h1, j1*h2, 1.0);//U[j1][i2];
            double f3 = mu[j2][i1] - u(i1*h1, j2*h2, 1.0);//U[j2][i1];
            double f4 = mu[j2][i2] - u(i2*h1, j2*h2, 1.0);//U[j2][i2];

            sum = sum + (0.25*(h1*h2))*(f1*f1 + f2*f2 + f3*f3 + f4*f4);
        }
    }

    //    printf("Sum: %.10f\n", sum);

    //    puts("Caculated integral.");
    //    puts("*********************************");
    return sum;
}

void HeadControl2D1::gradient(double step, const DoubleVector& f, DoubleVector& g)
{
    //    puts("*********************************");
    //    puts("Calculating gradient...");

    //    printf("E 1: %.10f %.10f\n", E[0], E[1]);
    //    printf("E 2: %.10f %.10f\n", E[2], E[3]);
    //    printf("E 3: %.10f %.10f\n", E[4], E[5]);
    //    puts("---------------------------------");

    calculateU(f, mu);
    calculateP(f, g);

    //    printf("g 1: %.10f %.10f\n", g[0], g[1]);
    //    printf("g 2: %.10f %.10f\n", g[2], g[3]);
    //    printf("g 3: %.10f %.10f\n", g[4], g[5]);
    //    puts("Calculated gradient...");
    //    puts("*********************************");
}

double HeadControl2D1::u(double x1, double x2, double t)
{
    return x1*x1 + x2*x2 + t*t;
}


double HeadControl2D1::fi(double x1, double x2)
{
    return u(x1, x2, 0.0);
}

double HeadControl2D1::m1(double x2, double t)
{
    return u(0.0, x2, t);
}

double HeadControl2D1::m2(double x2, double t)
{
    return u(1.0, x2, t);
}

double HeadControl2D1::m3(double x1, double t)
{
    return u(x1, 0.0, t);
}

double HeadControl2D1::m4(double x1, double t)
{
    return u(x1, 1.0, t);
}

double HeadControl2D1::psi_fi(int i, int j)
{
    return -2.0*(mu[j][i] - U[j][i]);
}

double HeadControl2D1::psi_m1(double x2, double t)
{
    return 0.0;
}

double HeadControl2D1::psi_m2(double x2, double t)
{
    return 0.0;
}

double HeadControl2D1::psi_m3(double x1, double t)
{
    return 0.0;
}

double HeadControl2D1::psi_m4(double x1, double t)
{
    return 0.0;
}

double HeadControl2D1::f(double x1, double x2, double t)
{
    return 2.0*t - 4.0;
}

void HeadControl2D1::initialize()
{
    puts("Initializing...");
    DoubleVector f;
    //    f.resize(C);

    //    for (unsigned int i=0; i<C; i++)
    //    {
    //        f[i] = 0.1;
    //    }

    //    for (unsigned int j=0; j<=N2; j++)
    //    {
    //        for (unsigned int i=0; i<=N1; i++)
    //        {
    //            U[j][i] = u(i*h1, j*h2, 1.0);
    //        }
    //    }

    calculateU(f, U);
    Printer::printMatrix(U, N2/10, N1/10);

    //    puts("print sum");
    //    fx(f);

    puts("Initialized.\n\n\n");
}

void HeadControl2D1::main()
{
    /* Function */
    HeadControl2D1 hc;
    hc.initialize();

    DoubleVector f0;
    f0.resize(hc.C);

    for (unsigned int i=0; i<f0.size(); i++) f0[i] = 0.0;

    /* Minimization */
    //    SteepestDescentGradient g1;
    //    g1.setFunction(&hc);
    //    g1.setEpsilon1(0.0000001);
    //    g1.setEpsilon2(0.0000001);
    //    g1.setGradientStep(0.000001);
    //    g1.setR1MinimizeEpsilon(0.01, 0.0000001);
    //    g1.setPrinter(new HeadControl2DPrinter);
    //    g1.calculate(E0);

    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setEpsilon1(0.0000001);
    g2.setEpsilon2(0.0000001);
    g2.setGradientStep(0.000001);
    g2.setR1MinimizeEpsilon(0.1, 0.0000001);
    g2.setPrinter(new HeadControl2D1Printer);
    g2.calculate(f0);

    //    printf("E0 1: %.10f %.10f\n", E0[0], E0[1]);
    //    printf("E0 2: %.10f %.10f\n", E0[2], E0[3]);
    //    printf("E0 3: %.10f %.10f\n", E0[4], E0[5]);
}

void HeadControl2D1Printer::print(unsigned int c, const DoubleVector &e, const DoubleVector &s, double alpha, RnFunction *f) const
{
    printf("J: %.16f\n", f->fx(e));
}

