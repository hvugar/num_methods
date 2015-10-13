#include "heatcontrol2delta.h"
#include <tomasmethod.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <math.h>
#include <stdio.h>

HeatControl2Delta::HeatControl2Delta(unsigned int m, unsigned int n2, unsigned int n1) : RnFunction()
{
    M = m;
    N2 = n2;
    N1 = n1;
    C = (M+1)*(N2+1)*(N1+1);
    L = 3;

    t0 = 0.0;
    t1 = 1.0;
    x10 = 0.0;
    x11 = 1.0;
    x20 = 0.0;
    x21 = 1.0;

    a1 = a2 = 1.0;

    ht = (t1-t0)/M;
    h1 = (x11-x10)/N1;
    h2 = (x21-x20)/N2;

    initialize();
}

HeatControl2Delta::~HeatControl2Delta()
{
}

double HeatControl2Delta::fx(const DoubleVector &e)
{
//    puts("--------------------");
//    printf("E[0]: %12.6f E[1]: %12.6f\n", e[0], e[1]);
//    printf("E[2]: %12.6f E[3]: %12.6f\n", e[2], e[3]);
//    printf("E[4]: %12.6f E[5]: %12.6f\n", e[4], e[5]);
//    puts("--------------------");

    calculateU(e, uT);

    double sum = 0.0;

    for (unsigned int j=0; j<N2; j++)
    {
        for (unsigned int i=0; i<N1; i++)
        {
            int j1 = j;
            int j2 = j+1;
            int i1 = i;
            int i2 = i+1;

            double f1 = uT[j1][i1] - U[j1][i1];
            double f2 = uT[j1][i2] - U[j1][i2];
            double f3 = uT[j2][i1] - U[j2][i1];
            double f4 = uT[j2][i2] - U[j2][i2];

            sum += (f1*f1 + f2*f2 + f3*f3 + f4*f4);
        }
    }

    sum = (0.25*(h1*h2))*sum;

    double norm = 0.0;

//    for (unsigned int k=0; k<M; k++)
//    {
//        for (unsigned int j=0; j<N2; j++)
//        {
//            for (unsigned int i=0; i<N1; i++)
//            {
//                unsigned int i1 = i;
//                unsigned int i2 = i+1;
//                unsigned int j1 = j;
//                unsigned int j2 = j+1;
//                unsigned int k1 = k;
//                unsigned int k2 = k+1;

//                double f1 = f[k1*(N2+1)*(N1+1) + j1*(N1+1) + i1] - fxt(i1*h1, j1*h2, k1*ht);
//                double f2 = f[k1*(N2+1)*(N1+1) + j1*(N1+1) + i2] - fxt(i2*h1, j1*h2, k1*ht);
//                double f3 = f[k1*(N2+1)*(N1+1) + j2*(N1+1) + i1] - fxt(i1*h1, j2*h2, k1*ht);
//                double f4 = f[k1*(N2+1)*(N1+1) + j2*(N1+1) + i2] - fxt(i2*h1, j2*h2, k1*ht);
//                double f5 = f[k2*(N2+1)*(N1+1) + j1*(N1+1) + i1] - fxt(i1*h1, j1*h2, k2*ht);
//                double f6 = f[k2*(N2+1)*(N1+1) + j1*(N1+1) + i2] - fxt(i2*h1, j1*h2, k2*ht);
//                double f7 = f[k2*(N2+1)*(N1+1) + j2*(N1+1) + i1] - fxt(i1*h1, j2*h2, k2*ht);
//                double f8 = f[k2*(N2+1)*(N1+1) + j2*(N1+1) + i2] - fxt(i2*h1, j2*h2, k2*ht);

//                norm += f1*f1 + f2*f2 + f3*f3 + f4*f4 + f5*f5 + f6*f6 + f7*f7 + f8*f8;
//            }
//        }
//    }
//    norm = norm * (h1*h2*ht)*0.125;

    return sum + norm;
}

void HeatControl2Delta::gradient(const DoubleVector &e, DoubleVector &g, double gradient_step)
{
    calculateU(e, uT);
    calculateP(e, g);
}

void HeatControl2Delta::calculateU(const DoubleVector &e, DoubleMatrix& u)
{
    DoubleMatrix u0;
    DoubleMatrix u1;

    u0.resize(N2+1); for (unsigned int j=0; j<=N2; j++) u0[j].resize(N1+1);
    u1.resize(N2+1); for (unsigned int j=0; j<=N2; j++) u1[j].resize(N1+1);

    double alpha1 = 0;
    double beta1  = 0;
    double alpha2 = 0;
    double beta2  = 0;

    DoubleVector a;
    DoubleVector b;
    DoubleVector c;
    DoubleVector d;
    DoubleVector x;

    for (unsigned int k=0; k<=M; k++)
    {
        if (k==0)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned i=0; i<=N1; i++)
                {
                    u0[j][i] = fi(i*h1, j*h2);
                }
            }
        }
        else
        {
            // Approximation to x1 direction
            alpha1 = -(a2*ht)/(2.0*h2*h2);
            beta1  = 1.0 + (a2*ht)/(h2*h2);
            alpha2 = (a1*ht)/(2.0*h1*h1);
            beta2  = 1.0 - (a1*ht)/(h1*h1);

            a.resize(N1-1);
            b.resize(N1-1);
            c.resize(N1-1);
            d.resize(N1-1);
            x.resize(N1-1);

            for (unsigned int i=1; i<N1; i++)
            {
                for (unsigned int j=1; j<N2; j++)
                {
                    a[j-1] = alpha1;
                    b[j-1] = beta1;
                    c[j-1] = alpha1;
                    d[j-1] = alpha2*u0[j][i-1] + beta2*u0[j][i] + alpha2*u0[j][i+1] + (ht/2.0) * /*f[k*(N1+1)*(N2+1)+i*(N2+1)+j]*/fxt(i*h1, j*h2, (k-1.0)*ht, e);
                }

                a[0]     = 0.0;
                c[N2-2]  = 0.0;
                d[0]    -= alpha1 * m3(h1*i, ht*(k));
                d[N2-2] -= alpha1 * m4(h1*i, ht*(k));

                TomasAlgorithm(a, b, c, d, x);

                u1[0][i]  = m3(h1*i, ht*(k-0.5));
                for (unsigned int j=1; j<N2; j++)
                {
                    u1[j][i] = x[j-1];
                }
                u1[N2][i] = m4(h1*i, ht*(k-0.5));
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                u1[j][0]  = m1(h2*j, ht*(k-0.5));
                u1[j][N2] = m2(h2*j, ht*(k-0.5));
            }

            a.clear();
            b.clear();
            c.clear();
            d.clear();
            x.clear();

            // Approximation to x2 direction
            alpha1 = -(a1*ht)/(2.0*h1*h1);
            beta1  = 1.0 + (a1*ht)/(h1*h1);
            alpha2 = (a2*ht)/(2.0*h2*h2);
            beta2  = 1.0 - (a2*ht)/(h2*h2);

            a.resize(N2-1);
            b.resize(N2-1);
            c.resize(N2-1);
            d.resize(N2-1);
            x.resize(N2-1);

            for (unsigned int j=1; j<N2; j++)
            {
                for (unsigned int i=1; i<N1; i++)
                {
                    a[i-1] = alpha1;
                    b[i-1] = beta1;
                    c[i-1] = alpha1;
                    d[i-1] = alpha2*u1[j-1][i] + beta2*u1[j][i] + alpha2*u1[j+1][i] + (ht/2.0) * /*f[k*(N1+1)*(N2+1)+j*(N1+1)+i]*/fxt(i*h1, j*h2, (k-1)*ht, e);
                }
                a[0]     = 0.0;
                c[N1-2]  = 0.0;
                d[0]    -= alpha1 * m1(h2*j, ht*(k));
                d[N1-2] -= alpha1 * m2(h2*j, ht*(k));
                TomasAlgorithm(a, b, c, d, x);

                u0[j][0]  = m1(h2*j, ht*(k));
                for (unsigned int i=1; i<N1; i++)
                {
                    u0[j][i] = x[i-1];
                }
                u0[j][N1] = m2(h2*j, ht*(k));
            }

            for (unsigned int i=0; i<=N1; i++)
            {
                u0[0][i]  = m3(h1*i, ht*(k));
                u0[N1][i] = m4(h1*i, ht*(k));
            }

            a.clear();
            b.clear();
            c.clear();
            d.clear();
            x.clear();
        }
    }

    u = u0;
}

void HeatControl2Delta::calculateP(const DoubleVector &e, DoubleVector &g)
{
    DoubleMatrix psi0;
    DoubleMatrix psi1;

    psi0.resize(N2+1); for (unsigned int j=0; j<=N2; j++) psi0[j].resize(N1+1);
    psi1.resize(N2+1); for (unsigned int j=0; j<=N2; j++) psi1[j].resize(N1+1);

    double alpha1 = 0;
    double beta1  = 0;
    double alpha2 = 0;
    double beta2  = 0;

    DoubleVector a;
    DoubleVector b;
    DoubleVector c;
    DoubleVector d;
    DoubleVector x;

    g[0] = g[1] = g[2] = g[3] = g[4] = g[5] = 0.0;

    for (unsigned int k1=0; k1<=M; k1++)
    {
        unsigned int k = M-k1;

        if (k==M)
        {
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned i=0; i<=N1; i++)
                {
                    psi0[j][i] = -2.0*(uT[j][i] - U[j][i]);
                }
            }
        }
        else
        {
            // Approximation to x1 direction
            alpha1 = -(a1*ht)/(2.0*h1*h1);
            beta1  = 1.0 + (a1*ht)/(h1*h1);
            alpha2 = +(a2*ht)/(2.0*h2*h2);
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
                    d[i-1] = alpha2*psi0[j-1][i] + beta2*psi0[j][i] + alpha2*psi0[j+1][i];
                }

                a[0]     = 0.0;
                c[N1-2]  = 0.0;
                d[0]    -= alpha1 * pm1(h2*j, ht*(k+0.5));
                d[N1-2] -= alpha1 * pm2(h2*j, ht*(k+0.5));

                TomasAlgorithm(a, b, c, d, x);

                psi1[j][0]  = pm1(h2*j, ht*(k+0.5));
                for (unsigned int i=1; i<N1; i++)
                {
                    psi1[j][i] = x[i-1];
                }
                psi1[j][N1] = pm2(h2*j, ht*(k+0.5));
            }

            for (unsigned int i=0; i<=N1; i++)
            {
                psi1[0][i]  = pm3(h1*i, ht*(k+0.5));
                psi1[N2][i] = pm4(h1*i, ht*(k+0.5));
            }

            a.clear();
            b.clear();
            c.clear();
            d.clear();
            x.clear();

            // Approximation to x2 direction
            alpha1 = -(a2*ht)/(2.0*h2*h2);
            beta1  = 1.0 + (a2*ht)/(h2*h2);
            alpha2 = +(a1*ht)/(2.0*h1*h1);
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
                    d[j-1] = alpha2*psi1[j][i-1] + beta2*psi1[j][i] + alpha2*psi1[j][i+1];
                }
                a[0]     = 0.0;
                c[N2-2]  = 0.0;
                d[0]    -= alpha1 * pm3(h1*i, ht*k);
                d[N2-2] -= alpha1 * pm4(h1*i, ht*k);
                TomasAlgorithm(a, b, c, d, x);

                psi0[0][i]  = pm3(h1*i, ht*k);
                for (unsigned int j=1; j<N2; j++)
                {
                    psi0[j][i] = x[j-1];
                }
                psi0[N2][i] = pm4(h1*i, ht*k);
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                psi0[j][0]  = pm1(h2*j, ht*k);
                psi0[j][N1] = pm2(h2*j, ht*k);
            }

            a.clear();
            b.clear();
            c.clear();
            d.clear();
            x.clear();
        }

        // calculating gradient
        if (k==M || k==0)
        {
            unsigned int i;
            unsigned int j;

            i = (unsigned int)round(e[0]/h1);
            j = (unsigned int)round(e[1]/h2);
            g[0] += f1(k*ht) * ((psi0[j][i+1]-psi0[j][i-1])/(2.0*h1));
            g[1] += f1(k*ht) * ((psi0[j+1][i]-psi0[j-1][i])/(2.0*h2));

            i = (unsigned int)round(e[2]/h1);
            j = (unsigned int)round(e[3]/h2);
            g[2] += f2(k*ht) * ((psi0[j][i+1]-psi0[j][i-1])/(2.0*h1));
            g[3] += f2(k*ht) * ((psi0[j+1][i]-psi0[j-1][i])/(2.0*h2));

            i = (unsigned int)round(e[4]/h1);
            j = (unsigned int)round(e[5]/h2);
            g[4] += f3(k*ht) * ((psi0[j][i+1]-psi0[j][i-1])/(2.0*h1));
            g[5] += f3(k*ht) * ((psi0[j+1][i]-psi0[j-1][i])/(2.0*h2));
        }
        else
        {
            unsigned int i;
            unsigned int j;

            i = (unsigned int)round(e[0]/h1);
            j = (unsigned int)round(e[1]/h2);
            g[0] += 2.0*(f1(k*ht) * ((psi0[j][i+1]-psi0[j][i-1])/(2.0*h1)));
            g[1] += 2.0*(f1(k*ht) * ((psi0[j+1][i]-psi0[j-1][i])/(2.0*h2)));

            i = (unsigned int)round(e[2]/h1);
            j = (unsigned int)round(e[3]/h2);
            g[2] += 2.0*(f2(k*ht) * ((psi0[j][i+1]-psi0[j][i-1])/(2.0*h1)));
            g[3] += 2.0*(f2(k*ht) * ((psi0[j+1][i]-psi0[j-1][i])/(2.0*h2)));

            i = (unsigned int)round(e[4]/h1);
            j = (unsigned int)round(e[5]/h2);
            g[4] += 2.0*(f3(k*ht) * ((psi0[j][i+1]-psi0[j][i-1])/(2.0*h1)));
            g[5] += 2.0*(f3(k*ht) * ((psi0[j+1][i]-psi0[j-1][i])/(2.0*h2)));
        }
    }

    g[0] = -(ht/2.0) * g[0];
    g[1] = -(ht/2.0) * g[1];
    g[2] = -(ht/2.0) * g[2];
    g[3] = -(ht/2.0) * g[3];
    g[4] = -(ht/2.0) * g[4];
    g[5] = -(ht/2.0) * g[5];
}

double HeatControl2Delta::fxt(double x1, double x2, double t, const DoubleVector &e)
{
    double sum = 0.0;
    double epsilon = 0.00000001;
    if (fabs(x1-e[0]) < epsilon && fabs(x2 - e[1]) < epsilon) { sum += f1(t); }
    if (fabs(x1-e[2]) < epsilon && fabs(x2 - e[3]) < epsilon) { sum += f2(t); }
    if (fabs(x1-e[4]) < epsilon && fabs(x2 - e[5]) < epsilon) { sum += f3(t); }
    return sum;
}

void HeatControl2Delta::main()
{
    HeatControl2Delta hc(10, 100, 100);

    DoubleVector e;
    e.resize(2*hc.L);

//    DoubleVector g;
//    g.resize(e.size());
//    for (unsigned int i=0; i<e.size(); i++) g[i] = 0.5;

//    e[0] = 0.2; e[1] = 0.5;
//    e[2] = 0.5; e[3] = 0.6;
//    e[4] = 0.8; e[5] = 0.7;
//    printf("Result: %12.8f\n", hc.fx(e));

    for (unsigned int i=0; i<e.size(); i++) e[i] = 0.511;

    printf("Result: %.16f\n", hc.fx(e));

//    hc.gradient(e, g, 0.0);

//    puts("+++++++++++++++++++++++++++++++++++++++++");
//    printf("e[0]: %16.8f e[1]: %16.8f\n", e[0], e[1]);
//    printf("e[2]: %16.8f e[3]: %16.8f\n", e[2], e[3]);
//    printf("e[4]: %16.8f e[5]: %16.8f\n", e[4], e[5]);
//    printf("g[0]: %16.8f g[1]: %16.8f\n", g[0], g[1]);
//    printf("g[2]: %16.8f g[3]: %16.8f\n", g[2], g[3]);
//    printf("g[4]: %16.8f g[5]: %16.8f\n", g[4], g[5]);
//    puts("+++++++++++++++++++++++++++++++++++++++++");

    /* Minimization */
//        SteepestDescentGradient g1;
//        g1.setFunction(&hc);
//        g1.setEpsilon1(0.0000001);
//        g1.setEpsilon2(0.0000001);
//        g1.setGradientStep(0.000001);
//        g1.setR1MinimizeEpsilon(0.0001, 0.0000001);
//        g1.setNormalize(true);
//        g1.setPrinter(new HeatControl2DPrinter);
//        g1.calculate(f);

    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setEpsilon1(0.000000001);
    g2.setEpsilon2(0.000000001);
    g2.setGradientStep(0.000001);
    g2.setR1MinimizeEpsilon(0.01, 0.0000001);
    g2.setPrinter(new HeatControl2DeltaPrinter);
    g2.setNormalize(false);
    g2.calculate(e);

    printf("e[0]: %16.8f e[1]: %16.8f\n", e[0], e[1]);
    printf("e[2]: %16.8f e[3]: %16.8f\n", e[2], e[3]);
    printf("e[4]: %16.8f e[5]: %16.8f\n", e[4], e[5]);

//    Printer::printMatrix(hc.U, hc.N2/10, hc.N1/10);
}

void HeatControl2Delta::initialize()
{
    DoubleVector E;
    E.resize(2*L);
    E[0] = 0.2; E[1] = 0.5;
    E[2] = 0.5; E[3] = 0.6;
    E[4] = 0.8; E[5] = 0.7;

    U.resize(N2+1);
    for (unsigned int j=0; j<=N2; j++)
    {
        U[j].resize(N1+1);
    }

    calculateU(E, U);

    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
    Printer::printMatrix(U, N2/10, N1/10);
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");

//    FILE *f = fopen("e:\\U.txt", "w");
//    for (unsigned int j=0; j<=N2; j++)
//    {
//        for (unsigned int i=0; i<=N1; i++)
//        {
//            fprintf(f, "%12.8f ", U[j][i]);
//        }
//        fputs("\n", f);
//        fflush(f);
//    }
//    fclose(f);
}

void HeatControl2DeltaPrinter::print(unsigned int i, const DoubleVector &e, const DoubleVector &s, double a, RnFunction *f) const
{
    printf("e[0]: %16.8f e[1]: %16.8f\n", e[0], e[1]);
    printf("e[2]: %16.8f e[3]: %16.8f\n", e[2], e[3]);
    printf("e[4]: %16.8f e[5]: %16.8f\n", e[4], e[5]);
    printf("J: %.16f\n", f->fx(e));
}
