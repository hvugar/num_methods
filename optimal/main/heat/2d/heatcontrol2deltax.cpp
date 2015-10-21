#include "heatcontrol2deltax.h"
#include <tomasmethod.h>
#include <gradient_cjt.h>
#include <math.h>
#include <stdlib.h>

void HeatControl2DeltaX::main()
{
    HeatControl2DeltaX hc(100, 100, 100);

    DoubleVector e;
    e.resize(2*hc.L);

    //Optimal
    //e[0] = 0.7; e[1] = 0.2;
    //e[2] = 0.5; e[3] = 0.8;
    //e[4] = 0.2; e[5] = 0.3;

    e[0] = 0.75; e[1] = 0.25;
    e[2] = 0.55; e[3] = 0.85;
    e[4] = 0.25; e[5] = 0.35;

    //e[0] = 0.65; e[1] = 0.15;
    //e[2] = 0.45; e[3] = 0.75;
    //e[4] = 0.15; e[5] = 0.25;

    //e[0] = 0.5; e[1] = 0.6;
    //e[2] = 0.8; e[3] = 0.2;
    //e[4] = 0.4; e[5] = 0.3;

    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setEpsilon1(0.000000001);
    g2.setEpsilon2(0.000000001);
    g2.setGradientStep(0.000001);
    g2.setR1MinimizeEpsilon(0.1, 0.001);
    g2.setPrinter(&hc);
    g2.setProjection(&hc);
    g2.setNormalize(true);
    g2.calculate(e);
}

HeatControl2DeltaX::HeatControl2DeltaX(unsigned int M, unsigned int N2, unsigned int N1)
{
    alpha = 1.0;

    this->M  = M;
    this->N2 = N2;
    this->N1 = N1;
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

    U.resize(N2+1);
    for (unsigned int j=0; j<=N2; j++) U[j].resize(N1+1);

    initialize();
}

double HeatControl2DeltaX::fx(const DoubleVector& e)
{
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

            double f11 = uT[j1][i1] - U[j1][i1];
            double f12 = uT[j1][i2] - U[j1][i2];
            double f21 = uT[j2][i1] - U[j2][i1];
            double f22 = uT[j2][i2] - U[j2][i2];

            sum += (f11*f11 + f12*f12 + f21*f21 + f22*f22);
        }
    }
    sum = (0.25*(h1*h2))*sum;

    //    double p;
    //    double norm = 0.0;
    //    for (unsigned int l=0; l<L; l++)
    //    {
    //        for (unsigned int k=0; k<=M; k++)
    //        {
    //            if (k==0 || k==M) p = 1.0; else p = 2.0;

    //            switch(l)
    //            {
    //            case 0: { norm += p*(ht/2.0)*(f[0*(M+1)+k] - f1(k*ht))*(f[0*(M+1)+k] - f1(k*ht)); } break;
    //            case 1: { norm += p*(ht/2.0)*(f[1*(M+1)+k] - f2(k*ht))*(f[1*(M+1)+k] - f2(k*ht)); } break;
    //            case 2: { norm += p*(ht/2.0)*(f[2*(M+1)+k] - f3(k*ht))*(f[2*(M+1)+k] - f3(k*ht)); } break;
    //            }
    //        }
    //    }

    return sum;
}

void HeatControl2DeltaX::gradient(const DoubleVector& e, DoubleVector& g, double gradient_step)
{
    calculateU(e, uT);
    calculateP(e, g);
    puts("-----------------------------------------------------------");
    printf("e1: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", e[0], e[1], e[2], e[3], e[4], e[5]);
    printf("g1: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", g[0], g[1], g[2], g[3], g[4], g[5]);
    //calculateG2(e, g);
}

void HeatControl2DeltaX::calculateU(const DoubleVector &e, DoubleMatrix& u)
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
                    d[j-1] = alpha2*u0[j][i-1] + beta2*u0[j][i] + alpha2*u0[j][i+1] + (ht/2.0) * fxt(i, j, k, e);
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
                u1[j][N1] = m2(h2*j, ht*(k-0.5));
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
                    d[i-1] = alpha2*u1[j-1][i] + beta2*u1[j][i] + alpha2*u1[j+1][i] + (ht/2.0) * fxt(i, j, k, e);
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
                u0[N2][i] = m4(h1*i, ht*(k));
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

void HeatControl2DeltaX::calculateP(const DoubleVector &e, DoubleVector &g)
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

    for (unsigned int i=0; i<g.size(); i++) g[i] = 0.0;

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

        calculateGX(e, psi0, g, k);
        calculateGF(e, g);
    }
}

void HeatControl2DeltaX::calculateGX(const DoubleVector& e, const DoubleMatrix& psi, DoubleVector& g, unsigned int k)
{
    unsigned int i,j;
    if (k==0 || k==M)
    {
        i = (unsigned int)round(e[0]/h1);
        j = (unsigned int)round(e[1]/h2);
        g[0] = g[0] + f1((k)*ht) * (psi[j][i+1] - psi[j][i-1])/(2.0*h1);
        g[1] = g[1] + f1((k)*ht) * (psi[j+1][i] - psi[j-1][i])/(2.0*h2);
        i = (unsigned int)round(e[2]/h1);
        j = (unsigned int)round(e[3]/h2);
        g[2] = g[2] + f2((k)*ht) * (psi[j][i+1] - psi[j][i-1])/(2.0*h1);
        g[3] = g[3] + f2((k)*ht) * (psi[j+1][i] - psi[j-1][i])/(2.0*h2);
        i = (unsigned int)round(e[4]/h1);
        j = (unsigned int)round(e[5]/h2);
        g[4] = g[4] + f3((k)*ht) * (psi[j][i+1] - psi[j][i-1])/(2.0*h1);
        g[5] = g[5] + f3((k)*ht) * (psi[j+1][i] - psi[j-1][i])/(2.0*h2);
    }
    else
    {
        i = (unsigned int)round(e[0]/h1);
        j = (unsigned int)round(e[1]/h2);
        g[0] = g[0] + 2.0*f1((k)*ht) * (psi[j][i+1] - psi[j][i-1])/(2.0*h1);
        g[1] = g[1] + 2.0*f1((k)*ht) * (psi[j+1][i] - psi[j-1][i])/(2.0*h2);
        i = (unsigned int)round(e[2]/h1);
        j = (unsigned int)round(e[3]/h2);
        g[2] = g[2] + 2.0*f2((k)*ht) * (psi[j][i+1] - psi[j][i-1])/(2.0*h1);
        g[3] = g[3] + 2.0*f2((k)*ht) * (psi[j+1][i] - psi[j-1][i])/(2.0*h2);
        i = (unsigned int)round(e[4]/h1);
        j = (unsigned int)round(e[5]/h2);
        g[4] = g[4] + 2.0*f3((k)*ht) * (psi[j][i+1] - psi[j][i-1])/(2.0*h1);
        g[5] = g[5] + 2.0*f3((k)*ht) * (psi[j+1][i] - psi[j-1][i])/(2.0*h2);
    }

    if (k==0)
    {
        g[0] = -(ht/2.0)*g[0];
        g[1] = -(ht/2.0)*g[1];
        g[2] = -(ht/2.0)*g[2];
        g[3] = -(ht/2.0)*g[3];
        g[4] = -(ht/2.0)*g[4];
        g[5] = -(ht/2.0)*g[5];
    }
}

void HeatControl2DeltaX::calculateGF(const DoubleVector &e, DoubleVector& g)
{
}

void HeatControl2DeltaX::calculateG2(const DoubleVector &e, DoubleVector& g1)
{
    double h = 0.01;
    DoubleVector E(2*L);
    DoubleVector g(2*L);
    double f0 = fx(e);

    for (unsigned int i=0; i<e.size(); i++)
    {
        E = e;
        E[i] += h;
        double f1 = fx(E);
        g[i] = (f1-f0)/h;
    }

    printf("e2: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", e[0], e[1], e[2], e[3], e[4], e[5]);
    printf("g2: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", g[0], g[1], g[2], g[3], g[4], g[5]);
}

double HeatControl2DeltaX::fxt(unsigned int i, unsigned int j, unsigned k, const DoubleVector& e)
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t  = k*ht;
    double sum = 0.0;

    if (fabs(x1-e[0])<=h1 && fabs(x2-e[1])<=h2)
    {
        sum += f1(t) * ((h1-fabs(x1-e[0]))/(h1*h1))*((h2-fabs(x2-e[1]))/(h2*h2));
    }
    if (fabs(x1-e[2])<=h1 && fabs(x2-e[3])<=h2)
    {
        sum += f2(t) * ((h1-fabs(x1-e[2]))/(h1*h1))*((h2-fabs(x2-e[3]))/(h2*h2));
    }
    if (fabs(x1-e[4])<=h1 && fabs(x2-e[5])<=h2)
    {
        sum += f3(t) * ((h1-fabs(x1-e[4]))/(h1*h1))*((h2-fabs(x2-e[5]))/(h2*h2));
    }
    return sum;
}

void HeatControl2DeltaX::initialize()
{
    DoubleVector E;
    E.resize(2*L);

    E[0] = 0.7; E[1] = 0.2;
    E[2] = 0.5; E[3] = 0.8;
    E[4] = 0.2; E[5] = 0.3;

    calculateU(E, U);

    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
    Printer::printMatrix(U, N2/10, N1/10);
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");

    write("optimal.txt", U);
}

//Printer
void HeatControl2DeltaX::print(unsigned int i, const DoubleVector& e, const DoubleVector &gradient, double alpha, RnFunction* fn) const
{
    HeatControl2DeltaX *hc = dynamic_cast<HeatControl2DeltaX*>(fn);
    printf("J[%d]: %.16f\n", i, hc->fx(e));
    printf("e2: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", e[0], e[1], e[2], e[3], e[4], e[5]);

//    hc->calculateU(e, hc->uT);
//    char buffer [12];
//    int n=sprintf (buffer, "file%d.txt", i);
//    hc->write(buffer, hc->uT);
}

//Projection
void HeatControl2DeltaX::HeatControl2DeltaX::project(DoubleVector &e, int index)
{
    for (unsigned int i=0; i<e.size(); i++)
    {
        if (e[i]>1.0) e[i]=1.0;
        if (e[i]<0.0) e[i]=0.0;
    }
}

void HeatControl2DeltaX::write(const char *fileName, const DoubleMatrix& m)
{
    FILE* f = fopen(fileName, "w");
    for (unsigned int j=0; j<m.size(); j++)
    {
        for (unsigned int i=0; i<m[j].size(); i++)
        {
            if (i==0)
                fprintf(f, "%.10f", m[j][i]);
            else
                fprintf(f, " %.10f", m[j][i]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}
