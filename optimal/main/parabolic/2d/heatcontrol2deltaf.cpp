#include "heatcontrol2deltaf.h"
#include <tomasmethod.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <math.h>
#include <stdio.h>

void HeatControl2DeltaF::main()
{
    HeatControl2DeltaF hc(100, 10, 10);

    unsigned int count = (hc.M+1)*hc.L;
    DoubleVector f;
    f.resize(count);

    for (unsigned int i=0; i<f.size(); i++) f[i] = 0.5;

    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setEpsilon1(0.000000001);
    g2.setEpsilon2(0.000000001);
    g2.setR1MinimizeEpsilon(0.1, 0.000000001);
    g2.setPrinter(&hc);
    g2.setProjection(&hc);
    g2.setNormalize(false);
    g2.calculate(f);

    for (unsigned int l=0; l<hc.L; l++)
    {
        for (unsigned int k=0; k<=hc.M; k++)
        {
            if (k%(hc.M/10)==0) printf("%12.4f", f[l*(hc.M+1)+k]);
        }
        puts("");
    }
    puts("+++++++++++++++++++++++++++++++++++++++++");
}

HeatControl2DeltaF::HeatControl2DeltaF(unsigned int m, unsigned int n2, unsigned int n1) : RnFunction()
{
    alpha = 1.0;

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

HeatControl2DeltaF::~HeatControl2DeltaF() {}

double HeatControl2DeltaF::fx(const DoubleVector &f)
{
    calculateU(f, uT);

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

    double p;
    double norm = 0.0;
    for (unsigned int l=0; l<L; l++)
    {
        for (unsigned int k=0; k<=M; k++)
        {
            if (k==0 || k==M) p = 1.0; else p = 2.0;

            switch(l)
            {
            case 0: { norm += p*(ht/2.0)*(f[0*(M+1)+k] - f1(k*ht))*(f[0*(M+1)+k] - f1(k*ht)); } break;
            case 1: { norm += p*(ht/2.0)*(f[1*(M+1)+k] - f2(k*ht))*(f[1*(M+1)+k] - f2(k*ht)); } break;
            case 2: { norm += p*(ht/2.0)*(f[2*(M+1)+k] - f3(k*ht))*(f[2*(M+1)+k] - f3(k*ht)); } break;
            }
        }
    }

    return sum + alpha*norm;
}

void HeatControl2DeltaF::gradient(const DoubleVector &f, DoubleVector &g, double gradient_step)
{
    calculateU(f, uT);
    calculateP(f, g);
}

void HeatControl2DeltaF::calculateU(const DoubleVector &f, DoubleMatrix& u)
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
                    d[j-1] = alpha2*u0[j][i-1] + beta2*u0[j][i] + alpha2*u0[j][i+1] + (ht/2.0) * fxt(i, j, k, f);
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
                    d[i-1] = alpha2*u1[j-1][i] + beta2*u1[j][i] + alpha2*u1[j+1][i] + (ht/2.0) * fxt(i, j, k, f);
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

void HeatControl2DeltaF::calculateP(const DoubleVector &f, DoubleVector &g)
{
    DoubleMatrix psi0;
    DoubleMatrix psi1;

    std::vector<DoubleMatrix> psi;
    psi.resize(M+1);
    for (unsigned int k=0; k<=M; k++)
    {
        psi[k].resize(N2+1);
        for (unsigned int j=0; j<=N2; j++)
            psi[k][j].resize(N1+1);
    }

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

        psi[k] = psi0;
    }

    calculateG(f, g, psi);
}

void HeatControl2DeltaF::psiDerivative(double& psiX1, double& psiX2, double e1, double e2, const DoubleMatrix& psi)
{
    unsigned int i = (unsigned int)round(e1/h1);
    unsigned int j = (unsigned int)round(e2/h2);

    if (i==0)
        psiX1  = (psi[j][i+1] - psi[j][i])/h1;
    else if (i==N1)
        psiX1 = (psi[j][i] - psi[j][i-1])/h1;
    else
        psiX1 = (psi[j][i+1] - psi[j][i-1])/(2.0*h1);

    if (j==0)
        psiX2 = (psi[j+1][i] - psi[j][i])/h2;
    else if (j==N2)
        psiX2 = (psi[j][i] - psi[j-1][i])/h2;
    else
        psiX2 = (psi[j+1][i] - psi[j-1][i])/(2.0*h2);
}

void HeatControl2DeltaF::calculateG(const DoubleVector &f, DoubleVector &g, const std::vector<DoubleMatrix>& psi)
{
//    unsigned int i,j;
//    /* gradient of power */
//    for (unsigned k=0; k<=M; k++)
//    {
//        i = (unsigned int)round(E[0]/h1);
//        j = (unsigned int)round(E[1]/h2);
//        g[0*(M+1)+k] = -psi[k][j][i] + 2.0*alpha*(f[0*(M+1)+k] - f1(k*ht));
//        i = (unsigned int)round(E[2]/h1);
//        j = (unsigned int)round(E[3]/h2);
//        g[1*(M+1)+k] = -psi[k][j][i] + 2.0*alpha*(f[1*(M+1)+k] - f2(k*ht));
//        i = (unsigned int)round(E[4]/h1);
//        j = (unsigned int)round(E[5]/h2);
//        g[2*(M+1)+k] = -psi[k][j][i] + 2.0*alpha*(f[2*(M+1)+k] - f3(k*ht));
//    }

    /* gradient of power */
    for (unsigned k=0; k<=M; k++)
    {
        double p[L];
        p[0]=p[1]=p[2]=0.0;
        for (unsigned int j=0; j<=N2; j++)
        {
            for (unsigned int i=0; i<=N1; i++)
            {
                if (fabs(E[0]-i*h1)<h1 && fabs(E[1]-j*h2)<h2) p[0] += psi[k][j][i]*((h1-fabs(E[0]-i*h1))/h1)*((h2-fabs(E[1]-j*h2))/h2);
                if (fabs(E[2]-i*h1)<h1 && fabs(E[3]-j*h2)<h2) p[1] += psi[k][j][i]*((h1-fabs(E[2]-i*h1))/h1)*((h2-fabs(E[3]-j*h2))/h2);
                if (fabs(E[4]-i*h1)<h1 && fabs(E[5]-j*h2)<h2) p[2] += psi[k][j][i]*((h1-fabs(E[4]-i*h1))/h1)*((h2-fabs(E[5]-j*h2))/h2);
            }
        }
        g[0*(M+1)+k] = -p[0] + 2.0*alpha*(f[0*(M+1)+k] - f1(k*ht));
        g[1*(M+1)+k] = -p[1] + 2.0*alpha*(f[1*(M+1)+k] - f2(k*ht));
        g[2*(M+1)+k] = -p[2] + 2.0*alpha*(f[2*(M+1)+k] - f3(k*ht));
    }
}

double HeatControl2DeltaF::fxt(unsigned int i, unsigned int j, unsigned k, const DoubleVector &f)
{
    double x1 = i*h1;
    double x2 = j*h2;
    //double t  = k*ht;
    double sum = 0.0;

    if (fabs(x1-E[0])<=h1 && fabs(x2-E[1])<=h2)
    {
        //sum += f1(t) * ((h1-fabs(x1-E[0]))/(h1*h1))*((h2-fabs(x2-E[1]))/(h2*h2));
        sum += f[0*(M+1)+k] * ((h1-fabs(x1-E[0]))/(h1*h1))*((h2-fabs(x2-E[1]))/(h2*h2));
        //sum += (e[6+0*(M+1)+k] - f1(t)) * ((h1-fabs(x1-e[0]))/(h1*h1))*((h2-fabs(x2-e[1]))/(h2*h2));
    }
    if (fabs(x1-E[2])<=h1 && fabs(x2-E[3])<=h2)
    {
        //sum += f2(t) * ((h1-fabs(x1-E[2]))/(h1*h1))*((h2-fabs(x2-E[3]))/(h2*h2));
        sum += f[1*(M+1)+k] * ((h1-fabs(x1-E[2]))/(h1*h1))*((h2-fabs(x2-E[3]))/(h2*h2));
        //sum += (e[6+1*(M+1)+k] - f2(t)) * ((h1-fabs(x1-e[2]))/(h1*h1))*((h2-fabs(x2-e[3]))/(h2*h2));
    }
    if (fabs(x1-E[4])<=h1 && fabs(x2-E[5])<=h2)
    {
        //sum += f3(t) * ((h1-fabs(x1-E[4]))/(h1*h1))*((h2-fabs(x2-E[5]))/(h2*h2));
        sum += f[2*(M+1)+k] * ((h1-fabs(x1-E[4]))/(h1*h1))*((h2-fabs(x2-E[5]))/(h2*h2));
        //sum += (e[6+2*(M+1)+k]-f3(t)) * ((h1-fabs(x1-e[4]))/(h1*h1))*((h2-fabs(x2-e[5]))/(h2*h2));
    }

    return sum;
}

void HeatControl2DeltaF::initialize()
{
    DoubleVector f;
    f.resize((M+1)*L);

    E.resize(2*L);
    E[0] = 0.2; E[1] = 0.5;
    E[2] = 0.5; E[3] = 0.6;
    E[4] = 0.8; E[5] = 0.7;

    for (unsigned int k=0; k<=M; k++)
    {
        f[0*(M+1)+k] = f1(k*ht);
        f[1*(M+1)+k] = f2(k*ht);
        f[2*(M+1)+k] = f3(k*ht);
    }

    U.resize(N2+1);
    for (unsigned int j=0; j<=N2; j++)
    {
        U[j].resize(N1+1);
    }

    calculateU(f, U);

    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
    Printer::printMatrix(U, N2/10, N1/10);
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
    //exit(-1);
}

void HeatControl2DeltaF::print(unsigned int i, const DoubleVector &e, const DoubleVector &s, double a, RnFunction *f) const
{
    printf("J: %.16f\n", f->fx(e));
}

void HeatControl2DeltaF::project(DoubleVector &e, int index)
{
//    for (unsigned int i=0; i<6; i++)
//    {
//        if (e[i] < 0.0) e[i] = 0.0;
//        if (e[i] > 1.0) e[i] = 1.0;
//    }
}
