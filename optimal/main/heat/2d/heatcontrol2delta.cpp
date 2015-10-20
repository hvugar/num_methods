#include "heatcontrol2delta.h"
#include <tomasmethod.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <math.h>
#include <stdio.h>

void HeatControl2Delta::main()
{
    HeatControl2Delta hc(100, 10, 10);

    unsigned int count = 2*hc.L + (hc.M+1)*hc.L;
    DoubleVector e;
    e.resize(count);
    for (unsigned int i=2*hc.L; i<e.size(); i++) e[i] = 10.5;
    for (unsigned int i=0; i<2*hc.L; i++) e[i] = 0.5;

    //    e[0] = 0.2; e[1] = 0.5;
    //    e[2] = 0.5; e[3] = 0.6;
    //    e[4] = 0.8; e[5] = 0.7;

    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setEpsilon1(0.000000001);
    g2.setEpsilon2(0.000000001);
    g2.setGradientStep(0.000001);
    g2.setR1MinimizeEpsilon(0.01, 0.000001);
    g2.setPrinter(&hc);
    g2.setProjection(&hc);
    g2.setNormalize(false);
    g2.calculate(e);

    puts("+++++++++++++++++++++++++++++++++++++++++");
    printf("e[0]: %12.8f e[1]: %12.8f\n", e[0], e[1]);
    printf("e[2]: %12.8f e[3]: %12.8f\n", e[2], e[3]);
    printf("e[4]: %12.8f e[5]: %12.8f\n", e[4], e[5]);

    for (unsigned int l=0; l<hc.L; l++)
    {
        for (unsigned int k=0; k<=hc.M; k++)
        {
            if (k%(hc.M/10)==0) printf("%12.4f", e[6+l*(hc.M+1)+k]);
        }
        puts("");
    }
    puts("+++++++++++++++++++++++++++++++++++++++++");
}

HeatControl2Delta::HeatControl2Delta(unsigned int m, unsigned int n2, unsigned int n1) : RnFunction()
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

HeatControl2Delta::~HeatControl2Delta()
{}

double HeatControl2Delta::fx(const DoubleVector &e)
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
            case 0: { norm += p*(ht/2.0)*(e[6+0*(M+1)+k] - f1(k*ht))*(e[6+0*(M+1)+k] - f1(k*ht)); } break;
            case 1: { norm += p*(ht/2.0)*(e[6+1*(M+1)+k] - f2(k*ht))*(e[6+1*(M+1)+k] - f2(k*ht)); } break;
            case 2: { norm += p*(ht/2.0)*(e[6+2*(M+1)+k] - f3(k*ht))*(e[6+2*(M+1)+k] - f3(k*ht)); } break;
            }
        }
    }

    return sum + alpha*norm;
}

void HeatControl2Delta::gradient(const DoubleVector &e, DoubleVector &g, double gradient_step)
{
    calculateU(e, uT);
    calculateP(e, g);
}

void HeatControl2Delta::calculateU(const DoubleVector &E1, DoubleMatrix& u)
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
                    d[j-1] = alpha2*u0[j][i-1] + beta2*u0[j][i] + alpha2*u0[j][i+1] + (ht/2.0) * fxt(i, j, k, E1);
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
                    d[i-1] = alpha2*u1[j-1][i] + beta2*u1[j][i] + alpha2*u1[j+1][i] + (ht/2.0) * fxt(i, j, k, E1);
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

void HeatControl2Delta::calculateP(const DoubleVector &e, DoubleVector &g)
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

    calculateG(e, g, psi);
}

void HeatControl2Delta::calculateG(const DoubleVector &e, DoubleVector &g, const std::vector<DoubleMatrix>& psi)
{
    //unsigned int i1,j1;
    double psiX1 = 0.0;
    double psiX2 = 0.0;

    /* gradient of power */
    for (unsigned k=0; k<=M; k++)
    {
        double p[L];
        p[0]=p[1]=p[2]=0.0;
        for (unsigned int j=0; j<=N2; j++)
        {
            for (unsigned int i=0; i<=N1; i++)
            {
                if (fabs(e[0]-i*h1)<h1 && fabs(e[1]-j*h2)<h2) p[0] += psi[k][j][i]*((h1-fabs(e[0]-i*h1))/h1)*((h2-fabs(e[1]-j*h2))/h2);
                if (fabs(e[2]-i*h1)<h1 && fabs(e[3]-j*h2)<h2) p[1] += psi[k][j][i]*((h1-fabs(e[2]-i*h1))/h1)*((h2-fabs(e[3]-j*h2))/h2);
                if (fabs(e[4]-i*h1)<h1 && fabs(e[5]-j*h2)<h2) p[2] += psi[k][j][i]*((h1-fabs(e[4]-i*h1))/h1)*((h2-fabs(e[5]-j*h2))/h2);
            }
        }
        g[6+0*(M+1)+k] = -p[0] + 2.0*alpha*(e[6+0*(M+1)+k] - f1(k*ht));
        g[6+1*(M+1)+k] = -p[1] + 2.0*alpha*(e[6+1*(M+1)+k] - f2(k*ht));
        g[6+2*(M+1)+k] = -p[2] + 2.0*alpha*(e[6+2*(M+1)+k] - f3(k*ht));
    }

    /* gradient of replacement */
    for (unsigned k=0; k<=M; k++)
    {
        double p = 1.0;
        if (k==M || k==0)
            p = 1.0 * -(ht/2.0);
        else
            p = 2.0 * -(ht/2.0);
        //double t = k*ht;

        psiDerivative(psiX1, psiX2, e[0], e[1], psi[k]);
        g[0] = g[0] + p * e[6+0*(M+1)+k] * psiX1;
        g[1] = g[1] + p * e[6+0*(M+1)+k] * psiX2;
        //g[0] = g[0] + p * f1(t) * psiX1;
        //g[1] = g[1] + p * f1(t) * psiX2;

        psiDerivative(psiX1, psiX2, e[2], e[3], psi[k]);
        g[2] = g[2] + p * e[6+1*(M+1)+k] * psiX1;
        g[3] = g[3] + p * e[6+1*(M+1)+k] * psiX2;
        //g[2] = g[2] + p * f2(t) * psiX1;
        //g[3] = g[3] + p * f2(t) * psiX2;

        psiDerivative(psiX1, psiX2, e[4], e[5], psi[k]);
        g[4] = g[4] + p * e[6+2*(M+1)+k] * psiX1;
        g[5] = g[5] + p * e[6+2*(M+1)+k] * psiX2;
        //g[4] = g[4] + p * f3(t) * psiX1;
        //g[5] = g[5] + p * f3(t) * psiX2;
    }
}

double HeatControl2Delta::fxt(unsigned int i, unsigned int j, unsigned k, const DoubleVector& e)
{
    double x1 = i*h1;
    double x2 = j*h2;
    //double t  = k*ht;
    double sum = 0.0;

    if (fabs(x1-e[0])<h1 && fabs(x2-e[1])<h2)
    {
        //sum += f1(t) * ((h1-fabs(x1-e[0]))/(h1*h1))*((h2-fabs(x2-e[1]))/(h2*h2));
        sum += e[6+0*(M+1)+k] * ((h1-fabs(x1-e[0]))/(h1*h1))*((h2-fabs(x2-e[1]))/(h2*h2));
        //sum += (e[6+0*(M+1)+k] - f1(t)) * ((h1-fabs(x1-e[0]))/(h1*h1))*((h2-fabs(x2-e[1]))/(h2*h2));
    }
    if (fabs(x1-e[2])<h1 && fabs(x2-e[3])<h2)
    {
        //sum += f2(t) * ((h1-fabs(x1-e[2]))/(h1*h1))*((h2-fabs(x2-e[3]))/(h2*h2));
        sum += e[6+1*(M+1)+k] * ((h1-fabs(x1-e[2]))/(h1*h1))*((h2-fabs(x2-e[3]))/(h2*h2));
        //sum += (e[6+1*(M+1)+k] - f2(t)) * ((h1-fabs(x1-e[2]))/(h1*h1))*((h2-fabs(x2-e[3]))/(h2*h2));
    }
    if (fabs(x1-e[4])<h1 && fabs(x2-e[5])<h2)
    {
        //sum += f3(t) * ((h1-fabs(x1-e[4]))/(h1*h1))*((h2-fabs(x2-e[5]))/(h2*h2));
        sum += e[6+2*(M+1)+k] * ((h1-fabs(x1-e[4]))/(h1*h1))*((h2-fabs(x2-e[5]))/(h2*h2));
        //sum += (e[6+2*(M+1)+k]-f3(t)) * ((h1-fabs(x1-e[4]))/(h1*h1))*((h2-fabs(x2-e[5]))/(h2*h2));
    }

    return sum;
}

void HeatControl2Delta::psiDerivative(double& psiX1, double& psiX2, double e1, double e2, const DoubleMatrix& psi)
{
    for (unsigned int j=0; j<=N2-1; j++)
    {
        for (unsigned int i=0; i<=N1-1; i++)
        {
            if (i*h1<=e1 && e1<(i+1)*h1 && j*h2<=e2 && e2<(j+1)*h2)
            {
                psiX1 = (psi[j][i+1]*((fabs(e1-i*h1))/h1) - psi[j][i]*((h1-fabs(e1-i*h1)/h1)))/h1;
                psiX2 = (psi[j+1][i]*((fabs(e2-j*h2))/h2) - psi[j][i]*((h2-fabs(e2-i*h2)/h2)))/h2;
                //psiX1  = (psi[j][i+1] - psi[j][i])/h1;
                //psiX2  = (psi[j+1][i] - psi[j][i])/h2;
            }
        }
    }

//    unsigned int i = (unsigned int)round(e1/h1);
//    unsigned int j = (unsigned int)round(e2/h2);

//    if (i==0)
//        psiX1  = (psi[j][i+1] - psi[j][i])/h1;
//    else if (i==N1)
//        psiX1 = (psi[j][i] - psi[j][i-1])/h1;
//    else
//        psiX1 = (psi[j][i+1] - psi[j][i-1])/(2.0*h1);

//    if (j==0)
//        psiX2 = (psi[j+1][i] - psi[j][i])/h2;
//    else if (j==N2)
//        psiX2 = (psi[j][i] - psi[j-1][i])/h2;
//    else
//        psiX2 = (psi[j+1][i] - psi[j-1][i])/(2.0*h2);
}

void HeatControl2Delta::initialize()
{
    DoubleVector E;
    E.resize(2*L + (M+1)*L);
    E[0] = 0.2; E[1] = 0.5;
    E[2] = 0.5; E[3] = 0.6;
    E[4] = 0.8; E[5] = 0.7;
    for (unsigned int k=0; k<=M; k++)
    {
        E[6+0*(M+1)+k] = f1(k*ht);
        E[6+1*(M+1)+k] = f2(k*ht);
        E[6+2*(M+1)+k] = f3(k*ht);
    }

    U.resize(N2+1);
    for (unsigned int j=0; j<=N2; j++)
    {
        U[j].resize(N1+1);
    }

    calculateU(E, U);

    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
    Printer::printMatrix(U, N2/10, N1/10);
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
    //exit(-1);
}

void HeatControl2Delta::print(unsigned int i, const DoubleVector &e, const DoubleVector &g, double a, RnFunction *f) const
{
    //HeatControl2Delta *hc = dynamic_cast<HeatControl2Delta*>(f);
    printf("J[%d]: %.16f %.16f\n", i, f->fx(e), a);
    printf("e[0]:%10.6f e[1]:%10.6f e[0]:%10.6f e[1]:%10.6f e[0]:%10.6f e[1]:%10.6f\n", e[0], e[1], e[2], e[3], e[4], e[5]);
    printf("g[0]:%10.6f g[1]:%10.6f g[0]:%10.6f g[1]:%10.6f g[0]:%10.6f g[1]:%10.6f\n", g[0], g[1], g[2], g[3], g[4], g[5]);
//    puts("-------------------------------------------------------------------------");
//    for (unsigned int l=0; l<hc->L; l++)
//    {
//        for (unsigned int k=0; k<=hc->M; k++)
//        {
//            if (k%(hc->M/10)==0) printf("%12.4f", e[6+l*(hc->M+1)+k]);
//        }
//        puts("");
//    }
    puts("+++++++++++++++++++++++++++++++++++++++++");
}

void HeatControl2Delta::project(DoubleVector &e, int index)
{
    for (unsigned int i=0; i<6; i++)
    {
        if (e[i] < 0.0) e[i] = 0.0;
        if (e[i] > 1.0) e[i] = 1.0;
    }
}
