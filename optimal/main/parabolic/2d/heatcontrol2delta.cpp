#include "heatcontrol2delta.h"
#include <tomasmethod.h>
#include <gradient_cjt.h>
#include <math.h>
#include <stdlib.h>
#include <windows.h>

//#define PLACE_OPTIMIZE
//#define POWER_OPTIMIZE

void HeatControl2Delta::main()
{
    HeatControl2Delta hc(100, 100, 100);
    hc.optimal.resize(2*hc.L);
    hc.optimal[0] = 0.50;
    hc.optimal[1] = 0.80;
    hc.optimal[2] = 0.70;
    hc.optimal[3] = 0.20;
    hc.optimal[4] = 0.20;
    hc.optimal[5] = 0.30;

    hc.initialize();

    DoubleVector x;
#ifdef POWER_OPTIMIZE
    x.resize( 2*hc.L + (hc.M+1)*hc.L );
#else
    x.resize( 2*hc.L );
#endif

    //x[0] = 0.50; x[1] = 0.80; x[2] = 0.70; x[3] = 0.20; x[4] = 0.20; x[5] = 0.30; optimal
    x[0] = 0.60; x[1] = 0.70; x[2] = 0.65; x[3] = 0.25; x[4] = 0.25; x[5] = 0.35; //ishleyir eps1:0.0001 eps2:0.0001 min:1.0, 0.0001
    //x[0] = 0.60; x[1] = 0.70; x[2] = 0.60; x[3] = 0.30; x[4] = 0.30; x[5] = 0.40; //ishleyir eps1:0.0001 eps2:0.0001 min:1.0, 0.0001

//    x[0] = 0.5; x[1] = 0.8;
//    x[2] = 0.7; x[3] = 0.2;
//    x[4] = 0.2; x[5] = 0.3;
//    FILE* file = fopen("e.txt", "w");
//    for (unsigned int i=0; i<=hc.N1; i++)
//    {
//        x[5] = i*hc.h1;
//        double result = hc.fx(x);
//        fprintf(file, "%.16f\n", result);
//        fflush(file);
//    }
//    fclose(file);
//    return;



#ifdef POWER_OPTIMIZE
    for (unsigned int k=0; k<=hc.M; k++)
    {
        x[2*hc.L + 0*(hc.M+1) + k] = 0.0;//hc.g1(k*hc.ht);
        x[2*hc.L + 1*(hc.M+1) + k] = 0.0;//hc.g2(k*hc.ht);
        x[2*hc.L + 2*(hc.M+1) + k] = 0.0;//hc.g3(k*hc.ht);
    }
#endif

    hc.initial.resize(2*hc.L);
    hc.initial[0] = x[0];
    hc.initial[1] = x[1];
    hc.initial[2] = x[2];
    hc.initial[3] = x[3];
    hc.initial[4] = x[4];
    hc.initial[5] = x[5];

    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setGradient(&hc);
    g2.setEpsilon1(0.00001);
    g2.setEpsilon2(0.00001);
    g2.setEpsilon3(0.00001);
    g2.setR1MinimizeEpsilon(1.0, 0.0001);
    g2.setPrinter(&hc);
    g2.setProjection(&hc);
    g2.setNormalize(true);
    g2.calculate(x);

    printf("J[%d]: %.16f\n", 0, hc.fx(x));
    printf("eo: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", 0.50, 0.80, 0.70, 0.20, 0.20, 0.30);
    printf("e1: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", x[0], x[1], x[2], x[3], x[4], x[5]);
}

HeatControl2Delta::HeatControl2Delta(unsigned int M, unsigned int N2, unsigned int N1)
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
}

double HeatControl2Delta::fx(const DoubleVector& x)
{
    DoubleMatrix u;
    calculateU(x, u);

    double sum = 0.0;
    for (unsigned int j=0; j<=N2-1; j++)
    {
        for (unsigned int i=0; i<=N1-1; i++)
        {
            int j1 = j;
            int j2 = j+1;
            int i1 = i;
            int i2 = i+1;

            double f11 = u[j1][i1] - U[j1][i1];
            double f12 = u[j1][i2] - U[j1][i2];
            double f21 = u[j2][i1] - U[j2][i1];
            double f22 = u[j2][i2] - U[j2][i2];

            sum += (f11*f11 + f12*f12 + f21*f21 + f22*f22);
        }
    }
    sum = (0.25*(h1*h2))*sum;

#ifdef POWER_OPTIMIZE
    sum = sum + alpha*norm(x);
#endif
    return sum;
}

double HeatControl2Delta::norm(const DoubleVector& x) const
{
    double p;
    double nrm = 0.0;
    for (unsigned int l=0; l<L; l++)
    {
        for (unsigned int k=0; k<=M; k++)
        {
            if (k==0 || k==M) p = 1.0; else p = 2.0;

            switch(l)
            {
            case 0: { nrm += p*(ht/2.0)*(x[2*L+0*(M+1)+k] - g1(k*ht))*(x[2*L+0*(M+1)+k] - g1(k*ht)); } break;
            case 1: { nrm += p*(ht/2.0)*(x[2*L+1*(M+1)+k] - g2(k*ht))*(x[2*L+1*(M+1)+k] - g2(k*ht)); } break;
            case 2: { nrm += p*(ht/2.0)*(x[2*L+2*(M+1)+k] - g3(k*ht))*(x[2*L+2*(M+1)+k] - g3(k*ht)); } break;
            }
        }
    }

    return nrm;
}

void HeatControl2Delta::gradient(const DoubleVector& x, DoubleVector& g)
{
    //    double nrm = norm(x);
    //    if (nrm < 0.0000201) alpha = 0.0;

    //static int i=1;
    DoubleMatrix u;
    calculateU(x, u);
    calculateP(x, g, u);

    //    char filename[20];
    //    int n = sprintf(filename, "heat%d.txt", i++);
    //    filename[n] = '\0';
    //    write(filename, u);
    //puts("-----------------------------------------------------------");
    //printf("e1: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", e[0], e[1], e[2], e[3], e[4], e[5]);
    //printf("g1: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", g[0], g[1], g[2], g[3], g[4], g[5]);
    //calculateG2(e, g);
}

void HeatControl2Delta::calculateU(const DoubleVector &x, DoubleMatrix& u)
{
    px = &x;

    DoubleMatrix u0(N2+1, N1+1);
    DoubleMatrix u1(N2+1, N1+1);

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
                    u0[j][i] = fi(i*h1, j*h2);
                }
            }
        }
        else
        {
            // Approximation to x1 direction
            //            da1.resize(N1-1);
            //            db1.resize(N1-1);
            //            dc1.resize(N1-1);
            //            dd1.resize(N1-1);
            //            rx1.resize(N1-1);


            for (unsigned int i=1; i<N1; i++)
            {
                for (unsigned int j=1; j<N2; j++)
                {
                    da1[j-1] = x1_alpha1;
                    db1[j-1] = x1_beta1;
                    dc1[j-1] = x1_alpha1;
                    dd1[j-1] = x1_alpha2*u0[j][i-1] + x1_beta2*u0[j][i] + x1_alpha2*u0[j][i+1] + (ht/2.0) * f(i, j, k);
                }

                da1[0]     = 0.0;
                dc1[N2-2]  = 0.0;
                dd1[0]    -= x1_alpha1 * m3(h1*i, ht*(k));
                dd1[N2-2] -= x1_alpha1 * m4(h1*i, ht*(k));

                TomasAlgorithm(da1, db1, dc1, dd1, rx1);

                u1[0][i]  = m3(h1*i, ht*k - ht*0.5);
                for (unsigned int j=1; j<N2; j++)
                {
                    u1[j][i] = rx1[j-1];
                }
                u1[N2][i] = m4(h1*i, ht*k - ht*0.5);
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                u1[j][0]  = m1(h2*j, ht*k - ht*0.5);
                u1[j][N1] = m2(h2*j, ht*k - ht*0.5);
            }

            //            da1.clear();
            //            db1.clear();
            //            dc1.clear();
            //            dd1.clear();
            //            rx1.clear();

            // Approximation to x2 direction
            //            da2.resize(N2-1);
            //            db2.resize(N2-1);
            //            dc2.resize(N2-1);
            //            dd2.resize(N2-1);
            //            rx2.resize(N2-1);

            gause1 = gause2 = gause3 = 0.0;
            //file = fopen("gause.txt", "w");
            for (unsigned int j=1; j<N2; j++)
            {
                for (unsigned int i=1; i<N1; i++)
                {
                    da2[i-1] = x2_alpha1;
                    db2[i-1] = x2_beta1;
                    dc2[i-1] = x2_alpha1;
                    dd2[i-1] = x2_alpha2*u1[j-1][i] + x2_beta2*u1[j][i] + x2_alpha2*u1[j+1][i] + (ht/2.0) * f(i, j, k);
                }
                da2[0]     = 0.0;
                dc2[N1-2]  = 0.0;
                dd2[0]    -= x2_alpha1 * m1(h2*j, ht*(k));
                dd2[N1-2] -= x2_alpha1 * m2(h2*j, ht*(k));

                TomasAlgorithm(da2, db2, dc2, dd2, rx2);

                u0[j][0]  = m1(h2*j, ht*(k));
                for (unsigned int i=1; i<N1; i++)
                {
                    u0[j][i] = rx2[i-1];
                }
                u0[j][N1] = m2(h2*j, ht*(k));
            }

            for (unsigned int i=0; i<=N1; i++)
            {
                u0[0][i]  = m3(h1*i, ht*(k));
                u0[N2][i] = m4(h1*i, ht*(k));
            }
            //fclose(file);
//            puts("------------------------------------------------------------------------");
//            printf("gause1: %.8f %.8f [%.8f %.8f]\n", gause1, gause1 * h1 * h2, x[0], x[1]);
//            printf("gause2: %.8f %.8f [%.8f %.8f]\n", gause2, gause2 * h1 * h2, x[2], x[3]);
//            printf("gause3: %.8f %.8f [%.8f %.8f]\n", gause3, gause3 * h1 * h2, x[4], x[5]);

            //            da2.clear();
            //            db2.clear();
            //            dc2.clear();
            //            dd2.clear();
            //            rx2.clear();
        }
    }

    u = u0;
}

void HeatControl2Delta::calculateP(const DoubleVector &x, DoubleVector &g, const DoubleMatrix &u)
{
    DoubleMatrix psi0(N2+1, N1+1);
    DoubleMatrix psi1(N2+1, N1+1);

    //    psi0.resize(N2+1); for (unsigned int j=0; j<=N2; j++) psi0[j].resize(N1+1);
    //    psi1.resize(N2+1); for (unsigned int j=0; j<=N2; j++) psi1[j].resize(N1+1);

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

    //    DoubleVector da;
    //    DoubleVector db;
    //    DoubleVector dc;
    //    DoubleVector dd;
    //    DoubleVector rx;

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
                    psi0[j][i] = -2.0*(u[j][i] - U[j][i]);
                }
            }
        }
        else
        {
            // Approximation to x1 direction
            //            da1.resize(N1-1);
            //            db1.resize(N1-1);
            //            dc1.resize(N1-1);
            //            dd1.resize(N1-1);
            //            rx1.resize(N1-1);

            for (unsigned int j=1; j<N2; j++)
            {
                for (unsigned int i=1; i<N1; i++)
                {
                    da1[i-1] = x1_alpha1;
                    db1[i-1] = x1_beta1;
                    dc1[i-1] = x1_alpha1;
                    dd1[i-1] = x1_alpha2*psi0[j-1][i] + x1_beta2*psi0[j][i] + x1_alpha2*psi0[j+1][i];
                }

                da1[0]     = 0.0;
                dc1[N1-2]  = 0.0;
                dd1[0]    -= x1_alpha1 * pm1(h2*j, ht*(k+0.5));
                dd1[N1-2] -= x1_alpha1 * pm2(h2*j, ht*(k+0.5));

                TomasAlgorithm(da1, db1, dc1, dd1, rx1);

                psi1[j][0]  = pm1(h2*j, ht*(k+0.5));
                for (unsigned int i=1; i<N1; i++)
                {
                    psi1[j][i] = rx1[i-1];
                }
                psi1[j][N1] = pm2(h2*j, ht*(k+0.5));
            }

            for (unsigned int i=0; i<=N1; i++)
            {
                psi1[0][i]  = pm3(h1*i, ht*(k+0.5));
                psi1[N2][i] = pm4(h1*i, ht*(k+0.5));
            }

            //            da1.clear();
            //            db1.clear();
            //            dc1.clear();
            //            dd1.clear();
            //            rx1.clear();

            // Approximation to x2 direction
            //            da2.resize(N2-1);
            //            db2.resize(N2-1);
            //            dc2.resize(N2-1);
            //            dd2.resize(N2-1);
            //            rx2.resize(N2-1);

            for (unsigned int i=1; i<N1; i++)
            {
                for (unsigned int j=1; j<N2; j++)
                {
                    da2[j-1] = x2_alpha1;
                    db2[j-1] = x2_beta1;
                    dc2[j-1] = x2_alpha1;
                    dd2[j-1] = x2_alpha2*psi1[j][i-1] + x2_beta2*psi1[j][i] + x2_alpha2*psi1[j][i+1];
                }
                da2[0]     = 0.0;
                dc2[N2-2]  = 0.0;
                dd2[0]    -= x2_alpha1 * pm3(h1*i, ht*k);
                dd2[N2-2] -= x2_alpha1 * pm4(h1*i, ht*k);
                TomasAlgorithm(da2, db2, dc2, dd2, rx2);

                psi0[0][i]  = pm3(h1*i, ht*k);
                for (unsigned int j=1; j<N2; j++)
                {
                    psi0[j][i] = rx2[j-1];
                }
                psi0[N2][i] = pm4(h1*i, ht*k);
            }

            for (unsigned int j=0; j<=N2; j++)
            {
                psi0[j][0]  = pm1(h2*j, ht*k);
                psi0[j][N1] = pm2(h2*j, ht*k);
            }

            //            da2.clear();
            //            db2.clear();
            //            dc2.clear();
            //            dd2.clear();
            //            rx2.clear();
        }

        calculateGX(x, psi0, g, k);
#ifdef POWER_OPTIMIZE
        calculateGF(x, psi0, g, k);
#endif
    }
}

void HeatControl2Delta::calculateGX(const DoubleVector& x, const DoubleMatrix& psi, DoubleVector& g, unsigned int k)
{
    double psiX1;
    double psiX2;
    if (k==0 || k==M)
    {
        psiDerivative(psiX1, psiX2, x[0], x[1], psi);
        g[0] = g[0] + g1((k)*ht) * psiX1;
        g[1] = g[1] + g1((k)*ht) * psiX2;
        psiDerivative(psiX1, psiX2, x[2], x[3], psi);
        g[2] = g[2] + g2((k)*ht) * psiX1;
        g[3] = g[3] + g2((k)*ht) * psiX2;
        psiDerivative(psiX1, psiX2, x[4], x[5], psi);
        g[4] = g[4] + g3((k)*ht) * psiX1;
        g[5] = g[5] + g3((k)*ht) * psiX2;
    }
    else
    {
        psiDerivative(psiX1, psiX2, x[0], x[1], psi);
        g[0] = g[0] + 2.0*g1((k)*ht) * psiX1;
        g[1] = g[1] + 2.0*g1((k)*ht) * psiX2;
        psiDerivative(psiX1, psiX2, x[2], x[3], psi);
        g[2] = g[2] + 2.0*g2((k)*ht) * psiX1;
        g[3] = g[3] + 2.0*g2((k)*ht) * psiX2;
        psiDerivative(psiX1, psiX2, x[4], x[5], psi);
        g[4] = g[4] + 2.0*g3((k)*ht) * psiX1;
        g[5] = g[5] + 2.0*g3((k)*ht) * psiX2;
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

void HeatControl2Delta::psiDerivative(double &psiX1, double &psiX2, double x1, double x2, const DoubleMatrix &psi)
{
    unsigned int i = (unsigned int)round(x1/h1);
    unsigned int j = (unsigned int)round(x2/h2);

    if (i==0) psiX1  = (psi[j][i+1] - psi[j][i])/h1;
    else if (i==N1) psiX1 = (psi[j][i] - psi[j][i-1])/h1;
    else psiX1 = (psi[j][i+1] - psi[j][i-1])/(2.0*h1);

    if (j==0) psiX2 = (psi[j+1][i] - psi[j][i])/h2;
    else if (j==N2) psiX2 = (psi[j][i] - psi[j-1][i])/h2;
    else psiX2 = (psi[j+1][i] - psi[j-1][i])/(2.0*h2);
}

void HeatControl2Delta::calculateGF(const DoubleVector &x, const DoubleMatrix& psi, DoubleVector& g, unsigned int k)
{
//    for (unsigned k=0; k<=M; k++)
//    {
//        if (alpha < 1.0)
//        {
//            g[2*L+0*(M+1)+k] = 0.0;
//            g[2*L+1*(M+1)+k] = 0.0;
//            g[2*L+2*(M+1)+k] = 0.0;
//        }
//        else
//        {
//            double p[L];
//            p[0]=p[1]=p[2]=0.0;
//            for (unsigned int j=0; j<=N2; j++)
//            {
//                for (unsigned int i=0; i<=N1; i++)
//                {
//                    if (fabs(x[0]-i*h1)<h1 && fabs(x[1]-j*h2)<h2) p[0] += psi[j][i]*((h1-fabs(x[0]-i*h1))/h1)*((h2-fabs(x[1]-j*h2))/h2);
//                    if (fabs(x[2]-i*h1)<h1 && fabs(x[3]-j*h2)<h2) p[1] += psi[j][i]*((h1-fabs(x[2]-i*h1))/h1)*((h2-fabs(x[3]-j*h2))/h2);
//                    if (fabs(x[4]-i*h1)<h1 && fabs(x[5]-j*h2)<h2) p[2] += psi[j][i]*((h1-fabs(x[4]-i*h1))/h1)*((h2-fabs(x[5]-j*h2))/h2);
//                }
//            }
//            g[2*L+0*(M+1)+k] = -p[0] + 2.0*alpha*(x[2*L+0*(M+1)+k] - g1(k*ht));
//            g[2*L+1*(M+1)+k] = -p[1] + 2.0*alpha*(x[2*L+1*(M+1)+k] - g2(k*ht));
//            g[2*L+2*(M+1)+k] = -p[2] + 2.0*alpha*(x[2*L+2*(M+1)+k] - g3(k*ht));
//        }
//    }
}

void HeatControl2Delta::calculateG2(const DoubleVector &x, DoubleVector& g1)
{
//    double h = 0.01;
//    DoubleVector E(2*L);
//    DoubleVector g(2*L);
//    double f0 = fx(x);

//    for (unsigned int i=0; i<x.size(); i++)
//    {
//        E = x;
//        E[i] += h;
//        double f1 = fx(E);
//        g[i] = (f1-f0)/h;
//    }

//    printf("e2: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", x[0], x[1], x[2], x[3], x[4], x[5]);
//    printf("g2: %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", g[0], g[1], g[2], g[3], g[4], g[5]);
}

double HeatControl2Delta::f(unsigned int i, unsigned int j, unsigned int k)
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t  = k*ht;
    double sum = 0.0;

    DoubleVector e(6);
    e[0] = (*px)[0];
    e[1] = (*px)[1];
    e[2] = (*px)[2];
    e[3] = (*px)[3];
    e[4] = (*px)[4];
    e[5] = (*px)[5];


//    DoubleVector e(2*L);
//    for (unsigned int l=0; l<L; l++)
//    {
//        e[2*l+0] = (*px)[2*l+0];
//        e[2*l+1] = (*px)[2*l+1];
//    }

    //    if (fabs(x1-e[0])<=h1 && fabs(x2-e[1])<=h2)
    //    {
    //        sum += g1(t) * ((h1-fabs(x1-e[0]))/(h1*h1))*((h2-fabs(x2-e[1]))/(h2*h2));
    //    }
    //    if (fabs(x1-e[2])<=h1 && fabs(x2-e[3])<=h2)
    //    {
    //        sum += g2(t) * ((h1-fabs(x1-e[2]))/(h1*h1))*((h2-fabs(x2-e[3]))/(h2*h2));
    //    }
    //    if (fabs(x1-e[4])<=h1 && fabs(x2-e[5])<=h2)
    //    {
    //        sum += g3(t) * ((h1-fabs(x1-e[4]))/(h1*h1))*((h2-fabs  v c  v (x2-e[5]))/(h2*h2));
    //    }

    double sgm1 = 3.0*h1;
    double sgm2 = 3.0*h2;
    double a = 1.0/(2.0*M_PI*sgm1*sgm2);\
    double b = 2.0*sgm1*sgm2;

#ifdef POWER_OPTIMIZE
    double _g1 = (*px)[2*L + 0*(M+1) + k];
    double _g2 = (*px)[2*L + 1*(M+1) + k];
    double _g3 = (*px)[2*L + 2*(M+1) + k];

    sum += (_g1) * a * exp(-((x1-e[0])*(x1-e[0]) + (x2-e[1])*(x2-e[1]))/b);
    sum += (_g2) * a * exp(-((x1-e[2])*(x1-e[2]) + (x2-e[3])*(x2-e[3]))/b);
    sum += (_g3) * a * exp(-((x1-e[4])*(x1-e[4]) + (x2-e[5])*(x2-e[5]))/b);
#else
    sum += g1(t) * a * exp(-((x1-e[0])*(x1-e[0]) + (x2-e[1])*(x2-e[1]))/b);// * h1*h2;
    sum += g2(t) * a * exp(-((x1-e[2])*(x1-e[2]) + (x2-e[3])*(x2-e[3]))/b);// * h1*h2;
    sum += g3(t) * a * exp(-((x1-e[4])*(x1-e[4]) + (x2-e[5])*(x2-e[5]))/b);// * h1*h2;

    //double pp = a * exp(-((x1-e[0])*(x1-e[0]) + (x2-e[1])*(x2-e[1]))/b);
    //if (i==50 && j==80)
    //printf("%3d %3d %3d %.16f %.16f %.16f %.16f %.16f %.16f %.16f\n", k, j, i, pp, e[0], e[1], x1, x2, exp(-((x1-e[0])*(x1-e[0]) + (x2-e[1])*(x2-e[1]))/b), a * exp(-((x1-e[0])*(x1-e[0]) + (x2-e[1])*(x2-e[1]))/b));
    gause1 += a * exp(-((x1-e[0])*(x1-e[0]) + (x2-e[1])*(x2-e[1]))/b);
    gause2 += a * exp(-((x1-e[2])*(x1-e[2]) + (x2-e[3])*(x2-e[3]))/b);
    gause3 += a * exp(-((x1-e[4])*(x1-e[4]) + (x2-e[5])*(x2-e[5]))/b);

    //    if (k==M)//
    //    //if (pp > 0.0000001)
    //    {
    //        fprintf(file, "%.16f ", pp);
    //        if (i == 99) fputs("\n", file);
    //        //printf("%d %d %d\n", i, j, k);
    //    }
#endif

    e.clear();
    return sum;
}

void HeatControl2Delta::initialize()
{
    // cleaning U
    for (unsigned int j=0; j<U.size(); j++) U[j].clear();
    U.clear();

    DoubleVector x = optimal;
    // initializing U
    U.resize(N2+1);
    for (unsigned int j=0; j<=N2; j++) U[j].resize(N1+1);

#ifdef POWER_OPTIMIZE
    x.resize( 2*L + (M+1)*L );
#else
    x.resize( 2*L );
#endif

#ifdef POWER_OPTIMIZE
    for (unsigned int k=0; k<=M; k++)
    {
        x[2*L + 0*(M+1) + k] = g1(k*ht);
        x[2*L + 1*(M+1) + k] = g2(k*ht);
        x[2*L + 2*(M+1) + k] = g3(k*ht);
    }
#endif

    calculateU(optimal, U);

    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
    IPrinter::printMatrix(U, 10, 10);
    printf("eo: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", x[0], x[1], x[2], x[3], x[4], x[5]);
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");

    write("optimal.txt", U);
}

void HeatControl2Delta::print(unsigned int i, const DoubleVector& x, const DoubleVector &g, double alpha, RnFunction* fn) const
{
    HeatControl2Delta *hc = dynamic_cast<HeatControl2Delta*>(fn);
    printf("J[%d]: %.16f\n", i, hc->fx(x));
    printf("Norm: %.16f Alpha: %.16f %.16f\n", hc->norm(x), hc->alpha, alpha);
    printf("eo: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", 0.50, 0.80, 0.70, 0.20, 0.20, 0.30);
    printf("e1: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", x[0], x[1], x[2], x[3], x[4], x[5]);
    printf("g1: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", g[0], g[1], g[2], g[3], g[4], g[5]);

    //    DoubleMatrix u;
    //    hc->calculateU(x, u);
    //    char filename1[100];
    //    char filename2[100];
    //    int count1 = sprintf(filename1, "optimal%d.txt", i);
    //    filename1[count1] = '\0';
    //    int count2 = sprintf(filename2, "optimal%d.png", i);
    //    filename2[count2] = '\0';
    //    hc->write(filename2, u);
    //    system("imager.exe -w 101 -h 101 -i optimal1.txt -o optimal3.png");

    //    DoubleVector f(hc->M+1);
    //    for (unsigned int k=0; k<=M; k++) f[k] = x[2*hc->L+0*(M+1)+k];
    //    Printer::printVector(f, 10, "g1");
    //    for (unsigned int k=0; k<=M; k++) f[k] = x[2*hc->L+1*(M+1)+k];
    //    Printer::printVector(f, 10, "g2");
    //    for (unsigned int k=0; k<=M; k++) f[k] = x[2*hc->L+2*(M+1)+k];
    //    Printer::printVector(f, 10, "g3");

    //    DoubleVector fg(hc->M+1);
    //    for (unsigned int k=0; k<=M; k++) fg[k] = g[2*hc->L+0*(M+1)+k];
    //    Printer::printVector(fg, 10, "fg1");
    //    for (unsigned int k=0; k<=M; k++) fg[k] = g[2*hc->L+1*(M+1)+k];
    //    Printer::printVector(fg, 10, "fg2");
    //    for (unsigned int k=0; k<=M; k++) fg[k] = g[2*hc->L+2*(M+1)+k];
    //    Printer::printVector(fg, 10, "fg3");
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");

    //    hc->calculateU(e, hc->uT);
    //    char buffer [12];
    //    int n=sprintf (buffer, "file%d.txt", i);
    //    hc->write(buffer, hc->uT);
}

void HeatControl2Delta::project(DoubleVector &e, int index)
{
    if (index<6)
    {
        if (e[index] > 1.0) e[index] = 1.0;
        if (e[index] < 0.0) e[index] = 0.0;
    }
    //    for (unsigned int i=0; i<e.size(); i++)
    //    {
    //        if (e[i]>1.0) e[i]=1.0;
    //        if (e[i]<0.0) e[i]=0.0;
    //    }
}

void HeatControl2Delta::write(const char *fileName, const DoubleMatrix& m)
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
