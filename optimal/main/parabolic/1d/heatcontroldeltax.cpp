#include "heatcontroldeltax.h"
#include <gradient_cjt.h>
#include <math.h>
#include <tomasmethod.h>

void HeatControlDeltaX::main()
{
    /* Function */
    HeatControlDeltaX hc(1000, 10, 1.0);

    DoubleVector E;
    E.resize(hc.C);
    for (unsigned int i=0; i<hc.C; i++)
    {
        E[i] = 0.4;
    }
//    E[0] = 0.2;
//    E[1] = 0.8;

    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setEpsilon1(0.0000001);
    g2.setEpsilon2(0.0000001);
    g2.setGradientStep(0.000001);
    g2.setR1MinimizeEpsilon(10.0, 1.0);
    g2.setPrinter(&hc);
    g2.setProjection(&hc);
    g2.setNormalize(false);
    g2.calculate(E);
}

HeatControlDeltaX::HeatControlDeltaX(unsigned int M, unsigned int N, double a1)
{
    t0 = 0.0;
    t1 = 1.0;
    x0 = 0.0;
    x1 = 1.0;
    this->M = M;
    this->N = N;
    this->a1 = a1;
    ht = (t1-t0)/M;
    hx = (x1-x0)/N;
    this->L = 2;
    this->C = L;

    initializeU();
}

double HeatControlDeltaX::fx(const DoubleVector &f)
{
    calculateU(f, uT);

    double sum = 0.0;
    for (unsigned int i=0; i<N; i++)
    {
        int j = i+1;

        double f1 = uT[j] - U[j];
        double f2 = uT[i] - U[i];

        sum += (f1*f1+f2*f2);
    }
    sum = 0.5*hx*sum;

    double norm = 0.0;
    //    for (unsigned int j=0; j<M; j++)
    //    {
    //        unsigned int j1 = j;
    //        unsigned int j2 = j+1;

    //        double f11 = f[j1] - f1(j1*ht);
    //        double f12 = f[j2] - f1(j1*ht);

    //        norm += f11*f11 + f12*f12;
    //    }
    //    norm = (ht)*0.5*norm;

    return sum + norm;
}

void HeatControlDeltaX::gradient(const DoubleVector &f, DoubleVector &g, double gradient_step)
{
    calculateU(f, uT);
    calculateP(f, g);
}

void HeatControlDeltaX::calculateU(const DoubleVector &E, DoubleVector &u)
{
    DoubleVector u1;
    u1.resize(N+1);

    DoubleVector a;
    DoubleVector b;
    DoubleVector c;
    DoubleVector d;
    DoubleVector x;

    a.resize(N-1);
    b.resize(N-1);
    c.resize(N-1);
    d.resize(N-1);
    x.resize(N-1);

    double alpha = -(a1*ht)/(hx*hx);
    double beta  = 1.0 + (2.0*a1*ht)/(hx*hx);

    for (unsigned int j=0; j<=M; j++)
    {
        if (j == 0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u1[i] = fi(i*hx);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                a[i-1] = alpha;
                b[i-1] = beta;
                c[i-1] = alpha;
                d[i-1] = u1[i] + ht * fxt(i, j, E);
            }

            a[0]   = 0.0;
            c[N-2] = 0.0;
            d[0]   -= alpha * m1(j*ht);
            d[N-2] -= alpha * m2(j*ht);
            TomasAlgorithm(a, b, c, d, x);

            u1[0] = m1(j*ht);
            for (unsigned int i=1; i<=N-1; i++)
            {
                u1[i] = x[i-1];
            }
            u1[N] = m2(j*ht);
        }
    }

    a.clear();
    b.clear();
    c.clear();
    d.clear();
    x.clear();

    u = u1;

    u1.clear();
}

void HeatControlDeltaX::calculateP(const DoubleVector &E, DoubleVector &g)
{
    DoubleVector psi;
    psi.resize(N+1);

    DoubleVector a;
    DoubleVector b;
    DoubleVector c;
    DoubleVector d;
    DoubleVector x;

    a.resize(N-1);
    b.resize(N-1);
    c.resize(N-1);
    d.resize(N-1);
    x.resize(N-1);

    double alpha = -(a1*ht)/(hx*hx);
    double beta  = 1.0 + (2.0*a1*ht)/(hx*hx);

    for (unsigned int k=0; k<=M; k++)
    {
        unsigned int j = M-k;

        if (j == M)
        {
            for (unsigned i=0; i<=N; i++)
            {
                psi[i] = -2.0 * (uT[i] - U[i]);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                a[i-1] = alpha;
                b[i-1] = beta;
                c[i-1] = alpha;
                d[i-1] = psi[i];
            }

            a[0]   = 0.0;
            c[N-2] = 0.0;
            d[0]   -= alpha * pm1(j*ht);
            d[N-2] -= alpha * pm2(j*ht);

            TomasAlgorithm(a, b, c, d, x);

            psi[0] = pm1(j*ht);
            for (unsigned int i=1; i<=N-1; i++)
            {
                psi[i] = x[i-1];
            }
            psi[N] = pm2(j*ht);
        }

        double p = 1.0;
        if (k==M || k==0)
            p = 1.0 * -(ht/2.0);
        else
            p = 2.0 * -(ht/2.0);
        double t = j*ht;

        double psiX1=0.0;

        for (unsigned int i=0; i<=N-1; i++)
        {
            if (i*hx<=E[0] && E[0]<(i+1)*hx)
            {
                psiX1 = (psi[i+1]*((fabs(E[0]-i*hx))/hx) - psi[i]*((hx-fabs(E[0]-i*hx)/hx)))/hx;
            }
        }
        g[0] = g[0] + p * f1(t) * psiX1;

        for (unsigned int i=0; i<=N-1; i++)
        {
            if (i*hx<=E[1] && E[1]<(i+1)*hx)
            {
                psiX1 = (psi[i+1]*((fabs(E[1]-i*hx))/hx) - psi[i]*((hx-fabs(E[1]-i*hx)/hx)))/hx;
            }
        }
        g[1] = g[1] + p * f2(t) * psiX1;
    }

    a.clear();
    b.clear();
    c.clear();
    d.clear();
    x.clear();

    psi.clear();
}

void HeatControlDeltaX::calculateG(const DoubleVector& f, const std::vector<DoubleVector>& psi, DoubleVector& g, unsigned int j)
{
}

double HeatControlDeltaX::fxt(unsigned int i, unsigned int j, const DoubleVector& E)
{
    double x = i*hx;
    //double t  = k*ht;
    double sum = 0.0;
    if (fabs(x-E[0])<hx)
    {
        sum += (1.0/hx) * f1(j*ht) * ((hx-fabs(x-E[0]))/hx);
    }
    if (fabs(x-E[1])<hx)
    {
        sum += (1.0/hx) * f2(j*ht) * ((hx-fabs(x-E[1]))/hx);
    }
    return sum;
}

void HeatControlDeltaX::initializeU()
{
    DoubleVector E;
    E.resize(L);
    E[0] = 0.2;
    E[1] = 0.8;

    U.resize(N+1);
    calculateU(E, U);

    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
    Printer::printVector(U, 0, N/10);
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
}


void HeatControlDeltaX::print(unsigned int i, const DoubleVector &e, const DoubleVector &s, double a, RnFunction *f) const
{
    HeatControlDeltaX *hc = dynamic_cast<HeatControlDeltaX*>(f);
    printf("J: %.16f\n", hc->fx(e));
    printf("Printing e[0]: %.12f e[1]: %.12f\n", e[0], e[1]);
    printf("Printing g[0]: %.12f g[1]: %.12f\n", s[0], s[1]);
}

void HeatControlDeltaX::project(DoubleVector &e, int index)
{
    for (unsigned int l=0; l<L; l++)
    {
        if (e[l] < x0) e[l] = x0;
        if (e[l] > x1) e[l] = x1;
    }
}
