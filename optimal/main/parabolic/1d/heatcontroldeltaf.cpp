#include "heatcontroldeltaf.h"
#include <tomasmethod.h>
#include <gradient_cjt.h>
#include <math.h>

void HeatControlDeltaF::main()
{
    /* Function */
    HeatControlDeltaF hc(1000, 10, 1.0);

    DoubleVector f;
    f.resize(hc.C);
    for (unsigned int i=0; i<hc.C; i++)
    {
//        unsigned int k = i/((hc.N1+1)*(hc.N2+1));
//        double t = k*hc.ht;
        f[i] = 20.0;//2.0*t - 4.0;
    }

    /* Minimization */
    ConjugateGradient g;
    g.setGradient(&hc);
    g.setFunction(&hc);
    g.setEpsilon1(0.0000001);
    g.setEpsilon2(0.0000001);
    g.setR1MinimizeEpsilon(0.1, 0.0000001);
    g.setPrinter(&hc);
    g.setNormalize(false);
    g.calculate(f);
}

HeatControlDeltaF::HeatControlDeltaF(unsigned int M, unsigned int N, double a1)
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
    this->L = 1;
    this->C =L * (M+1);

    initializeU();
}

double HeatControlDeltaF::fx(const DoubleVector &f)
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
    for (unsigned int j=0; j<M; j++)
    {
        unsigned int j1 = j;
        unsigned int j2 = j+1;

        double f11 = f[j1] - f1(j1*ht);
        double f12 = f[j2] - f1(j1*ht);

        norm += f11*f11 + f12*f12;
    }
    norm = (ht)*0.5*norm;

    return sum + norm;
}

void HeatControlDeltaF::gradient(const DoubleVector &f, DoubleVector &g)
{
    calculateU(f, uT);
    calculateP(f, g);
}

void HeatControlDeltaF::calculateU(const DoubleVector &f, DoubleVector &u)
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
                d[i-1] = u1[i] + ht * fxt(i, j, f);
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

void HeatControlDeltaF::calculateP(const DoubleVector &f, DoubleVector &g)
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

        //calculateG(f, psi, g, j);

        g[j] = 0.0;
        for (unsigned int i=0; i<=N; i++)
        {
            if (fabs(E[0]-i*hx)<hx)
                g[j] += psi[i] * ((hx-fabs(E[0]-i*hx))/hx);
        }
        g[j] = -g[j] + 2.0*(f[j]-f1(j*ht));
    }

    a.clear();
    b.clear();
    c.clear();
    d.clear();
    x.clear();

    psi.clear();
}

void HeatControlDeltaF::calculateG(const DoubleVector& f, const std::vector<DoubleVector>& psi, DoubleVector& g, unsigned int j)
{
}

double HeatControlDeltaF::fxt(unsigned int i, unsigned int j, const DoubleVector& f)
{
    double x = i*hx;
    //double t  = k*ht;
    double sum = 0.0;
    if (fabs(x-E[0])<hx)
    {
        sum += (1.0/hx) * f[j] * ((hx-fabs(x-E[0]))/hx);
    }
    return sum;
}

void HeatControlDeltaF::initializeU()
{
    DoubleVector f;
    f.resize((M+1)*L);

    E.resize(L);
    E[0] = 0.2;

    for (unsigned int j=0; j<=M; j++)
    {
        f[0*(M+1)+j] = f1(j*ht);
    }

    U.resize(N+1);

    calculateU(f, U);

    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
    IPrinter::printVector(U, 0, N/10);
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");

}

void HeatControlDeltaF::print(unsigned int i, const DoubleVector &f0, const DoubleVector &s, double a, RnFunction *f) const
{
    HeatControlDeltaF *hc = dynamic_cast<HeatControlDeltaF*>(f);
    printf("J: %.16f\n", hc->fx(f0));
    printf("Printing f: ");
    for (unsigned int j=0; j<=hc->M; j++)
    {
        if (j%(hc->M/10)==0)
        {
            printf("%12.6f", f0[j]);
        }
    }
    puts("");
    printf("Printing g: ");
    for (unsigned int j=0; j<=hc->M; j++)
    {
        if (j%(hc->M/10)==0)
        {
            printf("%12.6f", s[j]);
        }
    }
    puts("");
}

