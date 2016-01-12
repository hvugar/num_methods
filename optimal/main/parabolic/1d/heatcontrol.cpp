#include "heatcontrol.h"
#include <tomasmethod.h>

HeatControl::HeatControl()
{
    this->t0 = 0.0;
    this->t1 = 1.0;
    this->x0 = 0.0;
    this->x1 = 1.0;
    this->a1 = 1.0;

    N = 100;
    M = 100;
    C = (M+1)*(N+1);

    this->hx  = (x1-x0)/N;
    this->ht = (t1-t0)/M;
    initializeU();
}

double HeatControl::fx(const DoubleVector &f)
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
        for (unsigned int i=0; i<N; i++)
        {
            unsigned int i1 = i;
            unsigned int i2 = i+1;
            unsigned int j1 = j;
            unsigned int j2 = j+1;

            double f11 = f[j1*(N+1) + i1] - fxt(i1*hx, j1*ht);
            double f12 = f[j1*(N+1) + i2] - fxt(i2*hx, j1*ht);
            double f21 = f[j2*(N+1) + i1] - fxt(i1*hx, j2*ht);
            double f22 = f[j2*(N+1) + i2] - fxt(i2*hx, j2*ht);

            norm += f11*f11 + f12*f12 + f21*f21 + f22*f22;
        }
    }
    norm = (hx*ht)*0.25*norm;

    return sum + norm;
}

void HeatControl::gradient(const DoubleVector &f, DoubleVector &g)
{
    calculateU(f, uT);
    calculateP(f, g);
}

void HeatControl::calculateU(const DoubleVector &f, DoubleVector& u)
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
                //d[i-1] = u1[i] + ht * fxt(i*h, j*dt);
                d[i-1] = u1[i] + ht * f[j*(N+1)+i];
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

void HeatControl::calculateU1(const DoubleVector &f)
{
    DoubleVector u2;
    u2.resize(N+1);
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

    double alpha = ((a1*ht)/(hx*hx));
    double beta  = 1.0 - (2.0*a1*ht)/(hx*hx);

    for (unsigned int j1=0; j1<=M; j1++)
    {
        unsigned int j=M-j1;

        if (j == M)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u1[i] = u(i*hx, 1.0);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                a[i-1] = alpha;
                b[i-1] = beta;
                c[i-1] = alpha;
                d[i-1] = u1[i] + ht * fxt(i*hx, (j+1)*ht);
                //d[i-1] = u1[i] + ht * f[j*(N+1)+i];
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

    u1.clear();
}

void HeatControl::calculateP(const DoubleVector &f, DoubleVector &g)
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

        calculateG(f, psi, g, j);
    }

    a.clear();
    b.clear();
    c.clear();
    d.clear();
    x.clear();

    psi.clear();
}

void HeatControl::calculateG(const DoubleVector& f, const DoubleVector& psi, DoubleVector& g, unsigned int j)
{
    for (unsigned int i=0; i<=N; i++)
    {
        int k = j*(N+1)+i;
        g[k] = -psi[i] + 2.0*(f[k] - fxt(i*hx, j*ht));
    }
}

void HeatControl::initializeU()
{
    U.resize(N+1);
    for (unsigned int i=0; i<=N; i++) U[i] = u(i*hx, t1);
}

void HeatControl::main()
{
    /* Function */
    HeatControl hc;

    DoubleVector f0;
    f0.resize(hc.C);
    for (unsigned int i=0; i<hc.C; i++)
    {
        //int j = i/(hc.N+1);
        //f0[i] = 2.0*j*hc.ht - 2.0;
        f0[i] = 2.0;
    }
    //DoubleVector u;
    //hc.calculateU(f0, u);
    //hc.calculateU1(f0);
    //return;

    /* Minimization */
    //SteepestDescentGradient g2;
    ConjugateGradient g2;
    g2.setGradient(&hc);
    g2.setFunction(&hc);
    g2.setEpsilon1(0.0000001);
    g2.setEpsilon2(0.0000001);
    g2.setR1MinimizeEpsilon(0.1, 0.0000001);
    g2.setPrinter(&hc);
    g2.setNormalize(true);
    g2.calculate(f0);

    for (unsigned int j=0; j<=hc.M; j++)
    {
        if (j%(hc.M/10)==0)
        {
            printf("%6d|", j);
            for (unsigned int i=0; i<=hc.N; i++)
            {
                if (i%(hc.N/10)==0)
                    printf("%12.6f", f0[j*(hc.N+1)+i]);
            }
            puts("");
        }
    }
}

void HeatControl::print(unsigned int i, const DoubleVector &f0, const DoubleVector &s, double a, RnFunction *f) const
{
    HeatControl *hc = dynamic_cast<HeatControl*>(f);
    printf("J: %.16f\n", hc->fx(f0));
    printf("Printing f-------------------------------\n");
    for (unsigned int j=0; j<=hc->M; j++)
    {
        if (j%(hc->M/10)==0)
        {
            printf("%6d|", j);
            for (unsigned int i=0; i<=hc->N; i++)
            {
                if (i%(hc->N/10)==0)
                    printf("%12.6f", f0[j*(hc->N+1)+i]);
            }
            puts("");
        }
    }
    printf("Printing g-------------------------------\n");
    for (unsigned int j=0; j<=hc->M; j++)
    {
        if (j%(hc->M/10)==0)
        {
            printf("%6d|", j);
            for (unsigned int i=0; i<=hc->N; i++)
            {
                if (i%(hc->N/10)==0)
                    printf("%12.6f", s[j*(hc->N+1)+i]);
            }
            puts("");
        }
    }
}
