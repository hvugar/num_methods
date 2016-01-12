#include "heatcontrol1.h"

void HeatControl1::main()
{
    HeatControl1 hc;

    DoubleVector f0((hc.M+1)*(hc.N+1));
    for (unsigned int k=0; k<f0.size(); k++)
    {
        //unsigned int j = k / (hc.M+1);
        //unsigned int i = k % (dh.N+1);
        //double t = j * hc.ht;
        //f0[k] = 2.0*t - 2.0*hc.a;
        f0[k] = 10.0;
    }

    Printer::printAsMatrix(f0, hc.M, hc.N);
    printf("----\n");

    ConjugateGradient g;
    g.setGradient(&hc);
    g.setFunction(&hc);
    g.setEpsilon1(0.0000001);
    g.setEpsilon2(0.0000001);
    g.setR1MinimizeEpsilon(0.1, 0.0000001);
    g.setPrinter(&hc);
    g.calculate(f0);

    Printer::printAsMatrix(f0, hc.M, hc.N);
}

HeatControl1::HeatControl1()
{
    t0 = 0.0;
    t1 = 1.0;
    x0 = 0.0;
    x1 = 1.0;
    M = 100;
    N = 100;
    ht = (t1-t0)/M;
    hx = (x1-x0)/N;
    a = 1.0;

    U.resize(N+1);
    for (unsigned int i=0; i<U.size(); i++)
    {
        double x = i*hx;
        U[i] = x*x + 1.0;
    }
}

double HeatControl1::fx(const DoubleVector &f)
{
    double sum  = 0.0;

    DoubleVector u;
    pf = &f;
    calculateU(u, hx, ht, N, M, a);

    double alpha;
    for (unsigned int i=0; i<u.size(); i++)
    {
        alpha = 1.0;
        if (i==0 || i==u.size()-1) alpha = 0.5;
        sum += alpha*(u[i] - U[i])*(u[i] - U[i]);
    }
    sum = hx*sum;

    double norm = 0.0;
    for (unsigned int j=0; j<=M; j++)
    {
        for (unsigned int i=0; i<=N; i++)
        {
            double b = 1.0;
            if (i==0 || i==N || j==0 || j==M) b = 0.5;
            if ((i==0 && j==0) || (i==0 && j==M) || (i==N && j==0) || (i==N && j==M)) b = 0.25;
            double f1 = ((*pf)[j*(N+1)+i] - fxt(i*hx, j*ht));
            norm += b*f1*f1;
        }
    }
    norm = hx*ht*norm;

    return sum+norm;
}

void HeatControl1::gradient(const DoubleVector &f, DoubleVector &g)
{
    DoubleVector u;
    pf = &f;
    calculateU(u, hx, ht, N, M, a);
    calculateP(u, g);
}

void HeatControl1::print(unsigned int iteration, const DoubleVector &f, const DoubleVector &gradient, double alpha, RnFunction *fn) const
{
    //HeatControl1* dh =dynamic_cast<HeatControl1*>(fn);
    printf("J[%d]: %.16f\n", iteration, fn->fx(f));
    //Printer::printAsMatrix(f, dh->M, dh->N);
    //puts("-------");
}

double HeatControl1::fi(unsigned int i) const
{
    double x = hx*i;
    return x*x;
}

double HeatControl1::m1(unsigned int j) const
{
    double t = j*ht;
    return t*t;
}

double HeatControl1::m2(unsigned int j) const
{
    double t = j*ht;
    return t*t+1.0;
}

double HeatControl1::f(unsigned int i, unsigned int j) const
{
    return (*pf)[j*(N+1)+i];
}

void HeatControl1::calculateP(const DoubleVector &u, DoubleVector &g)
{
    DoubleVector psi;
    psi.resize(N+1);

    DoubleVector a1;
    DoubleVector b1;
    DoubleVector c1;
    DoubleVector d1;
    DoubleVector x1;

    a1.resize(N-1);
    b1.resize(N-1);
    c1.resize(N-1);
    d1.resize(N-1);
    x1.resize(N-1);

    double alpha = -(a*ht)/(hx*hx);
    double beta  = 1.0 + (2.0*a*ht)/(hx*hx);

    for (unsigned int k=0; k<=M; k++)
    {
        unsigned int j = M-k;

        if (j == M)
        {
            for (unsigned i=0; i<=N; i++)
            {
                psi[i] = -2.0 * (u[i] - U[i]);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                a1[i-1] = alpha;
                b1[i-1] = beta;
                c1[i-1] = alpha;
                d1[i-1] = psi[i];
            }

            a1[0]   = 0.0;
            c1[N-2] = 0.0;
            d1[0]   -= alpha * pm1(j*ht);
            d1[N-2] -= alpha * pm2(j*ht);

            TomasAlgorithm(a1, b1, c1, d1, x1);

            psi[0] = pm1(j*ht);
            for (unsigned int i=1; i<=N-1; i++)
            {
                psi[i] = x1[i-1];
            }
            psi[N] = pm2(j*ht);
        }

        calculateG((*pf), psi, g, j);
    }

    a1.clear();
    b1.clear();
    c1.clear();
    d1.clear();
    x1.clear();

    psi.clear();
}

void HeatControl1::calculateG(const DoubleVector &f, const DoubleVector &psi, DoubleVector &g, unsigned int j)
{
    for (unsigned int i=0; i<=N; i++)
    {
        int k = j*(N+1)+i;
        g[k] = -psi[i] + 2.0*(f[k] - fxt(i*hx, j*ht));
    }
}
