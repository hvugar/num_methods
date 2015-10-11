#include "heatcontrol.h"
#include <tomasmethod.h>

HeatControl::HeatControl()
{
    this->t0 = 0.0;
    this->t1 = 1.0;

    this->x0 = 0.0;
    this->x1 = 1.0;

    N = 100;
    M = 100;
    C = (M+1)*(N+1);

    a1 = 1.0;

    this->h = (x1-x0)/N;
    this->dt = (t1-t0)/M;

    U.resize(N+1);
    for (unsigned int i=0; i<=N; i++) U[i] = u(i*h, t1);
}

double HeatControl::fx(const DoubleVector &f)
{
    calculateU(f);

    double sum = 0.0;
    for (unsigned int i=0; i<N; i++)
    {
        int j = i+1;

        double f1 = uT[j] - U[j];//U(j*h);
        double f2 = uT[i] - U[i];//U(i*h);

        sum += 0.5*h*(f1*f1+f2*f2);
    }

    double norm = 0.0;

    for (unsigned int i=0; i<C; i++)
    {
        norm += (f[i] - fxt((i%(M+1))*h, (i/(N+1))*dt)) * (f[i] - fxt((i%(M+1))*h, (i/(N+1))*dt));
    }

    return sum + norm;
}

void HeatControl::gradient(double step, const DoubleVector &f, DoubleVector &g)
{
    calculateU(f);
    calculateP(f, g);
}

void HeatControl::calculateU(const DoubleVector &f0)
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

    double alpha = -(a1*dt)/(h*h);
    double beta  = 1.0 + (2.0*a1*dt)/(h*h);

    for (unsigned int j=0; j<=M; j++)
    {
        if (j == 0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u1[i] = fi(i*h);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                a[i-1] = alpha;
                b[i-1] = beta;
                c[i-1] = alpha;
                d[i-1] = u1[i] + dt * f0[j*(N+1)+i];//fxt(i*h, j*dt);
            }

            a[0]   = 0.0;
            c[N-2] = 0.0;
            d[0]   -= alpha * m1(j*dt);
            d[N-2] -= alpha * m2(j*dt);
            TomasAlgorithm(a, b, c, d, x);

            u1[0] = m1(j*dt);
            for (unsigned int i=1; i<=N-1; i++)
            {
                u1[i] = x[i-1];
            }
            u1[N] = m2(j*dt);
        }
    }

    a.clear();
    b.clear();
    c.clear();
    d.clear();
    x.clear();

    uT = u1;

    u1.clear();
}

void HeatControl::calculateP(const DoubleVector &f0, DoubleVector &g)
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

    double alpha = -(a1*dt)/(h*h);
    double beta  = 1.0 + (2.0*a1*dt)/(h*h);

    for (unsigned int k=0; k<=M; k++)
    {
        unsigned int j = M-k;

        if (j == M)
        {
            for (unsigned i=0; i<=N; i++)
            {
                psi[i] = -2.0 * (uT[i] - U[i]/*(i*h)*/);
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
            d[0]   -= alpha * pm1(j*dt);
            d[N-2] -= alpha * pm2(j*dt);

            TomasAlgorithm(a, b, c, d, x);

            psi[0] = pm1(j*dt);
            for (unsigned int i=1; i<=N-1; i++)
            {
                psi[i] = x[i-1];
            }
            psi[N] = pm2(j*dt);
        }

        for (unsigned int i=0; i<=N; i++)
        {
            int k = j*(N+1)+i;
            g[k] = -psi[i] + 2.0*(f0[k] - fxt(i*h, j*dt));
        }
    }

    a.clear();
    b.clear();
    c.clear();
    d.clear();
    x.clear();

    psi.clear();
}

double HeatControl::u(double x, double t) const
{
    return x*x + t*t + 2*x;
}

//double HeatControl::U(double x) const
//{
//    return u(x, t1);
//}

double HeatControl::fi(double x)
{
    return u(x, t0);
}

double HeatControl::m1(double t)
{
    return u(x0, t);
}

double HeatControl::m2(double t)
{
    return u(x1, t);
}

double HeatControl::fxt(double x, double t)
{
    return 2.0*t - 2.0;
}

double HeatControl::pfi(double x)
{
    return 0.0;
}

double HeatControl::pm1(double t)
{
    return 0.0;
}

double HeatControl::pm2(double t)
{
    return 0.0;
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
        f0[i] = 2.0;//2.0*j*hc.dt - 2.0;
    }

    /* Minimization */
//    SteepestDescentGradient g1;
//    g1.setFunction(&hc);
//    g1.setEpsilon1(0.0000001);
//    g1.setEpsilon2(0.0000001);
//    g1.setGradientStep(0.000001);
//    g1.setR1MinimizeEpsilon(0.01, 0.0000001);
//    g1.setPrinter(new HeatControlPrinter);
//    g1.calculate(f0);

    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setEpsilon1(0.0000001);
    g2.setEpsilon2(0.0000001);
    g2.setGradientStep(0.000001);
    g2.setR1MinimizeEpsilon(0.01, 0.0000001);
    g2.setPrinter(new HeatControlPrinter);
    g2.setNormalize(false);
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

void HeatControlPrinter::print(unsigned int i, const DoubleVector &f0, const DoubleVector &s, double a, RnFunction *f) const
{
    printf("J: %.16f\n", f->fx(f0));
}
