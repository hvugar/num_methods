#include "heatcontrol.h"

void HeatControl::main()
{
    /* Function */
    HeatControl hc;

    DoubleVector f0((hc.M+1)*(hc.N+1));
    for (unsigned int i=0; i<f0.size(); i++)
    {
        //int j = i/(hc.N+1);
        //f0[i] = 2.0*j*hc.ht - 2.0;
        f0[i] = 2.0;
    }

    /* Minimization */
    ConjugateGradient g2;
    g2.setGradient(&hc);
    g2.setFunction(&hc);
    g2.setEpsilon1(0.0000001);
    g2.setEpsilon2(0.0000001);
    g2.setEpsilon3(0.0000001);
    g2.setR1MinimizeEpsilon(0.1, 0.0000001);
    g2.setPrinter(&hc);
    g2.setNormalize(true);
    g2.calculate(f0);

    IPrinter::printAsMatrix(f0, hc.M, hc.N);
}

HeatControl::HeatControl()
{
    this->t0 = 0.0;
    this->t1 = 1.0;
    this->x0 = 0.0;
    this->x1 = 1.0;
    this->a  = 1.0;

    N = 100;
    M = 100;
    this->hx  = (x1-x0)/N;
    this->ht = (t1-t0)/M;
    // initialize U
    U.resize(N+1);
    for (unsigned int i=0; i<=N; i++) U[i] = u(i*hx, t1);
}

double HeatControl::fx(const DoubleVector &f)
{
    pf = &f;
    DoubleVector u;
    IParabolicEquation::calculateU(u, hx, ht, N, M, a);

    double sum = 0.0;
    double alpha;
    for (unsigned int i=0; i<=N; i++)
    {
        alpha = 1.0;
        if (i==0 || i==N) alpha = 0.5;
        sum += alpha*(u[i] - U[i])*(u[i] - U[i]);
    }
    sum = hx*sum;

    double norm = 0.0;
    for (unsigned int j=0; j<=M; j++)
    {
        for (unsigned int i=0; i<=N; i++)
        {
            double alpha = 1.0;
            if (i==0 || i==N || j==0 || j==M) alpha = 0.5;
            if (i==0 && j==0) alpha = 0.25;
            if (i==0 && j==M) alpha = 0.25;
            if (i==N && j==0) alpha = 0.25;
            if (i==N && j==M) alpha = 0.25;
            norm += alpha*(f[j*(N+1) + i] - fxt(i*hx, j*ht))*(f[j*(N+1) + i] - fxt(i*hx, j*ht));
        }
    }
    norm = hx*ht*norm;

    return sum + norm;
}

void HeatControl::gradient(const DoubleVector &f, DoubleVector &g)
{
    pf = &f;
    DoubleVector u;
    IParabolicEquation::calculateU(u, hx, ht, N, M, a);

    pu = &u;
    DoubleMatrix psi;
    IBackwardParabolicEquation::calculateU(psi, hx, ht, N, M, a);

    for (unsigned int j=0; j<=M; j++)
    {
        for (unsigned int i=0; i<=N; i++)
        {
            int k = j*(N+1)+i;
            g[k] = -psi[j][i] + 2.0*(f[k] - fxt(i*hx, j*ht));
        }
    }
}

double HeatControl::fi(unsigned int i) const
{
    double x = i*hx;
    return u(x, t0);
}

double HeatControl::m1(unsigned int j) const
{
    double t = j*ht;
    return u(x0, t);
}

double HeatControl::m2(unsigned int j) const
{
    double t = j*ht;
    return u(x1, t);
}

double HeatControl::f(unsigned int i, unsigned int j) const
{
    return (*pf)[j*(N+1)+i];
}

double HeatControl::bfi(unsigned int i) const
{
    return -2.0 * ((*pu)[i] - U[i]);;
}

double HeatControl::bm1(unsigned int j) const
{
    C_UNUSED(j);
    return 0.0;
}

double HeatControl::bm2(unsigned int j) const
{
    C_UNUSED(j);
    return 0.0;
}

double HeatControl::bf(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

void HeatControl::print(unsigned int i, const DoubleVector &f0, const DoubleVector &g, double a, RnFunction *f) const
{
    C_UNUSED(g);
    C_UNUSED(a);

    HeatControl *hc = dynamic_cast<HeatControl*>(f);
    printf("J[%d]: %.20f\n", i, hc->fx(f0));
    //    printf("Printing f-------------------------------\n");
    //    IPrinter::printAsMatrix(f0, M, N);
    //    printf("Printing g-------------------------------\n");
    //    IPrinter::printAsMatrix(g, M, N);
}

double HeatControl::u(double x, double t) const
{
    return x*x+t*t;
}

double HeatControl::fxt(double x, double t)
{
    C_UNUSED(x);
    return 2.0*t - 2.0*a;
}

