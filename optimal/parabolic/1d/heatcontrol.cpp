#include "heatcontrol.h"

void HeatControl::Main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    /* Function */
    HeatControl hc;

    DoubleVector f0((hc.M+1)*(hc.N+1));
    for (unsigned int i=0; i<f0.size(); i++)
    {
        f0[i] = sin(i*hc.hx);
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
    a  = 1.0;

    N = 100;
    M = 100;
    hx = 0.01;
    ht = 0.01;

    // initialize U
    U.resize(N+1);
    for (unsigned int i=0; i<=N; i++) U[i] = u(i*hx, 1.0);
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

double HeatControl::initial(unsigned int i) const
{
    double x = i*hx;
    return x*x;
}

double HeatControl::boundary(Boundary type, unsigned int j) const
{
    double t = j*ht;
    if (type == Left)  return t*t;
    if (type == Right) return t*t+1.0;
    return 0.0;
}

double HeatControl::f(unsigned int i, unsigned int j) const
{
    return (*pf)[j*(N+1)+i];
}

double HeatControl::binitial(unsigned int i) const
{
    return -2.0 * ((*pu)[i] - U[i]);
}

double HeatControl::bboundary(Boundary type, unsigned int j) const
{
    C_UNUSED(type);
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
    printf("J[%d]: %.20f\n", i, f->fx(f0));
}

double HeatControl::u(double x, double t) const
{
    return x*x+t*t;
}

double HeatControl::fxt(double x, double t)
{
    C_UNUSED(x);
    return 2.0*t - 2.0;//*a;
}
