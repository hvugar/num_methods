#include "heatcontrol.h"

void HeatControl::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    /* Function */
    HeatControl hc;

    DoubleVector f0((hc.M+1)*(hc.N+1));
    for (unsigned int i=0; i<f0.length(); i++)
    {
        f0[i] = sin(i*hc.hx);
    }

    /* Minimization */
    ConjugateGradient g;
    g.setGradient(&hc);
    g.setFunction(&hc);
    g.setEpsilon1(0.0000001);
    g.setEpsilon2(0.0000001);
    g.setEpsilon3(0.0000001);
    g.setR1MinimizeEpsilon(0.1, 0.0000001);
    g.setPrinter(&hc);
    g.setNormalize(true);
    g.calculate(f0);

    IPrinter::printAsMatrix(f0, hc.M, hc.N);
}

HeatControl::HeatControl()
{
    N = 100;
    M = 100;
    hx = 0.01;
    ht = 0.01;
    a  = 1.0;

    // initialize U
    U.resize(N+1);
    for (unsigned int i=0; i<=N; i++) U[i] = u(i*hx, 1.0);
}

double HeatControl::fx(const DoubleVector &f) const
{
    const_cast<HeatControl*>(this)->pf = &f;
    DoubleVector u;
    IParabolicEquation::calculateU(u, hx, ht, N, M, a);

    double sum = 0.0;
    sum += 0.5*(u[0] - U[0])*(u[0] - U[0]);
    for (unsigned int i=1; i<=N-1; i++)
    {
        sum += (u[i] - U[i])*(u[i] - U[i]);
    }
    sum += 0.5*(u[N] - U[N])*(u[N] - U[N]);
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

void HeatControl::print(unsigned int i, const DoubleVector&, const DoubleVector &, double fx, GradientMethod::MethodResult result) const
{
    C_UNUSED(result);
    printf("J[%d]: %.14f\n", i, fx);
}

double HeatControl::u(double x, double t) const
{
    return x*x+t*t;
}

double HeatControl::fxt(double x, double t) const
{
    C_UNUSED(x);
    return 2.0*t - 2.0;//*a;
}
