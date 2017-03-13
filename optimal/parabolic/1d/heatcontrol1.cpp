#include "heatcontrol1.h"

void HeatControl1::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    /* Function */
    HeatControl1 hc;

    DoubleVector f0((hc.M+1)*(hc.N+1));
    for (unsigned int i=0; i<f0.size(); i++)
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

HeatControl1::HeatControl1()
{
    N = 100;
    M = 100;
    hx = 0.01;
    ht = 0.01;

    // initialize U
    U.resize(N+1);
    for (unsigned int i=0; i<=N; i++) U[i] = u(i*hx, 1.0);

    forward.setTimeDimension(Dimension(ht,M,0));
    forward.addSpaceDimension(Dimension(hx,N,0));

    backward.setTimeDimension(Dimension(ht,M,0));
    backward.addSpaceDimension(Dimension(hx,N,0));
    backward.pU = &U;
}

double HeatControl1::fx(const DoubleVector &f) const
{
    const_cast<HeatControl1*>(this)->forward.pf = &f;
    DoubleVector u;
    forward.gridMethod(u);

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

void HeatControl1::gradient(const DoubleVector &f, DoubleVector &g)
{
    DoubleVector u;
    forward.pf = &f;
    forward.gridMethod(u);

    DoubleVector p;
    DoubleMatrix psi(M+1, N+1);
    backward.pp = &psi;
    backward.pu = &u;
    backward.gridMethod(p);

    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            int k = m*(N+1)+n;
            g[k] = -psi[m][n] + 2.0*(f[k] - fxt(n*hx, m*ht));
        }
    }
}

void HeatControl1::print(unsigned int i, const DoubleVector &, const DoubleVector &, double fx, GradientMethod::MethodResult result) const
{
    printf("J[%d]: %.14f\n", i, fx);
}

double HeatControl1::u(double x, double t) const
{
    return x*x+t*t;
}

double HeatControl1::fxt(double x, double t) const
{
    return 2.0*t - 2.0;//*a;
}

double HeatControl1::CParabolicIBVP::initial(const SpaceNode &sn) const
{
    return sn.x*sn.x;
}

double HeatControl1::CParabolicIBVP::boundary(const SpaceNode &, const TimeNode &tn, BoundaryType boundary) const
{
    double t = tn.t;
    if (boundary==Left)  return t*t;
    if (boundary==Right) return t*t + 1.0;
    return 0.0;
}

double HeatControl1::CParabolicIBVP::f(const SpaceNode &sn, const TimeNode &tn) const
{
    Dimension dim1 = spaceDimension(Dimension::Dim1);
    unsigned int N = dim1.sizeN();
    return (*pf)[tn.i*(N+1)+sn.i];
}

double HeatControl1::CParabolicIBVP::a(const SpaceNode&, const TimeNode&) const
{
    return +1.0;
}

double HeatControl1::CBackwardParabolicIBVP::initial(const SpaceNode &sn) const
{
    return -2.0 * ((*pu)[sn.i] - (*pU)[sn.i]);
}

double HeatControl1::CBackwardParabolicIBVP::boundary(const SpaceNode&, const TimeNode&, BoundaryType) const
{
    return 0.0;
}

double HeatControl1::CBackwardParabolicIBVP::f(const SpaceNode&, const TimeNode&) const
{
    return 0.0;
}

double HeatControl1::CBackwardParabolicIBVP::a(const SpaceNode&, const TimeNode&) const
{
    return +1.0;
}

void HeatControl1::CBackwardParabolicIBVP::layerInfo(const DoubleVector & p, unsigned int m) const
{
    Dimension dim1 = spaceDimension(Dimension::Dim1);
    unsigned int N = dim1.sizeN();
    for (unsigned int n=0; n<=N; n++)
    {
        (*pp)[m][n] = p[n];
    }
}
