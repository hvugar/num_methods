#include "heatcontrol1.h"

void HeatControl1::Main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    /* Function */
    HeatControl1 hc;

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

HeatControl1::HeatControl1()
{
    N = 100;
    M = 100;
    hx = 0.01;
    ht = 0.01;

    // initialize U
    U.resize(N+1);
    for (unsigned int i=0; i<=N; i++) U[i] = u(i*hx, 1.0);

    pibvp.setTimeDimension(Dimension(0.01,100,0));
    pibvp.addSpaceDimension(Dimension(0.01,100,0));

    bpibvp.setTimeDimension(Dimension(0.01,100,0));
    bpibvp.addSpaceDimension(Dimension(0.01,100,0));
    bpibvp.pU = &U;
}

double HeatControl1::fx(const DoubleVector &f) const
{
    DoubleVector u;
    const_cast<HeatControl1*>(this)->pibvp.pf = &f;
    pibvp.gridMethod(u);

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

void HeatControl1::gradient(const DoubleVector &f, DoubleVector &g)
{
    DoubleVector u;
    pibvp.pf = &f;
    pibvp.gridMethod(u);

    DoubleMatrix psi;
    bpibvp.pu = &u;
    bpibvp.pU = &U;
    bpibvp.gridMethod(psi);

    for (unsigned int j=0; j<=M; j++)
    {
        for (unsigned int i=0; i<=N; i++)
        {
            int k = j*(N+1)+i;
            g[k] = -psi[j][i] + 2.0*(f[k] - fxt(i*hx, j*ht));
        }
    }
}

void HeatControl1::print(unsigned int i, const DoubleVector& f0, const DoubleVector &g, double fx) const
{
    C_UNUSED(g);
    printf("J[%d]: %.20f\n", i, const_cast<HeatControl1*>(this)->fx(f0));
}

double HeatControl1::u(double x, double t) const
{
    return x*x+t*t;
}

double HeatControl1::fxt(double x, double t) const
{
    C_UNUSED(x);
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
    return (*pf)[tn.i*(dim1.sizeN()+1)+sn.i];
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
