#include "hyperboliccontrol1d.h"
#include <gradient_cjt.h>

void HyperbolicControl1D::main()
{
    HyperbolicControl1D hc;

    DoubleVector v(2*(hc.M+1));
    for (unsigned int j=0; j<=hc.M; j++)
    {
        double t = j*hc.ht;
        v[0*(hc.M+1)+j] = 2.0*t;
        v[1*(hc.M+1)+j] = 2.0*t + 2.0;
    }

    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setGradient(&hc);
    g2.setEpsilon1(0.0001);
    g2.setEpsilon2(0.0001);
    g2.setEpsilon3(0.0001);
    g2.setR1MinimizeEpsilon(0.1, 0.000001);
    g2.setPrinter(&hc);
    g2.setProjection(&hc);
    g2.setNormalize(false);
    g2.calculate(v);

    DoubleVector v1(hc.M+1); for (unsigned j=0; j<=hc.M; j++) v1[j] = v[j];
    DoubleVector v2(hc.M+1); for (unsigned j=0; j<=hc.M; j++) v2[j] = v[hc.M+1+j];
    IPrinter::printVector(v1);
    IPrinter::printVector(v2);
}

HyperbolicControl1D::HyperbolicControl1D() : RnFunction(), IPrinter()
{
    t0 = 0.0; t1 = 1.0;
    x0 = 0.0; x1 = 1.0;
    a = 1.0;
    M  = 100;
    N  = 100;
    ht = (t1-t0)/M;
    hx = (x1-x0)/N;
    lamda = 0.25;

    DoubleVector v(2*(M+1));
    for (unsigned int j=0; j<=M; j++)
    {
        double t = j*ht;
        v[0*(M+1)+j] = v1(t);
        v[1*(M+1)+j] = v2(t);
    }
    pv = &v;
    IHyperbolicEquation::calculateU(U, hx, ht, M, N);
    puts("-------------------------");
    IPrinter::printVector(U);
    puts("-------------------------");
}

double HyperbolicControl1D::fx(const DoubleVector& v) const
{
    const_cast<HyperbolicControl1D*>(this)->pv = &v;
    DoubleVector u;
    IHyperbolicEquation::calculateU(u, hx, ht, M, N);

    double sum = 0.0;
    for (unsigned int j=0; j<=M; j++)
    {
        for (unsigned int i=0; i<=N; i++)
        {
            double alpha = 1.0;
            if (i==0 || i==N || j==M || j==0) alpha = 0.5;
            if (i==0 && j==M) alpha = 0.25;
            if (i==N && j==M) alpha = 0.25;
            sum += alpha*(u[i]-U[i])*(u[i]-U[i]);
        }
    }
    sum = hx*ht*sum;

    double norm1 = 0.0;
    for (unsigned int j=0; j<=M; j++)
    {
        double alpha = 1.0;
        if (j==0 || j==M) alpha = 0.5;
        norm1 += alpha*(v[0*(M+1)+j] - v1(j*ht))*(v[0*(M+1)+j] - v1(j*ht));
    }
    norm1 = ht*norm1;

    double norm2 = 0.0;
    for (unsigned int j=0; j<=M; j++)
    {
        double alpha = 1.0;
        if (j==0 || j==M) alpha = 0.5;
        norm2 += alpha*(v[1*(M+1)+j] - v2(j*ht))*(v[1*(M+1)+j] - v2(j*ht));
    }
    norm2 = ht*norm2;

    return sum + norm1 + norm2;
}

void HyperbolicControl1D::gradient(const DoubleVector& v, DoubleVector& g)
{
    pv = &v;
    DoubleVector u;
    IHyperbolicEquation::calculateU(u, hx, ht, M, N);

    pu = &u;
    DoubleMatrix psi;
    IBackwardHyperbolicEquation::calculateU(psi, hx, ht, M, N);

    for (unsigned int j=0; j<=M; j++)
    {
        double t = j*ht;
        g[0*(M+1)+j] = -(a*a)*(psi[j][1]  -psi[j][0])/ht + 2.0*(v[0*(M+1)+j] - v1(t));
        g[1*(M+1)+j] = -(a*a)*(psi[j][N-1]-psi[j][N])/ht + 2.0*(v[1*(M+1)+j] - v2(t));
    }
}

double HyperbolicControl1D::initial1(unsigned int i) const
{
    double x = i*hx;
    return x*x;
}

double HyperbolicControl1D::initial2(unsigned int i) const
{
    C_UNUSED(i);
    return 0.0;
}

double HyperbolicControl1D::boundary(Boundary type, unsigned int j) const
{
    const DoubleVector &v = *pv;
    if (type==Left)  return v[0*(M+1)+j];
    if (type==Right) return v[1*(M+1)+j];
    return 0.0;
}

double HyperbolicControl1D::f(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl1D::binitial1(unsigned int i) const
{
    C_UNUSED(i);
    return 0.0;
}

double HyperbolicControl1D::binitial2(unsigned int i) const
{
    const DoubleVector &u = *pu;
    return 2.0*(u[i] - U[i]);
}

double HyperbolicControl1D::bboundary(Boundary type, unsigned int j) const
{
    C_UNUSED(type);
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl1D::bf(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

void HyperbolicControl1D::print(unsigned int i, const DoubleVector& v, const DoubleVector &g, double fx) const
{
    C_UNUSED(g);
    HyperbolicControl1D *hc = const_cast<HyperbolicControl1D*>(this);
    printf("J[%d]: %.16f\n", i, hc->fx(v));
}

void HyperbolicControl1D::HyperbolicControl1D::project(DoubleVector &x, int index)
{
    C_UNUSED(x);
    C_UNUSED(index);
}
