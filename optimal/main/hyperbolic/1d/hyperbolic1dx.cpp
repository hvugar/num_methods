#include "hyperbolic1dx.h"

void Hyperbolic1DX::main()
{
    Hyperbolic1DX hc(100, 100);

    DoubleVector e(hc.L);
    e[0] = 0.30; e[1] = 0.50; e[2] = 0.60;
    //    e[0] = 0.10; e[1] = 0.50; e[2] = 0.80;

    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setGradient(&hc);
    g2.setEpsilon1(0.0000001);
    g2.setEpsilon2(0.0000001);
    g2.setEpsilon3(0.0000001);
    g2.setR1MinimizeEpsilon(1.0, 0.0000001);
    g2.setPrinter(&hc);
    g2.setProjection(&hc);
    g2.setNormalize(true);
    g2.calculate(e);
}

Hyperbolic1DX::Hyperbolic1DX(unsigned int m, unsigned int n)
{
    t0 = 0.0;
    t1 = 1.0;
    x0 = 0.0;
    x1 = 1.0;
    this->M = m;
    this->N = n;
    hx = (x1-x0)/N;
    ht = (t1-t0)/M;
    L = 3;
    a = 1.0;
    lamda = 0.25;

    //initialize
    DoubleVector E(L);
    E[0] = 0.20; E[1] = 0.40; E[2] = 0.70;
    pe = &E;
    IHyperbolicEquation::calculateU(U, hx, ht, M, N);
    //calculateU2(E, U);
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
    //FILE *file = fopen("file1.txt", "w");
    IPrinter::printVector(U, NULL, 10, 0, 0);
    //fclose(file);
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
    printf("eo: %12.8f %12.8f %12.8f\n", E[0], E[1], E[2]);
}

double Hyperbolic1DX::fx(const DoubleVector &e)
{
    pe = &e;
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
    return sum;
}

void Hyperbolic1DX::gradient(const DoubleVector &e, DoubleVector &g)
{
    pe = &e;
    DoubleVector u;
    IHyperbolicEquation::calculateU(u, hx, ht, M, N);

    pu = &u;
    DoubleMatrix p;
    IBackwardHyperbolicEquation::calculateU(p, hx, ht, M, N);

    for (unsigned int j=0; j<=M; j++)
    {
        calculateG(e, p[j], g, j);
    }
}

void Hyperbolic1DX::project(DoubleVector &e, int i)
{
    if (e[i]>1.0) e[i]=1.0;
    if (e[i]<0.0) e[i]=0.0;
}

void Hyperbolic1DX::print(unsigned int i, const DoubleVector &e, const DoubleVector &g, double alpha, RnFunction *fn) const
{
    C_UNUSED(g);
    C_UNUSED(alpha);
    Hyperbolic1DX *hc = dynamic_cast<Hyperbolic1DX*>(fn);
    printf("J[%d]: %.16f\n", i, hc->fx(e));
    printf("e1: %12.8f %12.8f %12.8f\n", e[0], e[1], e[2]);
    printf("g1: %12.8f %12.8f %12.8f\n", g[0], g[1], g[2]);
}

double Hyperbolic1DX::initial1(unsigned int i) const
{
    double x = i*hx;
    return x*x + t0*t0;
}

double Hyperbolic1DX::initial2(unsigned int i) const
{
    C_UNUSED(i);
    return 0.0;
}

double Hyperbolic1DX::boundary(Boundary type, unsigned int j) const
{
    double t = j*ht;
    if (type==Left) return x0*x0 + t*t;
    if (type==Right)   return x1*x1 + t*t;
    return 0.0;
}

//double Hyperbolic1DX::m1(unsigned int j) const
//{
//    double t = j*ht;
//    return x0*x0 + t*t;
//}

//double Hyperbolic1DX::m2(unsigned int j) const
//{
//    double t = j*ht;
//    return x1*x1 + t*t;
//}

double Hyperbolic1DX::f(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    double t = j*ht;
    double sum = 0.0;
    const DoubleVector &e = *pe;

    //    if (fabs(x-e[0])<=hx)
    //    {
    //        sum += f1(t) * ((hx-fabs(x-e[0]))/(hx*hx));
    //    }
    //    if (fabs(x-e[1])<=hx)
    //    {
    //        sum += f2(t) * ((hx-fabs(x-e[1]))/(hx*hx));
    //    }
    //    if (fabs(x-e[2])<=hx)
    //    {
    //        sum += f3(t) * ((hx-fabs(x-e[2]))/(hx*hx));
    //    }

    double sgm = 3.0*hx;
    double a1 = 1.0 / (sgm*sqrt(2.0*M_PI));
    double b1 = 2.0*sgm*sgm;
    sum += v1(t) * a1*exp(-((x-e[0])*(x-e[0]))/b1);
    sum += v2(t) * a1*exp(-((x-e[1])*(x-e[1]))/b1);
    sum += v3(t) * a1*exp(-((x-e[2])*(x-e[2]))/b1);

    return sum;
}

double Hyperbolic1DX::binitial1(unsigned int i) const
{
    C_UNUSED(i);
    return 0.0;
}

double Hyperbolic1DX::binitial2(unsigned int i) const
{
    const DoubleVector &u = *pu;
    return 2.0*(u[i] - U[i]);
}

double Hyperbolic1DX::bboundary(Boundary type, unsigned int j) const
{
    C_UNUSED(j);
    C_UNUSED(type);
    return 0.0;
}

double Hyperbolic1DX::bf(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

void Hyperbolic1DX::calculateG(const DoubleVector &e, const DoubleVector& psi, DoubleVector& g, unsigned int j)
{
    double psiX;
    if (j==0 || j==M)
    {
        psiDerivative(psiX, e[0], psi);
        g[0] = g[0] + v1(j*ht) * psiX;
        psiDerivative(psiX, e[1], psi);
        g[1] = g[1] + v2(j*ht) * psiX;
        psiDerivative(psiX, e[2], psi);
        g[2] = g[2] + v3(j*ht) * psiX;
    }
    else
    {
        psiDerivative(psiX, e[0], psi);
        g[0] = g[0] + 2.0*v1(j*ht) * psiX;
        psiDerivative(psiX, e[1], psi);
        g[1] = g[1] + 2.0*v2(j*ht) * psiX;
        psiDerivative(psiX, e[2], psi);
        g[2] = g[2] + 2.0*v3(j*ht) * psiX;
    }

    if (j==M)
    {
        g[0] = -(ht/2.0)*g[0];
        g[1] = -(ht/2.0)*g[1];
        g[2] = -(ht/2.0)*g[2];
    }
}

void Hyperbolic1DX::psiDerivative(double &psiX, double e, const DoubleVector& psi)
{
    unsigned int i = (unsigned int)round(e/hx);

    if (i==0)
        psiX  = (psi[i+1] - psi[i])/hx;
    else if (i==N)
        psiX = (psi[i] - psi[i-1])/hx;
    else
        psiX = (psi[i+1] - psi[i-1])/(2.0*hx);
}
