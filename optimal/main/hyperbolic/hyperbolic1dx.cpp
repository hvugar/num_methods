#include "hyperbolic1dx.h"
#include <tomasmethod.h>
#include <gradient_cjt.h>
#include <stdlib.h>
#include <math.h>

Hyperbolic1DX::Hyperbolic1DX(unsigned int M, unsigned int N)
{
    t0 = 0.0;
    t1 = 1.0;
    x0 = 0.0;
    x1 = 1.0;
    this->M = M;
    this->N = N;
    L = 3;
    hx = (x1-x0)/N;
    ht = (t1-t0)/M;
    a = 1.0;
    lamda = 0.25;

    initialize();
}

void Hyperbolic1DX::initialize()
{
    DoubleVector E(L);
    E[0] = 0.20; E[1] = 0.40; E[2] = 0.70;
    calculateU(E, U);
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
    IPrinter::printVector(U);
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
    printf("eo: %12.8f %12.8f %12.8f\n", E[0], E[1], E[2]);

}

double Hyperbolic1DX::fx(const DoubleVector &e)
{
    calculateU(e, uT);
    return calculateIntegral(e); //calculateNorm(e);
}

double Hyperbolic1DX::calculateIntegral(const DoubleVector& e)
{
    double sum = 0.0;
    for (unsigned int i=0; i<N; i++)
    {
        int j = i+1;

        double f1 = uT[j] - U[j];
        double f2 = uT[i] - U[i];

        sum += (f1*f1+f2*f2);
    }
    return 0.5*hx*sum;
}

double Hyperbolic1DX::calculateNorm(const DoubleVector& e)
{
    double norm = 0.0;
//    for (unsigned int j=0; j<M; j++)
//    {
//        for (unsigned int i=0; i<N; i++)
//        {
//            unsigned int i1 = i;
//            unsigned int i2 = i+1;
//            unsigned int j1 = j;
//            unsigned int j2 = j+1;

//            double f11 = f[j1*(N+1) + i1] - fxt(i1*hx, j1*ht);
//            double f12 = f[j1*(N+1) + i2] - fxt(i2*hx, j1*ht);
//            double f21 = f[j2*(N+1) + i1] - fxt(i1*hx, j2*ht);
//            double f22 = f[j2*(N+1) + i2] - fxt(i2*hx, j2*ht);

//            norm += f11*f11 + f12*f12 + f21*f21 + f22*f22;
//        }
//    }
    return (hx*ht)*0.25*norm;
}

void Hyperbolic1DX::gradient(const DoubleVector &e, DoubleVector &g)
{
    DoubleVector u;
    calculateU(e, uT);
    calculateP(e, g);
//    printf("e1: %12.8f %12.8f %12.8f\n", e[0], e[1], e[2]);
//    printf("g1: %12.8f %12.8f %12.8f\n", g[0], g[1], g[2]);
//    calculateG1(e, g);
}

void Hyperbolic1DX::project(DoubleVector &e, int k)
{
    for (unsigned int i=0; i<e.size(); i++)
    {
        if (e[i]>1.0) e[i]=1.0;
        if (e[i]<0.0) e[i]=0.0;
    }
}

void Hyperbolic1DX::print(unsigned int i, const DoubleVector &e, const DoubleVector &g, double a, RnFunction *fn) const
{
    Hyperbolic1DX *hc = dynamic_cast<Hyperbolic1DX*>(fn);
    printf("J[%d]: %.16f\n", i, hc->fx(e));
    printf("e1: %12.8f %12.8f %12.8f\n", e[0], e[1], e[2]);
    printf("g1: %12.8f %12.8f %12.8f\n", g[0], g[1], g[2]);
}

void Hyperbolic1DX::calculateU(const DoubleVector &e, DoubleVector &u)
{
    u.clear();
    u.resize(N+1);

    DoubleVector u0(N+1);
    DoubleVector u1(N+1);

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

    double alpha1 = -(lamda*a*a)*((ht*ht)/(hx*hx));
    double beta1  = 1.0 + (2.0*lamda*a*a*(ht*ht))/(hx*hx);
    double alpha2 = (1.0-2.0*lamda)*a*a*((ht*ht)/(hx*hx));
    double alpha3 = +(lamda*a*a)*((ht*ht)/(hx*hx));

    for (unsigned int j=0; j<=M-1; j++)
    {
        if (j==0)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                u0[i] = fi1(i*hx);
                u1[i] = fi1(i*hx) + ht*fi2(i*hx);
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                a1[i-1] = alpha1;
                b1[i-1] = beta1;
                c1[i-1] = alpha1;
                d1[i-1] = alpha2*(u1[i-1]-2.0*u1[i]+u1[i+1]) + 2.0*u1[i] - u0[i] + alpha3*(u0[i+1] - 2.0*u0[i] + u0[i-1]) + (ht*ht)*fxt(i, j, e);
            }

            a1[0]   = 0.0;
            c1[N-2] = 0.0;
            d1[0]   -= alpha1 * mu1((j+1)*ht);
            d1[N-2] -= alpha1 * mu2((j+1)*ht);
            TomasAlgorithm(a1, b1, c1, d1, x1);

            u[0] = mu1((j+1)*ht);
            for (unsigned int i=1; i<=N-1; i++)
            {
                u[i] = x1[i-1];
            }
            u[N] = mu2((j+1)*ht);

            for (unsigned int i=0; i<=N; i++)
            {
                u0[i] = u1[i];
                u1[i] = u[i];
            }
        }
    }

    a1.clear();
    b1.clear();
    c1.clear();
    d1.clear();
    x1.clear();

    //Printer::printVector(u, N/10);
}

void Hyperbolic1DX::calculateP(const DoubleVector &e, DoubleVector &g)
{
    DoubleVector p(N+1);
    DoubleVector p0(N+1);
    DoubleVector p1(N+1);

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

    double alpha1 = -(lamda*a*a)*((ht*ht)/(hx*hx));
    double beta1  = 1.0 + (2.0*lamda*a*a*(ht*ht))/(hx*hx);
    double alpha2 = (1.0-2.0*lamda)*a*a*((ht*ht)/(hx*hx));
    double alpha3 = +(lamda*a*a)*((ht*ht)/(hx*hx));

    for (unsigned int j1=0; j1<=M-1; j1++)
    {
        unsigned int j = M - j1;

        if (j==M)
        {
            for (unsigned int i=0; i<=N; i++)
            {
                p0[i] = pfi1(i*hx);
                p1[i] = pfi1(i*hx) - ht*pfi2(i);
            }
            calculateG(e, p0, g, M);
            calculateG(e, p1, g, M-1);
            //Printer::printVector(p0, N/10);
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                a1[i-1] = alpha1;
                b1[i-1] = beta1;
                c1[i-1] = alpha1;
                d1[i-1] = alpha2*(p1[i-1]-2.0*p1[i]+p1[i+1]) + 2.0*p1[i] - p0[i] + alpha3*(p0[i+1] - 2.0*p0[i] + p0[i-1]);
            }

            a1[0]   = 0.0;
            c1[N-2] = 0.0;
            d1[0]   -= alpha1 * pmu1((j-1)*ht);
            d1[N-2] -= alpha1 * pmu2((j-1)*ht);
            TomasAlgorithm(a1, b1, c1, d1, x1);

            p[0] = pmu1((j-1)*ht);
            for (unsigned int i=1; i<=N-1; i++)
            {
                p[i] = x1[i-1];
            }
            p[N] = pmu2((j-1)*ht);

            for (unsigned int i=0; i<=N; i++)
            {
                p0[i] = p1[i];
                p1[i] = p[i];
            }
            calculateG(e, p, g, j-1);
        }
    }

    a1.clear();
    b1.clear();
    c1.clear();
    d1.clear();
    x1.clear();

    p0.clear();
    p1.clear();

    //Printer::printVector(p, N/10);
}

void Hyperbolic1DX::calculateG(const DoubleVector& e, const DoubleVector& psi, DoubleVector& g, unsigned int j)
{
    double psiX;
    if (j==0 || j==M)
    {
        psiDerivative(psiX, e[0], psi);
        g[0] = g[0] + f1(j*ht) * psiX;
        psiDerivative(psiX, e[1], psi);
        g[1] = g[1] + f2(j*ht) * psiX;
        psiDerivative(psiX, e[2], psi);
        g[2] = g[2] + f3(j*ht) * psiX;
    }
    else
    {
        psiDerivative(psiX, e[0], psi);
        g[0] = g[0] + 2.0*f1(j*ht) * psiX;
        psiDerivative(psiX, e[1], psi);
        g[1] = g[1] + 2.0*f2(j*ht) * psiX;
        psiDerivative(psiX, e[2], psi);
        g[2] = g[2] + 2.0*f3(j*ht) * psiX;
    }

    if (j==0)
    {
        g[0] = -(ht/2.0)*g[0];
        g[1] = -(ht/2.0)*g[1];
        g[2] = -(ht/2.0)*g[2];
    }
}

void Hyperbolic1DX::calculateG1(const DoubleVector &e, DoubleVector &g)
{
    double h = 0.000001;
    DoubleVector E(L);
    double f0 = fx(e);

    for (unsigned int i=0; i<e.size(); i++)
    {
        E = e;
        E[i] += h;
        double f1 = fx(E);
        g[i] = (f1-f0)/h;
    }

    printf("e2: %12.8f, %12.8f %12.8f\n", e[0], e[1], e[2]);
    printf("g2: %12.8f, %12.8f %12.8f\n", g[0], g[1], g[2]);
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

double Hyperbolic1DX::u(double x, double t)
{
    return x*x + t*t;
}

double Hyperbolic1DX::fi1(double x)
{
    return x*x;
}

double Hyperbolic1DX::fi2(double x)
{
    return 0.0;
}

double Hyperbolic1DX::mu1(double t)
{
    return t*t;
}

double Hyperbolic1DX::mu2(double t)
{
    return t*t+1.0;
}

double Hyperbolic1DX::f(double x, double t)
{
    return 2.0-2.0*a*a;
}

double Hyperbolic1DX::fxt(unsigned int i, unsigned int j, const DoubleVector &e)
{
    double x = i*hx;
    double t = j*ht;

//    return 2.0-2.0*a*a;

    double sum = 0.0;

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

    double sgm = 10.0*hx;
    double a1 = 1.0 / (sgm*sqrt(2.0*M_PI));

    sum += f1(t) * a1*exp(-((x-e[0])*(x-e[0]))/(2.0*sgm*sgm));
    sum += f2(t) * a1*exp(-((x-e[1])*(x-e[1]))/(2.0*sgm*sgm));
    sum += f3(t) * a1*exp(-((x-e[2])*(x-e[2]))/(2.0*sgm*sgm));

    return sum;
}

double Hyperbolic1DX::pfi1(double x) const { return 0.0; }
double Hyperbolic1DX::pfi2(unsigned int i) const { return 2.0*(uT[i] - U[i]); }
double Hyperbolic1DX::pmu1(double t) const { return 0.0; }
double Hyperbolic1DX::pmu2(double t) const { return 0.0; }

void Hyperbolic1DX::main()
{
    Hyperbolic1DX hc(100, 100);
//    hc.test(2);

    DoubleVector e(hc.L);
    e[0] = 0.6; e[1] = 0.9; e[2] = 0.3;
    //e[0] = 0.2; e[1] = 0.4; e[2] = 0.7;

    //DoubleVector u;
    //DoubleVector g(hc.L);
    //hc.calculateU(e, u);
    //hc.calculateP(e, g);
    //printf("%f %f %f\n", g[0], g[1], g[2]);
    //return;

    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setGradient(&hc);
    g2.setEpsilon1(0.000000001);
    g2.setEpsilon2(0.000000001);
    g2.setR1MinimizeEpsilon(1.0, 0.00001);
    g2.setPrinter(&hc);
    g2.setProjection(&hc);
    g2.setNormalize(true);
    g2.calculate(e);
}

void Hyperbolic1DX::test(int j)
{
    DoubleVector e(L);
    //Optimal
    e[0] = 0.2; e[1] = 0.4; e[2] = 0.7;

    for (unsigned int i=0; i<=N; i++)
    {
        e[j] = i*hx;
        double J = fx(e);
        printf("%.12f\n", J);
    }
    exit(-1);
}


