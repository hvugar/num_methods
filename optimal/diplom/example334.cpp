#include "example334.h"

// Optimal points [0.50,0.80] [0.70,0.20] [0.20,0.30]
// Working points [0.60,0.70] [0.65,0.25] [0.25,0.35] epsilon1 0.0001 epsilon2 0.0001 epsilon3: 0.0001. min:1.0 0.0001
// Working points [0.60,0.70] [0.60,0.30] [0.30,0.40] epsilon1 0.0001 epsilon2 0.0001 epsilon3: 0.0001. min:1.0 0.0001

void Parabolic1DControl334::main()
{
    Parabolic1DControl334 hc;

    DoubleVector x(2*hc.L);
    x[0] = 0.60; x[1] = 0.70; x[2] = 0.65; x[3] = 0.25; x[4] = 0.25; x[5] = 0.35;

    /* Minimization */
    ConjugateGradient g2;
    g2.setFunction(&hc);
    g2.setGradient(&hc);
    g2.setEpsilon1(0.0001);
    g2.setEpsilon2(0.0001);
    g2.setEpsilon3(0.0001);
    g2.setR1MinimizeEpsilon(1.0, 0.0001);
    g2.setPrinter(&hc);
    g2.setProjection(&hc);
    g2.setNormalize(true);
    g2.calculate(x);

    DoubleVector gr1(x.size());
    hc.gradient(x, gr1);
    gr1.L2Normalize();

    DoubleVector gr2(x.size());
    IGradient::Gradient(&hc, 0.00001, x, gr2);
    gr2.L2Normalize();

    printf("J[%d]: %.16f\n", 0, hc.fx(x));
    printf("eo: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", 0.50, 0.80, 0.70, 0.20, 0.20, 0.30);
    printf("e1: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", x[0], x[1], x[2], x[3], x[4], x[5]);
    printf("gr1: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", gr1[0], gr1[1], gr1[2], gr1[3], gr1[4], gr1[5]);
    printf("gr2: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", gr2[0], gr2[1], gr2[2], gr2[3], gr2[4], gr2[5]);
}

Parabolic1DControl334::Parabolic1DControl334()
{
    alpha = 1.0;

    x10 = 0.0;
    x11 = 1.0;
    x20 = 0.0;
    x21 = 1.0;
    t0 = 0.0;
    t1 = 1.0;

    h1 = 0.01;
    h2 = 0.01;
    ht = 0.01;

    N1 = (unsigned int)(ceil(x11-x10)/h1);
    N2 = (unsigned int)(ceil(x21-x20)/h2);
    M  = (unsigned int)(ceil(t1-t0)/ht);

    this->M  = M;
    this->N2 = N2;
    this->N1 = N1;
    this->L  = 3;

    a1 = a2 = 1.0;

    double sgm1 = 3.0*h1;
    double sgm2 = 3.0*h2;
    gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    gause_b = 2.0*sgm1*sgm2;

    // initialize
    DoubleVector e(2*L);
    e[0] = 0.50;
    e[1] = 0.80;
    e[2] = 0.70;
    e[3] = 0.20;
    e[4] = 0.20;
    e[5] = 0.30;

    px = &e;
    IParabolicEquation2D::caluclateMVD(U, h1, h2, ht, N1, N2, M, a1, a2);

    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
    IPrinter::printMatrix(U, 10, 10);
    printf("eo: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", e[0], e[1], e[2], e[3], e[4], e[5]);
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");

    FILE* f = fopen("heat_optimal_e.txt", "w");
    IPrinter::printMatrix(U, N1, N2, NULL, f);
    fclose(f);
}

double Parabolic1DControl334::fx(const DoubleVector& x)
{
    px = &x;
    DoubleMatrix u;
    IParabolicEquation2D::caluclateMVD(u, h1, h2, ht, N1, N2, M, a1, a2);

    double sum = 0.0;
    for (unsigned int j=0; j<=N2; j++)
    {
        for (unsigned int i=0; i<=N1; i++)
        {
            double betta = 1.0;
            if (i==0 || i==N1 || j==0 || j==N2) betta = 0.5;
            if (i==0  && j==0)  betta = 0.25;
            if (i==0  && j==N2) betta = 0.25;
            if (i==N1 && j==0)  betta = 0.25;
            if (i==N1 && j==N2) betta = 0.25;
            sum += betta*(u[j][i] - U[j][i])*(u[j][i] - U[j][i]);
        }
    }
    sum = (h1*h2)*sum;

    double nrm = 0.0;
    //nrm = norm(x);
    return sum + alpha*nrm;
}

double Parabolic1DControl334::norm(const DoubleVector& v) const
{
    double nrm = 0.0;
    for (unsigned int k=0; k<=M; k++)
    {
        double betta = 1.0;
        if (k==0 || k==M) betta = 0.5;
        nrm += betta*(v[2*L+0*(M+1)+k] - v1(k*ht))*(v[2*L+0*(M+1)+k] - v1(k*ht));
        nrm += betta*(v[2*L+1*(M+1)+k] - v2(k*ht))*(v[2*L+1*(M+1)+k] - v2(k*ht));
        nrm += betta*(v[2*L+2*(M+1)+k] - v3(k*ht))*(v[2*L+2*(M+1)+k] - v3(k*ht));
    }
    nrm = ht * nrm;
    return nrm;
}

void Parabolic1DControl334::gradient(const DoubleVector& x, DoubleVector& g)
{
    px = &x;
    DoubleMatrix u;
    IParabolicEquation2D::caluclateMVD(u, h1, h2, ht, N1, N2, M, a1, a2);

    pu = &u;
    DoubleCube psi;
    IBackwardParabolicEquation2D::caluclateMVD(psi, h1, h2, ht, N1, N2, M, a1, a2);

    for (unsigned int i=0; i<g.size(); i++) g[i] = 0.0;

    for (unsigned int k=M; k>=1; k--)
    {
        calculateGX(x, psi[k], g, k);
    }

    psi.clear();
}

void Parabolic1DControl334::calculateGX(const DoubleVector& x, const DoubleMatrix& psi, DoubleVector& g, unsigned int k)
{
    double psiX1;
    double psiX2;
    if (k==1 || k==M)
    {
        psiDerivative(psiX1, psiX2, x[0], x[1], psi);
        g[0] = g[0] + v1(k*ht) * psiX1;
        g[1] = g[1] + v1(k*ht) * psiX2;
        psiDerivative(psiX1, psiX2, x[2], x[3], psi);
        g[2] = g[2] + v2(k*ht) * psiX1;
        g[3] = g[3] + v2(k*ht) * psiX2;
        psiDerivative(psiX1, psiX2, x[4], x[5], psi);
        g[4] = g[4] + v3(k*ht) * psiX1;
        g[5] = g[5] + v3(k*ht) * psiX2;
    }
    else
    {
        psiDerivative(psiX1, psiX2, x[0], x[1], psi);
        g[0] = g[0] + 2.0*v1(k*ht) * psiX1;
        g[1] = g[1] + 2.0*v1(k*ht) * psiX2;
        psiDerivative(psiX1, psiX2, x[2], x[3], psi);
        g[2] = g[2] + 2.0*v2(k*ht) * psiX1;
        g[3] = g[3] + 2.0*v2(k*ht) * psiX2;
        psiDerivative(psiX1, psiX2, x[4], x[5], psi);
        g[4] = g[4] + 2.0*v3(k*ht) * psiX1;
        g[5] = g[5] + 2.0*v3(k*ht) * psiX2;
    }

    if (k==1)
    {
        g[0] = -(ht/2.0)*g[0];
        g[1] = -(ht/2.0)*g[1];
        g[2] = -(ht/2.0)*g[2];
        g[3] = -(ht/2.0)*g[3];
        g[4] = -(ht/2.0)*g[4];
        g[5] = -(ht/2.0)*g[5];
    }
}

void Parabolic1DControl334::psiDerivative(double &psiX1, double &psiX2, double x1, double x2, const DoubleMatrix &psi)
{
    unsigned int i = (unsigned int)round(x1/h1);
    unsigned int j = (unsigned int)round(x2/h2);

    if (i==0) psiX1  = (psi[j][i+1] - psi[j][i])/h1;
    else if (i==N1) psiX1 = (psi[j][i] - psi[j][i-1])/h1;
    else psiX1 = (psi[j][i+1] - psi[j][i-1])/(2.0*h1);

    if (j==0) psiX2 = (psi[j+1][i] - psi[j][i])/h2;
    else if (j==N2) psiX2 = (psi[j][i] - psi[j-1][i])/h2;
    else psiX2 = (psi[j+1][i] - psi[j-1][i])/(2.0*h2);
}

double Parabolic1DControl334::fi(unsigned int i, unsigned int j) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    return u(x1, x2, t0);
}

double Parabolic1DControl334::m1(unsigned int j, unsigned int k) const
{
    double x2 = j*h2;
    double t  = 0.5*k*ht;
    return u(x10, x2, t);
}

double Parabolic1DControl334::m2(unsigned int j, unsigned int k) const
{
    double x2 = j*h2;
    double t  = 0.5*(k*ht);
    return u(x11, x2, t);
}

double Parabolic1DControl334::m3(unsigned int i, unsigned int k) const
{
    double x1 = i*h1;
    double t  = 0.5*k*ht;
    return u(x1, x20, t);
}

double Parabolic1DControl334::m4(unsigned int i, unsigned int k) const
{
    double x1 = i*h1;
    double t  = 0.5*k*ht;
    return u(x1, x21, t);
}

double Parabolic1DControl334::f(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t  = 0.5*k*ht;
    double sum = 0.0;

    //    double sgm1 = 3.0*h1;
    //    double sgm2 = 3.0*h2;
    //    double gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    //    double gause_b = 2.0*sgm1*sgm2;

    const DoubleVector &e = *px;
    sum += v1(t) * gause_a * exp(-((x1-e[0])*(x1-e[0]) + (x2-e[1])*(x2-e[1]))/gause_b);
    sum += v2(t) * gause_a * exp(-((x1-e[2])*(x1-e[2]) + (x2-e[3])*(x2-e[3]))/gause_b);
    sum += v3(t) * gause_a * exp(-((x1-e[4])*(x1-e[4]) + (x2-e[5])*(x2-e[5]))/gause_b);
    return sum;

    //    if (fabs(x1-e[0])<=h1 && fabs(x2-e[1])<=h2)
    //    {
    //        sum += v1(t) * ((h1-fabs(x1-e[0]))/(h1*h1))*((h2-fabs(x2-e[1]))/(h2*h2));
    //    }
    //    if (fabs(x1-e[2])<=h1 && fabs(x2-e[3])<=h2)
    //    {
    //        sum += v1(t) * ((h1-fabs(x1-e[2]))/(h1*h1))*((h2-fabs(x2-e[3]))/(h2*h2));
    //    }
    //    if (fabs(x1-e[4])<=h1 && fabs(x2-e[5])<=h2)
    //    {
    //        sum += v1(t) * ((h1-fabs(x1-e[4]))/(h1*h1))*((h2-fabs(x2-e[5]))/(h2*h2));
    //    }
}

double Parabolic1DControl334::bfi(unsigned int i, unsigned int j) const
{
    return -2.0*((*pu)[j][i] - U[j][i]);
}

double Parabolic1DControl334::bm1(unsigned int j, unsigned int k) const
{
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double Parabolic1DControl334::bm2(unsigned int j, unsigned int k) const
{
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double Parabolic1DControl334::bm3(unsigned int i, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(k);
    return 0.0;
}

double Parabolic1DControl334::bm4(unsigned int i, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(k);
    return 0.0;
}

double Parabolic1DControl334::bf(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

void Parabolic1DControl334::print(unsigned int i, const DoubleVector& x, const DoubleVector &g, double alpha, RnFunction* fn) const
{
    C_UNUSED(alpha);
    Parabolic1DControl334 *hc = dynamic_cast<Parabolic1DControl334*>(fn);
    printf("J[%d]: %.16f\n", i, hc->fx(x));
    DoubleVector g1 = g;
    g1.L2Normalize();
    printf("eo: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", 0.50, 0.80, 0.70, 0.20, 0.20, 0.30);
    printf("e1: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", x[0], x[1], x[2], x[3], x[4], x[5]);
    printf("g1: [%12.8f, %12.8f] [%12.8f, %12.8f] [%12.8f, %12.8f]\n", g1[0], g1[1], g1[2], g1[3], g1[4], g1[5]);
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
}

void Parabolic1DControl334::project(DoubleVector &e, int index)
{
    if (index<6)
    {
        if (e[index] > 1.0) e[index] = 1.0;
        if (e[index] < 0.0) e[index] = 0.0;
    }
}
