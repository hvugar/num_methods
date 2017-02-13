#include "hyperboliccontrol2dmx.h"

void HyperbolicControl2DMX::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    HyperbolicControl2DMX hc;
    hc.fx(1.0);
}

HyperbolicControl2DMX::HyperbolicControl2DMX()
{
}

HyperbolicControl2DMX::~HyperbolicControl2DMX()
{
}

double HyperbolicControl2DMX::fx(double t)
{
    h1 = 0.01;
    h2 = 0.01;
    ht = 0.005;

    N1 = 100;
    N2 = 100;
    M  = 200;
    L = 2;

    e.resize(2);
    e[0] = 0.2;
    e[1] = 0.2;

    alpha0 = 1.0;//4.0;
    alpha1 = 10.0;//5.0;
    alpha2 = 2.0;//6.0;
    alpha3 = 1.0;//4.0;
    qamma = 0.2;

    U0 = 0.0;
    U1 = 0.0;
    a = 1.0;

    DoubleVector x(2*L);
    x[0] = 0.3;
    x[1] = 0.4;
    x[2] = 0.7;
    x[3] = 0.7;

    //    double min_step = 1.0;
    //    double gold_eps = 0.001;

    //    ConjugateGradient cg;
    //    cg.setFunction(this);
    //    cg.setGradient(this);
    //    cg.setEpsilon1(0.001);
    //    cg.setEpsilon2(0.001);
    //    cg.setEpsilon3(0.001);
    //    cg.setR1MinimizeEpsilon(min_step, gold_eps);
    //    cg.setPrinter(this);
    //    cg.setProjection(this);
    //    cg.setNormalize(false);
    //    cg.showEndMessage(false);
    //    cg.calculate(x);

    double rf = fx(x);

    double h = 0.001;
    DoubleVector g1(x.size());
    gradient(x, g1);
    g1.L2Normalize();
    DoubleVector g2(x.size());
    IGradient::Gradient(this, h, x, g2);
    g2.L2Normalize();

    FILE *file = fopen("gradient_x.txt", "a");
    fprintf(file, "--------------------------------------------------------------------\n");
    IPrinter::printDateTime(file);
    fprintf(file, "T: %f L: %d h:%f Functional: %.20f N1: %d N2: %d M: %d h1: %f h2: %f ht: %f\n", t, L, h, rf, N1, N2, M, h1, h2, ht);
    fprintf(file, "x: %.8f %.8f %.8f %.8f\n", x[0], x[1], x[2], x[3]);
    fprintf(file, "Agx: %.8f %.8f %.8f %.8f\n", g1[0], g1[1], g1[2], g1[3]);
    fprintf(file, "Ngx: %.8f %.8f %.8f %.8f\n", g2[0], g2[1], g2[2], g2[3]);
    IPrinter::printDateTime(file);
    //    fprintf(file, "U\n");
    //    DoubleCube c;
    //    IHyperbolicEquation2D::calculateU1(c, h1, h2, ht, N1, N2, M, a, a, qamma);
    //    IPrinter::printMatrix(c[c.size()-1], N2, N1, NULL, file);
    fclose(file);
    x.clear();

    return rf;
}

double HyperbolicControl2DMX::fx(const DoubleVector &x)
{
    px = &x;
    DoubleCube c;
    IHyperbolicEquation2D::calculateU1(c, h1, h2, ht, N1, N2, M, a, a, qamma);

    const DoubleMatrix &u0 = c.matrix(M);
    const DoubleMatrix &u1 = c.matrix(M-2);

    double sum = 0.0;

    double sum1 = 0.0;
    for (unsigned int j=0; j<=N2; j++)
    {
        for (unsigned int i=0; i<=N1; i++)
        {
            double k = 1.0;
            if (i==0 || i==N1) k *= 0.5;
            if (j==0 || j==N2) k *= 0.5;
            sum1 = sum1 + k * (u0[j][i]-U0) * (u0[j][i]-U0);
        }
    }
    sum1 = h1*h2*sum1;

    double sum2 = 0.0;
    for (unsigned int j=0; j<=N2; j++)
    {
        for (unsigned int i=0; i<=N1; i++)
        {
            double k = 1.0;
            if (i==0 || i==N1) k *= 0.5;
            if (j==0 || j==N2) k *= 0.5;
            sum2 = sum2 + k * ( (u0[j][i]-u1[j][i])/(2.0*ht)-U1)*((u0[j][i]-u1[j][i])/(2.0*ht)-U1);
        }
    }
    sum2 = h1*h2*sum2;

    sum = sum1 + alpha0*sum2;
    return sum;// + norm(v);
}

double HyperbolicControl2DMX::norm(const DoubleVector& v) const
{
    C_UNUSED(v);
    double nrm = 0.0;
    for (unsigned int k=0; k<=M; k++)
    {
        double betta = 1.0;
        if (k==0 || k==M) betta = 0.5;
        double m1 = 0.0;//(v[0*(M+1)+k] - v1(k*ht));
        nrm += betta * m1*m1;
    }
    nrm = ht * nrm;
    return nrm;
}

void HyperbolicControl2DMX::gradient(const DoubleVector &x, DoubleVector &g)
{
    px = &x;
    DoubleCube u;
    IHyperbolicEquation2D::calculateU1(u, h1, h2, ht, N1, N2, M, a, a, qamma);

    pu = &u;
    DoubleCube p;
    IBackwardHyperbolicEquation2D::calculateU1(p, h1, h2, ht, N1, N2, M, a, a, qamma);

    for (unsigned int k=0; k<=M; k++)
    {
        double psiX1;
        double psiX2;
        double m = 1.0;
        if (k==0 || k==M) m *= 0.5;

        psiDerivative(psiX1, psiX2, x[0], x[1], p.matrix(k));
        g[0] = g[0] + m * v1(k*ht) * psiX1;
        g[1] = g[1] + m * v1(k*ht) * psiX2;

        psiDerivative(psiX1, psiX2, x[2], x[3], p.matrix(k));
        g[2] = g[2] + m * v2(k*ht) * psiX1;
        g[3] = g[3] + m * v2(k*ht) * psiX2;
    }
    g[0] = -ht*g[0];
    g[1] = -ht*g[1];
    g[2] = -ht*g[2];
    g[3] = -ht*g[3];

    printf("gx: %.8f %.8f %.8f %.8f\n", g[0], g[1], g[2], g[3]);
}

void HyperbolicControl2DMX::psiDerivative(double &psiX1, double &psiX2, double x1, double x2, const DoubleMatrix &psi)
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

double HyperbolicControl2DMX::initial1(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    return 0.0;
}

double HyperbolicControl2DMX::initial2(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    return 0.0;
}

double HyperbolicControl2DMX::boundary(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2DMX::f(unsigned int i, unsigned int j, unsigned int k) const
{
    double sum = 0.0;
    static double sgm1 = 3.0*h1;
    static double sgm2 = 3.0*h2;
    static double gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    static double gause_b = 2.0*sgm1*sgm2;

    double x1 = i*h1;
    double x2 = j*h2;
    double t = k*ht;
    const DoubleVector &x = *px;

    double _v1 = v1(t);
    double _v2 = v2(t);

    sum += _v1 * gause_a * exp(-((x1-x[0])*(x1-x[0]) + (x2-x[1])*(x2-x[1]))/gause_b);
    sum += _v2 * gause_a * exp(-((x1-x[2])*(x1-x[2]) + (x2-x[3])*(x2-x[3]))/gause_b);

    sum += fxt(i, j, k);

    return sum;
}

double HyperbolicControl2DMX::binitial1(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu).matrix(M);
    const DoubleMatrix &u1 = (*pu).matrix(M-2);
    return -2.0 * alpha1 * ( (u0[j][i]-u1[j][i])/(2.0*ht) - U1);
}

double HyperbolicControl2DMX::binitial2(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu).matrix(M);
    return +2.0 * (u0[j][i] - U0) + qamma*(binitial1(i,j));
}

double HyperbolicControl2DMX::bboundary(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2DMX::bf(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2DMX::fxt(unsigned int i, unsigned int j, unsigned int k) const
{
    double sum = 0.0;
    double x1 = i*h1;
    double x2 = j*h2;
    double t = k*ht;
    sum += alpha1*exp(-alpha2*((x1-e[0])*(x1-e[0])+(x2-e[1])*(x2-e[1]))-alpha3*t);
    return sum;
}

void HyperbolicControl2DMX::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double alpha, RnFunction *fn) const
{
    C_UNUSED(g);
    C_UNUSED(alpha);
    printf("J[%d]: %.16f\n", i, fn->fx(x));
}

void HyperbolicControl2DMX::project(DoubleVector &x, int i)
{
    if (i==0)
    {
        if (x[i] <= 0.0) x[i] = 0.0 + h1;
        if (x[i] >= 0.5) x[i] = 0.5 - h1;
    }
    if (i==1)
    {
        if (x[i] <= 0.0) x[i] = 0.0 + h1;
        if (x[i] >= 0.5) x[i] = 0.5 - h1;
    }
    if (i==2)
    {
        if (x[i] <= 0.5) x[i] = 0.5 + h1;
        if (x[i] >= 1.0) x[i] = 1.0 - h1;
    }
    if (i==3)
    {
        if (x[i] <= 0.5) x[i] = 0.5 + h1;
        if (x[i] >= 1.0) x[i] = 1.0 - h1;
    }
}

double HyperbolicControl2DMX::v1(double t) const
{
    C_UNUSED(t);
    return 2.0;
}

double HyperbolicControl2DMX::v2(double t) const
{
    C_UNUSED(t);
    return 2.0;
}
