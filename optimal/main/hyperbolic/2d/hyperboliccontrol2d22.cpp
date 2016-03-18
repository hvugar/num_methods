#include "hyperboliccontrol2d22.h"

void HyperbolicControl2D22::main()
{
    HyperbolicControl2D22 hc;
    hc.file = fopen("hyperboliccontrol2d22.txt", "w");
    //hc.file = stdout;
    for (double t=0.8; t<=10.1; t+=0.1)
    {
        double rf = hc.fx(t);
        fprintf(hc.file, "%f %.16f\n", t, rf);
        fflush(hc.file);
    }
    fclose(hc.file);
}

HyperbolicControl2D22::HyperbolicControl2D22()
{
    x10 = 0.0;
    x11 = 1.0;
    x20 = 0.0;
    x21 = 1.0;
    t0 = 0.0;
    t1 = 1.0;

    h1 = 0.01;
    h2 = 0.01;
    ht = 0.005;
    N1 = (unsigned)ceil((x11 - x10)/h1);
    N2 = (unsigned)ceil((x21 - x20)/h2);
    M  = (unsigned)ceil((t1 - t0)/ht);

    L = 2;

    e.resize(2);
    e[0] = 0.2;
    e[1] = 0.2;

    alpha0 = 1.0;
    alpha1 = 5.0;
    alpha2 = 6.0;
    alpha3 = 4.0;
    qamma = 0.2;

    a1 = 1.0;
    a2 = 1.0;

    e.resize(2);
    e[0] = 0.2;
    e[1] = 0.2;

    alpha0 = 4.0;
    alpha1 = 5.0;
    alpha2 = 6.0;
    alpha3 = 4.0;
    qamma = 0.2;

    U0 = 0.0;
    U1 = 0.0;
    a1 = a2 = 1.0;

    L = 2;
    x.resize(2*L);
    x[0] = 0.3;
    x[1] = 0.4;
    x[2] = 0.7;
    x[3] = 0.7;

    pv = NULL;
    pu = NULL;
}

double HyperbolicControl2D22::fx(double T)
{
    t1 = T;
    M  = (unsigned)ceil((t1 - t0)/ht);

    DoubleVector v0((M+1)*L);
    for (unsigned int k=0; k<=M; k++)
    {
        v0[0*(M+1)+k] = 0.0;
        v0[1*(M+1)+k] = 0.0;
    }

    double min_step = 1.0;
    double gold_eps = 0.001;
    ConjugateGradient cg;
    cg.setFunction(this);
    cg.setGradient(this);
    cg.setEpsilon1(0.001);
    cg.setEpsilon2(0.001);
    cg.setEpsilon3(0.001);
    cg.setR1MinimizeEpsilon(min_step, gold_eps);
    cg.setPrinter(this);
    //cg.setProjection(this);
    cg.setNormalize(true);
    //cg.showEndMessage(false);
    cg.calculate(v0);

    IPrinter::printVector(v0);
    IPrinter::printVector(v0, NULL, v0.size(), 0, 0, file);
    double rf = fx(v0);
    return rf;
}

double HyperbolicControl2D22::fx(const DoubleVector &v)
{
    pv = &v;
    DoubleCube c;
    IHyperbolicEquation2D::calculateU1(c, h1, h2, ht, N1, N2, M, a1, a2, qamma);

    DoubleMatrix &u0 = c[M];
    DoubleMatrix &u1 = c[M-2];
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
            sum2 = sum2 + k * ((u0[j][i]-u1[j][i])/(2.0*ht) - U1) * ((u0[j][i]-u1[j][i])/(2.0*ht) - U1);
        }
    }
    sum2 = h1*h2*sum2;

    sum = sum1 + alpha0*sum2;
    return sum;
}

void HyperbolicControl2D22::gradient(const DoubleVector &v, DoubleVector &g)
{
    pv = &v;
    DoubleCube u;
    IHyperbolicEquation2D::calculateU1(u, h1, h2, ht, N1, N2, M, a1, a2, qamma);

    //IPrinter::printMatrix(u[u.size()-1]);

    pu = &u;
    DoubleCube p;
    IBackwardHyperbolicEquation2D::calculateU1(p, h1, h2, ht, N1, N2, M, a1, a2, qamma);

    //puts("---");
    //IPrinter::printMatrix(p[p.size()-1]);

    unsigned int i,j;
    for (unsigned int k=0; k<=M; k++)
    {
        i = (unsigned int)round(x[0]/h1);
        j = (unsigned int)round(x[1]/h2);
        g[0*(M+1)+k] = -p[k][j][i];

        i = (unsigned int)round(x[2]/h1);
        j = (unsigned int)round(x[3]/h2);
        g[1*(M+1)+k] = -p[k][j][i];
    }

    //IPrinter::printVector(v);
    //IPrinter::printVector(g);
}

void HyperbolicControl2D22::print(unsigned int i, const DoubleVector &v, const DoubleVector &g, double alpha, RnFunction* fn) const
{
    C_UNUSED(g);
    C_UNUSED(alpha);
    //IPrinter::printVector(v);
    fprintf(file, "J[%d]: %.16f\n", i, fn->fx(v));
    fflush(file);
    printf("J[%d]: %.16f\n", i, fn->fx(v));
}


double HyperbolicControl2D22::fi1(unsigned int i, unsigned int j) const
{
    return 0.0;
}

double HyperbolicControl2D22::fi2(unsigned int i, unsigned int j) const
{
    return 0.0;
}

double HyperbolicControl2D22::m1(unsigned int j, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2D22::m2(unsigned int j, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2D22::m3(unsigned int i, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2D22::m4(unsigned int i, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2D22::f(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(k);

    double sum = 0.0;

    static double sgm1 = 3.0*h1;
    static double sgm2 = 3.0*h2;
    static double gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    static double gause_b = 2.0*sgm1*sgm2;

    double x1 = i*h1;
    double x2 = j*h2;
    const DoubleVector v = *pv;

    double _v1 = v[0*(M+1)+k];
    double _v2 = v[1*(M+1)+k];

    sum += _v1 * gause_a * exp(-((x1-x[0])*(x1-x[0]) + (x2-x[1])*(x2-x[1]))/gause_b);
    sum += _v2 * gause_a * exp(-((x1-x[2])*(x1-x[2]) + (x2-x[3])*(x2-x[3]))/gause_b);

    sum += fxt(i, j, k);
    return sum;
}

double HyperbolicControl2D22::bfi1(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu)[M];
    const DoubleMatrix &u1 = (*pu)[M-2];
    return -2.0 * alpha1 * ((u0[j][i]-u1[j][i])/(2.0*ht) - U1);
}

double HyperbolicControl2D22::bfi2(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu)[M];
    return +2.0 * (u0[j][i] - U0) + qamma*(bfi1(i,j));
}

double HyperbolicControl2D22::bm1(unsigned int j, unsigned int k) const
{
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2D22::bm2(unsigned int j, unsigned int k) const
{
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2D22::bm3(unsigned int i, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2D22::bm4(unsigned int i, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2D22::bf(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j );
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2D22::fxt(unsigned int i, unsigned int j, unsigned int k) const
{
    double sum = 0.0;
    double x1 = i*h1;
    double x2 = j*h2;
    double t = k*ht;
    sum += alpha1*exp(-alpha2*((x1-e[0])*(x1-e[0])+(x2-e[1])*(x2-e[1]))-alpha3*t);
    return sum;
}
