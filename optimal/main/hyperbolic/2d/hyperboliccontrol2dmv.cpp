#include "hyperboliccontrol2dmv.h"

void HyperbolicControl2DMV::main()
{
    HyperbolicControl2DMV hc;
    //    hc.fx(0.4);
    //    hc.fx(0.6);
    //    hc.fx(0.8);
    hc.fx(1.0);
    //    hc.fx(1.2);
    //    hc.fx(1.4);
    //    hc.fx(1.6);
    //    hc.fx(1.8);
    //    hc.fx(2.0);
}

HyperbolicControl2DMV::HyperbolicControl2DMV()
{
}

HyperbolicControl2DMV::~HyperbolicControl2DMV()
{
}

double HyperbolicControl2DMV::fx(double T)
{
    x10 = 0.0;
    x11 = 1.0;
    x20 = 0.0;
    x21 = 1.0;
    t0 = 0.0;
    t1 = T;

    h1 = 0.01;
    h2 = 0.01;
    ht = 0.005;

    N1 = (unsigned)ceil((x11 - x10)/h1);
    N2 = (unsigned)ceil((x21 - x20)/h2);
    M  = (unsigned)ceil((t1 - t0)/ht);

    printf("%d %d %d\n", N1, N2, M);

    L = 2;
    d.resize(2*L);
    d[0] = 0.3;
    d[1] = 0.4;
    d[2] = 0.7;
    d[3] = 0.7;

    e.resize(2);
    e[0] = 0.2;
    e[1] = 0.2;

    // limits of v
    vd = -2.0;
    vu = +2.0;

    alpha0 = 4.0;
    alpha1 = 5.0;
    alpha2 = 6.0;
    alpha3 = 4.0;
    qamma = 0.2;

    U0 = 0.0;
    U1 = 0.0;
    a = 1.0;

    DoubleVector v((M+1)*L);
    for (unsigned int k=0; k<=M; k++)
    {
        v[0*(M+1)+k] = 2.0;
        v[1*(M+1)+k] = 2.0;
    }

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
    //    cg.setNormalize(true);
    //    cg.showEndMessage(false);
    //    cg.calculate(v);

    double rf = fx(v);

    double h = 0.001;
    DoubleVector ag(v.size());
    DoubleVector ng(v.size());
    FILE *file = fopen("gradient_v.txt", "a");
    fprintf(file, "--------------------------------------------------------------------\n");
    IPrinter::printDateTime(file);

    gradient(v, ag);
    ag.L2Normalize();
    IGradient::Gradient(this, h, v, ng);
    ng.L2Normalize();

    fprintf(file, "T: %f L: %d h:%f Functional: %.20f N1: %d N2: %d M: %d h1: %f h2: %f ht: %f\n", t1, L, h, rf, N1, N2, M, h1, h2, ht);
    IPrinter::printVector(v,  "v1: ", (M+1), 0*(M+1), 0*(M+1)+M, file);
    IPrinter::printVector(ag, "AG1:", (M+1), 0*(M+1), 0*(M+1)+M, file);
    IPrinter::printVector(ng, "NG1:", (M+1), 0*(M+1), 0*(M+1)+M, file);

    IPrinter::printVector(v,  "v2: ", (M+1), 1*(M+1), 1*(M+1)+M, file);
    IPrinter::printVector(ag, "AG2:", (M+1), 1*(M+1), 1*(M+1)+M, file);
    IPrinter::printVector(ng, "NG2:", (M+1), 1*(M+1), 1*(M+1)+M, file);
    //fprintf(file, "Analytic Gradients: %.20f\n", g1.L2Norm());
    //fprintf(file, "Numerical Gradients: %.20f\n", g2.L2Norm());
    IPrinter::printDateTime(file);
    //    fprintf(file, "U\n");
    //    DoubleCube c;
    //    IHyperbolicEquation2D::calculateU1(c, h1, h2, ht, N1, N2, M, a, a, qamma);
    //    IPrinter::printMatrix(c[c.size()-1], N2, N1, NULL, file);
    fclose(file);

    v.clear();
    return rf;
}

double HyperbolicControl2DMV::fx(const DoubleVector &v)
{
    pv = &v;
    DoubleCube c;
    IHyperbolicEquation2D::calculateU1(c, h1, h2, ht, N1, N2, M, a, a, qamma);

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
            sum2 = sum2 + k * ( (u0[j][i]-u1[j][i])/(2.0*ht)-U1)*((u0[j][i]-u1[j][i])/(2.0*ht)-U1);
        }
    }
    sum2 = h1*h2*sum2;

    sum = sum1 + alpha0*sum2;
    return sum;// + norm(v);
}

double HyperbolicControl2DMV::norm(const DoubleVector& v) const
{
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

void HyperbolicControl2DMV::gradient(const DoubleVector &v, DoubleVector &g)
{
    pv = &v;
    DoubleCube u;
    IHyperbolicEquation2D::calculateU1(u, h1, h2, ht, N1, N2, M, a, a, qamma);

    IPrinter::printMatrix(u[u.size()-1]);

    pu = &u;
    DoubleCube p;
    IBackwardHyperbolicEquation2D::calculateU1(p, h1, h2, ht, N1, N2, M, a, a, qamma);

    puts("---");
    IPrinter::printMatrix(p[p.size()-1]);

    unsigned int i,j;
    for (unsigned int k=0; k<=M; k++)
    {
        i = (unsigned int)round(d[0]/h1);
        j = (unsigned int)round(d[1]/h2);
        g[0*(M+1)+k] = -p[k][j][i];

        i = (unsigned int)round(d[2]/h1);
        j = (unsigned int)round(d[3]/h2);
        g[1*(M+1)+k] = -p[k][j][i];
    }
}

double HyperbolicControl2DMV::fi1(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl2DMV::fi2(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl2DMV::m1(unsigned int j, unsigned int k) const
{
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2DMV::m2(unsigned int j, unsigned int k) const
{
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2DMV::m3(unsigned int i, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2DMV::m4(unsigned int i, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2DMV::f(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(k);
    double sum = 0.0;
    static double sgm1 = 3.0*h1;
    static double sgm2 = 3.0*h2;
    static double gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    static double gause_b = 2.0*sgm1*sgm2;

    double x1 = i*h1;
    double x2 = j*h2;
    const DoubleVector &v = *pv;

    double _v1 = v[0*(M+1)+k];
    double _v2 = v[1*(M+1)+k];

    sum += _v1 * gause_a * exp(-((x1-d[0])*(x1-d[0]) + (x2-d[1])*(x2-d[1]))/gause_b);
    sum += _v2 * gause_a * exp(-((x1-d[2])*(x1-d[2]) + (x2-d[3])*(x2-d[3]))/gause_b);

    sum += fxt(i, j, k);

    return sum;
}

double HyperbolicControl2DMV::bfi1(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu)[M];
    const DoubleMatrix &u1 = (*pu)[M-2];
    return -2.0 * alpha1 * ((u0[j][i]-u1[j][i])/(2.0*ht) - U1);
}

double HyperbolicControl2DMV::bfi2(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu)[M];
    return +2.0 * (u0[j][i] - U0) + qamma*(bfi1(i,j));
}

double HyperbolicControl2DMV::bm1(unsigned int j, unsigned int k) const
{
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2DMV::bm2(unsigned int j, unsigned int k) const
{
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2DMV::bm3(unsigned int i, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2DMV::bm4(unsigned int i, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2DMV::bf(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j );
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2DMV::fxt(unsigned int i, unsigned int j, unsigned int k) const
{
    double sum = 0.0;
    double x1 = i*h1;
    double x2 = j*h2;
    double t = k*ht;
    sum += alpha1*exp(-alpha2*((x1-e[0])*(x1-e[0])+(x2-e[1])*(x2-e[1]))-alpha3*t);
    return sum;
}

void HyperbolicControl2DMV::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double alpha, RnFunction *fn) const
{
    C_UNUSED(g);
    C_UNUSED(alpha);
    printf("J[%d]: %.16f\n", i, fn->fx(x));
}

void HyperbolicControl2DMV::project(DoubleVector &x, int i)
{
    if (x[i] < vd) x[i] = vd;
    if (x[i] > vu) x[i] = vu;
}
