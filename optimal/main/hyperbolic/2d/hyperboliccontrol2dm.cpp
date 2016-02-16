#include "hyperboliccontrol2dm.h"

void HyperbolicControl2DM::main()
{
    HyperbolicControl2DM hc;
    DoubleVector x(/*2*hc.L + */(hc.M+1)*hc.L);

    for (unsigned int k=0; k<=hc.M; k++)
    {
        x[/*2*hc.L + */0*(hc.M+1)+k] = 2.0;
        x[/*2*hc.L + */1*(hc.M+1)+k] = 2.0;
    }

    double min_step = 1.0;
    double gold_eps = 0.001;

    ConjugateGradient cg;
    cg.setFunction(&hc);
    cg.setGradient(&hc);
    cg.setEpsilon1(0.001);
    cg.setEpsilon2(0.001);
    cg.setEpsilon3(0.001);
    cg.setR1MinimizeEpsilon(min_step, gold_eps);
    cg.setPrinter(&hc);
    cg.setProjection(&hc);
    cg.setNormalize(true);
    //cg.showEndMessage(false);
    cg.calculate(x);

    DoubleVector g1(x.size());
    DoubleVector g2(x.size());
    FILE *file = fopen("gradients1.txt", "a");
    IPrinter::printDateTime(file);
    double h = 0.01;
    fprintf(file, "L: %d h:%f\n", hc.L, h);
    IPrinter::printVector(x, "v1:", (hc.M+1), 0*(hc.M+1), 0*(hc.M+1)+hc.M, file);
    IPrinter::printVector(x, "v2:", (hc.M+1), 1*(hc.M+1), 1*(hc.M+1)+hc.M, file);

    hc.gradient(x, g1);
    fprintf(file, "Analytic Gradients: %.20f\n", g1.L2Norm());
    g1.L2Normalize();
    IPrinter::printVector(g1, "g1:", (hc.M+1), 0*(hc.M+1), 0*(hc.M+1)+hc.M, file);
    IPrinter::printVector(g1, "g2:", (hc.M+1), 1*(hc.M+1), 1*(hc.M+1)+hc.M, file);

    IGradient::Gradient(&hc, h, x, g2);
    fprintf(file, "Numerical Gradients: %.20f\n", g2.L2Norm());
    g2.L2Normalize();
    IPrinter::printVector(g2, "g1:", (hc.M+1), 0*(hc.M+1), 0*(hc.M+1)+hc.M, file);
    IPrinter::printVector(g2, "g2:", (hc.M+1), 1*(hc.M+1), 1*(hc.M+1)+hc.M, file);
    IPrinter::printDateTime(file);
    fprintf(file, "--------------------------------------------------------------------\n");
    fclose(file);
}

HyperbolicControl2DM::HyperbolicControl2DM()
{
    x10 = 0.0;
    x11 = 1.0;
    x20 = 0.0;
    x21 = 1.0;
    t0 = 0.0;
    t1 = 1.0;
    N1 = 100;
    N2 = 100;
    M = 200;
    h1 = (x11 - x10) / N1;
    h2 = (x21 - x20) / N2;
    ht = (t1 - t0) / M;

    L = 2;
    g.resize(2*L);
    g[0] = 0.3;
    g[1] = 0.4;
    g[2] = 0.7;
    g[3] = 0.7;

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
}

HyperbolicControl2DM::~HyperbolicControl2DM()
{

}

double HyperbolicControl2DM::fx(double x)
{
    return 0.0;
}

double HyperbolicControl2DM::fx(const DoubleVector &v)
{
    px = &v;
    DoubleCube c;
    IHyperbolicEquation2D::calculateU1(c, h1, h2, ht, N1, N2, M);

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

    sum = alpha0*sum1 + alpha1*sum2;
    return sum;// + norm(v);
}

double HyperbolicControl2DM::norm(const DoubleVector& v) const
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

void HyperbolicControl2DM::gradient(const DoubleVector &x, DoubleVector &g)
{
    px = &x;
    DoubleCube u;
    IHyperbolicEquation2D::calculateU1(u, h1, h2, ht, N1, N2, M);

    pu = &u;
    DoubleCube p;
    IBackwardHyperbolicEquation2D::calculateU1(p, h1, h2, ht, N1, N2, M);

    unsigned int i,j;
    for (unsigned int k=0; k<=M; k++)
    {
        i = (unsigned int)round(g[0]/h1);
        j = (unsigned int)round(g[1]/h2);
        g[/*2*L+*/0*(M+1)+k] = -p[k][j][i];

        i = (unsigned int)round(g[0]/h1);
        j = (unsigned int)round(g[1]/h2);
        g[/*2*L+*/1*(M+1)+k] = -p[k][j][i];
    }
}

double HyperbolicControl2DM::fi1(unsigned int i, unsigned int j) const
{
    return 0.0;
}

double HyperbolicControl2DM::fi2(unsigned int i, unsigned int j) const
{
    return 0.0;
}

double HyperbolicControl2DM::m1(unsigned int j, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2DM::m2(unsigned int j, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2DM::m3(unsigned int j, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2DM::m4(unsigned int j, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2DM::f(unsigned int i, unsigned int j, unsigned int k) const
{
    double sum = 0.0;
    static double sgm1 = 3.0*h1;
    static double sgm2 = 3.0*h2;
    static double gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    static double gause_b = 2.0*sgm1*sgm2;

    double x1 = i*h1;
    double x2 = j*h2;
    const DoubleVector &x = *px;

    double _v1 = x[/*2*L+*/0*(M+1)+k];
    double _v2 = x[/*2*L+*/1*(M+1)+k];

    sum += _v1 * gause_a * exp(-((x1-g[0])*(x1-g[0]) + (x2-g[1])*(x2-g[1]))/gause_b);
    sum += _v2 * gause_a * exp(-((x1-g[2])*(x1-g[2]) + (x2-g[3])*(x2-g[3]))/gause_b);

    sum += fxt(i, j, k);

    return sum;
}

double HyperbolicControl2DM::bfi1(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu)[M];
    const DoubleMatrix &u1 = (*pu)[M-2];
    return -2.0 * alpha1 * ( (u0[j][i]-u1[j][i])/(2.0*ht) - U1);
}

double HyperbolicControl2DM::bfi2(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu)[M];
    return +2.0 * (u0[j][i] - U0) + qamma*(bfi1(i,j));
}

double HyperbolicControl2DM::bm1(unsigned int j, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2DM::bm2(unsigned int j, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2DM::bm3(unsigned int i, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2DM::bm4(unsigned int i, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2DM::bf(unsigned int i, unsigned int j, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2DM::fxt(unsigned int i, unsigned int j, unsigned int k) const
{
    double sum = 0.0;
    double x1 = i*h1;
    double x2 = j*h2;
    double t = k*ht;
    sum += alpha1*exp(-alpha2*((x1-e[0])*(x1-e[0])+(x2-e[1])*(x2-e[1]))-alpha3*t);
    return sum;
}

void HyperbolicControl2DM::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double alpha, RnFunction *fn) const
{
    C_UNUSED(g);
    C_UNUSED(alpha);
    printf("J[%d]: %.16f\n", i, fn->fx(x));
}

void HyperbolicControl2DM::project(DoubleVector &x, int i)
{
    //    if (i==0)
    //    {
    //        if (x[i] <= 0.0) x[i] = 0.0 + h1;
    //        if (x[i] >= 0.5) x[i] = 0.5 - h1;
    //    }
    //    if (i==1)
    //    {
    //        if (x[i] <= 0.0) x[i] = 0.0 + h1;
    //        if (x[i] >= 0.5) x[i] = 0.5 - h1;
    //    }
    //    if (i==2)
    //    {
    //        if (x[i] <= 0.5) x[i] = 0.5 + h1;
    //        if (x[i] >= 1.0) x[i] = 1.0 - h1;
    //    }
    //    if (i==3)
    //    {
    //        if (x[i] <= 0.5) x[i] = 0.5 + h1;
    //        if (x[i] >= 1.0) x[i] = 1.0 - h1;
    //    }
    //    if (i>3)
    {
        if (x[i] <= vd) x[i] = vd;
        if (x[i] >= vu) x[i] = vu;
    }
}
