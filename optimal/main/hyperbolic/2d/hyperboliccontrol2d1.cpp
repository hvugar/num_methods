#include "hyperboliccontrol2d1.h"

FILE *file;
void HyperbolicControl2D1::main()
{
    file = fopen("result.txt", "a");
    HyperbolicControl2D1 hc;
    //    for (double t=0.5; t<=2.1; t+=0.1)
    {
        hc.fx(1.4);
        fputs("-----------------------------------------------------------------------------------------------------------\n", file);
    }
    fclose(file);
}

HyperbolicControl2D1::HyperbolicControl2D1()
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
    h = 0.001;
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

    //double sigma = 0.0;

    a1 = 1.0;
    a2 = 1.0;
    count = 0;

    U0.resize(N2+1); for (unsigned int j=0; j<=N2; j++) U0[j].resize(N1+1);
    U1.resize(N2+1); for (unsigned int j=0; j<=N2; j++) U1[j].resize(N1+1);
    c.resize(2*L);

    // Initializing........................................
#ifdef ONLY_POWER
    DoubleVector x0((M+1)*L);
    for (unsigned int k=0; k<=M; k++)
    {
        x0[0*(M+1)+k] = v1(k*ht);
        x0[1*(M+1)+k] = v2(k*ht);
    }
    c[0] = 0.2;
    c[1] = 0.3;
    c[2] = 0.6;
    c[3] = 0.6;
#endif
#if defined(ONLY_COORDINATE)
    DoubleVector x0(2*L);
    x0[0] = 0.2;
    x0[1] = 0.3;
    x0[2] = 0.6;
    x0[3] = 0.6;
#endif

#if defined(POWER_COORDINATE)
    DoubleVector x0(2*L + (M+1)*L);
    for (unsigned int k=0; k<=M; k++)
    {
        x0[2*L + 0*(M+1)+k] = v1(k*ht);
        x0[2*L + 1*(M+1)+k] = v2(k*ht);
    }
    x0[0] = 0.2;
    x0[1] = 0.3;
    x0[2] = 0.6;
    x0[3] = 0.6;
#endif

    px = &x0;
    DoubleCube c;
    IHyperbolicEquation2D::calculateU1(c, h1, h2, ht, N1, N2, M, a1, a2, qamma);
    for (unsigned int j=0; j<=N2; j++)
    {
        for (unsigned int i=0; i<=N1; i++)
        {
            U0[j][i] = c[M][j][i];
            U1[j][i] = (c[M][j][i]-c[M-2][j][i])/(2.0*ht);
        }
    }

    IPrinter::printMatrix(U0, 10, 10, NULL, file);
    fprintf(file, "------------------\n");
    IPrinter::printMatrix(U1, 10, 10, NULL, file);
    fprintf(file, "------------------\n");
    //..................................................................
}

double HyperbolicControl2D1::fx(double T)
{
    t1 = T;
    M  = (unsigned)ceil((t1 - t0)/ht);

#ifdef ONLY_POWER
    DoubleVector x((M+1)*L);
    for (unsigned int k=0; k<=M; k++)
    {
        x[0*(M+1)+k] = 0.0;//v1(k*ht);
        x[1*(M+1)+k] = 0.0;//v2(k*ht);
    }
    c[0] = 0.2;
    c[1] = 0.3;
    c[2] = 0.6;
    c[3] = 0.6;
#endif

#if defined(ONLY_COORDINATE)
    DoubleVector x(2*L);
    x[0] = 0.3;
    x[1] = 0.4;
    x[2] = 0.7;
    x[3] = 0.7;
#endif

#ifdef POWER_COORDINATE
    DoubleVector x(2*L + (M+1)*L);
    for (unsigned int k=0; k<=M; k++)
    {
        x[2*L + 0*(M+1)+k] = 5.0;
        x[2*L + 1*(M+1)+k] = 5.0;
    }
    x[0] = 0.30;
    x[1] = 0.40;
    x[2] = 0.70;
    x[3] = 0.70;
#endif

    print(0, x, x, 0.0, this);
    // limits of v
    //vd = -2.0;
    //vu = +2.0;

    double min_step = 1.0;
    double gold_eps = 0.001;

    ConjugateGradient cg;
    cg.setFunction(this);
    cg.setGradient(this);
    cg.setEpsilon1(0.0001);
    cg.setEpsilon2(0.0001);
    cg.setEpsilon3(0.0001);
    cg.setR1MinimizeEpsilon(min_step, gold_eps);
    cg.setPrinter(this);
    cg.setProjection(this);
#ifdef POWER_COORDINATE
    // cg.setNormalize(false);
#endif
    cg.showEndMessage(true);
    cg.calculate(x);

    double rf = fx(x);
    count--;

    //    DoubleVector ag(x.size());
    //    gradient(x, ag);
    //    DoubleVector agx = ag.mid(0, 3);
    //    DoubleVector agv = ag.mid(4, ag.size()-1);
    //    agx.L2Normalize();
    //    agv.L2Normalize();

    //    DoubleVector ng(x.size());
    //    IGradient::Gradient(this, h, x, ng);
    //    DoubleVector ngx = ng.mid(0, 3);
    //    DoubleVector ngv = ng.mid(4, ng.size()-1);
    //    ngx.L2Normalize();
    //    ngv.L2Normalize();

    //FILE *file = stdout;//fopen("gradients_xv.txt", "a");
    //    fprintf(file, "--------------------------------------------------------------------\n");
    //    IPrinter::printDateTime(file);
    //    fprintf(file, "T: %f L: %d h:%f Functional: %.20f N1: %d N2: %d M: %d h1: %f h2: %f ht: %f\n", t1, L, h, rf, N1, N2, M, h1, h2, ht);
    //    fprintf(file, "x: %.8f %.8f %.8f %.8f\n", x[0], x[1], x[2], x[3]);
    //    fprintf(file, "AGx: %.8f %.8f %.8f %.8f\n", agx[0], agx[1], agx[2], agx[3]);
    //    fprintf(file, "NGx: %.8f %.8f %.8f %.8f\n", ngx[0], ngx[1], ngx[2], ngx[3]);
    //    unsigned int part = 10;//(M+1)
    //    IPrinter::printVector(x,   "v1: ", part, 0*(M+1)+2*L, 0*(M+1)+2*L+M, file);
    //    IPrinter::printVector(agv, "AG1:", part, 0*(M+1),     0*(M+1)+M,     file);
    //    IPrinter::printVector(ngv, "NG1:", part, 0*(M+1),     0*(M+1)+M,     file);
    //    IPrinter::printVector(x,   "v2: ", part, 1*(M+1)+2*L, 1*(M+1)+2*L+M, file);
    //    IPrinter::printVector(agv, "AG2:", part, 1*(M+1),     1*(M+1)+M,     file);
    //    IPrinter::printVector(ngv, "NG2:", part, 1*(M+1),     1*(M+1)+M,     file);
    //    IPrinter::printDateTime(file);
    //    fprintf(file, "U\n");
    //    DoubleCube c;
    //    IHyperbolicEquation2D::calculateU1(c, h1, h2, ht, N1, N2, M, a, a, qamma);
    //    IPrinter::printMatrix(c[c.size()-1], N2, N1, NULL, file);
    //    fclose(file);
    x.clear();

    return rf;
}

double HyperbolicControl2D1::fx(const DoubleVector &x)
{
    count++;
    px = &x;
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
            sum1 = sum1 + k * (u0[j][i]-U0[j][i]) * (u0[j][i]-U0[j][i]);

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
            sum2 = sum2 + k * ( (u0[j][i]-u1[j][i])/(2.0*ht) - U1[j][i]) * ( (u0[j][i]-u1[j][i])/(2.0*ht) - U1[j][i] );
        }
    }
    sum2 = h1*h2*sum2;

    sum = sum1 + alpha0*sum2;

#if defined(ONLY_POWER) || defined(POWER_COORDINATE)
    sum = sum + norm(x);
#endif
    return sum;
}

double HyperbolicControl2D1::norm(const DoubleVector& x) const
{
    double nrm = 0.0;
    for (unsigned int k=0; k<=M; k++)
    {
        double betta = 1.0;
        if (k==0 || k==M) betta = 0.5;
#if defined(ONLY_POWER)
        double m1 = (x[0*(M+1)+k] - v1(k*ht));
        double m2 = (x[1*(M+1)+k] - v2(k*ht));
        nrm += (betta * m1*m1 + betta * m2*m2);
#endif
#if defined(POWER_COORDINATE)
        double m1 = (x[2*L + 0*(M+1)+k] - v1(k*ht));
        double m2 = (x[2*L + 1*(M+1)+k] - v2(k*ht));
        nrm += (betta * m1*m1 + betta * m2*m2);
#endif
    }
    //    printf("%.f\n", nrm);
    nrm = ht * nrm;
    return nrm;
}

void HyperbolicControl2D1::gradient(const DoubleVector &x, DoubleVector &g)
{
    px = &x;
    DoubleCube u;
    IHyperbolicEquation2D::calculateU1(u, h1, h2, ht, N1, N2, M, a1, a2, qamma);

    pu = &u;
    DoubleCube p;
    IBackwardHyperbolicEquation2D::calculateU1(p, h1, h2, ht, N1, N2, M, a1, a2, qamma);

#if defined(ONLY_POWER)
    for (unsigned int k=0; k<=M; k++)
    {
        unsigned int i,j;
        i = (unsigned int)round(c[0]/h1);
        j = (unsigned int)round(c[1]/h2);
        g[0*(M+1)+k] = -p[k][j][i] + 2.0 * (x[0*(M+1)+k]-v1(k*ht));

        i = (unsigned int)round(c[2]/h1);
        j = (unsigned int)round(c[3]/h2);
        g[1*(M+1)+k] = -p[k][j][i] + 2.0 * (x[1*(M+1)+k]-v2(k*ht));
    }
#endif

#if defined(ONLY_COORDINATE)
    for (unsigned int k=0; k<=M; k++)
    {
        double psiX1;
        double psiX2;
        double m = 1.0;
        if (k==0 || k==M) m *= 0.5;
        psiDerivative(psiX1, psiX2, x[0], x[1], p[k]);
        g[0] = g[0] + m * v1(k*ht) * psiX1;
        g[1] = g[1] + m * v1(k*ht) * psiX2;
        psiDerivative(psiX1, psiX2, x[2], x[3], p[k]);
        g[2] = g[2] + m * v2(k*ht) * psiX1;
        g[3] = g[3] + m * v2(k*ht) * psiX2;
    }
    g[0] = -ht*g[0];
    g[1] = -ht*g[1];
    g[2] = -ht*g[2];
    g[3] = -ht*g[3];
#endif

#if defined(POWER_COORDINATE)

    for (unsigned int k=0; k<=M; k++)
    {
        double psiX1;
        double psiX2;
        double m = 1.0;
        if (k==0 || k==M) m *= 0.5;
        psiDerivative(psiX1, psiX2, x[0], x[1], p[k]);
        double _v1 = x[2*L+0*(M+1)+k];
        g[0] = g[0] + m * _v1 * psiX1;
        g[1] = g[1] + m * _v1 * psiX2;
        psiDerivative(psiX1, psiX2, x[2], x[3], p[k]);
        double _v2 = x[2*L+1*(M+1)+k];
        g[2] = g[2] + m * _v2 * psiX1;
        g[3] = g[3] + m * _v2 * psiX2;
    }
    g[0] = -ht*g[0];
    g[1] = -ht*g[1];
    g[2] = -ht*g[2];
    g[3] = -ht*g[3];

    for (unsigned int k=0; k<=M; k++)
    {
        unsigned int i,j;
        i = (unsigned int)round(x[0]/h1);
        j = (unsigned int)round(x[1]/h2);
        g[2*L+0*(M+1)+k] = -p[k][j][i] + 2.0 * (x[2*L+0*(M+1)+k]-v1(k*ht));

        i = (unsigned int)round(x[2]/h1);
        j = (unsigned int)round(x[3]/h2);
        g[2*L+1*(M+1)+k] = -p[k][j][i] + 2.0 * (x[2*L+1*(M+1)+k]-v2(k*ht));
    }
#endif
}

void HyperbolicControl2D1::psiDerivative(double &psiX1, double &psiX2, double x1, double x2, const DoubleMatrix &psi)
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

double HyperbolicControl2D1::fi1(unsigned int i, unsigned int j) const
{
    return u(i, j, 0);
}

double HyperbolicControl2D1::fi2(unsigned int i, unsigned int j) const
{
    return 0.0;
}

double HyperbolicControl2D1::m1(unsigned int j, unsigned int k) const
{
    return u(0, j, k);
}

double HyperbolicControl2D1::m2(unsigned int j, unsigned int k) const
{
    return u(N1, j, k);
}

double HyperbolicControl2D1::m3(unsigned int i, unsigned int k) const
{
    return u(i, 0, k);
}

double HyperbolicControl2D1::m4(unsigned int i, unsigned int k) const
{
    return u(i, N2, k);
}

double HyperbolicControl2D1::f(unsigned int i, unsigned int j, unsigned int k) const
{
    double sum = 0.0;
    static double sgm1 = 3.0*h1;
    static double sgm2 = 3.0*h2;
    static double gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    static double gause_b = 2.0*sgm1*sgm2;

    double x1 = i*h1;
    double x2 = j*h2;
    //double t = k*ht;
    const DoubleVector &x = *px;

#if defined(ONLY_POWER)
    double _v1 = x[0*(M+1)+k];
    double _v2 = x[1*(M+1)+k];
    sum += _v1 * gause_a * exp(-((x1-c[0])*(x1-c[0]) + (x2-c[1])*(x2-c[1]))/gause_b);
    sum += _v2 * gause_a * exp(-((x1-c[2])*(x1-c[2]) + (x2-c[3])*(x2-c[3]))/gause_b);
#endif

#if defined(ONLY_COORDINATE)
    double _v1 = v1(k*ht);
    double _v2 = v2(k*ht);;
    sum += _v1 * gause_a * exp(-((x1-x[0])*(x1-x[0]) + (x2-x[1])*(x2-x[1]))/gause_b);
    sum += _v2 * gause_a * exp(-((x1-x[2])*(x1-x[2]) + (x2-x[3])*(x2-x[3]))/gause_b);
#endif

#if defined(POWER_COORDINATE)
    double _v1 = x[2*L+0*(M+1)+k];
    double _v2 = x[2*L+1*(M+1)+k];
    sum += _v1 * gause_a * exp(-((x1-x[0])*(x1-x[0]) + (x2-x[1])*(x2-x[1]))/gause_b);
    sum += _v2 * gause_a * exp(-((x1-x[2])*(x1-x[2]) + (x2-x[3])*(x2-x[3]))/gause_b);
#endif

    sum += fxt(i, j, k);

    return sum;
}

double HyperbolicControl2D1::bfi1(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu)[M];
    const DoubleMatrix &u1 = (*pu)[M-2];
    return -2.0 * alpha1 * ( (u0[j][i]-u1[j][i])/(2.0*ht) - U1[j][i]);
}

double HyperbolicControl2D1::bfi2(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu)[M];
    return +2.0 * (u0[j][i] - U0[j][i]) + qamma*(bfi1(i,j));
}

double HyperbolicControl2D1::bm1(unsigned int j, unsigned int k) const
{
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2D1::bm2(unsigned int j, unsigned int k) const
{
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2D1::bm3(unsigned int i, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2D1::bm4(unsigned int i, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2D1::bf(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(k);
    return 0.0;
}

void HyperbolicControl2D1::print(unsigned int i, const DoubleVector& x, const DoubleVector &g, double alpha, RnFunction* fn) const
{
    C_UNUSED(g);
    C_UNUSED(alpha);

#ifdef ONLY_POWER
    printf("J[%d]: %18.16f   ", i, fn->fx(x));
    printf("c: %12.8f %12.8f %12.8f %12.8f\n", c[0], c[1], c[2], c[3]);
    IPrinter::printVector(x, "v1: ", 10, 0*(M+1), 0*(M+1)+M);
    IPrinter::printVector(x, "v2: ", 10, 1*(M+1), 1*(M+1)+M);
#endif

#ifdef ONLY_COORDINATE
    printf("J[%d]: %20.16f ", i, fn->fx(x));
    printf("x: %12.8f %12.8f %12.8f %12.8f\n", x[0], x[1], x[2], x[3]);
#endif

#ifdef POWER_COORDINATE
    double res = fn->fx(x);
    const_cast<HyperbolicControl2D1*>(this)->count--;
    printf("J[%d]: %20.16f T: %f Count: %d \n", i, res, t1, count);
    fprintf(file, "J[%d]: %20.16f T: %f Count: %d ", i, res, t1, count);
    fprintf(file, "x: %12.8f %12.8f %12.8f %12.8f\n", x[0], x[1], x[2], x[3]);
    IPrinter::printVector(x, "v1: ", 10, 2*L+0*(M+1), 2*L+0*(M+1)+M, file);
    IPrinter::printVector(x, "v2: ", 10, 2*L+1*(M+1), 2*L+1*(M+1)+M, file);
    IPrinter::printVector(x, "v1: ", M+1, 2*L+0*(M+1), 2*L+0*(M+1)+M, file);
    IPrinter::printVector(x, "v2: ", M+1, 2*L+1*(M+1), 2*L+1*(M+1)+M, file);
#endif
}

void HyperbolicControl2D1::project(DoubleVector &x, int i)
{
#if defined(ONLY_COORDINATE) || defined(POWER_COORDINATE)
    if (i<4)
    {
        if (x[i] <= 0.0) { x[i] = 0.0; }
        if (x[i] >= 1.0) { x[i] = 1.0; }
    }
#endif
}

double HyperbolicControl2D1::u(double i, double j, double k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t = k*ht;
    return x1*x1 + x2*x2 + t*t;
}

double HyperbolicControl2D1::fxt(unsigned int i, unsigned int j, unsigned int k) const
{
    double sum = 0.0;
    double x1 = i*h1;
    double x2 = j*h2;
    double t = k*ht;
    sum += alpha1*exp(-alpha2*((x1-e[0])*(x1-e[0])+(x2-e[1])*(x2-e[1]))-alpha3*t);
    return sum;
}

double HyperbolicControl2D1::v1(double t) const
{
    return 10.0*t;
}

double HyperbolicControl2D1::v2(double t) const
{
    return 12.0*t;
}
