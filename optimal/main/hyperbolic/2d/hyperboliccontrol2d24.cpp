#include "hyperboliccontrol2d24.h"

void HyperbolicControl2D24::main(int argc, char ** argv)
{
//    for (int i=0; i<argc; i++)
//    {

//    }


    C_UNUSED(argc);
    C_UNUSED(argv);

    HyperbolicControl2D24 hc;
    hc.file = fopen("20160409_1.txt", "w");
    //for (double t=0.1; t<=10.1; t+=0.1)
    {
        hc.fx(1.0);
    }
    fclose(hc.file);
}

HyperbolicControl2D24::HyperbolicControl2D24()
{
    x10 = x20 = t0 = 0.0;
    x11 = x21 = t1 = 1.0;

    h1 = h2 = 0.01;
    ht = 0.005;

    N1 = (unsigned)round((x11 - x10)/h1);
    N2 = (unsigned)round((x21 - x20)/h2);
    M  = (unsigned)round((t1 - t0)/ht);
    L = 2;

    e.resize(2);
    e[0] = 0.2;
    e[1] = 0.2;

    alpha0 = 1.0;

    alpha1 = 10.0;
    alpha2 = 2.0;
    alpha3 = 1.0;

    qamma = 0.2;

    a1 = 1.0;
    a2 = 1.0;

    e.resize(2);
    e[0] = 0.2;
    e[1] = 0.2;

    U0 = 0.0;
    U1 = 0.0;
    a1 = a2 = 1.0;

    px = NULL;
    pu = NULL;
    file = stdout;
}

double HyperbolicControl2D24::fx(double t)
{
    t1 = t;
    M  = (unsigned)ceil((t1 - t0)/ht);

    printf("t: %f M: %d\n--------------------\n", t1, M);

    DoubleVector x0(2*L + (M+1)*L);
    for (unsigned int k=0; k<=M; k++)
    {
        x0[2*L+0*(M+1)+k] = 0.0;
        x0[2*L+1*(M+1)+k] = 0.0;
    }
    x0[0] = 0.25;//0.3
    x0[1] = 0.25;//0.4
    x0[2] = 0.85;//0.7
    x0[3] = 0.85;//0.7

    printGradients(x0, 0, file);

    double min_step = 10.0;
    double gold_eps = 0.001;
    ConjugateGradient cg;
    cg.setFunction(this);
    cg.setGradient(this);
    cg.setEpsilon1(0.001);
    cg.setEpsilon2(0.001);
    cg.setEpsilon3(0.001);
    cg.setR1MinimizeEpsilon(min_step, gold_eps);
    cg.setPrinter(this);
    cg.setProjection(this);
    cg.setNormalize(true);
    //cg.showEndMessage(false);
    cg.calculate(x0);

    printGradients(x0, 0, file);

    double rf = fx(x0);
    fprintf(file, "%f %.16f\n", t, rf);
    fprintf(file, "%f %f %f %f\n", x0[0], x0[1], x0[2], x0[3]);
    IPrinter::printVector(x0, "v1", M+1, 2*L,       2*L+M,       file);
    IPrinter::printVector(x0, "v2", M+1, 2*L+(M+1), 2*L+(2*M+1), file);
    fflush(file);
    return rf;
}

void HyperbolicControl2D24::printGradients(const DoubleVector &x, unsigned int i, FILE* f) const
{
    printf("J[%d]: %.16f\n", i, const_cast<HyperbolicControl2D24*>(this)->fx(x));
    fprintf(f, "J[%d]: %.16f\n", i, const_cast<HyperbolicControl2D24*>(this)->fx(x));
    DoubleVector g(x.size());
    const_cast<HyperbolicControl2D24*>(this)->gradient(x, g);
    DoubleVector g1 = g.mid(0, 1);
    DoubleVector g2 = g.mid(2, 3);
    DoubleVector v1 = x.mid(4,   M+4);
    DoubleVector v2 = x.mid(M+5, 2*M+5);

    fprintf(f, "N: %10.6f G: %10.6f %10.6f P: %10.6f %10.6f\n", g1.L2Norm(), g1[0], g1[1], x[0], x[1]);
    g1.L2Normalize();
    fprintf(f, "N: %10.6f G: %10.6f %10.6f P: %10.6f %10.6f\n", g1.L2Norm(), g1[0], g1[1], x[0], x[1]);
    fprintf(f, "---\n");
    fprintf(f, "N: %10.6f G: %10.6f %10.6f P: %10.6f %10.6f\n", g2.L2Norm(), g2[0], g2[1], x[2], x[3]);
    g2.L2Normalize();
    fprintf(f, "N: %10.6f G: %10.6f %10.6f P: %10.6f %10.6f\n", g2.L2Norm(), g2[0], g2[1], x[2], x[3]);
    fprintf(f, "---\n");
    fprintf(f, "N:  %10.6f\n", v1.L2Norm());
    IPrinter::printVector(v1, "vg1:", 10, 0, M, f);
    v1.L2Normalize();
    IPrinter::printVector(v1, "ng1:", 10, 0, M, f);
    fprintf(f, "---\n");
    fprintf(f, "N:  %10.6f\n", v2.L2Norm());
    IPrinter::printVector(v2, "vg2:", 10, 0, M, f);
    v2.L2Normalize();
    IPrinter::printVector(v2, "ng2:", 10, 0, M, f);
    fprintf(f, "--------------------------------------------------------------------------------------------------\n");
    fflush(f);
}

void HyperbolicControl2D24::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double alpha, RnFunction* fn) const
{
    C_UNUSED(g);
    C_UNUSED(alpha);
    printGradients(x, i, file);
//    fprintf(file, "J[%d]: %.16f\n", i, fn->fx(x));
//    fprintf(file, "%f %f %f %f\n", x[0], x[1], x[2], x[3]);
//    IPrinter::printVector(x, "v1", M+1, 2*L,       2*L+M,       file);
//    IPrinter::printVector(x, "v2", M+1, 2*L+(M+1), 2*L+(2*M+1), file);
//    fflush(file);

//    printf("J[%d]: %.16f\n", i, fn->fx(x));
//    printf("%.6f %.6f %.6f %.6f\n", x[0], x[1], x[2], x[3]);
//    IPrinter::printVector(x, "v1", 10, 2*L,       2*L+M);
//    IPrinter::printVector(x, "v2", 10, 2*L+(M+1), 2*L+(2*M+1));
//    fprintf(file, "-----------------------\n");
//    fflush(file);

//    DoubleVector g1(x.size());
//    const_cast<HyperbolicControl2D24*>(this)->gradient(x, g1);
//    DoubleVector ge = g1.mid(0, 3);
//    DoubleVector v1 = g1.mid(4,   M+4);
//    DoubleVector v2 = g1.mid(M+5, 2*M+5);
//    printf("%.6f %.6f %.6f %.6f %.6f\n", ge[0], ge[1], ge[2], ge[3], ge.L2Norm());
//    ge.L2Normalize();
//    printf("%.6f %.6f %.6f %.6f\n", ge[0], ge[1], ge[2], ge[3]);
//    v1.L2Normalize();
//    IPrinter::printVector(v1, "gv1", 10, 0, M);
//    v2.L2Normalize();
//    IPrinter::printVector(v2, "gv2", 10, 0, M);
//    printf("-----------------------\n");
}

double HyperbolicControl2D24::fx(const DoubleVector &x)
{
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

void HyperbolicControl2D24::gradient(const DoubleVector &x, DoubleVector &g)
{
    px = &x;
    DoubleCube u;
    IHyperbolicEquation2D::calculateU1(u, h1, h2, ht, N1, N2, M, a1, a2, qamma);

    pu = &u;
    DoubleCube p;
    IBackwardHyperbolicEquation2D::calculateU1(p, h1, h2, ht, N1, N2, M, a1, a2, qamma);

    // placement gradients
    g[0] = g[1] = g[2] = g[3] = 0.0;
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

    // power gradients
    unsigned int i,j;
    for (unsigned int k=0; k<=M; k++)
    {
        i = (unsigned int)round(x[0]/h1);
        j = (unsigned int)round(x[1]/h2);
        g[2*L+0*(M+1)+k] = -p[k][j][i];

        i = (unsigned int)round(x[2]/h1);
        j = (unsigned int)round(x[3]/h2);
        g[2*L+1*(M+1)+k] = -p[k][j][i];
    }
}

void HyperbolicControl2D24::psiDerivative(double &psiX1, double &psiX2, double x1, double x2, const DoubleMatrix &psi)
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

double HyperbolicControl2D24::initial1(unsigned int i, unsigned int j) const
{
    return 0.0;
}

double HyperbolicControl2D24::initial2(unsigned int i, unsigned int j) const
{
    return 0.0;
}

double HyperbolicControl2D24::boundary(unsigned int i, unsigned int j, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2D24::f(unsigned int i, unsigned int j, unsigned int k) const
{
    double sum = 0.0;

    static double sgm1 = 3.0*h1;
    static double sgm2 = 3.0*h2;
    static double gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    static double gause_b = 2.0*sgm1*sgm2;

    double x1 = i*h1;
    double x2 = j*h2;
    const DoubleVector x = *px;

    double _v1 = x[2*L+0*(M+1)+k];
    double _v2 = x[2*L+1*(M+1)+k];

    sum += _v1 * gause_a * exp(-((x1-x[0])*(x1-x[0]) + (x2-x[1])*(x2-x[1]))/gause_b);
    sum += _v2 * gause_a * exp(-((x1-x[2])*(x1-x[2]) + (x2-x[3])*(x2-x[3]))/gause_b);

    sum += fxt(i, j, k);
    return sum;
}

double HyperbolicControl2D24::binitial1(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu)[M];
    const DoubleMatrix &u1 = (*pu)[M-2];
    return -2.0 * alpha0 * ((u0[j][i]-u1[j][i])/(2.0*ht) - U1);
}

double HyperbolicControl2D24::binitial2(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu)[M];
    return +2.0 * (u0[j][i] - U0) + qamma*(binitial1(i,j));
}

double HyperbolicControl2D24::bboundary(unsigned int i, unsigned int j, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2D24::bf(unsigned int i, unsigned int j, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2D24::fxt(unsigned int i, unsigned int j, unsigned int k) const
{
    double sum = 0.0;
    double x1 = i*h1;
    double x2 = j*h2;
    double t = k*ht;
    sum += alpha1*exp(-alpha2*((x1-e[0])*(x1-e[0])+(x2-e[1])*(x2-e[1]))-alpha3*t);
    return sum;
}

void HyperbolicControl2D24::project(DoubleVector &x, int i)
{  
    if (i==0)
    {
        if (x[i] <= 0.0 + 5.0*h1) { x[i] = 0.0 + 5.0*h1; }
        if (x[i] >= 0.5 - 5.0*h1) { x[i] = 0.5 - 5.0*h1; }
    }
    if (i==1)
    {
        if (x[i] <= 0.0 + 5.0*h2) { x[i] = 0.0 + 5.0*h2; }
        if (x[i] >= 0.5 - 5.0*h2) { x[i] = 0.5 - 5.0*h2; }
    }
    if (i==2)
    {
        if (x[i] <= 0.5 + 5.0*h1) { x[i] = 0.5 + 5.0*h1; }
        if (x[i] >= 1.0 - 5.0*h1) { x[i] = 1.0 - 5.0*h1; }
    }
    if (i==3)
    {
        if (x[i] <= 0.5 + 5.0*h2) { x[i] = 0.5 + 5.0*h2; }
        if (x[i] >= 1.0 - 5.0*h2) { x[i] = 1.0 - 5.0*h2; }
    }

//    if (i==0)
//    {
//        if (x[i] <= 0.0) { x[i] = 0.0; }
//        if (x[i] >= 0.5) { x[i] = 0.5; }

//    }
//    if (i==1)
//    {
//        if (x[i] <= 0.0) { x[i] = 0.0; }
//        if (x[i] >= 0.5) { x[i] = 0.5; }
//    }
//    if (i==2)
//    {
//        if (x[i] <= 0.5) { x[i] = 0.5; }
//        if (x[i] >= 1.0) { x[i] = 1.0; }
//    }
//    if (i==3)
//    {
//        if (x[i] <= 0.5) { x[i] = 0.5; }
//        if (x[i] >= 1.0) { x[i] = 1.0; }
//    }

//    if (i>3)
//    {
//        if (x[i] < -2.0) x[i] = -2.0;
//        if (x[i] > +2.0) x[i] = +2.0;
//    }
}
