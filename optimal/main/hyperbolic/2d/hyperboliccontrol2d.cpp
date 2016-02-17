#include "hyperboliccontrol2d.h"

void HyperbolicControl2D::main()
{
    HyperbolicControl2D hc;
    DoubleVector v(hc.L*(hc.M+1));
    for (unsigned int k=0; k<=hc.M; k++)
    {
        v[0*(hc.M+1)+k] = 1.0;//hc.v1(k*hc.ht);
        //v[1*(hc.M+1)+k] = 0.0;//hc.v2(k*hc.ht);
        //v[2*(hc.M+1)+k] = 0.0;//hc.v3(k*hc.ht);
    }
    DoubleVector g1(v.size());
    DoubleVector g2(v.size());

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
    cg.setNormalize(true);
    cg.showEndMessage(false);
    //cg.calculate(v);

    FILE *file = fopen("gradients.txt", "a");
    IPrinter::printDateTime(file);
    IPrinter::printVector(v, "v1:", (hc.M+1), 0*(hc.M+1), 0*(hc.M+1)+hc.M, file);
    //IPrinter::printVector(v, "v2:", (hc.M+1), 1*(hc.M+1), 1*(hc.M+1)+hc.M, file);
    //IPrinter::printVector(v, "v3:", 10, 2*(hc.M+1), 2*(hc.M+1)+hc.M);
    double h = 0.01;
    hc.gradient(v, g1);
    fprintf(file, "L: %d h:%f\n", hc.L, h);

    fprintf(file, "Analytic Gradients: %.20f\n", g1.L2Norm());
    g1.L2Normalize();
    IPrinter::printVector(g1, "g1:", (hc.M+1), 0*(hc.M+1), 0*(hc.M+1)+hc.M, file);
    //IPrinter::printVector(g1, "g2:", (hc.M+1), 1*(hc.M+1), 1*(hc.M+1)+hc.M, file);
    //IPrinter::printVector(g1, "g3:", 10, 2*(hc.M+1), 2*(hc.M+1)+hc.M);

    IGradient::Gradient(&hc, h, v, g2);
    fprintf(file, "Numerical Gradients: %.20f\n", g2.L2Norm());
    g2.L2Normalize();
    IPrinter::printVector(g2, "g1:", (hc.M+1), 0*(hc.M+1), 0*(hc.M+1)+hc.M, file);
    //IPrinter::printVector(g2, "g2:", (hc.M+1), 1*(hc.M+1), 1*(hc.M+1)+hc.M, file);
    //IPrinter::printVector(g2, "g3:", 10, 2*(hc.M+1), 2*(hc.M+1)+hc.M);
    IPrinter::printDateTime(file);
    fprintf(file, "--------------------------------------------------------------------\n");
    fclose(file);
}

HyperbolicControl2D::HyperbolicControl2D()
{
    x10 = 0.0;
    x11 = 1.0;
    x20 = 0.0;
    x21 = 1.0;
    t0 = 0.0;
    t1 = 1.0;
    a1 = a2 = 1.0;
    N1 = 100;
    N2 = 100;
    M = 200;
    h1 = (x11 - x10) / N1;
    h2 = (x21 - x20) / N2;
    ht = (t1 - t0) / M;
    alpha0 = 1.0;
    alpha1 = 1.0;
    qamma = 1.0;

    U0.resize(N2+1); for (unsigned int j=0; j<=N2; j++) U0[j].resize(N1+1);
    U1.resize(N2+1); for (unsigned int j=0; j<=N2; j++) U1[j].resize(N1+1);

    L = 1;
    E.resize(2*L);
    E[0] = 0.3;//30
    E[1] = 0.4;//40
    //E[2] = 0.8;//80
    //E[3] = 0.7;//70
    //E[4] = 0.2;//20
    //E[5] = 0.8;//80

    //initialize
    U0.resize(N2+1); for (unsigned int j=0; j<=N2; j++) U0[j].resize(N1+1);
    U1.resize(N2+1); for (unsigned int j=0; j<=N2; j++) U1[j].resize(N1+1);

    DoubleVector v(L*(M+1));
    for (unsigned int k=0; k<=M; k++)
    {
        v[0*(M+1)+k] = v1(k*ht);
        //v[1*(M+1)+k] = v2(k*ht);
        //v[2*(M+1)+k] = v3(k*ht);
    }
    //IPrinter::printVector(v, "v1:", 10, 0*(M+1), 0*(M+1)+M);
    //IPrinter::printVector(v, "v2:", 10, 1*(M+1), 1*(M+1)+M);

    pv = &v;
    DoubleCube c;
    IHyperbolicEquation2D::calculateU1(c, h1, h2, ht, N1, N2, M);
    DoubleMatrix &u0 = c[M];
    DoubleMatrix &u1 = c[M-2];
    for (unsigned int j=0; j<=N2; j++)
    {
        for (unsigned int i=0; i<=N1; i++)
        {
            U0[j][i] = u0[j][i];
            U1[j][i] = (u0[j][i] - u1[j][i])/(2.0*ht);
        }
    }

    //    printf("U0:\n");
    //    IPrinter::printMatrix(U0);
    //    printf("U1:\n");
    //    IPrinter::printMatrix(U1);

    FILE *file1 = fopen("U0.txt", "w");
    IPrinter::printMatrix(U0, N2, N1, NULL, file1);
    fclose(file1);

    FILE *file2 = fopen("U1.txt", "w");
    IPrinter::printMatrix(U1, N2, N1, NULL, file2);
    fclose(file2);
}

HyperbolicControl2D::~HyperbolicControl2D()
{
}

double HyperbolicControl2D::fx(double x)
{
    return 0.0;
}

double HyperbolicControl2D::fx(const DoubleVector &v)
{
    pv = &v;
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
            sum2 = sum2 + k * ( (u0[j][i]-u1[j][i])/(2.0*ht)-U1[j][i])*((u0[j][i]-u1[j][i])/(2.0*ht)-U1[j][i]);
        }
    }
    sum2 = h1*h2*sum2;

    sum = alpha0*sum1 + alpha1*sum2;
    return sum;// + norm(v);
}

double HyperbolicControl2D::norm(const DoubleVector& v) const
{
    double nrm = 0.0;
    for (unsigned int k=0; k<=M; k++)
    {
        double betta = 1.0;
        if (k==0 || k==M) betta = 0.5;
        double m1 = (v[0*(M+1)+k] - v1(k*ht));
        nrm += betta * m1*m1;
        //nrm += betta*(v[1*(M+1)+k] - v2(k*ht))*(v[1*(M+1)+k] - v2(k*ht));
        //nrm += betta*(v[2*(M+1)+k] - v3(k*ht))*(v[2*(M+1)+k] - v3(k*ht));
    }
    nrm = ht * nrm;
    return nrm;
}

void HyperbolicControl2D::gradient(const DoubleVector &v, DoubleVector &g)
{
    pv = &v;
    DoubleCube u;
    IHyperbolicEquation2D::calculateU1(u, h1, h2, ht, N1, N2, M);

    pu = &u;
    DoubleCube p;
    IBackwardHyperbolicEquation2D::calculateU2(p, h1, h2, ht, N1, N2, M);

    unsigned int i,j;
    for (unsigned int k=0; k<=M; k++)
    {
        i = (unsigned int)round(E[0]/h1);
        j = (unsigned int)round(E[1]/h2);
        g[0*(M+1)+k] = -p[k][j][i];// + 2.0*(v[0*(M+1)+k] - v1(k*ht));
        //i = (unsigned int)round(E[2]/h1);
        //j = (unsigned int)round(E[3]/h2);
        //g[1*(M+1)+k] = -p[k][j][i] + 2.0*(v[1*(M+1)+k] - v2(k*ht));
        //i = (unsigned int)round(E[4]/h1);
        //j = (unsigned int)round(E[5]/h2);
        //g[2*(M+1)+k] = -p[k][j][i] + 2.0*(v[2*(M+1)+k] - v3(k*ht));
    }
}

double HyperbolicControl2D::fi1(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return u(i, j, 0);
}

double HyperbolicControl2D::fi2(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl2D::m1(unsigned int j, unsigned int k) const
{
    return u(0, j, k);
}

double HyperbolicControl2D::m2(unsigned int j, unsigned int k) const
{
    return u(N1, j, k);
}

double HyperbolicControl2D::m3(unsigned int i, unsigned int k) const
{
    return u(i, 0, k);
}

double HyperbolicControl2D::m4(unsigned int i, unsigned int k) const
{
    return u(i, N2, k);
}

double HyperbolicControl2D::f(unsigned int i, unsigned int j, unsigned int k) const
{
    double sum = 0.0;
    //return 2.0 - 2.0*(a1*a1) - 2.0*(a2*a2);

    static double sgm1 = 3.0*h1;
    static double sgm2 = 3.0*h2;
    static double gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    static double gause_b = 2.0*sgm1*sgm2;

    double x1 = i*h1;
    double x2 = j*h2;
    //double t  = 0.5*k*ht;
    const DoubleVector &v = *pv;

    double _v1 = v[0*(M+1)+k];
    //double _v2 = v[1*(M+1)+k];
    //double _v3 = v[2*(M+1)+k];

    sum += _v1 * gause_a * exp(-((x1-E[0])*(x1-E[0]) + (x2-E[1])*(x2-E[1]))/gause_b);
    //sum += _v2 * gause_a * exp(-((x1-E[2])*(x1-E[2]) + (x2-E[3])*(x2-E[3]))/gause_b);
    //sum += _v3 * gause_a * exp(-((x1-E[4])*(x1-E[4]) + (x2-E[5])*(x2-E[5]))/gause_b);

    sum += fxt(i, j, k);

    return sum;
}

double HyperbolicControl2D::fxt(unsigned int i, unsigned int j, unsigned int k) const
{
    double sum = 0.0;
    double x1 = i*h1;
    double x2 = j*h2;
    double t = k*ht;

    double alpha_1 = 1.0;
    double alpha_2 = 1.0;
    double alpha_3 = 1.0;

    sum += alpha_1*exp(-alpha_2*((x1-0.2)*(x1-0.2)+(x2-0.2)*(x2-0.2))-alpha_3*t);
    //sum += alpha_1*exp(-alpha_2*((x1-E[2])*(x1-E[2])+(x2-E[3])*(x2-E[3]))-alpha_3*t);
    //sum += alpha_1*exp(-alpha_2*((x1-E[4])*(x1-E[4])+(x2-E[5])*(x2-E[5]))-alpha_3*t);

    return sum;
}

double HyperbolicControl2D::bfi1(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu)[M];
    const DoubleMatrix &u1 = (*pu)[M-2];
    return -2.0 * alpha1 * ( (u0[j][i]-u1[j][i])/(2.0*ht) - U1[j][i]);
}

double HyperbolicControl2D::bfi2(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu)[M];
    return +2.0 * alpha0 * (u0[j][i] - U0[j][i]) + qamma*(bfi1(i,j));
}

double HyperbolicControl2D::bm1(unsigned int j, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2D::bm2(unsigned int j, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2D::bm3(unsigned int i, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2D::bm4(unsigned int i, unsigned int k) const
{
    return 0.0;
}

double HyperbolicControl2D::bf(unsigned int i, unsigned int j, unsigned int k) const
{
    return 0.0;
}

void HyperbolicControl2D::print(unsigned int i, const DoubleVector &v, const DoubleVector &g, double alpha, RnFunction *fn) const
{
    C_UNUSED(g);
    C_UNUSED(alpha);
    printf("J[%d]: %.16f\n", i, fn->fx(v));
}

double HyperbolicControl2D::v1(double t) const
{
    return 20.0*t;
}

double HyperbolicControl2D::v2(double t) const
{
    return 10.0*t;
}

double HyperbolicControl2D::v3(double t) const
{
    return 15.0*t;
}

double HyperbolicControl2D::u(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t = k*ht;
    return x1*x1 + x2*x2 + t*t;
}
