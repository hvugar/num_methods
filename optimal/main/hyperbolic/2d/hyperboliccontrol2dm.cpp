#include "hyperboliccontrol2dm.h"

void HyperbolicControl2DM::main()
{
    HyperbolicControl2DM hc;
    //    hc.fx(0.4);
    //    hc.fx(0.6);
    //    hc.fx(0.8);
//    hc.fx(1.0);
        hc.fx(1.2);
    //    hc.fx(1.4);
    //    hc.fx(1.6);
    //    hc.fx(1.8);
    //    hc.fx(2.0);
}

HyperbolicControl2DM::HyperbolicControl2DM()
{
}

HyperbolicControl2DM::~HyperbolicControl2DM()
{
}

double HyperbolicControl2DM::fx(double T)
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
    L = 2;
    printf("T: %f %d %d %d\n", T, N1, N2, M);

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
    a = 1.0;

    DoubleVector x( 2*L + (M+1)*L);
    for (unsigned int k=0; k<=M; k++)
    {
        x[2*L + 0*(M+1)+k] = 2.0;
        x[2*L + 1*(M+1)+k] = 2.0;
    }
    x[0] = 0.3;
    x[1] = 0.4;
    x[2] = 0.7;
    x[3] = 0.7;

    h = 0.001;

    // limits of v
    vd = -2.0;
    vu = +2.0;

    double min_step = 0.1;
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
    //cg.setNormalize(false);
    cg.showEndMessage(true);
    cg.calculate(x);

    double rf = fx(x);

    DoubleVector ag(x.size());
    gradient(x, ag);
    DoubleVector agx = ag.mid(0, 3);
    DoubleVector agv = ag.mid(4, ag.size()-1);
    agx.L2Normalize();
    agv.L2Normalize();

    DoubleVector ng(x.size());
    IGradient::Gradient(this, h, x, ng);
    DoubleVector ngx = ng.mid(0, 3);
    DoubleVector ngv = ng.mid(4, ng.size()-1);
    ngx.L2Normalize();
    ngv.L2Normalize();

    FILE *file = fopen("gradients_xv.txt", "a");
    fprintf(file, "--------------------------------------------------------------------\n");
    IPrinter::printDateTime(file);
    fprintf(file, "T: %f L: %d h:%f Functional: %.20f N1: %d N2: %d M: %d h1: %f h2: %f ht: %f\n", t1, L, h, rf, N1, N2, M, h1, h2, ht);
    fprintf(file, "x: %.8f %.8f %.8f %.8f\n", x[0], x[1], x[2], x[3]);
    fprintf(file, "AGx: %.8f %.8f %.8f %.8f\n", agx[0], agx[1], agx[2], agx[3]);
    fprintf(file, "NGx: %.8f %.8f %.8f %.8f\n", ngx[0], ngx[1], ngx[2], ngx[3]);
    unsigned int part = 10;//(M+1)
    IPrinter::printVector(x,   "v1: ", part, 0*(M+1)+2*L, 0*(M+1)+2*L+M, file);
    IPrinter::printVector(agv, "AG1:", part, 0*(M+1),     0*(M+1)+M,     file);
    IPrinter::printVector(ngv, "NG1:", part, 0*(M+1),     0*(M+1)+M,     file);
    IPrinter::printVector(x,   "v2: ", part, 1*(M+1)+2*L, 1*(M+1)+2*L+M, file);
    IPrinter::printVector(agv, "AG2:", part, 1*(M+1),     1*(M+1)+M,     file);
    IPrinter::printVector(ngv, "NG2:", part, 1*(M+1),     1*(M+1)+M,     file);
    IPrinter::printDateTime(file);
    //    fprintf(file, "U\n");
    //    DoubleCube c;
    //    IHyperbolicEquation2D::calculateU1(c, h1, h2, ht, N1, N2, M, a, a, qamma);
    //    IPrinter::printMatrix(c[c.size()-1], N2, N1, NULL, file);
    fclose(file);
    x.clear();

    return rf;
}

double HyperbolicControl2DM::fx(const DoubleVector &x)
{
    px = &x;
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
    IHyperbolicEquation2D::calculateU1(u, h1, h2, ht, N1, N2, M, a, a, qamma);

    pu = &u;
    DoubleCube p;
    IBackwardHyperbolicEquation2D::calculateU1(p, h1, h2, ht, N1, N2, M, a, a, qamma);

    for (unsigned int k=0; k<=M; k++)
    {
        double psiX1;
        double psiX2;        double m = 1.0;
        if (k==0 || k==M) m *= 0.5;

        psiDerivative(psiX1, psiX2, x[0], x[1], p[k]);
        g[0] = g[0] + m * x[2*L+0*(M+1)+k] * psiX1;
        g[1] = g[1] + m * x[2*L+0*(M+1)+k] * psiX2;

        psiDerivative(psiX1, psiX2, x[2], x[3], p[k]);
        g[2] = g[2] + m * x[2*L+1*(M+1)+k] * psiX1;
        g[3] = g[3] + m * x[2*L+1*(M+1)+k] * psiX2;
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
        g[2*L+0*(M+1)+k] = -p[k][j][i];

        i = (unsigned int)round(x[2]/h1);
        j = (unsigned int)round(x[3]/h2);
        g[2*L+1*(M+1)+k] = -p[k][j][i];
    }
}

void HyperbolicControl2DM::psiDerivative(double &psiX1, double &psiX2, double x1, double x2, const DoubleMatrix &psi)
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

double HyperbolicControl2DM::initial1(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl2DM::initial2(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl2DM::boundary(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
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
    //double t = k*ht;
    const DoubleVector &x = *px;

    double _v1 = x[2*L+0*(M+1)+k];
    double _v2 = x[2*L+1*(M+1)+k];

    sum += _v1 * gause_a * exp(-((x1-x[0])*(x1-x[0]) + (x2-x[1])*(x2-x[1]))/gause_b);
    sum += _v2 * gause_a * exp(-((x1-x[2])*(x1-x[2]) + (x2-x[3])*(x2-x[3]))/gause_b);

    sum += fxt(i, j, k);

    return sum;
}

double HyperbolicControl2DM::binitial1(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu)[M];
    const DoubleMatrix &u1 = (*pu)[M-2];
    return -2.0 * alpha1 * ( (u0[j][i]-u1[j][i])/(2.0*ht) - U1);
}

double HyperbolicControl2DM::binitial2(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu)[M];
    return +2.0 * (u0[j][i] - U0) + qamma*(binitial1(i,j));
}

double HyperbolicControl2DM::bboundary(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2DM::bf(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
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
    printf("\nJ[%d]: %.16f ALPHA: %.8f\n", i, fn->fx(x), alpha);

    DoubleVector ng(x.size());
//    IGradient::Gradient(fn, h, x, ng);
    DoubleVector ngx = ng.mid(0, 3);
    DoubleVector ngv = ng.mid(4, ng.size()-1);
    ngx.L2Normalize();
    ngv.L2Normalize();

    DoubleVector agx = g.mid(0, 3);
    DoubleVector agv = g.mid(4, g.size()-1);
    agx.L2Normalize();
    agv.L2Normalize();

    FILE *file = fopen("gradients_xv.txt", "a");
    fprintf(file, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    IPrinter::printDateTime(file);
    fprintf(file, "T: %f L: %d Functional: %.20f N1: %d N2: %d M: %d h1: %f h2: %f ht: %f\n", t1, L, fn->fx(x), N1, N2, M, h1, h2, ht);
    fprintf(file, "x: %.8f %.8f %.8f %.8f\n", x[0], x[1], x[2], x[3]);
    fprintf(file, "AGx: %.8f %.8f %.8f %.8f\n", agx[0], agx[1], agx[2], agx[3]);
    fprintf(file, "NGx: %.8f %.8f %.8f %.8f\n", ngx[0], ngx[1], ngx[2], ngx[3]);
    IPrinter::printVector(x,   "v1: ", (M+1), 0*(M+1)+2*L, 0*(M+1)+2*L+M, file);
    IPrinter::printVector(agv, "AG1:", (M+1), 0*(M+1),     0*(M+1)+M,     file);
    IPrinter::printVector(ngv, "NG1:", (M+1), 0*(M+1),     0*(M+1)+M,     file);
    IPrinter::printVector(x,   "v2: ", (M+1), 1*(M+1)+2*L, 1*(M+1)+2*L+M, file);
    IPrinter::printVector(agv, "AG2:", (M+1), 1*(M+1),     1*(M+1)+M,     file);
    IPrinter::printVector(ngv, "NG2:", (M+1), 1*(M+1),     1*(M+1)+M,     file);
    IPrinter::printDateTime(file);
    //    fprintf(file, "U\n");
    //    DoubleCube c;
    //    IHyperbolicEquation2D::calculateU1(c, h1, h2, ht, N1, N2, M, a, a, qamma);
    //    IPrinter::printMatrix(c[c.size()-1], N2, N1, NULL, file);
    fprintf(file, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    fclose(file);
}

void HyperbolicControl2DM::project(DoubleVector &x, int i)
{
    if (i==0)
    {
        if (x[i] <= 0.0) { x[i] = 0.0; }
        if (x[i] >= 0.5) { x[i] = 0.5; }

    }
    if (i==1)
    {
        if (x[i] <= 0.0) { x[i] = 0.0; }
        if (x[i] >= 0.5) { x[i] = 0.5; }
    }
    if (i==2)
    {
        if (x[i] <= 0.5) { x[i] = 0.5; }
        if (x[i] >= 1.0) { x[i] = 1.0; }
    }
    if (i==3)
    {
        if (x[i] <= 0.5) { x[i] = 0.5; }
        if (x[i] >= 1.0) { x[i] = 1.0; }
    }
    if (i>3)
    {
        if (x[i] < vd) x[i] = vd;
        if (x[i] > vu) x[i] = vu;
    }
}
