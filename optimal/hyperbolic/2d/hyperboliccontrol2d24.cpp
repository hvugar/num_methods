#include "hyperboliccontrol2d24.h"

void HyperbolicControl2D24::main(int argc, char ** argv)
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    if (argc < 2)
    {
        puts("usage: -f filename -x1 x1_1 x_12 -x2 x_21 x_22 -v1 v_1 -v2 v_2 -T T");
        return;
    }

    HyperbolicControl2D24 hc;

    hc.file = stdout;
    if (strcmp(argv[1], "-f") == 0)hc.file = fopen(argv[2], "w");

    hc.x0[0] = atof(argv[4]);
    hc.x0[1] = atof(argv[5]);

    hc.x0[2] = atof(argv[7]);
    hc.x0[3] = atof(argv[8]);

    hc.v0[0] = atof(argv[10]);
    hc.v0[1] = atof(argv[12]);

    double T = atof(argv[14]);
    hc.fx(T);
    fclose(hc.file);
}

HyperbolicControl2D24::HyperbolicControl2D24()
{
    x10 = x20 = t0 = 0.0;
    x11 = x21 = t1 = 1.0;
    a1  = a2  = 1.0;

    // zerbe noqtesi
    e.resize(2);
    e[0] = 0.2;
    e[1] = 0.2;

    alpha0 = 1.0;

    // parametrler
    alpha1 = 10.0;
    alpha2 = 2.0;
    alpha3 = 1.0;
    qamma  = 0.2;

    h1 = h2 = 0.01;
    ht = 0.0025;
    N1 = (unsigned)round((x11 - x10)/h1);
    N2 = (unsigned)round((x21 - x20)/h2);
    M  = (unsigned)round((t1 - t0)/ht);

    L = 2;
    x0.resize(2*L);
    v0.resize(L);

    U0 = 0.0;
    U1 = 0.0;

    pw = NULL;
    pu = NULL;
}

double HyperbolicControl2D24::fx(double t) const
{
    HyperbolicControl2D24 *p = const_cast<HyperbolicControl2D24*>(this);
    p->t1 = t;
    p->M  = (unsigned)ceil((t1 - t0)/ht);

    printf("T: %f M: %d\n--------------------\n", t1, M);

    DoubleVector w(2*L + (M+1)*L, 0.0);
    for (unsigned int k=0; k<=M; k++)
    {
        w[2*L+0*(M+1)+k] = this->v0[0];
        w[2*L+1*(M+1)+k] = this->v0[1];
    }
    w[0] = this->x0[0];
    w[1] = this->x0[1];
    w[2] = this->x0[2];
    w[3] = this->x0[3];

    //    {
    //        DoubleVector ga1(x0.size());
    //        gradient(x0, ga1);
    //        DoubleVector ga = ga1.mid(0, 3);
    //        printf("%f %f %f %f\n", ga[0], ga[1], ga[2], ga[3]);
    //        ga.L2Normalize();
    //        printf("%f %f %f %f\n", ga[0], ga[1], ga[2], ga[3]);
    //        DoubleVector ga_v1 = ga1.mid(4+0*(M+1), 4+0*(M+1)+M);
    //        DoubleVector ga_v2 = ga1.mid(4+1*(M+1), 4+1*(M+1)+M);
    //        IPrinter::printVector(ga_v1, "ga_v1");
    //        IPrinter::printVector(ga_v2, "ga_v2");
    //        ga_v1.L2Normalize();
    //        ga_v2.L2Normalize();
    //        IPrinter::printVector(ga_v1, "ga_v1");
    //        IPrinter::printVector(ga_v2, "ga_v2");
    //    }

    //    {
    //        DoubleVector gn1(x0.size());
    //        IGradient::Gradient(this, 0.0001, x0, gn1);
    //        DoubleVector gn = gn1.mid(0, 3);
    //        printf("%f %f %f %f\n", gn[0], gn[1], gn[2], gn[3]);
    //        gn.L2Normalize();
    //        printf("%f %f %f %f\n", gn[0], gn[1], gn[2], gn[3]);
    //        DoubleVector gn_v1 = gn1.mid(4+0*(M+1), 4+0*(M+1)+M);
    //        DoubleVector gn_v2 = gn1.mid(4+1*(M+1), 4+1*(M+1)+M);
    //        IPrinter::printVector(gn_v1, "gn_v1");
    //        IPrinter::printVector(gn_v2, "gn_v2");
    //        gn_v1.L2Normalize();
    //        gn_v2.L2Normalize();
    //        IPrinter::printVector(gn_v1, "gn_v1");
    //        IPrinter::printVector(gn_v2, "gn_v2");

    //        return 0.0;
    //    }

    printGradients(w, 0, 0.0, file);

    double min_step = 1.0;
    double gold_eps = 0.001;
    ConjugateGradient cg;
    cg.setFunction(p);
    cg.setGradient(p);
    cg.setEpsilon1(0.001);
    cg.setEpsilon2(0.001);
    cg.setEpsilon3(0.001);
    cg.setR1MinimizeEpsilon(min_step, gold_eps);
    cg.setPrinter(p);
    cg.setProjection(p);
    cg.setNormalize(p);
    //cg.showEndMessage(false);
    cg.calculate(w);

    printGradients(w, 0, 0.0, file);

    double rf = fx(w);
    fprintf(file, "%f %.16f\n", t, rf);
    fprintf(file, "%f %f %f %f\n", w[0], w[1], w[2], w[3]);
    IPrinter::printVector(w, "v1", M+1, 2*L,       2*L+M,       file);
    IPrinter::printVector(w, "v2", M+1, 2*L+(M+1), 2*L+(2*M+1), file);
    fflush(file);
    return rf;
}

void HyperbolicControl2D24::printGradients(const DoubleVector &x, unsigned int i, double alpha, FILE* f) const
{
    double _fx = const_cast<HyperbolicControl2D24*>(this)->fx(x);
    printf("J[%d]: %.16f %.16f\n", i, _fx, alpha);
    fprintf(f, "J[%d]: %.16f %.16f\n", i, _fx, alpha);

    DoubleVector g(x.size());
    const_cast<HyperbolicControl2D24*>(this)->gradient(x, g);
    DoubleVector eg1 = g.mid(0, 1);
    DoubleVector eg2 = g.mid(2, 3);
    DoubleVector vg1 = g.mid(4,   M+4);
    DoubleVector vg2 = g.mid(M+5, 2*M+5);

    DoubleVector v1 = x.mid(4,   M+4);
    DoubleVector v2 = x.mid(M+5, 2*M+5);

    fprintf(f, "P1: %10.6f %10.6f Gr: %10.6f %10.6f Nr: %10.6f ", x[0], x[1], eg1[0], eg1[1], eg1.L2Norm());
    eg1.L2Normalize();
    fprintf(f, "Gr: %10.6f %10.6f Nr: %10.6f\n", eg1[0], eg1[1], eg1.L2Norm());
    fprintf(f, "P2: %10.6f %10.6f Gr: %10.6f %10.6f Nr: %10.6f ", x[2], x[3], eg2[0], eg2[1], eg2.L2Norm());
    eg2.L2Normalize();
    fprintf(f, "Gr: %10.6f %10.6f Nr: %10.6f\n", eg2[0], eg2[1], eg2.L2Norm());
    fprintf(f, "---\n");

    IPrinter::printVector(v1, "v1:", 10, 0, M, f);
    fprintf(f, "N:  %10.6f\n", vg1.L2Norm());
    IPrinter::printVector(vg1, "vg1:", 10, 0, M, f);
    vg1.L2Normalize();
    IPrinter::printVector(vg1, "ng1:", 10, 0, M, f);
    fprintf(f, "---\n");

    IPrinter::printVector(v2, "v2:", 10, 0, M, f);
    fprintf(f, "N:  %10.6f\n", vg2.L2Norm());
    IPrinter::printVector(vg2, "vg2:", 10, 0, M, f);
    vg2.L2Normalize();
    IPrinter::printVector(vg2, "ng2:", 10, 0, M, f);
    fprintf(f, "--------------------------------------------------------------------------------------------------\n");
    fflush(f);
}

void HyperbolicControl2D24::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double fx) const
{
    C_UNUSED(g);
    printGradients(x, i, 0.0, file);
}

double HyperbolicControl2D24::fx(const DoubleVector &w) const
{
    HyperbolicControl2D24 *p = const_cast<HyperbolicControl2D24*>(this);
    p->pw = &w;
    DoubleCube c;
    IHyperbolicEquation2D::calculateU1(c, h1, h2, ht, N1, N2, M, a1, a2, qamma);

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
            sum2 = sum2 + k * ((u0[j][i]-u1[j][i])/(2.0*ht) - U1) * ((u0[j][i]-u1[j][i])/(2.0*ht) - U1);
        }
    }
    sum2 = h1*h2*sum2;

    sum = sum1 + alpha0*sum2;
    return sum;
}

void HyperbolicControl2D24::gradient(const DoubleVector &x, DoubleVector &g)
{
    pw = &x;
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
        psiDerivative(psiX1, psiX2, x[0], x[1], p.matrix(k));
        double _v1 = x[2*L+0*(M+1)+k];
        g[0] = g[0] + m * _v1 * psiX1;
        g[1] = g[1] + m * _v1 * psiX2;
        psiDerivative(psiX1, psiX2, x[2], x[3], p.matrix(k));
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
        g[2*L+0*(M+1)+k] = -p.at(k,j,i);

        i = (unsigned int)round(x[2]/h1);
        j = (unsigned int)round(x[3]/h2);
        g[2*L+1*(M+1)+k] = -p.at(k,j,i);
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

double HyperbolicControl2D24::initial1(unsigned int, unsigned int) const
{
    return 0.0;
}

double HyperbolicControl2D24::initial2(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl2D24::boundary(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
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
    const DoubleVector x = *pw;

    double _v1 = x[2*L+0*(M+1)+k];
    double _v2 = x[2*L+1*(M+1)+k];

    sum += _v1 * gause_a * exp(-((x1-x[0])*(x1-x[0]) + (x2-x[1])*(x2-x[1]))/gause_b);
    sum += _v2 * gause_a * exp(-((x1-x[2])*(x1-x[2]) + (x2-x[3])*(x2-x[3]))/gause_b);

    sum += fxt(i, j, k);
    return sum;
}

double HyperbolicControl2D24::binitial1(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu).matrix(M);
    const DoubleMatrix &u1 = (*pu).matrix(M-2);
    return -2.0 * alpha0 * ((u0[j][i]-u1[j][i])/(2.0*ht) - U1);
}

double HyperbolicControl2D24::binitial2(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu).matrix(M);
    return +2.0 * (u0[j][i] - U0) + qamma*(binitial1(i,j));
}

double HyperbolicControl2D24::bboundary(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2D24::bf(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
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
//    if (i==0)
//    {
//        if (x[i] <= 0.0 + 5.0*h1) { x[i] = 0.0 + 5.0*h1; }
//        if (x[i] >= 1.0 - 5.0*h1) { x[i] = 1.0 - 5.0*h1; }
//    }
//    if (i==1)
//    {
//        if (x[i] <= 0.0 + 5.0*h2) { x[i] = 0.0 + 5.0*h2; }
//        if (x[i] >= 1.0 - 5.0*h2) { x[i] = 1.0 - 5.0*h2; }
//    }
//    if (i==2)
//    {
//        if (x[i] <= 0.0 + 5.0*h1) { x[i] = 0.0 + 5.0*h1; }
//        if (x[i] >= 1.0 - 5.0*h1) { x[i] = 1.0 - 5.0*h1; }
//    }
//    if (i==3)
//    {
//        if (x[i] <= 0.0 + 5.0*h2) { x[i] = 0.0 + 5.0*h2; }
//        if (x[i] >= 1.0 - 5.0*h2) { x[i] = 1.0 - 5.0*h2; }
//    }

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

    if (i>3)
    {
        if (x[i] < -2.0) x[i] = -2.0;
        if (x[i] > +2.0) x[i] = +2.0;
    }
}
