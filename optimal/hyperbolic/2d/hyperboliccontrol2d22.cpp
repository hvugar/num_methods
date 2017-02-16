#include "hyperboliccontrol2d22.h"

void HyperbolicControl2D22::Main(int argc, char* argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    HyperbolicControl2D22 hc;
    hc.file = stdout;
//    hc.file = fopen("20160322.txt", "w");
//    for (double t=0.1; t<=10.1; t+=0.1)
//    {
//        hc.fx(t);
//    }
//    fclose(hc.file);
    hc.fx(1.0);
}

HyperbolicControl2D22::HyperbolicControl2D22()
{
    h1 = 0.01;
    h2 = 0.01;
    ht = 0.005;

    N1 = 100;
    N2 = 100;
    M  = 200;
    L  = 2;

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

    x.resize(2*L);
    x[0] = 0.3;
    x[1] = 0.4;
    x[2] = 0.7;
    x[3] = 0.7;

    v0.resize((M+1)*L, 0.0);

    pv = NULL;
    pu = NULL;
}

double HyperbolicControl2D22::fx(double T) const
{
    HyperbolicControl2D22 *p = const_cast<HyperbolicControl2D22*>(this);
    p->M  = (unsigned int)(T*200);
    unsigned int m = M;

    printf("t: %f M: %d m %d\n", T, M, m);

    p->v0.resize((M+1)*L);

    for (unsigned int k=0; k<=m; k++)
    {
        p->v0[0*(M+1)+k] = 2.0;
        p->v0[1*(M+1)+k] = 2.0;
    }

    for (unsigned int k=m+1; k<=M; k++)
    {
        p->v0[0*(M+1)+k] = v0[0*(M+1)+m];
        p->v0[1*(M+1)+k] = v0[1*(M+1)+m];
    }

    double min_step = 10.0;
    double gold_eps = 0.001;
    ConjugateGradient cg;
    cg.setFunction(p);
    cg.setGradient(p);
    cg.setEpsilon1(0.001);
    cg.setEpsilon2(0.001);
    cg.setEpsilon3(0.001);
    cg.setR1MinimizeEpsilon(min_step, gold_eps);
    cg.setPrinter(p);
    //cg.setProjection(p);
    cg.setNormalize(p);
    //cg.showEndMessage(false);
    cg.calculate(p->v0);

    double rf = fx(v0);
    fprintf(file, "%f %.16f\n", T, rf);
    IPrinter::printVector(v0, "v21", v0.size()/2, 0, v0.size()/2-1, file);
    IPrinter::printVector(v0, "v22", v0.size()/2, v0.size()/2, v0.size()-1, file);
    fflush(file);
    return rf;
}

void HyperbolicControl2D22::print(unsigned int i, const DoubleVector &v, const DoubleVector &g, double fx) const
{
    C_UNUSED(g);
    fprintf(file, "J[%d]: %.16f\n", i, const_cast<HyperbolicControl2D22*>(this)->fx(v));
    //IPrinter::printVector(v, "v1 ", v.size()/2, 0, v.size()/2-1, file);
    //IPrinter::printVector(v, "v2 ", v.size()/2, v.size()/2, v.size()-1, file);
    fflush(file);
    printf("J[%d]: %.16f\n", i, const_cast<HyperbolicControl2D22*>(this)->fx(v));
}

double HyperbolicControl2D22::fx(const DoubleVector &v) const
{
    HyperbolicControl2D22 *p = const_cast<HyperbolicControl2D22*>(this);
    p->pv = &v;
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

void HyperbolicControl2D22::gradient(const DoubleVector &v, DoubleVector &g)
{
    pv = &v;
    DoubleCube u;
    IHyperbolicEquation2D::calculateU1(u, h1, h2, ht, N1, N2, M, a1, a2, qamma);

    pu = &u;
    DoubleCube p;
    IBackwardHyperbolicEquation2D::calculateU1(p, h1, h2, ht, N1, N2, M, a1, a2, qamma);

    unsigned int i,j;
    for (unsigned int k=0; k<=M; k++)
    {
        i = (unsigned int)round(x[0]/h1);
        j = (unsigned int)round(x[1]/h2);
        g[0*(M+1)+k] = -p.at(k,j,i);

        i = (unsigned int)round(x[2]/h1);
        j = (unsigned int)round(x[3]/h2);
        g[1*(M+1)+k] = -p.at(k,j,i);
    }
}

double HyperbolicControl2D22::initial1(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl2D22::initial2(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl2D22::boundary(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
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

double HyperbolicControl2D22::binitial1(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu).matrix(M);
    const DoubleMatrix &u1 = (*pu).matrix(M-2);
    return -2.0 * alpha1 * ((u0[j][i]-u1[j][i])/(2.0*ht) - U1);
}

double HyperbolicControl2D22::binitial2(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu).matrix(M);
    return +2.0 * (u0[j][i] - U0) + qamma*(binitial1(i,j));
}

double HyperbolicControl2D22::bboundary(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
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
