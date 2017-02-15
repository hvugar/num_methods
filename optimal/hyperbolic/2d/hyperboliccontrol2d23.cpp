#include "hyperboliccontrol2d23.h"

void HyperbolicControl2D23::main(int argc, char ** argv)
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    HyperbolicControl2D23 hc;
    hc.file = fopen("20160326.txt", "w");
    //hc.file = stdout;
    for (double t=0.1; t<=10.1; t+=0.1)
    {
        hc.fx(t);
    }
    fclose(hc.file);
}

HyperbolicControl2D23::HyperbolicControl2D23()
{
    h1 = 0.01;
    h2 = 0.01;
    ht = 0.005;

    N1 = 100;
    N2 = 100;
    M  = 2;
    D = 10;
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

    x.resize(2*L);
    x[0] = 0.3;
    x[1] = 0.4;
    x[2] = 0.7;
    x[3] = 0.7;

    v0.resize((M+D+1)*L, 0.0);

    pv = NULL;
    pu = NULL;
}

double HyperbolicControl2D23::fx(double t)
{
    DoubleVector vt = v0;
    unsigned int m = M+D;

    t1 = t;
    M  = (unsigned)round((t1 - t0)/ht);

    //size = v0.size();
    //fprintf(file, "size2: %d\n", size);
    //IPrinter::printVector(v0, "v01", size/2, 0, size/2-1, file);
    //IPrinter::printVector(v0, "v02", size/2, size/2, size-1, file);

    v0.clear();
    v0.resize((M+D+1)*L, 0.0);
    unsigned size = v0.size();
    printf("t: %f M: %d M+D: %d m %d size: %u\n", t1, M, M+D, m, size);

    //for (unsigned int k=0; k<=m; k++)
    //{
    //    v0[0*(M+D+1)+k] = vt[0*(m+1)+k];
    //    v0[1*(M+D+1)+k] = vt[1*(m+1)+k];
    //}
    //vt.clear();

    //for (unsigned int k=m+1; k<=M+D; k++)
    //{
    //    v0[0*(M+D+1)+k] = v0[0*(M+D+1)+m];
    //    v0[1*(M+D+1)+k] = v0[1*(M+D+1)+m];
    //}

//    size = v0.size();
//    fprintf(file, "size2: %d\n", size);
//    IPrinter::printVector(v0, "v11", size/2, 0, size/2-1, file);
//    IPrinter::printVector(v0, "v12", size/2, size/2, size-1, file);

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
    //cg.setProjection(this);
    cg.setNormalize(true);
    //cg.showEndMessage(false);
    cg.calculate(v0);

    double rf = fx(v0);
    fprintf(file, "%f %.16f\n", t, rf);
//    unsigned int size = v0.size();
//    fprintf(file, "size2: %d\n", size);
//    IPrinter::printVector(v0, "v21", size/2, 0, size/2-1, file);
//    IPrinter::printVector(v0, "v22", size/2, size/2, size-1, file);
//    fputs("---------------\n", file);
    fflush(file);
    return rf;
}

void HyperbolicControl2D23::print(unsigned int i, const DoubleVector &v, const DoubleVector &g, double fx) const
{
    C_UNUSED(g);
    fprintf(file, "J[%d]: %.16f\n", i, const_cast<HyperbolicControl2D23*>(this)->fx(v));
    //IPrinter::printVector(v, "v1 ", v.size()/2, 0, v.size()/2-1, file);
    //IPrinter::printVector(v, "v2 ", v.size()/2, v.size()/2, v.size()-1, file);
    fflush(file);
    printf("J[%d]: %.16f\n", i, const_cast<HyperbolicControl2D23*>(this)->fx(v));
}

double HyperbolicControl2D23::fx(const DoubleVector &v)
{
    pv = &v;
    DoubleCube c;
    IHyperbolicEquation2D::calculateU1(c, h1, h2, ht, N1, N2, M+D, a1, a2, qamma);

    //DoubleMatrix &u0 = c[M];
    //DoubleMatrix &u1 = c[M-2];
    double sum = 0.0;

    double sum1 = 0.0;
    for (unsigned int k=M; k<=M+D; k++)
    {
        for (unsigned int j=0; j<=N2; j++)
        {
            for (unsigned int i=0; i<=N1; i++)
            {
                double k1 = 1.0;
                if (k==M || k==M+D) k1 *= 0.5;
                if (i==0 || i==N1)  k1 *= 0.5;
                if (j==0 || j==N2)  k1 *= 0.5;
                sum1 = sum1 + k1 * (c.at(k,j,i)-U0) * (c.at(k,j,i)-U0);
            }
        }
    }
    sum1 = ht*h1*h2*sum1;

    double sum2 = 0.0;
    //    for (unsigned int j=0; j<=N2; j++)
    //    {
    //        for (unsigned int i=0; i<=N1; i++)
    //        {
    //            double k = 1.0;
    //            if (i==0 || i==N1) k *= 0.5;
    //            if (j==0 || j==N2) k *= 0.5;
    //            sum2 = sum2 + k * ((u0[j][i]-u1[j][i])/(2.0*ht) - U1) * ((u0[j][i]-u1[j][i])/(2.0*ht) - U1);
    //        }
    //    }
    //    sum2 = h1*h2*sum2;

    sum = sum1 + alpha0*sum2;
    return sum;
}

void HyperbolicControl2D23::gradient(const DoubleVector &v, DoubleVector &g)
{
    pv = &v;
    DoubleCube u;
    IHyperbolicEquation2D::calculateU1(u, h1, h2, ht, N1, N2, M+D, a1, a2, qamma);

    pu = &u;
    DoubleCube p;
    IBackwardHyperbolicEquation2D::calculateU1(p, h1, h2, ht, N1, N2, M+D, a1, a2, qamma);

    unsigned int i,j;
    for (unsigned int k=0; k<=M+D; k++)
    {
        i = (unsigned int)round(x[0]/h1);
        j = (unsigned int)round(x[1]/h2);
        g[0*(M+1)+k] = -p.at(k,j,i);

        i = (unsigned int)round(x[2]/h1);
        j = (unsigned int)round(x[3]/h2);
        g[1*(M+1)+k] = -p.at(k,j,i);
    }

//    unsigned int size = g.size();
//    size = v.size();
//    fprintf(file, "v: %d\n", size);
//    IPrinter::printVector(v, "v11", size/2, 0, size/2-1, file);
//    IPrinter::printVector(v, "v12", size/2, size/2, size-1, file);

//    fprintf(file, "g: %d\n", size);
//    IPrinter::printVector(g, "g1", size/2, 0, size/2-1, file);
//    IPrinter::printVector(g, "g2", size/2, size/2, size-1, file);
}

double HyperbolicControl2D23::initial1(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl2D23::initial2(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl2D23::boundary(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2D23::f(unsigned int i, unsigned int j, unsigned int k) const
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

    double _v1 = v[0*(M+D+1)+k];
    double _v2 = v[1*(M+D+1)+k];

    sum += _v1 * gause_a * exp(-((x1-x[0])*(x1-x[0]) + (x2-x[1])*(x2-x[1]))/gause_b);
    sum += _v2 * gause_a * exp(-((x1-x[2])*(x1-x[2]) + (x2-x[3])*(x2-x[3]))/gause_b);

    sum += fxt(i, j, k);

    return sum;
}

double HyperbolicControl2D23::binitial1(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double HyperbolicControl2D23::binitial2(unsigned int i, unsigned int j) const
{
    return qamma * binitial1(i,j);
}

double HyperbolicControl2D23::bboundary(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HyperbolicControl2D23::bf(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j );
    C_UNUSED(k);
    if (k>=M)
    {
        //const DoubleMatrix &u0 = (*pu)[M+D];
        return -2.0 * ((*pu).at(k,j,i) - U0);
    }
    return 0.0;
}

double HyperbolicControl2D23::fxt(unsigned int i, unsigned int j, unsigned int k) const
{
    double sum = 0.0;
    double x1 = i*h1;
    double x2 = j*h2;
    double t = k*ht;
    sum += alpha1*exp(-alpha2*((x1-e[0])*(x1-e[0])+(x2-e[1])*(x2-e[1]))-alpha3*t);
    return sum;
}
