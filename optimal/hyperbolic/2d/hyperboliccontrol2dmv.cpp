#include "hyperboliccontrol2dmv.h"

void HyperbolicControl2DMV::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    HyperbolicControl2DMV hc;
    hc.fx(1.0);
}

double HyperbolicControl2DMV::fx(double t UNUSED_PARAM) const
{
    HyperbolicControl2DMV *p = const_cast<HyperbolicControl2DMV*>(this);
    p->h1 = 0.01;
    p->h2 = 0.01;
    p->ht = 0.005;

    p->N1 = 100;
    p->N2 = 100;
    p->M  = 200;
    p->L = 2;

    p->e.resize(2);
    p->e[0] = 0.2;
    p->e[1] = 0.2;

    p->alpha0 = 1.0;//4.0;
    p->alpha1 = 10.0;//5.0;
    p->alpha2 = 2.0;//6.0;
    p->alpha3 = 1.0;//4.0;
    p->qamma = 0.2;

    p->U0 = 0.0;
    p->U1 = 0.0;
    p->a = 1.0;

    p->d.resize(2*L);
    p->d[0] = 0.3;
    p->d[1] = 0.4;
    p->d[2] = 0.7;
    p->d[3] = 0.7;

    // limits of v
    p->vd = -2.0;
    p->vu = +2.0;

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

    p->gradient(v, ag);
    //ag.L2Normalize();
    IPrinter::printVector(14,10, ag,NULL, ag.size(),0,0,file);

    DoubleVector agv1 = ag.mid(0,M);
    DoubleVector agv2 = ag.mid(M+1,2*M+1);
    agv1.L2Normalize();
    agv2.L2Normalize();
    //    ag.L2Normalize();
    //IGradient::Gradient(this, h, v, ng);
    ng.L2Normalize();

    fprintf(file, "T: %f L: %d h:%f Functional: %.20f N1: %d N2: %d M: %d h1: %f h2: %f ht: %f\n", t, L, h, rf, N1, N2, M, h1, h2, ht);
    IPrinter::printVector(v,  "v1: ", 10, 0*(M+1), 0*(M+1)+M, file);
    IPrinter::printVector(agv1, "AG1:", 10, 0, M, file);
    IPrinter::printVector(ng, "NG1:", 10, 0*(M+1), 0*(M+1)+M, file);

    IPrinter::printVector(v,  "v2: ", 10, 1*(M+1), 1*(M+1)+M, file);
    IPrinter::printVector(agv2, "AG2:", 10, 0, M, file);
    IPrinter::printVector(ng, "NG2:", 10, 1*(M+1), 1*(M+1)+M, file);
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

double HyperbolicControl2DMV::fx(const DoubleVector &v) const
{
    HyperbolicControl2DMV *p = const_cast<HyperbolicControl2DMV*>(this);
    p->pv = &v;
    DoubleCube c;
    IHyperbolicEquation2D::calculateU1(c, h1, h2, ht, N1, N2, M, a, a, qamma);

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
            sum2 = sum2 + k * ( (u0[j][i]-u1[j][i])/(2.0*ht)-U1)*((u0[j][i]-u1[j][i])/(2.0*ht)-U1);
        }
    }
    sum2 = h1*h2*sum2;

    sum = sum1 + alpha0*sum2;
    return sum;// + norm(v);
}

double HyperbolicControl2DMV::norm(const DoubleVector& v) const
{
    C_UNUSED(v);
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
    printf("%f %f %f %d %d %d %f %f %f\n", h1, h2, ht, N1, N2, M, a, a, qamma);
    HyperbolicControl2DMV *pm = const_cast<HyperbolicControl2DMV*>(this);
    pm->pv = &v;
    DoubleCube u;
    IHyperbolicEquation2D::calculateU1(u, h1, h2, ht, N1, N2, M, a, a, qamma);

    pm->pu = &u;
    DoubleCube p;
    IBackwardHyperbolicEquation2D::calculateU1(p, h1, h2, ht, N1, N2, M, a, a, qamma);

    for (unsigned int k=0; k<=M; k++)
    {
        unsigned int i,j;
        i = (unsigned int)round(d[0]/h1);
        j = (unsigned int)round(d[1]/h2);
        g[0*(M+1)+k] = -p.at(k,j,i);

        i = (unsigned int)round(d[2]/h1);
        j = (unsigned int)round(d[3]/h2);
        g[1*(M+1)+k] = -p.at(k,j,i);
    }

    FILE *file = fopen("xv.txt", "a");
    IPrinter::printVector(14,10,g,NULL,g.size(),0,0,file);
    fprintf(file, "--\n");
    //    IPrinter::printMatrix(15,10,u.matrix(M),10,10,NULL,file);
    //    fprintf(file, "--\n");
    //    IPrinter::printMatrix(15,10,u.matrix(M-2),10,10,NULL,file);
    //    fprintf(file, "--\n");
    IPrinter::printMatrix(15,10,p.matrix(0),10,10,NULL,file);
    fprintf(file, "--\n");
    IPrinter::printMatrix(15,10,p.matrix(M),10,10,NULL,file);
    fprintf(file, "--\n");
    fclose(file);


}

double HyperbolicControl2DMV::initial1(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    return 0.0;
}

double HyperbolicControl2DMV::initial2(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    return 0.0;
}

double HyperbolicControl2DMV::boundary(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
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

double HyperbolicControl2DMV::binitial1(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu).matrix(M);
    const DoubleMatrix &u1 = (*pu).matrix(M-2);
    return -2.0 * alpha0 * ((u0[j][i]-u1[j][i])/(2.0*ht) - U1);
}

double HyperbolicControl2DMV::binitial2(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu).matrix(M);
    return +2.0 * (u0[j][i] - U0) + qamma*(binitial1(i,j));
}

double HyperbolicControl2DMV::bboundary(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
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

void HyperbolicControl2DMV::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double fx) const
{
    C_UNUSED(g);
    printf("J[%d]: %.16f\n", i, const_cast<HyperbolicControl2DMV*>(this)->fx(x));
}

void HyperbolicControl2DMV::project(DoubleVector &x, int i)
{
    if (x[i] < vd) x[i] = vd;
    if (x[i] > vu) x[i] = vu;
}
