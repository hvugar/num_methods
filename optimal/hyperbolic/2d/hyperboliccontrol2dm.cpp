#include "hyperboliccontrol2dm.h"

void HyperbolicControl2DM::Main(int argc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    HyperbolicControl2DM hc;
    hc.fx(1.0);
}

HyperbolicControl2DM::HyperbolicControl2DM()
{
    h1 = 0.01;  N1 = 100;
    h2 = 0.01;  N2 = 100;
    ht = 0.005; M  = 200;
    h = 0.001;
    L = 2;

    e << 0.2 << 0.2;

    alpha0 = 1.0;
    alpha1 = 10.0;
    alpha2 = 2.0;
    alpha3 = 1.0;
    qamma = 0.2;

    U0 = U1 = 0.0;
    a = 1.0;

    // limits of v
    vd = -2.0;
    vu = +2.0;
}

double HyperbolicControl2DM::fx(double t UNUSED_PARAM) const
{
    HyperbolicControl2DM *p = const_cast<HyperbolicControl2DM*>(this);
    p->M  = 200;

    DoubleVector x0( 2*L + (M+1)*L);
    for (unsigned int k=0; k<=M; k++)
    {
        x0[2*L + 0*(M+1)+k] = 2.0;
        x0[2*L + 1*(M+1)+k] = 2.0;
    }
    x0[0] = 0.3;
    x0[1] = 0.4;
    x0[2] = 0.7;
    x0[3] = 0.7;

    double min_step = 1.0;
    double gold_eps = 0.0001;
    ConjugateGradient cg;
    cg.setFunction(p);
    cg.setGradient(p);
    cg.setEpsilon1(0.0001);
    cg.setEpsilon2(0.0001);
    cg.setEpsilon3(0.0001);
    cg.setR1MinimizeEpsilon(min_step, gold_eps);
    cg.setPrinter(p);
    cg.setProjection(p);
    cg.setNormalize(true);
    cg.showEndMessage(true);
    cg.calculate(x0);

    p->px = &x0;
    DoubleCube c;
    IHyperbolicEquation2D::calculateU1(c, h1, h2, ht, N1, N2, M, a, a, qamma);

    double MIN = +10000.0;
    double MAX = -10000.0;
    for (unsigned int i=0; i<=M; i++)
    {
        DoubleMatrix mx = c.matrix(i);

        double min = mx.min();
        double max = mx.max();
        if (MIN > min) MIN = min;
        if (MAX < max) MAX = max;

        char buffer[20];
        int n = 0;
        if (i<10) n = sprintf(buffer, "data/0000000%d.txt", i);
        if (i<100 && i>=10) n = sprintf(buffer, "data/000000%d.txt", i);
        if (i<1000 && i>=100) n = sprintf(buffer, "data/00000%d.txt", i);
        if (i<10000 && i>=1000) n = sprintf(buffer, "data/0000%d.txt", i);
        buffer[n] = '\0';
        FILE *file = fopen(buffer, "w");
        IPrinter::printMatrix(mx, N2, N1, NULL, file);
        fclose(file);

        printf("File: %s min: %.16f max: %.16f\n", buffer, MIN, MAX);
    }

    //-0.3892016224094016 max: 0.3951541064046917

    double rf = fx(x0);

    return rf;
}

double HyperbolicControl2DM::fx(const DoubleVector &x) const
{
    HyperbolicControl2DM *p = const_cast<HyperbolicControl2DM*>(this);
    p->px = &x;
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

//    sum1 = 0.0;
//    sum1 += 0.25 * (u0[0][0]-U0)   * (u0[0][0]-U0);
//    sum1 += 0.25 * (u0[0][N1]-U0)  * (u0[0][N1]-U0);
//    sum1 += 0.25 * (u0[N2][0]-U0)  * (u0[N2][0]-U0);
//    sum1 += 0.25 * (u0[N2][N1]-U0) * (u0[N2][N1]-U0);
//    for (unsigned int i=1; i<=N1-1; i++)
//    {
//        sum1 += 0.5*(u0[0][i]-U0)  * (u0[0][i]-U0);
//        sum1 += 0.5*(u0[N2][i]-U0) * (u0[N2][i]-U0);
//    }
//    for (unsigned int j=1; j<=N2-1; j++)
//    {
//        sum1 += 0.5*(u0[j][0]-U0)  * (u0[j][0]-U0);
//        for (unsigned int i=1; i<=N1-1; i++)
//        {
//            sum1 += (u0[j][i]-U0) * (u0[j][i]-U0);
//        }
//        sum1 += 0.5*(u0[j][N1]-U0) * (u0[j][N1]-U0);
//    }

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
    return sum;
}

double HyperbolicControl2DM::norm(const DoubleVector& v) const
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
        double psiX2;
        double m = 1.0;
        if (k==0 || k==M) m *= 0.5;
        psiDerivative(psiX1, psiX2, x[0], x[1], p.matrix(k));
        g[0] = g[0] + m * x[2*L+0*(M+1)+k] * psiX1;
        g[1] = g[1] + m * x[2*L+0*(M+1)+k] * psiX2;
        psiDerivative(psiX1, psiX2, x[2], x[3], p.matrix(k));
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
        g[2*L+0*(M+1)+k] = -p.at(k,j,i);

        i = (unsigned int)round(x[2]/h1);
        j = (unsigned int)round(x[3]/h2);
        g[2*L+1*(M+1)+k] = -p.at(k,j,i);
    }

    //    FILE *file = fopen("xv.txt", "w");
    //    IPrinter::printVector(14,10,g,NULL,g.size(),0,0,file);
    //    fprintf(file, "--\n");
    ////    IPrinter::printMatrix(15,10,u.matrix(M),10,10,NULL,file);
    ////    fprintf(file, "--\n");
    ////    IPrinter::printMatrix(15,10,u.matrix(M-2),10,10,NULL,file);
    ////    fprintf(file, "--\n");
    //    IPrinter::printMatrix(15,10,p.matrix(0),10,10,NULL,file);
    //    fprintf(file, "--\n");
    //    IPrinter::printMatrix(15,10,p.matrix(M),10,10,NULL,file);
    //    fprintf(file, "--\n");
    //    fclose(file);
}

void HyperbolicControl2DM::psiDerivative(double &psiX1, double &psiX2, double x1, double x2, const DoubleMatrix &psi)
{
    unsigned int i = (unsigned int)round(x1/h1);
    unsigned int j = (unsigned int)round(x2/h2);
    if (i==0)
    {
        psiX1  = (psi[j][i+1] - psi[j][i])/h1;
    }
    else if (i==N1)
    {
        psiX1 = (psi[j][i] - psi[j][i-1])/h1;
    }
    else
    {
        psiX1 = (psi[j][i+1] - psi[j][i-1])/(2.0*h1);
    }

    if (j==0)
    {
        psiX2 = (psi[j+1][i] - psi[j][i])/h2;
    }
    else if (j==N2)
    {
        psiX2 = (psi[j][i] - psi[j-1][i])/h2;
    }
    else
    {
        psiX2 = (psi[j+1][i] - psi[j-1][i])/(2.0*h2);
    }
}

double HyperbolicControl2DM::initial1(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
    return 0.0;
}

double HyperbolicControl2DM::initial2(unsigned int i UNUSED_PARAM, unsigned int j UNUSED_PARAM) const
{
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
    const DoubleMatrix &u0 = (*pu).matrix(M);
    const DoubleMatrix &u1 = (*pu).matrix(M-2);
    return -2.0 * alpha0 * ((u0[j][i]-u1[j][i])/(2.0*ht) - U1);
}

double HyperbolicControl2DM::binitial2(unsigned int i, unsigned int j) const
{
    const DoubleMatrix &u0 = (*pu).matrix(M);
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

void HyperbolicControl2DM::print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double r) const
{
    C_UNUSED(g);

    printf("J[%d]: %.16f\n", i, r);

//    if (i==0)
//    {
//        DoubleVector ag=g;
//        DoubleVector ge  = ag.mid(0,3);       ge.L2Normalize();
//        DoubleVector gv1 = ag.mid(4,M+4);     gv1.L2Normalize();
//        DoubleVector gv2 = ag.mid(M+5,2*M+5); gv2.L2Normalize();
//        printf("a:%14.10f %14.10f %14.10f %14.10f |%14.10f %14.10f %14.10f %14.10f |%14.10f %14.10f %14.10f %14.10f\n",
//               ge[0], ge[1], ge[2], ge[3],
//               gv1[40], gv1[80], gv1[120], gv1[160],
//               gv2[40], gv2[80], gv2[120], gv2[160]);

//        IGradient::Gradient(this, 0.01, x, ag);
//        ge  = ag.mid(0,3);       ge.L2Normalize();
//        gv1 = ag.mid(4,M+4);     gv1.L2Normalize();
//        gv2 = ag.mid(M+5,2*M+5); gv2.L2Normalize();
//        printf("n:%14.10f %14.10f %14.10f %14.10f |%14.10f %14.10f %14.10f %14.10f |%14.10f %14.10f %14.10f %14.10f\n",
//               ge[0], ge[1], ge[2], ge[3],
//               gv1[40], gv1[80], gv1[120], gv1[160],
//               gv2[40], gv2[80], gv2[120], gv2[160]);
//        IPrinter::printSeperatorLine();
//    }
}

void HyperbolicControl2DM::project(DoubleVector &x, int i)
{
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

    if (i>3)
    {
        if (x[i] < vd) x[i] = vd;
        if (x[i] > vu) x[i] = vu;
    }
}
