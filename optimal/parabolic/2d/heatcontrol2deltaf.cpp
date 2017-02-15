#include "heatcontrol2deltaf.h"

void HeatControl2DeltaF::main(int argc, char *argv[])
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    HeatControl2DeltaF hc(100, 100, 100);

    DoubleVector v((hc.M+1)*hc.L);
    for (unsigned int k=0; k<=hc.M; k++)
    {
        v[0*(hc.M+1) + k] = 1.0;//hc.g1(k*hc.ht);
        v[1*(hc.M+1) + k] = 1.0;//hc.g2(k*hc.ht);
        v[2*(hc.M+1) + k] = 1.0;//hc.g3(k*hc.ht);
    }

    /* Minimization */
    SteepestDescentGradient g2;
    g2.setFunction(&hc);
    g2.setGradient(&hc);
    g2.setEpsilon1(0.0001);
    g2.setEpsilon2(0.0001);
    g2.setEpsilon3(0.0001);
    g2.setR1MinimizeEpsilon(1.0, 0.0001);
    g2.setPrinter(&hc);
    g2.setProjection(&hc);
    g2.setNormalize(true);
    g2.calculate(v);

    DoubleVector gr1(v.size());
    IGradient::Gradient(&hc, 0.0001, v, gr1);
    gr1.L2Normalize();

    DoubleVector gr2(v.size());
    hc.gradient(v, gr2);
    gr2.L2Normalize();

    printf("J[%d]: %.16f\n", 0, hc.fx(v));
}

HeatControl2DeltaF::HeatControl2DeltaF(unsigned int m, unsigned int n2, unsigned int n1)
{
    alpha = 1.0;

    this->M  = m;
    this->N2 = n2;
    this->N1 = n1;
    this->L  = 3;

    t0 = 0.0;
    t1 = 1.0;

    x10 = 0.0;
    x11 = 1.0;

    x20 = 0.0;
    x21 = 1.0;

    a1 = a2 = 1.0;

    ht = (t1-t0)/m;
    h1 = (x11-x10)/n1;
    h2 = (x21-x20)/n2;

    double sgm1 = 3.0*h1;
    double sgm2 = 3.0*h2;
    gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    gause_b = 2.0*sgm1*sgm2;

    //initilize
    E.resize(2*L);
    E[0] = 0.50;
    E[1] = 0.80;
    E[2] = 0.70;
    E[3] = 0.20;
    E[4] = 0.20;
    E[5] = 0.30;

    DoubleVector v((M+1)*L);
    for (unsigned int k=0; k<=M; k++)
    {
        v[0*(m+1)+k] = v1(k*ht);
        v[1*(m+1)+k] = v2(k*ht);
        v[2*(m+1)+k] = v3(k*ht);
    }

    pv = &v;
    IParabolicEquation2D::calculateMVD(U, h1, h2, ht, N1, N2, M, a1, a2);

    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
    IPrinter::printMatrix(U, 10, 10);
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");

    FILE* f = fopen("heat_optimal_v.txt", "w");
    IPrinter::printMatrix(U, N1, N2, NULL, f);
    fclose(f);
}

double HeatControl2DeltaF::fx(const DoubleVector& v) const
{
    const_cast<HeatControl2DeltaF*>(this)->pv = &v;
    DoubleMatrix u;
    IParabolicEquation2D::calculateMVD(u, h1, h2, ht, N1, N2, M, a1, a2);

    double sum = 0.0;
    for (unsigned int j=0; j<=N2; j++)
    {
        for (unsigned int i=0; i<=N1; i++)
        {
            double betta = 1.0;
            if (i==0 || i==N1 || j==0 || j==N2) betta = 0.5;
            if (i==0  && j==0)  betta = 0.25;
            if (i==0  && j==N2) betta = 0.25;
            if (i==N1 && j==0)  betta = 0.25;
            if (i==N1 && j==N2) betta = 0.25;
            sum += betta*(u[j][i] - U[j][i])*(u[j][i] - U[j][i]);
        }
    }
    sum = (h1*h2)*sum;

    double nrm = 0.0;
    nrm = norm(v);
    return sum + alpha*nrm;
}

double HeatControl2DeltaF::norm(const DoubleVector& v) const
{
    double nrm = 0.0;
    for (unsigned int k=0; k<=M; k++)
    {
        double betta = 1.0;
        if (k==0 || k==M) betta = 0.5;
        nrm += betta*(v[0*(M+1)+k] - v1(k*ht))*(v[0*(M+1)+k] - v1(k*ht));
        nrm += betta*(v[1*(M+1)+k] - v2(k*ht))*(v[1*(M+1)+k] - v2(k*ht));
        nrm += betta*(v[2*(M+1)+k] - v3(k*ht))*(v[2*(M+1)+k] - v3(k*ht));
    }
    nrm = ht * nrm;
    return nrm;
}

void HeatControl2DeltaF::gradient(const DoubleVector& v, DoubleVector& g)
{
    pv = &v;
    DoubleMatrix u;
    IParabolicEquation2D::calculateMVD(u, h1, h2, ht, N1, N2, M, a1, a2);

    pu = &u;
    DoubleCube psi;
    IBackwardParabolicEquation2D::calculateMVD(psi, h1, h2, ht, N1, N2, M, a1, a2);

    //    for (unsigned int i=0; i<g.size(); i++) g[i] = 0.0;
    //    for (unsigned int k=M; k!=(unsigned int)0-1; k--)
    //    {
    //        calculateGF(v, psi[k], g, k);
    //    }

    for (unsigned int k=0; k<=M; k++)
    {
        unsigned int i1 = (unsigned int)round(E[0]/h1);
        unsigned int j1 = (unsigned int)round(E[1]/h2);
        g[0*(M+1)+k] = -psi.at(k,j1,i1) + 2.0*alpha*(v[0*(M+1)+k] - v1(k*ht));

        unsigned int i2 = (unsigned int)round(E[2]/h1);
        unsigned int j2 = (unsigned int)round(E[3]/h2);
        g[1*(M+1)+k] = -psi.at(k,j2,i2) + 2.0*alpha*(v[1*(M+1)+k] - v2(k*ht));

        unsigned int i3 = (unsigned int)round(E[4]/h1);
        unsigned int j3 = (unsigned int)round(E[5]/h2);
        g[2*(M+1)+k] = -psi.at(k,j3,i3) + 2.0*alpha*(v[2*(M+1)+k] - v3(k*ht));
    }

    psi.clear();
    //    IGradient::Gradient(this, 0.0001, x, g);
}

void HeatControl2DeltaF::psiDerivative(double &psiX1, double &psiX2, double x1, double x2, const DoubleMatrix &psi)
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

void HeatControl2DeltaF::calculateGF(const DoubleVector &v, const DoubleMatrix& psi, DoubleVector& g, unsigned int k)
{
    C_UNUSED(k);
    for (unsigned k=0; k<=M; k++)
    {
        if (alpha < 1.0)
        {
            g[2*L+0*(M+1)+k] = 0.0;
            g[2*L+1*(M+1)+k] = 0.0;
            g[2*L+2*(M+1)+k] = 0.0;
        }
        else
        {
            double p[L];
            p[0]=p[1]=p[2]=0.0;
            for (unsigned int j=0; j<=N2; j++)
            {
                for (unsigned int i=0; i<=N1; i++)
                {
                    if (fabs(E[0]-i*h1)<h1 && fabs(E[1]-j*h2)<h2) p[0] += psi[j][i]*((h1-fabs(E[0]-i*h1))/h1)*((h2-fabs(E[1]-j*h2))/h2);
                    if (fabs(E[2]-i*h1)<h1 && fabs(E[3]-j*h2)<h2) p[1] += psi[j][i]*((h1-fabs(E[2]-i*h1))/h1)*((h2-fabs(E[3]-j*h2))/h2);
                    if (fabs(E[4]-i*h1)<h1 && fabs(E[5]-j*h2)<h2) p[2] += psi[j][i]*((h1-fabs(E[4]-i*h1))/h1)*((h2-fabs(E[5]-j*h2))/h2);
                }
            }
            g[2*L+0*(M+1)+k] = -p[0] + 2.0*alpha*(v[0*(M+1)+k] - v1(k*ht));
            g[2*L+1*(M+1)+k] = -p[1] + 2.0*alpha*(v[1*(M+1)+k] - v2(k*ht));
            g[2*L+2*(M+1)+k] = -p[2] + 2.0*alpha*(v[2*(M+1)+k] - v3(k*ht));
        }
    }
}

double HeatControl2DeltaF::initial(unsigned int i, unsigned int j) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    return u(x1, x2, t0);
}

double HeatControl2DeltaF::boundary(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t  = 0.5*k*ht;
    return u(x1, x2, t);
}

double HeatControl2DeltaF::f(unsigned int i, unsigned int j, unsigned int k) const
{
    //double x1 = i*h1;
    //double x2 = j*h2;
    //double t  = 0.5*k*ht;

    double sum = 0.0;
    unsigned int k1 = (k%2==0 ? k/2 : (k+1)/2);

    //printf("k: %d k1: %d j: %d i: %d\n", k, k1, j, i);

    double _v1 = (*pv)[0*(M+1)+k1];
    double _v2 = (*pv)[1*(M+1)+k1];
    double _v3 = (*pv)[2*(M+1)+k1];

    unsigned int i1 = (unsigned int)round(E[0]/h1);
    unsigned int j1 = (unsigned int)round(E[1]/h2);
    if (i==i1 && j==j1) sum += (1.0/(h1*h2)) * _v1;

    unsigned int i2 = (unsigned int)round(E[2]/h1);
    unsigned int j2 = (unsigned int)round(E[3]/h2);
    if (i==i2 && j==j2) sum += (1.0/(h1*h2)) * _v2;

    unsigned int i3 = (unsigned int)round(E[4]/h1);
    unsigned int j3 = (unsigned int)round(E[5]/h2);
    if (i==i3 && j==j3) sum += (1.0/(h1*h2)) * _v3;

    //    double sgm1 = 3.0*h1;
    //    double sgm2 = 3.0*h2;
    //    double gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    //    double gause_b = 2.0*sgm1*sgm2;

    //    sum += _v1 * gause_a * exp(-((x1-E[0])*(x1-E[0]) + (x2-E[1])*(x2-E[1]))/gause_b);
    //    sum += _v2 * gause_a * exp(-((x1-E[2])*(x1-E[2]) + (x2-E[3])*(x2-E[3]))/gause_b);
    //    sum += _v3 * gause_a * exp(-((x1-E[4])*(x1-E[4]) + (x2-E[5])*(x2-E[5]))/gause_b);

    //    sum += g1(t) * gause_a * exp(-((x1-e[0])*(x1-e[0]) + (x2-e[1])*(x2-e[1]))/gause_b);// * h1*h2;
    //    sum += g2(t) * gause_a * exp(-((x1-e[2])*(x1-e[2]) + (x2-e[3])*(x2-e[3]))/gause_b);// * h1*h2;
    //    sum += g3(t) * gause_a * exp(-((x1-e[4])*(x1-e[4]) + (x2-e[5])*(x2-e[5]))/gause_b);// * h1*h2;

    //    DoubleVector e(2*L);
    //    for (unsigned int l=0; l<L; l++)
    //    {
    //        e[2*l+0] = (*px)[2*l+0];
    //        e[2*l+1] = (*px)[2*l+1];
    //    }

    //    if (fabs(x1-e[0])<=h1 && fabs(x2-e[1])<=h2)
    //    {
    //        sum += g1(t) * ((h1-fabs(x1-e[0]))/(h1*h1))*((h2-fabs(x2-e[1]))/(h2*h2));
    //    }
    //    if (fabs(x1-e[2])<=h1 && fabs(x2-e[3])<=h2)
    //    {
    //        sum += g2(t) * ((h1-fabs(x1-e[2]))/(h1*h1))*((h2-fabs(x2-e[3]))/(h2*h2));
    //    }
    //    if (fabs(x1-e[4])<=h1 && fabs(x2-e[5])<=h2)
    //    {
    //        sum += g3(t) * ((h1-fabs(x1-e[4]))/(h1*h1))*((h2-fabs  v c  v (x2-e[5]))/(h2*h2));
    //    }

    return sum;
}

double HeatControl2DeltaF::binitial(unsigned int i, unsigned int j) const
{
    return -2.0*((*pu)[j][i] - U[j][i]);
}

double HeatControl2DeltaF::bboundary(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

double HeatControl2DeltaF::bf(unsigned int i, unsigned int j, unsigned int k) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    C_UNUSED(k);
    return 0.0;
}

void HeatControl2DeltaF::print(unsigned int i, const DoubleVector& x, const DoubleVector &g, double fx) const
{
    C_UNUSED(g);
    C_UNUSED(alpha);
    HeatControl2DeltaF *hc = const_cast<HeatControl2DeltaF*>(this);
    printf("J[%d]: %.16f\n", i, hc->fx(x));
    IPrinter::printVector(x, "v1:", 10, 0*(M+1), 0*(M+1) + M);
    IPrinter::printVector(x, "v2:", 10, 1*(M+1), 1*(M+1) + M);
    IPrinter::printVector(x, "v3:", 10, 2*(M+1), 2*(M+1) + M);
}

void HeatControl2DeltaF::project(DoubleVector &e, int index)
{
    if (index<6)
    {
        if (e[index] > 1.0) e[index] = 1.0;
        if (e[index] < 0.0) e[index] = 0.0;
    }
}
