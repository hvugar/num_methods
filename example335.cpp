#include "example335.h"

void Parabolic1DControl335::main(int argc, char *argv[])
{
    Parabolic1DControl335 hc;

    DoubleVector v((hc.M+1)*hc.L);
    for (unsigned int k=0; k<=hc.M; k++)
    {
        v[0*(hc.M+1) + k] = 1.0;
        v[1*(hc.M+1) + k] = 1.0;
        v[2*(hc.M+1) + k] = 1.0;
        v[3*(hc.M+1) + k] = 1.0;
    }

    /* Minimization */
    ConjugateGradient g2;
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

Parabolic1DControl335::Parabolic1DControl335()
{
    alpha = 1.0;
    t0 = 0.0;
    t1 = 1.0;
    x10 = 0.0;
    x11 = 1.0;
    x20 = 0.0;
    x21 = 1.0;
    h1 = 0.01;
    h2 = 0.01;
    ht = 0.01;
    a1 = a2 = 1.0;

    N1 = (unsigned int)(ceil(x11-x10)/h1);
    N2 = (unsigned int)(ceil(x21-x20)/h2);
    M  = (unsigned int)(ceil(t1-t0)/ht);
    L  = 4;

    double sgm1 = 5.0*h1;
    double sgm2 = 5.0*h2;
    gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    gause_b = 2.0*sgm1*sgm2;

    //initilize
    E.resize(2*L);
    E[0] = 0.20;
    E[1] = 0.20;
    E[2] = 0.20;
    E[3] = 0.80;
    E[4] = 0.80;
    E[5] = 0.20;
    E[6] = 0.80;
    E[7] = 0.80;

    U.resize(N2+1); for (unsigned int j=0; j<=N2; j++) U[j].resize(N1+1);
    for (unsigned int j=0; j<=N2; j++)
    {
        for (unsigned int i=0; i<=N1; i++)
        {
            U[j][i] = 3.0;
        }
    }

    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
    IPrinter::printMatrix(U, 10, 10);
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");

    FILE* file = fopen("example335.txt", "w");
    IPrinter::printMatrix(U, N1, N2, NULL, file);
    fclose(file);
}

double Parabolic1DControl335::fx(const DoubleVector& v)
{
    pv = &v;
    DoubleMatrix u;
    IParabolicEquation2D::caluclateMVD(u, h1, h2, ht, N1, N2, M, a1, a2);

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

double Parabolic1DControl335::norm(const DoubleVector& v) const
{
    double nrm = 0.0;
    for (unsigned int k=0; k<=M; k++)
    {
        double betta = 1.0;
        if (k==0 || k==M) betta = 0.5;
        nrm += betta*(v[0*(M+1)+k])*(v[0*(M+1)+k]);
        nrm += betta*(v[1*(M+1)+k])*(v[1*(M+1)+k]);
        nrm += betta*(v[2*(M+1)+k])*(v[2*(M+1)+k]);
        nrm += betta*(v[3*(M+1)+k])*(v[3*(M+1)+k]);
    }
    nrm = ht * nrm;
    return nrm;
}

void Parabolic1DControl335::gradient(const DoubleVector& v, DoubleVector& g)
{
    pv = &v;
    DoubleMatrix u;
    IParabolicEquation2D::caluclateMVD(u, h1, h2, ht, N1, N2, M, a1, a2);

    pu = &u;
    DoubleCube psi;
    IBackwardParabolicEquation2D::caluclateMVD(psi, h1, h2, ht, N1, N2, M, a1, a2);

    for (unsigned int k=0; k<=M; k++)
    {
        unsigned int i1 = (unsigned int)round(E[0]/h1);
        unsigned int j1 = (unsigned int)round(E[1]/h2);
        g[0*(M+1)+k] = -psi[k][j1][i1] + 2.0*alpha*(v[0*(M+1)+k]);

        unsigned int i2 = (unsigned int)round(E[2]/h1);
        unsigned int j2 = (unsigned int)round(E[3]/h2);
        g[1*(M+1)+k] = -psi[k][j2][i2] + 2.0*alpha*(v[1*(M+1)+k]);

        unsigned int i3 = (unsigned int)round(E[4]/h1);
        unsigned int j3 = (unsigned int)round(E[5]/h2);
        g[2*(M+1)+k] = -psi[k][j3][i3] + 2.0*alpha*(v[2*(M+1)+k]);

        unsigned int i4 = (unsigned int)round(E[4]/h1);
        unsigned int j4 = (unsigned int)round(E[5]/h2);
        g[3*(M+1)+k] = -psi[k][j4][i4] + 2.0*alpha*(v[3*(M+1)+k]);
    }

    psi.clear();
}

double Parabolic1DControl335::initial(unsigned int i, unsigned int j) const
{
    return 2.0;
}

double Parabolic1DControl335::boundary(unsigned int i, unsigned int j, unsigned int k) const
{
    return 2.0;
}

double Parabolic1DControl335::f(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    //    double t  = 0.5*k*ht;

    double sum = 0.0;
    unsigned int k1 = (k%2==0 ? k/2 : (k+1)/2);

    //printf("k: %d k1: %d j: %d i: %d\n", k, k1, j, i);

    double _v1 = (*pv)[0*(M+1)+k1];
    double _v2 = (*pv)[1*(M+1)+k1];
    double _v3 = (*pv)[2*(M+1)+k1];
    double _v4 = (*pv)[3*(M+1)+k1];

    //    unsigned int i1 = (unsigned int)round(E[0]/h1);
    //    unsigned int j1 = (unsigned int)round(E[1]/h2);
    //    if (i==i1 && j==j1) sum += (1.0/(h1*h2)) * _v1;

    //    unsigned int i2 = (unsigned int)round(E[2]/h1);
    //    unsigned int j2 = (unsigned int)round(E[3]/h2);
    //    if (i==i2 && j==j2) sum += (1.0/(h1*h2)) * _v2;

    //    unsigned int i3 = (unsigned int)round(E[4]/h1);
    //    unsigned int j3 = (unsigned int)round(E[5]/h2);
    //    if (i==i3 && j==j3) sum += (1.0/(h1*h2)) * _v3;

    double sgm1 = 5.0*h1;
    double sgm2 = 5.0*h2;
    double gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    double gause_b = 2.0*sgm1*sgm2;

    sum += _v1 * gause_a * exp(-((x1-E[0])*(x1-E[0]) + (x2-E[1])*(x2-E[1]))/gause_b);
    sum += _v2 * gause_a * exp(-((x1-E[2])*(x1-E[2]) + (x2-E[3])*(x2-E[3]))/gause_b);
    sum += _v3 * gause_a * exp(-((x1-E[4])*(x1-E[4]) + (x2-E[5])*(x2-E[5]))/gause_b);
    sum += _v4 * gause_a * exp(-((x1-E[6])*(x1-E[6]) + (x2-E[7])*(x2-E[7]))/gause_b);

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

double Parabolic1DControl335::binitial(unsigned int i, unsigned int j) const
{
    return -2.0*((*pu)[j][i] - U[j][i]);
}

double Parabolic1DControl335::bboundary(unsigned int i, unsigned int j, unsigned int k) const
{
    return 0.0;
}

double Parabolic1DControl335::bf(unsigned int i, unsigned int j, unsigned int k) const
{
    return 0.0;
}

void Parabolic1DControl335::print(unsigned int i, const DoubleVector &v, const DoubleVector &g, double alpha, RnFunction* fn) const
{
    Parabolic1DControl335 *hc = dynamic_cast<Parabolic1DControl335*>(fn);
    printf("J[%d]: %.16f\n", i, hc->fx(v));
    IPrinter::printVector(v, "v1:", 10, 0*(M+1), 0*(M+1) + M);
    IPrinter::printVector(v, "v2:", 10, 1*(M+1), 1*(M+1) + M);
    IPrinter::printVector(v, "v3:", 10, 2*(M+1), 2*(M+1) + M);

    const_cast<Parabolic1DControl335*>(this)->pv = &v;
    DoubleMatrix u;
    IParabolicEquation2D::caluclateMVD(u, h1, h2, ht, N1, N2, M, a1, a2);

    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
    IPrinter::printMatrix(u, 10, 10);
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");

    FILE* file = fopen("example335.txt", "w");
    IPrinter::printMatrix(u, N1, N2, NULL, file);
    fclose(file);
}

void Parabolic1DControl335::project(DoubleVector &e, int index)
{
    //    if (index<6)
    //    {
    //        if (e[index] > 1.0) e[index] = 1.0;
    //        if (e[index] < 0.0) e[index] = 0.0;
    //    }

    // if (e[index] <= 10.0) e[index] = 10.0;
    //  if (e[index] >= 103.5) e[index] = 103.5;

}
