#include "example332.h"

void Parabolic1DControl332::main()
{
    /* Function */
    Parabolic1DControl332 pc;

    DoubleVector e0(pc.L);
    e0[0] = 0.15;
    e0[1] = 0.55;

    printf("J[%d]: %.16f ", 0, pc.fx(e0));
    printf("e:  %10.6f %10.6f\n", e0[0], e0[1]);

    /* Minimization */
    ConjugateGradient g2;
    g2.setGradient(&pc);
    g2.setFunction(&pc);
    g2.setEpsilon1(0.0001);
    g2.setEpsilon2(0.0001);
    g2.setEpsilon3(0.0001);
    g2.setR1MinimizeEpsilon(0.1, 0.00001);
    g2.setPrinter(&pc);
    g2.setProjection(&pc);
    g2.setNormalize(true);
    g2.calculate(e0);

    printf("J[%d]: %.16f ", 0, pc.fx(e0));
    printf("e:  %10.6f %10.6f\n", e0[0], e0[1]);
//    IGradient::Gradient(&pc, 0.0001, e0, gn);
//    //gn.L2Normalize();
//    printf("gn: %10.8f %10.8f\n", gn[0], gn[1]);
//    pc.gradient(e0, ga);
//    ga.L2Normalize();
//    printf("ga: %10.8f %10.8f\n", ga[0], ga[1]);
}

Parabolic1DControl332::Parabolic1DControl332()
{
    x0 = 0.0;
    x1 = 1.0;
    t0 = 0.0;
    t1 = 1.0;
    hx = 0.001;
    ht = 0.001;
    N = (unsigned int)(ceil(x1-x0)/hx);
    M = (unsigned int)(ceil(t1-t0)/ht);
    L = 2;
    a = 1.0;

    DoubleVector e;
    e.resize(L);
    e[0] = 0.25;
    e[1] = 0.65;
    pe = &e;
    IParabolicEquation::calculateU(U, hx, ht, N, M, a);
    IPrinter::printVector(U);
    puts("-----------------------");
    FILE *file = fopen("example332.txt", "w");
    IPrinter::printVector(U, NULL, N, 0, 0, file);
    fclose(file);
}

double Parabolic1DControl332::fx(const DoubleVector &e)
{
    pe = &e;
    DoubleVector u;
    IParabolicEquation::calculateU(u, hx, ht, N, M, a);

    double sum = 0.0;
    double alpha;
    for (unsigned int i=0; i<=N; i++)
    {
        alpha = 1.0;
        if (i==0 || i==N) alpha = 0.5;
        sum += alpha*(u[i] - U[i])*(u[i] - U[i]);
    }
    sum = hx*sum;

    return sum;
}


void Parabolic1DControl332::gradient(const DoubleVector &e, DoubleVector &g)
{
    pe = &e;
    DoubleVector u;
    IParabolicEquation::calculateU(u, hx, ht, N, M, a);

    pu = &u;
    DoubleMatrix psi;
    IBackwardParabolicEquation::calculateU(psi, hx, ht, N, M, a);

    g[0] = 0.0;
    g[1] = 0.0;
    for (unsigned int j=0; j<=M; j++)
    {
        double t = j*ht;
        double alpha = 1.0;
        if (j==0 || j==M) alpha = 0.5;

        unsigned int i1 = (unsigned int)round(e[0]/hx);
        unsigned int i2 = (unsigned int)round(e[1]/hx);

        if (j==0)
        {
            g[0] += alpha * v1(t) * (psi[j][1] - psi[j][0])/hx;
            g[1] += alpha * v2(t) * (psi[j][1] - psi[j][0])/hx;
        }
        else if (j==M)
        {
            g[0] += alpha * v1(t) * (psi[j][N] - psi[j][N-1])/hx;
            g[1] += alpha * v2(t) * (psi[j][N] - psi[j][N-1])/hx;
        }
        else
        {
            g[0] += alpha * v1(t) * (psi[j][i1+1] - psi[j][i1-1])/(2.0*hx);
            g[1] += alpha * v2(t) * (psi[j][i2+1] - psi[j][i2-1])/(2.0*hx);
        }
    }
    g[0] = -ht*g[0];
    g[1] = -ht*g[1];
}

void Parabolic1DControl332::print(unsigned int i, const DoubleVector &e, const DoubleVector &gr, double alpha, RnFunction *fn) const
{
    C_UNUSED(alpha);
    Parabolic1DControl332 *hc = dynamic_cast<Parabolic1DControl332*>(fn);
    printf("J[%d]: %.16f ", i, hc->fx(e));
    printf("e:  %10.6f %10.6f\n", e[0], e[1]);
}

void Parabolic1DControl332::project(DoubleVector &e, int i)
{
    if (e[i] < 0.0) e[i] = 0.0;
    if (e[i] > 1.0) e[i] = 1.0;
}

double Parabolic1DControl332::fi(unsigned int i) const
{
    double x = i*hx;
    return 0.5*x;
}

double Parabolic1DControl332::m1(unsigned int j) const
{
    double t = j*ht;
    return t*t;
}

double Parabolic1DControl332::m2(unsigned int j) const
{
    double t = j*ht;
    return t*t+0.5;
}

double Parabolic1DControl332::f(unsigned int i, unsigned int j) const
{
    const DoubleVector &e = *pe;
    double x = i*hx;
    double t = j*ht;
    double sgm = 3.0*hx;
    double gause_a = 1.0/(sqrt(2.0*M_PI)*sgm);
    double gause_b = 2.0*sgm*sgm;
    double sum = 0.0;
    sum += v1(t) * gause_a * exp(-((x-e[0])*(x-e[0]))/gause_b);
    sum += v2(t) * gause_a * exp(-((x-e[1])*(x-e[1]))/gause_b);
    return sum;
}

double Parabolic1DControl332::bfi(unsigned int i) const
{
    return -2.0 * ((*pu)[i] - U[i]);
}

double Parabolic1DControl332::bm1(unsigned int j) const
{
    C_UNUSED(j);
    return 0.0;
}

double Parabolic1DControl332::bm2(unsigned int j) const
{
    C_UNUSED(j);
    return 0.0;
}

double Parabolic1DControl332::bf(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

double Parabolic1DControl332::v1(double t) const
{
    return 16.0*t;
}

double Parabolic1DControl332::v2(double t) const
{
    return 18.0*t;
}

