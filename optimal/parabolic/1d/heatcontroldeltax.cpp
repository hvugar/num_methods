#include "heatcontroldeltax.h"

void HeatControlDeltaX::main(int argc, char ** argv)
{
    C_UNUSED(argc);
    C_UNUSED(argv);
    /* Function */
    HeatControlDeltaX hc;

    DoubleVector e(hc.L);
    e[0] = 0.40;
    e[1] = 0.70;
    e[2] = 0.30;

    /* Minimization */
    ConjugateGradient g2;
    g2.setGradient(&hc);
    g2.setFunction(&hc);
    g2.setEpsilon1(0.0001);
    g2.setEpsilon2(0.0001);
    g2.setEpsilon3(0.0001);
    g2.setR1MinimizeEpsilon(1.0, 0.0001);
    g2.setPrinter(&hc);
    g2.setProjection(&hc);
    g2.setNormalize(true);
    g2.calculate(e);

    DoubleVector gr1(e.size());
    hc.gradient(e, gr1);
    gr1.L2Normalize();

    DoubleVector gr2(e.size());
    IGradient::Gradient(&hc, 0.00001, e, gr2);
    gr2.L2Normalize();

    printf("J[%d]: %.16f\n", 0, hc.fx(e));
    printf("eo:  %12.8f %12.8f %12.8f\n", 0.40, 0.80, 0.20);
    printf("e1:  %12.8f %12.8f %12.8f\n", e[0], e[1], e[2]);
    printf("gr1: %12.8f %12.8f %12.8f\n", gr1[0], gr1[1], gr1[2]);
    printf("gr2: %12.8f %12.8f %12.8f\n", gr2[0], gr2[1], gr2[2]);
}

HeatControlDeltaX::HeatControlDeltaX()
{
    this->t0 = 0.0;
    this->t1 = 1.0;
    this->x0 = 0.0;
    this->x1 = 1.0;
    this->a  = 1.0;

    this->N  = 100;
    this->M  = 100;
    this->hx = 0.01;
    this->ht = 0.01;
    this->L  = 3;

    // initialize U
    DoubleVector e(L);
    e[0] = 0.4;
    e[1] = 0.8;
    e[2] = 0.2;
    pe = &e;
    IParabolicEquation::calculateU(U, hx, ht, N, M);
    //FILE *file = fopen("heat_e.txt", "w");
    //IPrinter::printVector(U, NULL, N, 0, 0, file);
    //fclose(file);
}

double HeatControlDeltaX::fx(const DoubleVector &e) const
{
    const_cast<HeatControlDeltaX*>(this)->pe = &e;
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

    double nrm = 0.0;
    return sum + nrm;
}

double HeatControlDeltaX::norm(const DoubleVector &v) const
{
    C_UNUSED(v);
    double norm = 0.0;
    //    for (unsigned int j=0; j<=M; j++)
    //    {
    //        for (unsigned int i=0; i<=N; i++)
    //        {
    //            double alpha = 1.0;
    //            if (i==0 || i==N || j==0 || j==M) alpha = 0.5;
    //            if (i==0 && j==0) alpha = 0.25;
    //            if (i==0 && j==M) alpha = 0.25;
    //            if (i==N && j==0) alpha = 0.25;
    //            if (i==N && j==M) alpha = 0.25;
    //            norm += alpha*(v[j] - f1(j*ht))*(v[j] - f1(j*ht));
    //        }
    //    }
    //    norm = hx*ht*norm;
    return norm;
}

void HeatControlDeltaX::gradient(const DoubleVector &e, DoubleVector &g)
{
    pe = &e;
    DoubleVector u;
    IParabolicEquation::calculateU(u, hx, ht, N, M, a);

    pu = &u;
    DoubleMatrix psi;
    IBackwardParabolicEquation::calculateU(psi, hx, ht, N, M, a);

    for (unsigned int j=0; j<=M; j++)
    {
        double alpha = 1.0;
        if (j==0 || j==M) alpha = 0.5;
        unsigned int i = (unsigned int)round(e[0]/hx);
        if (i==0)
            g[0] += alpha * v1(j*ht) * (psi[j][1] - psi[j][0])/(hx);
        else if (i==N)
            g[0] += alpha * v1(j*ht) * (psi[j][N] - psi[j][N-1])/(hx);
        else
            g[0] += alpha * v1(j*ht) * (psi[j][i+1] - psi[j][i-1])/(2.0*hx);
    }
    g[0] = -ht*g[0];

    for (unsigned int j=0; j<=M; j++)
    {
        double alpha = 1.0;
        if (j==0 || j==M) alpha = 0.5;
        unsigned int i = (unsigned int)round(e[1]/hx);
        if (i==0)
            g[1] += alpha * v2(j*ht) * (psi[j][1] - psi[j][0])/(hx);
        else if (i==N)
            g[1] += alpha * v2(j*ht) * (psi[j][N] - psi[j][N-1])/(hx);
        else
            g[1] += alpha * v2(j*ht) * (psi[j][i+1] - psi[j][i-1])/(2.0*hx);
    }
    g[1] = -ht*g[1];

    for (unsigned int j=0; j<=M; j++)
    {
        double alpha = 1.0;
        if (j==0 || j==M) alpha = 0.5;
        unsigned int i = (unsigned int)round(e[2]/hx);
        if (i==0)
            g[2] += alpha * v3(j*ht) * (psi[j][1] - psi[j][0])/(hx);
        else if (i==N)
            g[2] += alpha * v3(j*ht) * (psi[j][N] - psi[j][N-1])/(hx);
        else
            g[2] += alpha * v3(j*ht) * (psi[j][i+1] - psi[j][i-1])/(2.0*hx);
    }
    g[2] = -ht*g[2];
}

double HeatControlDeltaX::initial(unsigned int i) const
{
    double x = i*hx;
    return u(x, t0);
}

double HeatControlDeltaX::boundary(Boundary type, unsigned int j) const
{
    double t = j*ht;
    if (type == Left) return u(x0, t);
    if (type == Right) return u(x1, t);
    return 0.0;
}

double HeatControlDeltaX::f(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    double t = j*ht;
    const DoubleVector &e = *pe;
    double sum = 0.0;

    double sgm = 3.0*hx;
    double gause_a = 1.0/(sqrt(2.0*M_PI)*sgm);
    double gause_b = 2.0*sgm*sgm;

    sum += v1(t) * gause_a * exp(-((x-e[0])*(x-e[0]))/gause_b);// * h1*h2;
    sum += v2(t) * gause_a * exp(-((x-e[1])*(x-e[1]))/gause_b);// * h1*h2;
    sum += v3(t) * gause_a * exp(-((x-e[2])*(x-e[2]))/gause_b);// * h1*h2;

    //    if (fabs(x-e[0]) < (hx+0.000001))
    //    {
    //        sum += (1.0/hx) * g1(j*ht) * ((hx-fabs(x-e[0]))/hx);
    //    }

    //    if (fabs(x-e[1]) < (hx+0.000001))
    //    {
    //        sum += (1.0/hx) * g2(j*ht) * ((hx-fabs(x-e[1]))/hx);
    //    }

    return sum;
}

double HeatControlDeltaX::binitial(unsigned int i) const
{
    return -2.0 * ((*pu)[i] - U[i]);;
}

double HeatControlDeltaX::bboundary(Boundary type, unsigned int j) const
{
    C_UNUSED(type);
    C_UNUSED(j);
    return 0.0;
}

double HeatControlDeltaX::bf(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

void HeatControlDeltaX::print(unsigned int i, const DoubleVector& e, const DoubleVector &g, double fx, GradientMethod::MethodResult result) const
{
    C_UNUSED(g);
    printf("J[%d]: %.20f\n", i, const_cast<HeatControlDeltaX*>(this)->fx(e));
    DoubleVector g1 = g;
    g1.L2Normalize();
    printf("eo: %12.8f %12.8f %12.8f\n", 0.40, 0.80, 0.20);
    printf("e1: %12.8f %12.8f %12.8f\n", e[0], e[1], e[2]);
    printf("g1: %12.8f %12.8f %12.8f\n", g1[0], g1[1], g1[2]);
    puts("+------------------------------------------------------------------------------------------------------------------------------------------------------------------+");
}


void HeatControlDeltaX::project(DoubleVector &e, int index)
{
    C_UNUSED(index)
            for (unsigned int l=0; l<L; l++)
    {
        if (e[l] < x0) e[l] = x0+hx;
        if (e[l] > x1) e[l] = x1-hx;
    }
}
