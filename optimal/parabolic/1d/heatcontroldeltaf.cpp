#include "heatcontroldeltaf.h"

void HeatControlDeltaF::main(int argc, char ** argv)
{
    C_UNUSED(argc);
    C_UNUSED(argv);

    /* Function */
    HeatControlDeltaF hc;

    DoubleVector f0((hc.M+1)*hc.L);
    for (unsigned int i=0; i<f0.size(); i++)
    {
        //int j = i/(hc.N+1);
        //f0[i] = 2.0*j*hc.ht - 2.0;
        f0[i] = 20.0;
    }

    /* Minimization */
    ConjugateGradient g2;
    g2.setGradient(&hc);
    g2.setFunction(&hc);
    g2.setEpsilon1(0.00001);
    g2.setEpsilon2(0.00001);
    g2.setEpsilon3(0.00001);
    g2.setR1MinimizeEpsilon(0.1, 0.00001);
    g2.setPrinter(&hc);
    g2.setNormalize(true);
    g2.calculate(f0);

    IPrinter::printVector(f0, "f:");
}

HeatControlDeltaF::HeatControlDeltaF()
{
    this->t0 = 0.0;
    this->t1 = 1.0;
    this->x0 = 0.0;
    this->x1 = 1.0;
    this->a  = 1.0;

    this->N = 100;
    this->M = 100;
    this->hx = 0.01;
    this->ht = 0.01;

    this->L = 1;
    this->e = 0.2;
    this->E = 20;

    // initialize U
    DoubleVector f((M+1)*L);
    for (unsigned int j=0; j<=M; j++) f[j] = f1(j*ht);
    U.resize(N+1);
    pf = &f;
    IParabolicEquation::calculateU(U, hx, ht, N, M);
    //FILE *file = fopen("heat.txt", "w");
    //IPrinter::printVector(U,NULL,N,0,0,file);
    //fclose(file);
}

double HeatControlDeltaF::fx(const DoubleVector &f) const
{
    const_cast<HeatControlDeltaF*>(this)->pf = &f;
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

    double norm = 0.0;
    for (unsigned int j=0; j<=M; j++)
    {
        for (unsigned int i=0; i<=N; i++)
        {
            double alpha = 1.0;
            if (i==0 || i==N || j==0 || j==M) alpha = 0.5;
            if (i==0 && j==0) alpha = 0.25;
            if (i==0 && j==M) alpha = 0.25;
            if (i==N && j==0) alpha = 0.25;
            if (i==N && j==M) alpha = 0.25;
            norm += alpha*(f[j] - f1(j*ht))*(f[j] - f1(j*ht));
        }
    }
    norm = hx*ht*norm;

    return sum + norm;
}

void HeatControlDeltaF::gradient(const DoubleVector &f, DoubleVector &g)
{
    pf = &f;
    DoubleVector u;
    IParabolicEquation::calculateU(u, hx, ht, N, M, a);

    pu = &u;
    DoubleMatrix psi;
    IBackwardParabolicEquation::calculateU(psi, hx, ht, N, M, a);

    for (unsigned int j=0; j<=M; j++)
    {
//        g[j] = 0.0;
//        for (unsigned int i=0; i<=N; i++)
//        {
//            if (fabs(E[0]-i*hx)<hx) g[j] += psi[j][i] * ((hx-fabs(E[0]-i*hx))/hx);
//        }
//        g[j] = -g[j] + 2.0*(f[j]-f1(j*ht));
        g[j] = -psi[j][E]+2.0*(f[j]-f1(j*ht));
    }
//    IGradient::Gradient(this, 0.00001, f, g);
}

double HeatControlDeltaF::initial(unsigned int i) const
{
    double x = i*hx;
    return u(x, t0);
}

double HeatControlDeltaF::boundary(Boundary type, unsigned int j) const
{
    double t = j*ht;
    if (type == Left) return u(x0, t);
    if (type == Right) return u(x1, t);
    return 0.0;
}

double HeatControlDeltaF::f(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);

//    double x = i*hx;
    double sum = 0.0;
    if (i==E) sum = (1.0/hx) * (*pf)[j];
//    if (fabs(x-E[0])<hx)
//    {
//        sum += (1.0/hx) * (*pf)[j] * ((hx-fabs(x-E[0]))/hx);
//    }
    return sum;
}

double HeatControlDeltaF::binitial(unsigned int i) const
{
    return -2.0 * ((*pu)[i] - U[i]);
}

double HeatControlDeltaF::bboundary(Boundary type, unsigned int j) const
{
    C_UNUSED(type);
    C_UNUSED(j);
    return 0.0;
}

double HeatControlDeltaF::bf(unsigned int i, unsigned int j) const
{
    C_UNUSED(i);
    C_UNUSED(j);
    return 0.0;
}

void HeatControlDeltaF::print(unsigned int i, const DoubleVector& f0, const DoubleVector &g, double fx) const
{
    C_UNUSED(g);
    printf("J[%d]: %.20f\n", i, const_cast<HeatControlDeltaF*>(this)->fx(f0));
}
