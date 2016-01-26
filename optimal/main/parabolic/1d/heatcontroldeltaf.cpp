#include "heatcontroldeltaf.h"

void HeatControlDeltaF::main()
{
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
    SteepestDescentGradient g2;
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
    this->hx = (x1-x0)/N;
    this->ht = (t1-t0)/M;

    this->L = 1;
    this->e = 0.2;
    this->E = (unsigned int)ceil(0.2/hx);

    // initialize U
    DoubleVector f;
    f.resize((M+1)*L);
    for (unsigned int j=0; j<=M; j++) f[j] = f1(j*ht);
    U.resize(N+1);
    pf = &f;
    IParabolicEquation::calculateU(U,hx,ht,N,M);
    FILE *file = fopen("heat.txt", "w");
    IPrinter::printVector(U,NULL,N,0,0,file);
    fclose(file);
}

double HeatControlDeltaF::fx(const DoubleVector &f)
{
    pf = &f;
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

double HeatControlDeltaF::fi(unsigned int i) const
{
    double x = i*hx;
    return u(x, t0);
}

double HeatControlDeltaF::m1(unsigned int j) const
{
    double t = j*ht;
    return u(x0, t);
}

double HeatControlDeltaF::m2(unsigned int j) const
{
    double t = j*ht;
    return u(x1, t);
}

double HeatControlDeltaF::f(unsigned int i, unsigned int j) const
{
//    double x = i*hx;
    double sum = 0.0;
    if (i==E) sum = (1.0/hx) * (*pf)[j];
//    if (fabs(x-E[0])<hx)
//    {
//        sum += (1.0/hx) * (*pf)[j] * ((hx-fabs(x-E[0]))/hx);
//    }
    return sum;
}

double HeatControlDeltaF::bfi(unsigned int i) const
{
    return -2.0 * ((*pu)[i] - U[i]);;
}

double HeatControlDeltaF::bm1(unsigned int j) const
{
    return 0.0;
}

double HeatControlDeltaF::bm2(unsigned int j) const
{
    return 0.0;
}

double HeatControlDeltaF::bf(unsigned int i, unsigned int j) const
{
    return 0.0;
}

void HeatControlDeltaF::print(unsigned int i, const DoubleVector &f0, const DoubleVector &g, double a, RnFunction *f) const
{
    HeatControlDeltaF *hc = dynamic_cast<HeatControlDeltaF*>(f);
    printf("J[%d]: %.20f\n", i, hc->fx(f0));
    //    printf("Printing f-------------------------------\n");
    //    IPrinter::printAsMatrix(f0, M, N);
    //    printf("Printing g-------------------------------\n");
    //    IPrinter::printAsMatrix(g, M, N);
}
