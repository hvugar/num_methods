#include "discreteheat.h"

void DiscreteHeat::main()
{
    DiscreteHeat dh;

    DoubleVector f0((dh.M+1)*(dh.N+1));
    for (unsigned int k=0; k<f0.size(); k++)
    {
        //unsigned int j = k / (dh.M+1);
        //unsigned int i = k % (dh.N+1);
        //double t = j * dh.ht;
        //f0[k] = 2.0*t - 2.0*dh.a;
        f0[k] = 0.0;
    }

    //Printer::printAsMatrix(f0, dh.M, dh.N);
    //printf("----\n");

    ConjugateGradient g2;
    g2.setGradient(&dh);
    g2.setFunction(&dh);
    g2.setEpsilon1(0.0001);
    g2.setEpsilon2(0.0001);
    g2.setEpsilon3(0.0001);
    g2.setR1MinimizeEpsilon(0.1, 0.001);
    g2.setPrinter(&dh);
    g2.calculate(f0);

    IPrinter::printAsMatrix(f0, dh.M, dh.N);
}

DiscreteHeat::DiscreteHeat()
{
    t0 = 0.0; t1 = 1.0;
    x0 = 0.0; x1 = 1.0;
    N = 100;
    M = 100;
    hx = (t1-t0)/N;
    ht = (x1-x0)/M;
    a = 1.0;
    U.resize(N+1);
    for (unsigned int i=0; i<U.size(); i++)
    {
        double x = i*hx;
        U[i] = x*x + 1.0;
    }
}

double DiscreteHeat::fx(const DoubleVector& f)
{
    double sum  = 0.0;

    pf = &f;
    DoubleVector u;
    IParabolicEquation::calculateU(u, hx, ht, N, M, a);

    double alpha;
    for (unsigned int i=0; i<u.size(); i++)
    {
        alpha = 1.0;
        if (i==0 || i==u.size()-1) alpha = 0.5;
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

            double f1 = (f[j*(N+1)+i] - fxt(i*hx, j*ht));
            norm += alpha*f1*f1;
        }
    }
    norm = hx*ht*norm;

    return sum+norm;
}

void DiscreteHeat::gradient(const DoubleVector& f, DoubleVector& g)
{
    pf = &f;
    DoubleVector u;
    IParabolicEquation::calculateU(u, hx, ht, N, M, a);

    DoubleMatrix psi;
    calculateP(f, u, psi, g);


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

            unsigned int k = j*(N+1)+i;
            g[k] = -psi[j][i] + 2.0*alpha*(f[k]-fxt(i*hx, j*ht));
        }
    }
}

double DiscreteHeat::initial(unsigned int i) const
{
    double x = hx*i;
    return x*x;
}

double DiscreteHeat::boundary(Boundary type, unsigned int j) const
{
    double t = j*ht;
    if (type == Left) return t*t;
    if (type == Right) return t*t+1.0;
    return 0.0;
}

double DiscreteHeat::f(unsigned int i, unsigned int j) const
{
    return (*pf)[j*(N+1)+i];
}

double DiscreteHeat::fxt(double x, double t)
{
    C_UNUSED(x);
    return 2.0*t - 2.0*a;
}

void DiscreteHeat::calculateP(const DoubleVector &f, const DoubleVector &u, DoubleMatrix &psi, DoubleVector &g)
{
    C_UNUSED(f);
    C_UNUSED(g);

    double lamda = -(a*ht)/(hx*hx);
    double k = 1.0-2.0*lamda;

    psi.clear();
    psi.resize(M+1, N+1);

    DoubleVector a1(N-1);
    DoubleVector b1(N-1);
    DoubleVector c1(N-1);
    DoubleVector d1(N-1);
    DoubleVector x1(N-1);

    for (unsigned int m=0; m<=M; m++)
    {
        unsigned int j = M-m;

        if (j==M)
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                a1[i-1] = lamda;
                b1[i-1] = k;
                c1[i-1] = lamda;
                d1[i-1] = -2.0*hx*(u[i]-U[i]);
            }
            a1[0]   = 0.0;
            c1[N-2] = 0.0;
            tomasAlgorithm(a1.data(), b1.data(), c1.data(), d1.data(), x1.data(), x1.size());

            for (unsigned int i=1; i<=N-1; i++)
            {
                psi[j][i] = x1[i-1];
            }
            psi[j][0] = -lamda*psi[j][1]   -2.0*hx*0.5*(u[0]-U[0]);
            psi[j][N] = -lamda*psi[j][N-1] -2.0*hx*0.5*(u[N]-U[N]);
        }
        else if (j==0)
        {
            psi[j][0] = -lamda*psi[j][1];
            psi[j][N] = -lamda*psi[j][N-1];

            for (unsigned int i=1; i<=N-1; i++)
            {
                psi[j][i] = psi[j+1][i];
            }
        }
        else
        {
            for (unsigned int i=1; i<=N-1; i++)
            {
                a1[i-1] = lamda;
                b1[i-1] = k;
                c1[i-1] = lamda;
                d1[i-1] = +psi[j+1][i];
            }
            a1[0]   = 0.0;
            c1[N-2] = 0.0;
            tomasAlgorithm(a1.data(), b1.data(), c1.data(), d1.data(), x1.data(), x1.size());

            for (unsigned int i=1; i<=N-1; i++)
            {
                psi[j][i] = x1[i-1];
            }
            psi[j][0] = -lamda*psi[j][1];
            psi[j][N] = -lamda*psi[j][N-1];
        }
    }
    a1.clear();
    b1.clear();
    c1.clear();
    d1.clear();
    x1.clear();
}

void DiscreteHeat::print(unsigned int i, const DoubleVector &f, const DoubleVector &g, double alpha, RnFunction *fn) const
{
    C_UNUSED(g);
    C_UNUSED(alpha);
    printf("J[%d]: %.12f\n", i, fn->fx(f));
}
